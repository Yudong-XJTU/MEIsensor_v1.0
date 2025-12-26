#!/usr/bin/env python3
#
# MEIsensor
import copy
import gc
import logging
import math
import multiprocessing
import os
import threading
import time
from argparse import Namespace
from collections import deque
from dataclasses import dataclass
from typing import Optional, Union, Callable

import pysam

import combine
import preprocess
import postprocess
import SV
from region import Region
from result import Result, ErrorResult, CallResult, GenotypeResult, CombineResult


@dataclass
class Task:
    """
    A task is a generic unit of work sent to a child process to be worked on in parallel. Must be pickleable.
    """
    id: int
    sv_id: int
    contig: str
    start: int
    end: int
    config: Namespace
    assigned_process_id: Optional[int] = None
    lead_provider: preprocess.LeadProvider = None
    bam: object = None
    genotype_svs: list = None
    regions: list[Region] = None
    _logger = None
    result: Result = None

    def __str__(self):
        return f'Task #{self.id}'

    @property
    def logger(self) -> logging.Logger:
        if self._logger is None:
            self._logger = logging.getLogger(f'MEIsensor.progress')

        return self._logger

    @property
    def done(self) -> bool:
        return self.result is not None

    @property
    def success(self) -> bool:
        return self.done and not self.result.error

    def add_result(self, result: Result) -> None:
        self.result = result

    def execute(self, worker: 'MEIsensorWorker' = None) -> Optional[Result]:
        """
        Execute this Task, returning a Result object
        :param worker is the worker executing this task
        """
        raise NotImplemented

    def build_leadtab(self):
        assert (self.lead_provider is None)

        config = self.config

        if config.input_is_cram and config.reference is not None:
            self.bam = pysam.AlignmentFile(config.input, config.input_mode, require_index=True, reference_filename=config.reference)
        else:
            self.bam = pysam.AlignmentFile(config.input, config.input_mode, require_index=True)
        self.lead_provider = preprocess.LeadProvider(config, self.id * config.task_read_id_offset_mult)
        externals = self.lead_provider.build_leadtab(self.regions if self.regions else [Region(self.contig, self.start, self.end)], self.bam)
        return externals, self.lead_provider.read_count

    def call_candidates(self, keep_qc_fails, config):
        candidates = []
        for svtype in SV.TYPES:
            for svcluster in combine.resolve(svtype, self.lead_provider, config):
                for svcall in SV.call_from(svcluster, config, keep_qc_fails, self):
                    if config.dev_trace_read is not False:
                        cluster_has_read = False
                        for ld in svcluster.leads:
                            if ld.read_qname == config.dev_trace_read:
                                cluster_has_read = True
                        if cluster_has_read:
                            import copy
                            svcall_copy = copy.deepcopy(svcall)
                            svcall_copy.postprocess = None
                            print(f"[DEV_TRACE_READ] [3/4] [Task.call_candidates] Read {config.dev_trace_read} -> Cluster {svcluster.id} -> preliminary SVCall {svcall_copy}")
                    candidates.append(svcall)

        self.coverage_average_fwd, self.coverage_average_rev = postprocess.coverage(candidates, self.lead_provider, config)
        self.coverage_average_total = self.coverage_average_fwd + self.coverage_average_rev
        return candidates

    def finalize_candidates(self, candidates: list['SVCall'], keep_qc_fails, config):
        passed = []
        for svcall in candidates:
            svcall.qc = svcall.qc and postprocess.qc_sv(svcall, config)
            if not keep_qc_fails and not svcall.qc:
                continue
            svcall.qc = svcall.qc and postprocess.qc_sv_support(svcall, self.coverage_average_total, config)
            if not keep_qc_fails and not svcall.qc:
                continue

            postprocess.annotate_sv(svcall, config)

            svcall.qc = svcall.qc and postprocess.qc_sv_post_annotate(svcall, config)

            if config.dev_trace_read:
                cluster_has_read = False
                for ld in svcall.postprocess.cluster.leads:
                    if ld.read_qname == config.dev_trace_read:
                        cluster_has_read = True
                if cluster_has_read:
                    import copy
                    svcall_copy = copy.deepcopy(svcall)
                    svcall_copy.postprocess = None
                    print(f"[DEV_TRACE_READ] [4/4] [Task.finalize_candidates] Read {config.dev_trace_read} -> Cluster {svcall.postprocess.cluster.id} -> finalized SVCall, QC={svcall_copy.qc}: {svcall_copy}")

            if not keep_qc_fails and not svcall.qc:
                continue

            svcall.finalize()  # Remove internal information (not written to output) before sending to mainthread for VCF writing
            passed.append(svcall)
        return passed


class CallTask(Task):
    """
    """

    def execute(self, worker: 'MEIsensorWorker' = None) -> CallResult:
        config = self.config

        # if config.snf is not None or config.no_qc:
        if config.no_qc:
            qc = False
        else:
            qc = True

        _, read_count = self.build_leadtab()
        svcandidates = self.call_candidates(qc, config)
        svcalls = self.finalize_candidates(svcandidates, not qc, config)
        if not config.no_qc:
            svcalls = [s for s in svcalls if s.qc]

        if config.sort:
            svcalls = sorted(svcalls, key=lambda svcall: svcall.pos)

        from result import CallResult
        result = CallResult(self, svcalls, read_count)
        '''
        if config.snf is not None:  # and len(svcandidates):
            snf_filename = f"{config.snf}.tmp_{self.id}.snf"

            with open(snf_filename, "wb") as handle:
                snf_out = snf.SNFile(config, handle)
                for cand in svcandidates:
                    snf_out.store(cand)
                snf_out.annotate_block_coverages(self.lead_provider)
                snf_out.write_and_index()
                handle.close()
            result.snf_filename = snf_filename
            result.snf_index = snf_out.get_index()
            result.snf_total_length = snf_out.get_total_length()
            result.snf_candidate_count = len(svcandidates)
            result.has_snf = True
        '''
        result.coverage_average_total = self.coverage_average_total

        return result


class GenotypeTask(Task):
    def execute(self, worker: 'MEIsensorWorker' = None) -> Optional[GenotypeResult]:
        config = self.config

        qc = False
        _, read_count = self.build_leadtab()
        svcandidates = self.call_candidates(qc, config=config)
        svcalls = self.finalize_candidates(svcandidates, not qc, config=config)

        binsize = 5000
        binedge = int(binsize / 10)
        genotype_svs_svtypes_bins = {svtype: {} for svtype in SV.TYPES}
        for genotype_sv in self.genotype_svs:
            genotype_sv.genotype_match_sv = None
            genotype_sv.genotype_match_dist = math.inf

            if genotype_sv.svtype not in genotype_svs_svtypes_bins:
                logging.getLogger('MEIsensor').warning(f'Unsupported SVTYPE: {genotype_sv.svtype}')
                continue

            bins = [int(genotype_sv.pos / binsize) * binsize]
            if genotype_sv.pos % binsize < binedge:
                bins.append((int(genotype_sv.pos / binsize) - 1) * binsize)
            if genotype_sv.pos % binsize > binsize - binedge:
                bins.append((int(genotype_sv.pos / binsize) + 1) * binsize)

            for bin in bins:
                if bin not in genotype_svs_svtypes_bins[genotype_sv.svtype]:
                    genotype_svs_svtypes_bins[genotype_sv.svtype][bin] = []
                genotype_svs_svtypes_bins[genotype_sv.svtype][bin].append(genotype_sv)

        for cand in svcandidates:
            bin = int(cand.pos / binsize) * binsize
            if bin not in genotype_svs_svtypes_bins[cand.svtype]:
                continue
            for genotype_sv in genotype_svs_svtypes_bins[cand.svtype][bin]:
                dist = abs(genotype_sv.pos - cand.pos) + abs(abs(genotype_sv.svlen) - abs(cand.svlen))
                minlen = float(min(abs(genotype_sv.svlen), abs(cand.svlen)))
                if minlen > 0 and dist < genotype_sv.genotype_match_dist and dist <= config.combine_match * math.sqrt(
                        minlen) and dist <= config.combine_match_max:
                    genotype_sv.genotype_match_sv = cand
                    genotype_sv.genotype_match_dist = dist

        postprocess.coverage(self.genotype_svs, self.lead_provider, config)

        # Determine genotypes for unmatched input SVs
        for svcall in self.genotype_svs:
            coverage_list = [svcall.coverage_start, svcall.coverage_center, svcall.coverage_end]
            coverage_list = [c for c in coverage_list if c is not None]
            if len(coverage_list) == 0:
                return
            coverage = round(sum(coverage_list) / len(coverage_list))
            svcall.genotypes = {}
            if coverage > 0:
                svcall.genotypes[0] = (0, 0, 0, coverage, 0, (None, None))
            else:
                svcall.genotypes[0] = config.genotype_none

        from result import GenotypeResult
        return GenotypeResult(self, self.genotype_svs, read_count)


class ShutdownTask:
    id = None

    def __str__(self):
        return 'Shutdown Request'

    def execute(self, *args, **kwargs) -> Result:
        raise MEIsensorWorker.Shutdown


class MEIsensorWorker:
    """
    Handle for a worker process. Since we're forking, this class will be available in
    both the parent and the worker processes.
    """
    id: int  # sequential ID of this worker, starting with 0 for the first
    externals: list = None
    recycle: bool = False
    running = True
    pid: int = None
    # Event to shut down heartbeat threads
    _shutdown: threading.Event
    _heartbeat: float = 0  # last heartbeat received
    HEARTBEAT_INTERVAL = 3  # in seconds
    HEARTBEAT_TIMEOUT = 10  # in seconds

    class Shutdown(Exception):
        """
        Indicates this worker process should shut down
        """

    def __init__(self, process_id: int, config: Namespace, tasks: deque[Task], recycle_hint: Union[bool, Callable] = None):
        self.id = process_id
        self.config = config
        self.tasks = tasks
        self.task = None
        self.finished_tasks = []
        self.recycle = recycle_hint

        self.pipe_main, self.pipe_worker = multiprocessing.Pipe()
        self.heartbeat_main, self.heartbeat_worker = multiprocessing.Pipe()

        self.process = multiprocessing.Process(
            target=self.run_worker,
            daemon=True
        )

        self._logger = logging.getLogger('MEIsensor.worker')

    def __str__(self):
        return f'Worker {self.id} @ process {self.pid}'

    def start(self) -> None:
        self._logger.info(f'Starting worker {self.id}')
        self.running = True
        self.process.start()
        self._heartbeat = time.monotonic()

    def maybe_recycle(self):
        """
        Recycle this worker if that has been requested
        """
        recycle = self.recycle(self.id, self.process.pid) if callable(self.recycle) else self.recycle

        if recycle:
            self._logger.info(f'Recycling worker {self.id}')
            # Shut down current worker process
            self.pipe_main.send(ShutdownTask())
            self.process.join(2)
            # Start new one
            self.process = multiprocessing.Process(
                target=self.run_worker,
                daemon=True
            )
            self.process.start()
            self._heartbeat = time.monotonic()

    def run_parent(self) -> bool:
        """
        Worker thread, running in parent process
        """
        try:
            if self.task is None:
                # we are not working on something...
                if len(self.tasks) > 0:
                    # ...but there is more work to be done
                    self.maybe_recycle()

                    try:
                        self.task = self.tasks.popleft()
                    except IndexError:
                        # another worker may have taken the last task
                        self._logger.debug(f'No more tasks to do for {self.id}')
                    else:
                        self.pipe_main.send(self.task)
                        self._logger.info(f'Dispatched task #{self.task.id} to worker {self.id} ({len(self.tasks)}  tasks left)')
                else:
                    # ...and no more work available, so we shut down this worker
                    self._logger.info(f'Worker {self.id} shutting down...')
                    self.pipe_main.send(ShutdownTask())
                    self.running = False
            else:
                if self.pipe_main.poll(0.01):
                    self._logger.debug(f'Worker {self.id} got result for task {self.task.id}...')
                    result: Result = self.pipe_main.recv()

                    if result.error:
                        self._logger.error(f'Worker {self.id} received error: {result}')
                    else:
                        self._logger.info(f'Worker {self.id} got result for task #{result.task_id}')

                    self.task.add_result(result)
                    self.finished_tasks.append(self.task)
                    self.task = None

                if self.heartbeat_main.poll():
                    hb = self.heartbeat_main.recv()
                    self._heartbeat = time.monotonic()
                    # self._logger.debug(f'Worker {self.id} got heartbeat #{hb}')

                if self._heartbeat < time.monotonic() - self.HEARTBEAT_TIMEOUT:
                    self._logger.debug(f'Worker {self.id} missed heartbeat!')
                    try:
                        self.process.join(0.2)  # try collecting process remains...
                    except:  # noqa
                        ...
                    if self.process.exitcode is not None:
                        # if we got an exitcode, the process really was killed
                        self._logger.warning(f'Worker {self.id} found dead!')
                        if self.task:  # if we were working on a task, requeue it to have it picked up by another worker...
                            self.tasks.appendleft(self.task)
                        self.running = False  # ...and shut down
        except:
            self._logger.exception(f'Unhandled error in worker {self.id}. This may result in an orphened worker process.')
            try:
                self.process.kill()
            except:
                ...

        return self.running

    def finalize(self):
        self.process.join(10)

        if self.process.exitcode is None:
            self._logger.warning(f'Worker {self.id} refused to shut down gracefully, killing it.')
            self.process.kill()
            self.process.join(2)
        self._logger.info(f'Worker {self.id} done (code {self.process.exitcode}).')

    def run_worker(self):
        """
        Entry point/main loop for the worker process
        """
        self.pid = os.getpid()
        self._shutdown = threading.Event()

        t = threading.Thread(target=self.run_worker_heartbeats, daemon=True)
        t.start()

        while self.running:
            self._logger.debug(f'Worker {self.id} ({self.pid}) waiting for tasks...')

            task = self.pipe_worker.recv()

            self._logger.debug(f'Worker {self.id} got task {task}')

            try:
                result = task.execute(self)
            except self.Shutdown:
                self.running = False
                self._shutdown.set()
            except Exception as e:
                self._logger.exception(msg := f'Error in worker process while executing {task}')
                self.pipe_worker.send(ErrorResult(msg))
            else:
                self._logger.debug(f'Worker {self.id} finished executing {task}, sending back result...')

                if result is not None:
                    self.pipe_worker.send(result)

            del task
            gc.collect()

        t.join(1.0)

    def run_worker_heartbeats(self):
        hb = 0
        while self.running:
            hb += 1
            self.heartbeat_worker.send(hb)
            self._shutdown.wait(self.HEARTBEAT_INTERVAL)


def execute_task(task: Task):
    logging.getLogger('MEIsensor.parallel').info(f'Working on {task}')
    return task.execute()


class MEIsensorParentWorker(MEIsensorWorker):
    """
    A worker class without multiprocessing, i.e. running in the main process. Used for profiling.
    """
    id: int = 0

    def __init__(self, config: Namespace, tasks: list[Task], **kwargs):  # noqa
        self.tasks = tasks
        self.task = None
        self.config = config
        self.finished_tasks: list[Task] = []
        self._log = logging.getLogger('MEIsensor.worker')
        self._log.info(f'Using parent worker')

    def start(self) -> None:
        ...

    def run_parent(self) -> bool:
        count = len(self.tasks)
        for i, task in enumerate(self.tasks):
            self._log.info(f'Executing {task} ({i+1}/{count})')
            result = task.execute(self)
            task.add_result(result)
            self.finished_tasks.append(task)
        self._log.info(f'All tasks done.')

        return False

    def finalize(self):
        ...