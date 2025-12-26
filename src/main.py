#!/usr/bin/env python3
#
# MEIsensor
import logging
import logging.config
import multiprocessing
from collections import deque
from typing import Optional

from UTILS.resmon import ResourceMonitor

import sys

if not sys.version_info >= (3, 10):
    print(f"Error: MEIsensor must be run with Python version 3.10 or above (detected Python version: {sys.version_info.major}.{sys.version_info.minor}). Exiting")
    exit(1)

import math
import time
import os

import pysam

from config import Config
import vcf
import parallel_task
import utils

# TODO: Dev/Debugging only - Remove for prod
DEV_MONITOR_MEM = False

if DEV_MONITOR_MEM:
    try:
        import psutil
    except ImportError:
        logging.getLogger('MEIsensor.memory').warning('psutil not available.')

        def dbg_get_total_memory_usage_MB():
            pass
    else:
        logging.getLogger('MEIsensor.memory').info('Watching memory')


        def dbg_get_total_memory_usage_MB():
            total = 0
            n = 0
            proc = psutil.Process(os.getpid())
            for child in proc.children(recursive=True):
                total += child.memory_info().rss
                n += 1
            total += proc.memory_info().rss
            return total / (1000.0 * 1000.0)


# """
# END:TODO


def Main(processes: list[parallel_task.MEIsensorWorker]):
    # Determine MEIsensor run mode
    '''
    config = Config(
        # input='sample/HG00513/alignment/HG00513.hifi.minimap2.sorted.bam',
        input='../sample/test_alignment/HG00512.hifi.minimap2.sorted.matched_reads.bam',
        reference='../sample/reference/GRCh38.fasta',
        vcf='../sample/test_alignment/test_result/test.MEI.vcf',
    )
    '''
    config = Config()
    input_ext = config.input.split(".")[-1].lower()
    if input_ext not in ["bam"]:
        sys.exit(f"Error: Input file must be a .bam file. Provided file: {config.input}")

    # needed for running on osx
    if sys.platform == "darwin":
        multiprocessing.set_start_method("fork")

    if "bam" in input_ext or "cram" in input_ext:
        if input_ext.count("bam") + input_ext.count("cram") > 1:
            utils.fatal_error_main(f"Please specify max 1 .bam//.cram file as input (got {input_ext.count('bam')})")
        config.input = config.input

        if config.genotype_vcf is not None:
            config.mode = "genotype_vcf"
        else:
            config.mode = "call_sample"

        config.input_is_cram = False
        if "bam" in input_ext:
            config.input_mode = r"rb"
        elif "cram" in input_ext:
            config.input_mode = r"rc"
            config.input_is_cram = True

    log = logging.getLogger('MEIsensor.main')
    if config.dev_debug_log:
        logging.getLogger().setLevel(logging.DEBUG)
    if config.dev_progress_log:
        logging.getLogger('MEIsensor.progress').setLevel(logging.INFO)

    if config.mode == "call_sample":
        if config.sample_id is None:
            # config.sample_id,_=os.path.splitext(os.path.basename(config.input))
            config.sample_ids_vcf = [(0, "SAMPLE")]
        else:
            config.sample_ids_vcf = [(0, config.sample_id)]
    log.info(f"  Run Mode: {config.mode}")
    log.info("==============================")

    rkwargs = {}  # result kwargs

    monitor = ResourceMonitor(config)
    if monitor and monitor.filename is not None:
        logging.getLogger('MEIsensor.resources').info(f'Logging memory usage to {monitor.filename}')


    #
    # call_sample/genotype_vcf: Open .bam file for single calling / .snf creation
    #
    contig_tandem_repeats = {}
    if config.mode == "call_sample" or config.mode == "genotype_vcf":
        log.info(f"Opening for reading: {config.input}")
        bam_in = pysam.AlignmentFile(config.input, config.input_mode)
        try:
            has_index = bam_in.check_index()
            if not has_index:
                raise ValueError
        except ValueError:
            utils.fatal_error_main(f"Unable to load index for input file '{config.input}'. Please verify that your input file is sorted + indexed and that the index .bai file is valid and in the right location.")

    #
    # genotype_vcf: Read SVs from VCF to be genotyped
    #
    if config.mode == "genotype_vcf":
        path, ext = os.path.splitext(config.genotype_vcf)
        ext = ext.lower()
        if ext == ".gz":
            vcf_in_handle = pysam.BGZFile(config.genotype_vcf, "rb")
        elif ext == ".vcf":
            vcf_in_handle = open(config.genotype_vcf, "r")
        else:
            utils.fatal_error_main("Expected a .vcf or .vcf.gz file for genotyping using --genotype-vcf")
        vcf_in = vcf.VCF(config, vcf_in_handle)

        genotype_lineindex_order = []
        genotype_lineindex_svs = {}
        genotype_contig_svs = {}
        for svcall in vcf_in.read_svs_iter():
            if svcall.contig not in genotype_contig_svs:
                genotype_contig_svs[svcall.contig] = []
            assert (svcall.raw_vcf_line_index not in genotype_lineindex_svs)
            genotype_lineindex_order.append(svcall.raw_vcf_line_index)
            genotype_lineindex_svs[svcall.raw_vcf_line_index] = svcall
            genotype_contig_svs[svcall.contig].append(svcall)
        rkwargs['genotype_lineindex_order'] = genotype_lineindex_order
        log.info(f"Opening for reading: {config.genotype_vcf} (read {len(genotype_lineindex_svs)} SVs to be genotyped)")

    #
    # Open output files
    #
    vcf_out = None
    if config.vcf is not None:

        vcf_output_info = []
        if config.sort:
            vcf_output_info.append("sorted")
        if config.vcf_output_bgz:
            vcf_output_info.append("bgzipped")
            vcf_output_info.append("tabix-indexed")

        if len(vcf_output_info) == 0:
            vcf_output_info_str = ""
        else:
            vcf_output_info_str = f"({', '.join(vcf_output_info)})"

        if os.path.exists(config.vcf) and not config.allow_overwrite:
            utils.fatal_error_main(f"Output file '{config.vcf}' already exists! Use --allow-overwrite to ignore this check and overwrite.")

        if config.vcf_output_bgz:
            if not config.sort:
                utils.fatal_error_main(".gz (bgzip) output is only supported with sorting enabled")
            vcf_handle = pysam.BGZFile(config.vcf, "w")
        else:
            vcf_handle = open(config.vcf, "w")

        vcf_out = vcf.VCF(config, vcf_handle)

        #if config.mode == "call_sample" or config.mode == "combine":
        if config.mode == "call_sample":
            if config.reference is not None:
                log.info(f"Opening for reading: {config.reference}")
            vcf_out.open_reference()

        log.info(f"Opening for writing: {config.vcf} {vcf_output_info_str}")
    #
    # Plan multiprocessing tasks
    #
    task_id = 0
    tasks = deque()
    contigs = []
    contig_tasks_intervals = {}

    if config.mode == "call_sample" or config.mode == "genotype_vcf":
        #
        # Process .bam header
        #
        task_classes = {
            'call_sample': parallel_task.CallTask,
            'genotype_vcf': parallel_task.GenotypeTask,
        }

        total_mapped = bam_in.mapped
        if (config.threads == 1 and not config.low_memory) or config.task_count_multiplier == 0:
            task_max_reads = total_mapped
        else:
            task_max_reads = max(1, math.floor(total_mapped / (config.threads * config.task_count_multiplier)))

        if total_mapped == 0:
            # Total mapped returns 0 for CRAM files
            config.task_read_id_offset_mult = 10 ** 9
        else:
            # BAM file
            config.task_read_id_offset_mult = 10 ** math.ceil(math.log(total_mapped) + 1)

        contig_lengths = []
        contigs_with_tr_annotations = 0
        for contig in bam_in.get_index_statistics():
            if task_max_reads == 0:
                task_count = 1
            else:
                task_count = max(1, math.ceil(contig.mapped / float(task_max_reads)))
            contig_str = str(contig.contig)

            if config.contig and contig_str not in config.contig:
                continue

            if config.regions_by_contig and contig_str not in config.regions_by_contig:
                continue

            contigs.append(contig_str)
            contig_length = bam_in.get_reference_length(contig_str)
            contig_lengths.append((contig_str, contig_length))
            task_length = math.floor(contig_length / float(task_count))
            contigs_with_tr_annotations += int(contig_str in contig_tandem_repeats)
            startpos = 0

            while startpos < contig_length - 1:
                endpos = min(contig_length - 1, startpos + task_length)
                if config.genotype_vcf is not None:
                    if contig_str in genotype_contig_svs:
                        genotype_svs = [target_sv for target_sv in genotype_contig_svs[contig_str] if target_sv.pos >= startpos and target_sv.pos < endpos]
                    else:
                        genotype_svs = []
                else:
                    genotype_svs = None

                task = task_classes[config.mode](
                    id=task_id,
                    contig=contig_str,
                    start=startpos,
                    end=endpos,
                    assigned_process_id=None,
                    # tandem_repeats=contig_tandem_repeats[contig_str] if contig_str in contig_tandem_repeats else None,
                    genotype_svs=genotype_svs,
                    sv_id=0,
                    config=config,
                    regions=config.regions_by_contig.get(contig_str),
                )
                tasks.append(task)
                if contig_str not in contig_tasks_intervals:
                    contig_tasks_intervals[contig_str] = []
                contig_tasks_intervals[contig_str].append((task.start, task.end, task))
                startpos += task_length
                task_id += 1
        config.contig_lengths = contig_lengths

    if config.mode != "genotype_vcf" and config.vcf is not None:
        vcf_out.write_header(contig_lengths)
    elif config.mode == "genotype_vcf":
        vcf_out.rewrite_header_genotype(vcf_in.header_str)

    #
    # Start workers
    #
    if config.threads:
        for pnum in range(config.threads):
            processes.append(parallel_task.MEIsensorWorker(process_id=pnum, config=config, tasks=tasks, recycle_hint=monitor))
    else:
        processes.append(parallel_task.MEIsensorParentWorker(config=config, tasks=tasks))

    if config.vcf is not None and config.sort:
        task_id_calls = {}

    log.info("")
    if config.mode == "call_sample" or config.mode == "genotype_vcf":
        if config.input_is_cram:
            # CRAM file
            log.info(f"Analyzing alignments... (progress display disabled for CRAM input)")
        else:
            log.info(f"Analyzing {total_mapped} alignments total...")
    log.info("")

    #
    # Distribute analysis tasks to workers and collect results
    #
    analysis_start_time = time.monotonic()

    for p in processes:
        p.start()

    finished_tasks: list[task.Task] = []

    while any([p.run_parent() for p in processes if p.running]):
        time.sleep(0.1)

    for p in processes:
        p.finalize()
        finished_tasks.extend(p.finished_tasks)

    finished_tasks.sort(key=lambda task: task.id)

    for t in finished_tasks:
        t.result.emit(vcf_out=vcf_out, **rkwargs)

    if config.vcf is not None:
        vcf_out.close()
        if config.vcf_output_bgz:
            vcf_index_start_time = time.time()
            log.info(f"Generating index for {config.vcf}...")
            try:
                pysam.tabix_index(config.vcf, preset="vcf", force=True)
            except:
                log.exception(f'Error indexing VCF.')
            else:
                log.info(f"Indexing VCF output took {time.time() - vcf_index_start_time:.2f}s.")

    log.info(f"Took {time.monotonic() - analysis_start_time:.2f}s.")
    log.info("")

    if (config.mode == "call_sample" or config.mode == "combine") and config.vcf is not None:
        log.info(f"Wrote {vcf_out.call_count} called SVs to {config.vcf} {vcf_output_info_str}")

    if monitor:
        log.debug(f'Stopping resource monitoring.')
        monitor.stop()


if __name__ == "__main__":
    processes = []

    try:
        logging.config.dictConfig({
            'version': 1,
            'formatters': {
                'default': {
                    'format': '%(asctime)s %(levelname)s %(name)s (%(process)d): %(message)s'
                }
            },
            'handlers': {
                'console': {
                    'class': 'logging.StreamHandler',
                    'formatter': 'default',
                    'stream': 'ext://sys.stdout',
                }
            },
            'loggers': {
                'MEIsensor.progress': {
                    'level': logging.WARNING,
                },
                'MEIsensor.vcf': {
                    'level': logging.INFO,
                }
            },
            'root': {
                'level': logging.INFO,
                'handlers': ['console'],
            },
            'disable_existing_loggers': False,
        })
    except (ValueError, TypeError, AttributeError, ImportError):
        logging.exception(f'Error configuring loggers.')

    try:
        Main(processes)
    except (utils.MEIsensorExit, SystemExit) as exit_code:
        if len(processes):
            # Allow time for child process error messages to propagate
            print("MEIsensorMain: Shutting down workers")
            time.sleep(10)
        for proc in processes:
            try:
                proc.process.terminate()
            except:
                pass

        for proc in processes:
            try:
                proc.process.join()
            except:
                pass
        exit(exit_code.code)
    except:
        logging.getLogger('MEIsensor.main').exception(f'Unhandled error while running MEIsensor.')
        exit(1)