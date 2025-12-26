#!/usr/bin/env python3
#
# MEIsensor
#

import os
import sys
import datetime
import argparse
from collections import defaultdict

from typing import Union, Optional

import utils
from region import Region

VERSION = "MEIsensor"
BUILD = "1.0.0"


class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


def tobool(v):
    if v is True or v is False:
        return v
    elif v.strip().lower() == "true" or v.strip() == "1":
        return True
    elif v.strip().lower() == "false" or v.strip() == "0":
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value (True | False) required for argument")




class Config(argparse.Namespace):
    '''
    input: str
    vcf: str
    reference: str
    phase: bool
    threads: int
    contig: Optional[str]
    run_id: str
    '''

    def add_main_args(self, parser):
        main_args = parser.add_argument_group("Common parameters")
        main_args.add_argument("-i", "--input", metavar="input.bam", type=str,
                               help="(Required)For single-sample calling: A coordinate-sorted and indexed .bam/.cram (BAM/CRAM format) file containing aligned reads.)",
                               required=True)
        main_args.add_argument("-v", "--vcf", metavar="output.vcf", type=str,
                               help="VCF output filename to write the called and refined SVs to. If the given filename ends with .gz, the VCF file will be automatically bgzipped and a .tbi index built for it.",
                               required=False)
        main_args.add_argument("--reference", metavar="reference.fasta", type=str,
                               help="(Optional) Reference sequence the reads were aligned against. To enable output of deletion SV sequences, this parameter must be set.",
                               default=None)
        main_args.add_argument("-t", "--threads", metavar="N", type=int,
                               help="Number of parallel threads to use (speed-up for multi-core CPUs)", default=4)
        main_args.add_argument("-m", "--model", metavar="model.ph", type=str,
                               help="(Required)Trained model")

    def __init__(self, *args, **kwargs):
        # super().__init__(**kwargs)

        parser = argparse.ArgumentParser(description="", epilog=' ',
                                         formatter_class=lambda prog: ArgFormatter(prog, max_help_position=100,
                                                                                   width=150), usage=' ')
        parser.add_argument("--version", action="version", version=f"MEIsensor, Version {BUILD}")

        self.add_main_args(parser)
        parser.parse_args(args=args or None, namespace=self)
        '''
        self.input = kwargs.get("input", None)
        self.vcf = kwargs.get("vcf", None)
        self.reference = kwargs.get("reference", None)
        self.threads = kwargs.get("threads", 4)
        '''
        self.contig = kwargs.get("contig", None)
        self.regions = kwargs.get("regions", None)

        # filter
        self.minsupport = kwargs.get("minsupport", "auto")
        self.minsupport_auto_mult = kwargs.get("minsupport_auto_mult", None)
        self.minsvlen = kwargs.get("minsvlen", 50)
        self.minsvlen_screen_ratio = kwargs.get("minsvlen_screen_ratio", 0.9)
        self.mapq = kwargs.get("mapq", 0)  # argparse.SUPPRESS if not provided
        self.no_qc = kwargs.get("no_qc", False)
        self.qc_stdev = kwargs.get("qc_stdev", True)
        self.qc_stdev_abs_max = kwargs.get("qc_stdev_abs_max", 500)
        self.qc_strand = kwargs.get("qc_strand", False)
        self.qc_coverage = kwargs.get("qc_coverage", 1)
        self.long_ins_length = kwargs.get("long_ins_length", 2500)
        self.max_splits_kb = kwargs.get("max_splits_kb", 0.1)
        self.max_splits_base = kwargs.get("max_splits_base", 3)
        self.min_alignment_length = kwargs.get("min_alignment_length", 0)
        self.detect_large_ins = kwargs.get("detect_large_ins", True)

        # SV Clustering parameters
        self.cluster_binsize = kwargs.get("cluster_binsize", 100)
        self.cluster_r = kwargs.get("cluster_r", 2.5)
        self.cluster_repeat_h = kwargs.get("cluster_repeat_h", 1.5)
        self.cluster_repeat_h_max = kwargs.get("cluster_repeat_h_max", 1000)
        self.cluster_merge_pos = kwargs.get("cluster_merge_pos", 150)
        self.cluster_merge_len = kwargs.get("cluster_merge_len", 0.33)
        self.cluster_merge_bnd = kwargs.get("cluster_merge_bnd", 1000)

        # genotype
        self.genotype_ploidy = kwargs.get("genotype_ploidy", 2)
        self.genotype_error = kwargs.get("genotype_error", 0.05)
        self.sample_id = kwargs.get("sample_id", None)
        self.genotype_vcf = kwargs.get("genotype_vcf", None)

        self.re_qc = kwargs.get("re_qc", "auto")

        # Additional parameter
        self.allow_overwrite = kwargs.get("allow_overwrite", False)

        # Postprocessing, QC and output parameters
        self.output_rnames = kwargs.get("output_rnames", False)
        self.no_consensus = kwargs.get("no_consensus", False)
        self.no_sort = kwargs.get("no_sort", False)
        self.no_progress = kwargs.get("no_progress", False)
        self.quiet = kwargs.get("quiet", False)
        self.max_del_seq_len = kwargs.get("max_del_seq_len", 50000)
        self.symbolic = kwargs.get("symbolic", False)
        self.allow_overwrite = kwargs.get("allow_overwrite", False)

        # Mosaic calling mode parameters
        self.mosaic = kwargs.get("mosaic", False)
        self.mosaic_af_max = kwargs.get("mosaic_af_max", 0.218)
        self.mosaic_af_min = kwargs.get("mosaic_af_min", 0.05)
        self.mosaic_qc_invdup_min_length = kwargs.get("mosaic_qc_invdup_min_length", 500)
        self.mosaic_qc_nm = kwargs.get("mosaic_qc_nm", True)
        self.mosaic_qc_nm_mult = kwargs.get("mosaic_qc_nm_mult", 1.66)
        self.mosaic_qc_coverage_max_change_frac = kwargs.get("mosaic_qc_coverage_max_change_frac", 0.1)
        self.mosaic_qc_strand = kwargs.get("mosaic_qc_strand", True)
        self.mosaic_include_germline = kwargs.get("mosaic_include_germline", False)

        # Developer parameters
        self.dev_emit_sv_lengths = kwargs.get("dev_emit_sv_lengths", False)
        self.dev_cache = kwargs.get("dev_cache", False)
        self.dev_cache_dir = kwargs.get("dev_cache_dir", None)
        self.dev_debug_svtyping = kwargs.get("dev_debug_svtyping", False)
        self.dev_keep_lowqual_splits = kwargs.get("dev_keep_lowqual_splits", False)
        self.dev_call_region = kwargs.get("dev_call_region", None)
        self.dev_dump_clusters = kwargs.get("dev_dump_clusters", False)
        self.dev_merge_inline = kwargs.get("dev_merge_inline", False)
        self.dev_seq_cache_maxlen = kwargs.get("dev_seq_cache_maxlen", 50000)
        self.consensus_max_reads = kwargs.get("consensus_max_reads", 20)
        self.consensus_max_reads_bin = kwargs.get("consensus_max_reads_bin", 10)
        self.combine_consensus = kwargs.get("combine_consensus", False)
        self.dev_dump_coverage = kwargs.get("dev_dump_coverage", False)
        self.dev_no_resplit = kwargs.get("dev_no_resplit", False)
        self.dev_no_resplit_repeat = kwargs.get("dev_no_resplit_repeat", False)
        self.low_memory = kwargs.get("low_memory", False)
        self.qc_nm = kwargs.get("qc_nm", False)
        self.qc_nm_mult = kwargs.get("qc_nm_mult", 1.66)
        self.qc_coverage_max_change_frac = kwargs.get("qc_coverage_max_change_frac", -1)
        self.coverage_updown_bins = kwargs.get("coverage_updown_bins", 5)
        self.coverage_shift_bins = kwargs.get("coverage_shift_bins", 3)
        self.coverage_shift_bins_min_aln_length = kwargs.get("coverage_shift_bins_min_aln_length", 1000)
        self.cluster_binsize_combine_mult = kwargs.get("cluster_binsize_combine_mult", 5)
        self.cluster_resplit_binsize = kwargs.get("cluster_resplit_binsize", 20)
        self.dev_trace_read = kwargs.get("dev_trace_read", False)
        self.dev_split_max_query_distance_mult = kwargs.get("dev_split_max_query_distance_mult", 5)
        self.dev_no_qc = kwargs.get("dev_no_qc", False)
        self.dev_disable_interblock_threads = kwargs.get("dev_disable_interblock_threads", False)
        self.dev_combine_medians = kwargs.get("dev_combine_medians", False)
        self.dev_monitor_memory = kwargs.get("dev_monitor_memory", 0)
        self.dev_monitor_filename = kwargs.get("dev_monitor_filename", None)
        self.dev_debug_log = kwargs.get("dev_debug_log", False)
        self.dev_progress_log = kwargs.get("dev_progress_log", False)

        # supplement parameter
        self.task_count_multiplier = 0
        self.coverage_binsize = self.cluster_binsize
        self.coverage_binsize_combine = self.cluster_binsize * self.cluster_binsize_combine_mult
        self.run_id = f'{os.environ.get("SLURM_JOB_ID") or os.getpid()}'

        self.qc_nm_measure = self.qc_nm

        if self.quiet:
            sys.stdout = open(os.devnull, "w")

        self.start_date = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        self.run_id = f'{os.environ.get("SLURM_JOB_ID") or os.getpid()}'

        self.task_count_multiplier = 0

        self.version = VERSION
        self.build = BUILD
        self.command = " ".join(sys.argv)

        if self.dev_call_region is not None:
            region_contig, region_startend = self.dev_call_region.replace(",", "").split(":")
            start, end = region_startend.split("-")
            self.dev_call_region = dict(contig=region_contig, start=int(start), end=int(end))

        if self.contig and self.regions:
            utils.fatal_error('Please provide either --contig or --regions, not both.')

        if self.regions is not None:
            regions = defaultdict(list)
            with open(self.regions, 'r') as f:
                for line in f.readlines():
                    if line.startswith('#') or line.strip() == '':
                        continue
                    else:
                        r = Region.from_bed_line(line)
                        if r is not None:
                            regions[r.contig].append(r)
            self.regions_by_contig = regions
        else:
            self.regions_by_contig = {}

        # "--minsvlen" parameter is for final output filtering
        # for intermediate steps, a lower threshold is used to account for sequencing, mapping imprecision
        self.minsvlen_screen = int(self.minsvlen_screen_ratio * self.minsvlen)
        # config.minsupport_screen=max(1,int(0.333*config.minsupport*(config.cluster_binsize/100.0)))

        if self.minsupport != "auto":
            self.minsupport = int(self.minsupport)

        if self.dev_no_qc:
            self.no_qc = True

        if self.re_qc == 'auto':
            self.reqc = 'auto'
        elif self.re_qc in ('0', '1'):
            self.reqc = bool(int(self.re_qc))
        else:
            utils.fatal_error('Invalid value for --re-qc, allowed values are: auto, 0, 1')

        if not hasattr(self, 'mapq'):
            self.mapq = 0 if self.dev_no_qc else 20
        if not hasattr(self, 'min_alignment_length'):
            self.min_alignment_length = 0 if self.dev_no_qc else 1000

        # --minsupport auto defaults
        self.minsupport_auto_base = 1.5
        self.minsupport_auto_regional_coverage_weight = 0.75

        if self.minsupport_auto_mult is None:
            self.minsupport_auto_mult = 0.1

        self.coverage_binsize = self.cluster_binsize
        self.coverage_binsize_combine = self.cluster_binsize * self.cluster_binsize_combine_mult

        # INS Consensus parameters
        # config.consensus_max_reads=20
        # config.consensus_max_reads_bin=10
        self.consensus_min_reads = 4
        self.consensus_kmer_len = 6
        self.consensus_kmer_skip_base = 3
        self.consensus_kmer_skip_seqlen_mult = 1.0 / 500.0
        self.consensus_low_threshold = 0.0  # 0.15

        # Large INS
        self.long_ins_rescale_base = 1.66
        self.long_ins_rescale_mult = 0.33

        # BND
        self.bnd_cluster_length = 1000

        # Genotyping
        self.genotype_format = "GT:GQ:DR:DV"
        self.genotype_none = (".", ".", 0, 0, 0, (None, None))
        self.genotype_null = (0, 0, 0, 0, 0, (None, None))
        self.genotype_min_z_score = 5
        if self.genotype_ploidy != 2:
            utils.fatal_error("Currently only --genotype-ploidy 2 is supported")

        # Combine
        self.combine_exhaustive = False
        self.combine_relabel_rare = False
        self.combine_overlap_abs = 2500
        self.combine_min_size = 100

        # Misc
        self.precise = 25  # Max. sum of pos and length stdev for SVs to be labelled PRECISE
        self.tandem_repeat_region_pad = 500
        self.id_prefix = "MEIsensor."
        self.phase_identifiers = ["1", "2"]

        self.dev_profile = False

        self.workdir = os.getcwd()

        # Mosaic
        if self.mosaic_include_germline:
            self.mosaic = True

        self.qc_nm_measure = self.qc_nm
        if self.mosaic:
            # config.qc_coverage_max_change_frac=config.mosaic_qc_coverage_max_change_frac
            self.qc_nm_measure = self.qc_nm_measure or self.mosaic_qc_nm
            # config.qc_nm_mult=config.mosaic_qc_nm_mult
            # config.qc_strand=config.mosaic_qc_strand

        Config.GLOBAL = self

    @property
    def sort(self) -> bool:
        """
        Output is sorted unless explicitly asked for otherwise, and only if the output is not bgzipped
        """
        return self.vcf_output_bgz or not self.no_sort

    @property
    def vcf_output_bgz(self) -> Optional[bool]:
        """
        Should the output vcf file be compressed?
        """
        if self.vcf:
            path, ext = os.path.splitext(self.vcf)
            return ext == ".gz" or ext == ".bgz"
        return None