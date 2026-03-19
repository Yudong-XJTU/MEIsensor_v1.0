#!/usr/bin/env python3
#
# MEIsensor

import logging

import pysam
import os

import SV
import utils
from config import Config

import torch
from trained_model.model import VariantSeqClassify


log = logging.getLogger(__name__)


def format_info(k, v):
    if isinstance(v, float):
        return f"{k}={v:.3f}"
    elif isinstance(v, list):
        return f"{k}={','.join(v)}"
    else:
        return f"{k}={v}"




def format_genotype(gt):
    """
    hp_i is the index of the haplotype in config.phase_identifiers:
    HP:1 => index 0 => phased genotype in the form of 1|0
    HP:2 => index 1 => phased genotype in the form of 0|1
    """
    if len(gt) == 5:
        a, b, qual, dr, dv = gt
        gt_sep = "/"
        return f"{a}{gt_sep}{b}:{qual}:{dr}:{dv}"
    else:
        a, b, qual, dr, dv, svid = gt
        gt_sep = "/"
        return f"{a}{gt_sep}{b}:{qual}:{dr}:{dv}:{svid}"

def one_hot_encoding(seq):
        base_to_one_hot = {
            'A': [1, 0, 0, 0],
            'C': [0, 1, 0, 0],
            'G': [0, 0, 1, 0],
            'T': [0, 0, 0, 1],
            'N': [0, 0, 0, 0]
        }
        one_hot_tensor = torch.zeros((4, len(seq)), dtype=torch.float32)
        for i, base in enumerate(seq.upper()):
            if base in base_to_one_hot:
                one_hot_tensor[:, i] = torch.tensor(base_to_one_hot[base], dtype=torch.float32)
        return one_hot_tensor


def classify_INS(variant_seq, Config):
    cuda_condition = torch.cuda.is_available()
    device = torch.device('cuda' if cuda_condition else 'cpu')
    model = VariantSeqClassify()
    model.load_state_dict(torch.load(Config.model))
    model.eval()
    model = model.to(device)
    target_length = 7000
    if len(variant_seq) < target_length:
        variant_seq += 'N' * (target_length - len(variant_seq))
    seq_tensor = one_hot_encoding(variant_seq).unsqueeze(0).to(device)
    model_prediction = torch.argmax(model(seq_tensor), dim=1).item()
    LINE_score = model(seq_tensor).squeeze().tolist()[2]
    Alu_INS = (model_prediction == 1)
    LINE1_INS = (model_prediction == 2)
    SVA_INS = (model_prediction == 3)
    other_INS = (model_prediction == 0)
    return Alu_INS, LINE1_INS, SVA_INS, other_INS, LINE_score



class VCF:
    def __init__(self, config: Config, handle):
        self.config = config
        self.handle = handle
        self.call_count = 0
        self.info_order = ["SVTYPE", "SVLEN", "END", "SUPPORT", "RNAMES", "COVERAGE", "STRAND"]
        if config.qc_nm_measure:
            self.info_order.append("NM")

        if config.dev_emit_sv_lengths:
            self.info_order.append("SVLENGTHS")

        self.default_genotype = config.genotype_none

        # Add phasing if needed
        self.genotype_format = config.genotype_format
        self.reference_handle = None
        self.header_str = ""

    def open_reference(self):
        if self.config.reference is None:
            return

        if not os.path.exists(self.config.reference + ".fai") and not os.path.exists(self.config.reference + ".gzi"):
            print(f"Info: Fasta index for {self.config.reference} not found. Generating with pysam.faidx "
                  f"(this may take a while)")
            pysam.faidx(self.config.reference)
        self.reference_handle = pysam.FastaFile(self.config.reference)

    def write_header(self, contigs_lengths):
        self.write_header_line("fileformat=VCFv4.2")
        for contig, contig_len in contigs_lengths:
            self.write_header_line(f"contig=<ID={contig},length={contig_len}>")

        self.write_header_line('ALT=<ID=INS,Description="Insertion">')
        self.write_header_line('FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        self.write_header_line('FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">')
        self.write_header_line('FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">')
        self.write_header_line('FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">')
        self.write_header_line('FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase-block, zero if none or not phased">')
        self.write_header_line('FORMAT=<ID=ID,Number=1,Type=String,Description="Individual sample SV ID for multi-sample output">')

        self.write_header_line('FILTER=<ID=PASS,Description="All filters passed">')
        self.write_header_line('FILTER=<ID=GT,Description="Genotype filter">')
        self.write_header_line('FILTER=<ID=SUPPORT_MIN,Description="Minimum read support filter">')
        self.write_header_line('FILTER=<ID=STDEV_POS,Description="SV Breakpoint standard deviation filter">')
        self.write_header_line('FILTER=<ID=STDEV_LEN,Description="SV length standard deviation filter">')
        self.write_header_line('FILTER=<ID=COV_MIN,Description="Minimum coverage filter">')
        self.write_header_line('FILTER=<ID=COV_MIN_GT,Description="Minimum coverage filter (missing genotype)">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_INS,Description="Coverage change filter for INS">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_US,Description="Coverage fractional change filter: upstream-start">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_SC,Description="Coverage fractional change filter: start-center">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_CE,Description="Coverage fractional change filter: center-end">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_ED,Description="Coverage fractional change filter: end-downstream">')
        self.write_header_line('FILTER=<ID=ALN_NM,Description="Length adjusted mismatch filter">')
        self.write_header_line('FILTER=<ID=STRAND,Description="Strand support filter for germline SVs">')
        self.write_header_line('FILTER=<ID=SVLEN_MIN,Description="SV length filter">')
        self.write_header_line('INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">')
        self.write_header_line('INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">')
        
        self.write_header_line('INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">')
        self.write_header_line('INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">')
        self.write_header_line('INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">')
        self.write_header_line('INFO=<ID=SUPPORT_INLINE,Number=1,Type=Integer,Description="Number of reads supporting an INS/DEL SV (non-split events only)">')
        self.write_header_line('INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="Number of soft-clipped reads putatively supporting the long insertion SV">')
        self.write_header_line('INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">')
        self.write_header_line('INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Standard deviation of structural variation start position">')
        self.write_header_line('INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Standard deviation of structural variation length">')
        self.write_header_line('INFO=<ID=COVERAGE,Number=.,Type=Float,Description="Coverages near upstream, start, center, end, downstream of structural variation">')
        self.write_header_line('INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of supporting reads for structural variant">')
        self.write_header_line('INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count, summed up over all samples">')
        self.write_header_line('INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="List of read support for all samples">')
        self.write_header_line('INFO=<ID=CONSENSUS_SUPPORT,Number=1,Type=Integer,Description="Number of reads that support the generated insertion (INS) consensus sequence">')
        self.write_header_line('INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads (if enabled with --output-rnames)">')
        self.write_header_line('INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Fraction">')
        self.write_header_line('INFO=<ID=NM,Number=.,Type=Float,Description="Mean number of query alignment length adjusted mismatches of supporting reads">')
        self.write_header_line('INFO=<ID=PHASE,Number=.,Type=String,Description="Phasing information derived from supporting reads, represented as list of: HAPLOTYPE,PHASESET,HAPLOTYPE_SUPPORT,PHASESET_SUPPORT,HAPLOTYPE_FILTER,PHASESET_FILTER">')
        samples_header = "\t".join(sample_id for _, sample_id in self.config.sample_ids_vcf)
        self.write_raw(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples_header}")
        
    def write_raw(self, text, endl="\n"):
        if self.config.vcf_output_bgz:
            self.handle.write(text.encode())
            self.handle.write(endl.encode())
        else:
            self.handle.write(text)
            self.handle.write(endl)

    def write_header_line(self, text):
        self.write_raw("##" + text)

    def write_call(self, call):
        if call.svtype != "INS":
            return
        if call.alt == "<INS>":
            return
        # pysam coordinates are 0-based, VCF 1-based
        # but VCF also requires the index of the base before the SV to be reported,
        # so we are fine without offsetting
        end = call.end
        pos = call.pos if call.pos > 0 else 1

        # Determine genotypes columns
        ac = 0
        supvec = []
        sample_genotypes = []
        for internal_id, _ in self.config.sample_ids_vcf:
            if internal_id in call.genotypes and call.genotypes[internal_id] is not None:
                gt_curr = call.genotypes[internal_id]
                sample_genotypes.append(format_genotype(gt_curr))
                if gt_curr[0] != "." and gt_curr[4] > 0:  # Not non-genotype and has supporting reads
                    ac += sum(call.genotypes[internal_id][:2])
                    supp = "1"
                else:
                    supp = "0"
            else:
                sample_genotypes.append(format_genotype(self.default_genotype))
                supp = "0"
            supvec.append(supp)

        if len(self.config.sample_ids_vcf) > 1:
            call.set_info("AC", ac)
            call.set_info("SUPP_VEC", svec := "".join(supvec))

            if int(svec) == 0:
                log.debug(f'Dropped {call} due to all zero support vector.')
                return

            if ac == 0:
                call.filter = "GT"

        # Output core SV attributes
        infos = {
            "SVTYPE": call.svtype,
            "SVLEN": call.svlen,
            "SVLENGTHS": ",".join(map(str, call.svlens)) if call.svlens else None,
            "END": end,
            "SUPPORT": call.support,
            "RNAMES": call.rnames if self.config.output_rnames else None,
            "COVERAGE": f"{call.coverage_upstream},{call.coverage_start},{call.coverage_center},{call.coverage_end},"
                        f"{call.coverage_downstream}",
            "STRAND": ("+" if call.fwd > 0 else "") + ("-" if call.rev > 0 else ""),
            "NM": call.nm
        }
        infos_ordered = ["PRECISE" if call.precise else "IMPRECISE"]
        af = call.get_info("VAF")
        af = af if af is not None else 0
        sv_is_mosaic = af <= self.config.mosaic_af_max
        if sv_is_mosaic and self.config.mosaic:
            infos_ordered.append("MOSAIC")
        infos_ordered.extend(format_info(k, infos[k]) for k in self.info_order if infos[k] is not None)
        info_str = ";".join(infos_ordered)

        # Output call specific additional information
        for k in sorted(call.info):
            if call.info[k] is None:
                continue
            info_str += ";" + format_info(k, call.info[k])

        # if call.id==None:
        #    call.id=f"{call.svtype}.{self.call_count+1:06}"

        if (not self.config.symbolic and call.svtype == "DEL" and self.reference_handle is not None
                and abs(call.svlen) <= self.config.max_del_seq_len):
            return

        if self.config.symbolic:
            call.ref = "N"
            call.alt = f"<{call.svtype}>"
        else:
            if self.reference_handle is not None and call.ref == 'N':
                try:
                    call.ref = self.reference_handle.fetch(call.contig, start := max(0, call.pos - 1), start + 1)
                except (KeyError, ValueError):
                    ...
                else:
                    if call.svtype == "INS" and call.alt != '<INS>':
                        INS_seq = call.ref + call.alt
                        Alu_INS, LINE1_INS, SVA_INS, other_INS, score = classify_INS(INS_seq, self.config)
                        if SVA_INS or Alu_INS or LINE1_INS:
                            call.alt = INS_seq
                            if SVA_INS:
                                call.id = 'SVA.INS'
                            elif Alu_INS:
                                call.id = "ALU.INS"
                            elif LINE1_INS:
                                call.id = f"LINE1.INS"
                        else:
                            log.debug(f"Dropped {call} due to non-specific INS type.")
                            return

        call.qual = max(0, min(60, call.qual)) if call.qual is not None else None

        self.write_raw("\t".join(str(v) for v in [call.contig, pos, self.config.id_prefix + call.id, call.ref,
                                                  call.alt, call.qual if call.qual is not None else '.', call.filter, info_str, self.genotype_format] +
                                 sample_genotypes))
        self.call_count += 1

    def read_svs_iter(self):
        self.header_str = ""
        line_index = 0
        for line in self.handle:
            try:
                if isinstance(line, bytes):
                    line = line.decode("utf-8")
                line_index += 1
                line_strip = line.strip()
                if line_strip == "" or line_strip[0] == "#":
                    if line_strip[0] == "#":
                        self.header_str += line_strip + "\n"
                    continue
                CHROM, POS, _, REF, ALT, QUAL, FILTER, INFO = line.split("\t")[:8]
                info_dict = {}
                for info_item in INFO.split(";"):
                    if "=" in info_item:
                        key, value = info_item.split("=")
                    else:
                        key, value = info_item, True
                    info_dict[key] = value
                call = SV.SVCall(contig=CHROM,
                                 pos=int(POS) - 1,
                                 id=line_index,
                                 ref=REF,
                                 alt=ALT,
                                 qual=int(QUAL) if QUAL != '.' else None,
                                 filter=FILTER,
                                 info=info_dict,
                                 svtype=None,
                                 svlen=None,
                                 end=None,
                                 rnames=None,
                                 qc=True,
                                 postprocess=None,
                                 genotypes=None,
                                 precise=None,
                                 support=0,
                                 fwd=0,
                                 rev=0,
                                 nm=-1)
                if len(call.alt) > len(call.ref):
                    call.svtype = "INS"
                    call.svlen = len(call.alt)
                    call.end = call.pos
                else:
                    call.svtype = "DEL"
                    call.svlen = -len(call.ref)
                    call.end = call.pos + call.svlen

                if "SVTYPE" in info_dict:
                    call.svtype = info_dict["SVTYPE"]
                    if call.svtype == "TRA":
                        call.svtype = "BND"

                if "SVLEN" in info_dict:
                    call.svlen = int(info_dict["SVLEN"])
                if "SVLENGTHS`" in info_dict:
                    call.svlens = info_dict["SVLENGTHS"]

                if "END" in info_dict:
                    call.end = int(info_dict["END"])

                call.raw_vcf_line = line_strip
                call.raw_vcf_line_index = line_index
                yield call
            except Exception as e:
                utils.fatal_error(f"Error parsing input VCF: Line {line_index}: {e}")

    def rewrite_genotype(self, svcall):
        parts_no_gt = svcall.raw_vcf_line.split("\t")[:8]
        gt_format = self.config.genotype_format
        if svcall.genotype_match_sv != None:
            if len(svcall.genotype_match_sv.genotypes) > 0:
                gt = svcall.genotype_match_sv.genotypes[0]
            else:
                gt = svcall.genotypes[0]
        else:
            gt = svcall.genotypes[0]
        parts = parts_no_gt + [gt_format, format_genotype(gt)]
        self.write_raw("\t".join(parts))

    def rewrite_header_genotype(self,orig_header):
        header_lines=orig_header.split("\n")
        header_lines.insert(1,'##genotypeFileDate="'+self.config.start_date+'"')
        header_lines.insert(1,'##genotypeCommand="'+self.config.command+'"')
        header_lines.insert(1,f"##genotypeSource={self.config.version}_{self.config.build}")

        has_gt_headers = {
            "GT": False,
            "GQ": False,
            "DR": False,
            "DV": False,
        }
        for header_line in header_lines:
            for gt in has_gt_headers.keys():
                if "##FORMAT=<ID="+gt+"," in header_line:
                    has_gt_headers[gt] = True

        if not has_gt_headers["GT"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        if not has_gt_headers["GQ"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">')
        if not has_gt_headers["DR"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">')
        if not has_gt_headers["DV"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">')

        self.write_raw("\n".join(header_lines), endl="")

    def close(self):
        self.handle.close()