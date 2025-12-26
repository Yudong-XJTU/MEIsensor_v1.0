# cluster.py

from dataclasses import dataclass
import statistics
import math
from typing import Optional

import SV
import preprocess


@dataclass
class Cluster:
    id: str
    svtype: str
    contig: str
    start: int
    end: int
    seed: int
    leads: list
    leads_long: Optional[list]

    @property
    def span(self) -> Optional[int]:
        if self.end is None or self.start is None:
            return None

        return self.end - self.start

    def compute_metrics(self, max_n=100):
        n = min(len(self.leads), max_n)
        if n == 0:
            self.mean_svlen = 0
            self.stdev_start = 0
            return

        step = int(len(self.leads) / n)
        if n > 1:
            self.mean_svlen = sum(self.leads[i].svlen for i in range(0, len(self.leads), step)) / float(n)
            self.stdev_start = statistics.stdev(self.leads[i].ref_start for i in range(0, len(self.leads), step))
        else:
            self.mean_svlen = self.leads[0].svlen
            self.stdev_start = 0


def merge_inner(cluster, threshold):
    read_seq = {}
    for ld in cluster.leads:
        if ld.read_qname not in read_seq:
            read_seq[ld.read_qname] = []
        read_seq[ld.read_qname].append(ld)

    cluster.leads = []
    for qname in read_seq:
        read_seq[qname].sort(key=lambda k: k.ref_start)
        to_merge = read_seq[qname][0]

        curr_lead = to_merge

        last_ref_end = to_merge.ref_end
        last_qry_end = to_merge.qry_end
        last_ref_start = to_merge.ref_start
        last_qry_start = to_merge.qry_start

        for to_merge in read_seq[qname][1:]:
            merge = (threshold == -1) or ((abs(to_merge.ref_start - last_ref_end) < threshold or abs(to_merge.ref_start - last_ref_start) < threshold) and (
                        abs(to_merge.qry_start - last_qry_end) < threshold or abs(to_merge.qry_start - last_qry_start) < threshold))
            if merge:
                curr_lead.svlen += to_merge.svlen
                if to_merge.seq is None or curr_lead.seq is None:
                    curr_lead.seq = None
                else:
                    curr_lead.seq += to_merge.seq
            else:
                cluster.leads.append(curr_lead)
                curr_lead = to_merge
            last_ref_end = to_merge.ref_end
            last_qry_end = to_merge.qry_end
            last_ref_start = to_merge.ref_start
            last_qry_start = to_merge.qry_start

        cluster.leads.append(curr_lead)
    return cluster


def resplit(cluster, prop, binsize, merge_threshold_min, merge_threshold_frac):
    bins_leads = {}
    for lead in cluster.leads:
        bin = int(abs(prop(lead)) / binsize) * binsize
        if not bin in bins_leads:
            bins_leads[bin] = [lead]
        else:
            bins_leads[bin].append(lead)

    new_clusters = list(sorted(bins_leads.keys()))
    i = 1
    while len(new_clusters) > 1 and i < len(new_clusters):
        last_cluster = new_clusters[i - 1]
        curr_cluster = new_clusters[i]
        merge_threshold = max(merge_threshold_min, min(curr_cluster, last_cluster) * merge_threshold_frac)
        merge = abs(curr_cluster - last_cluster) <= merge_threshold
        if merge:
            bins_leads[new_clusters[i]].extend(bins_leads[new_clusters[i - 1]])
            new_clusters.pop(i - 1)
            i = max(0, i - 2)
        else:
            i += 1

    for cluster_index in new_clusters:
        new_cluster = Cluster(id=cluster.id + f".{cluster_index}",
                              svtype=cluster.svtype,
                              contig=cluster.contig,
                              start=cluster.start,
                              end=cluster.end,
                              seed=cluster.seed,
                              leads=bins_leads[cluster_index],
                              leads_long=cluster.leads_long)
        yield new_cluster



def resolve(svtype, leadtab_provider, config):
    leadtab = leadtab_provider.leadtab[svtype]
    seeds = sorted(leadtab_provider.leadtab[svtype])

    if len(seeds) == 0:
        return []

    clusters = []
    for seed_index, seed in enumerate(seeds):
        if config.dev_call_region != None:
            if seed < config.dev_call_region["start"] or seed > config.dev_call_region["end"]:
                continue
        if svtype == "INS":
            leads = [lead for lead in leadtab[seed] if lead.svlen != None]
            leads_long = [lead for lead in leadtab[seed] if lead.svlen == None]

        else:
            leads = leadtab[seed]
            leads_long = None

        cluster = Cluster(id=f"CL.{svtype}.{leadtab_provider.contig}.{leadtab_provider.start}.{seed_index}",
                          svtype=svtype,
                          contig=leadtab_provider.contig,
                          start=seed,
                          end=seed + config.cluster_binsize,
                          seed=seed,
                          leads=leads,
                          leads_long=leads_long)

        cluster.compute_metrics()
        clusters.append(cluster)

    # Merge clusters
    cluster_count_initial = len(clusters)
    i = 0
    while i < len(clusters) - 1:
        curr_cluster = clusters[i]
        next_cluster = clusters[i + 1]

        inner_dist = (next_cluster.start - curr_cluster.end)
        outer_dist = (next_cluster.end - curr_cluster.start)
        merge = inner_dist <= min(curr_cluster.stdev_start, next_cluster.stdev_start) * config.cluster_r
        merge = merge

        if merge:
            clusters.pop(i + 1)
            curr_cluster.leads += next_cluster.leads
            if svtype == "INS":
                curr_cluster.leads_long += next_cluster.leads_long
            curr_cluster.end = next_cluster.end
            curr_cluster.compute_metrics()
            i = max(0, i - 2)
        i += 1

    if config.dev_trace_read:
        for c in clusters:
            for ld in c.leads:
                if ld.read_qname == config.dev_trace_read:
                    print(f"[DEV_TRACE_READ [2/4] [cluster.resolve] Read lead {ld} is in cluster {c.id}, containing a total of {len(c.leads)} leads")

    if config.dev_dump_clusters:
        filename = f"{config.input}.clusters.{svtype}.{leadtab_provider.contig}.{leadtab_provider.start}.{leadtab_provider.end}.bed"
        print(f"Dumping clusters to {filename}")
        with open(filename, "w") as h:
            for c in clusters:
                info = f"ID={c.id}, #LEADS={len(c.leads)}; "
                for ld in c.leads:
                    info += f"(ref_start={ld.ref_start},svlen={ld.svlen},source={ld.source}); "
                h.write(f"{c.contig}\t{c.start}\t{c.end}\t\"{info}\"\n")

    for cluster in clusters:
        if len(cluster.leads) == 0:
            continue
        if svtype == "INS":
            merge_inner_threshold = config.cluster_merge_pos
            merge_inner(cluster, merge_inner_threshold)

        if not config.dev_no_resplit:
            for new_cluster in resplit(cluster,
                                       prop=lambda lead: lead.svlen,
                                       binsize=config.cluster_resplit_binsize,
                                       merge_threshold_min=config.minsvlen,
                                       merge_threshold_frac=config.cluster_merge_len):
                yield new_cluster
        else:
            yield cluster

'''
def resolve_block_groups(svtype, svcands, groups_initial, config):

    # TODO: Remove sorting
    groups = groups_initial
    for svcand in sorted(svcands, key=lambda cand: cand.support, reverse=True):
        best_group = None
        best_dist = math.inf
        for group in groups:
            # TODO: Favor bigger groups in placement
            dist = abs(group.pos_mean - svcand.pos) + abs(abs(group.len_mean) - abs(
                svcand.svlen))  # check if group.pos_mean is updated or stays the same for the first SV starting the group
            minlen = float(min(abs(group.len_mean), abs(svcand.svlen)))
            if minlen > 0 and dist < best_dist and dist <= config.combine_match * math.sqrt(
                    minlen) and dist <= config.combine_match_max:
                if (
                        not config.combine_separate_intra or svcand.sample_internal_id not in group.included_samples) and group.align_call(
                        svcand, config.combine_pctseq):
                    best_group = group
                    best_dist = dist
        if best_group is None:
            groups.append(
                sv.SVGroup.from_candidate(svcand)
            )
        else:
            best_group.add_candidate(svcand)

    return groups
'''