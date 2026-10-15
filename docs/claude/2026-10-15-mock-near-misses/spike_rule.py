#!/usr/bin/env python3
"""The spike-gene rule protal applies (src/Profiling/SpikeGenes.h, spike_genes::SpikesOf), on a SAM: per taxon, the
genes whose records are far more than the taxon's other genes explain. The defaults are protal's (the first prototype,
before the rule was fixed, took genes one at a time with --min-records 20 --fold 10 and a Poisson bound):

- the background is the taxon's records on its genes other than the --max-genes with the most records per base, over
  those genes' lengths (all its genes, from reference.map, hit or not; protal takes only the genes with unique k-mers,
  whose records the profile counts, and so differs where a pile lies on a gene without them);
- those top genes are taken in order while each holds at least --min-records records and --fold times what the
  background gives a gene of its length (and, with --p below 1, beyond the Poisson bound at --p over the taxon's genes),
  and while the taxon keeps more than --min-rest-genes genes with records;
- the genes taken are spikes only if together they hold --min-share of the taxon's records.

Prints per taxon with a spike: its records, genes hit, the spike genes with their records and expectation, and the
records left; with --truth (one taxid per line) whether the taxon is present. Summary: taxa with spikes, present and
absent, and the records they lose.

The last line ("# sample:") is protal's sample check (spike_genes::Judge): taxa with spikes, the spikes' share of the
records, and how many of the ten taxa with the most records have spikes.

usage: zstdcat X.sam.zst | spike_rule.py REFERENCE_MAP [--truth FILE] [--min-records 50] [--fold 50] [--p 1]
       [--min-share 0.5] [--min-rest-genes 5] [--max-genes 3] [--quiet]
"""
import argparse
import collections
import math
import sys


def poisson_tail_at_least(k, mean):
    """P(X >= k) for X ~ Poisson(mean), by the regularized gamma function's series (k >= 1)."""
    if k <= 0:
        return 1.0
    if mean <= 0:
        return 0.0
    # P(X >= k) = P(k, mean), the lower regularized gamma; series sum mean^n / (k (k+1) ... (k+n))
    term = 1.0 / k
    total = term
    n = 1
    while term > total * 1e-15 and n < 100000:
        term *= mean / (k + n)
        total += term
        n += 1
    log = k * math.log(mean) - mean - math.lgamma(k) + math.log(total)
    return math.exp(min(log, 0.0))


def spikes_of(counts, lengths, args):
    """counts: {gene: records}; lengths: {gene: length} of the taxon's hittable genes. As protal's
    spike_genes::SpikesOf (src/Profiling/SpikeGenes.h): the background leaves out the --max-genes genes with the most
    records per base; those are taken in order while each holds --min-records and --fold times the background's
    expectation; the taken genes count only if they hold --min-share of the records. Returns (gene, records, fold)."""
    counts = {g: c for g, c in counts.items() if lengths.get(g, 0) > 0}
    total_length = sum(lengths.values())
    total = sum(counts.values())
    if total_length == 0 or total == 0 or len(counts) <= args.min_rest_genes + 1:
        return []
    order = sorted(counts, key=lambda g: (-counts[g] / lengths[g], int(g) if g.isdigit() else g))
    top = order[:min(args.max_genes, len(counts) - args.min_rest_genes - 1)]
    background = (total - sum(counts[g] for g in top)) / (total_length - sum(lengths[g] for g in top))
    taken = []
    for g in top:
        expected = background * lengths[g]
        if counts[g] < args.min_records or counts[g] < args.fold * expected:
            break
        if args.p < 1 and poisson_tail_at_least(counts[g], expected) >= args.p / max(len(lengths), 1):
            break
        taken.append((g, counts[g], counts[g] / expected if expected > 0 else math.inf))
    if not taken or sum(c for _, c, _ in taken) < args.min_share * total:
        return []
    return taken


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference_map")
    ap.add_argument("--truth")
    ap.add_argument("--min-records", type=int, default=50)
    ap.add_argument("--fold", type=float, default=50.0)
    ap.add_argument("--p", type=float, default=1.0, help="a Poisson bound on top of the fold (1: none, as protal)")
    ap.add_argument("--min-share", type=float, default=0.5)
    ap.add_argument("--min-rest-genes", type=int, default=5)
    ap.add_argument("--max-genes", type=int, default=3)
    ap.add_argument("--quiet", action="store_true", help="the summary only")
    args = ap.parse_args()

    lengths = collections.defaultdict(dict)
    with open(args.reference_map) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 4 and f[0].isdigit():
                lengths[f[0]][f[1]] = int(f[3]) - int(f[2]) - 1
    truth = set()
    if args.truth:
        with open(args.truth) as fh:
            truth = {line.strip() for line in fh if line.strip()}

    counts = collections.defaultdict(collections.Counter)
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.split("\t", 4)
        if int(f[1]) & 0x104 or f[2] == "*":  # as protal: supplementary records (a long read's other genes) count
            continue
        taxon, gene = f[2].split("_", 1)
        counts[taxon][gene] += 1

    summary = collections.Counter()
    sample_records = sum(sum(g.values()) for g in counts.values())
    spike_records = 0
    spike_taxa = 0
    biggest = []  # the taxa with the most records: how many of the top ten have spikes
    for taxon, genes in sorted(counts.items(), key=lambda kv: -sum(kv[1].values()))[:10]:
        biggest.append(bool(spikes_of(genes, lengths.get(taxon, {}), args)))
    for taxon, genes in sorted(counts.items(), key=lambda kv: -sum(kv[1].values())):
        spikes = spikes_of(genes, lengths.get(taxon, {}), args)
        present = taxon in truth
        if not spikes:
            continue
        n = sum(genes.values())
        lost = sum(c for _, c, _ in spikes)
        spike_records += lost
        spike_taxa += 1
        key = "present" if present else "absent"
        summary[key + "_taxa"] += 1
        summary[key + "_records_lost"] += lost
        if not args.quiet or present:
            print(f"{taxon}\t{'present' if present else 'absent'}\trecords {n}\tgenes {len(genes)}\t"
                  f"spikes {','.join(f'{g}:{c}:x{fold:.0f}' for g, c, fold in spikes)}\tleft {n - lost}")
    total_taxa = len(counts)
    present_taxa = sum(1 for t in counts if t in truth)
    print(f"# taxa {total_taxa} (present {present_taxa}); with spikes: present {summary['present_taxa']} "
          f"(records lost {summary['present_records_lost']}), absent {summary['absent_taxa']} "
          f"(records lost {summary['absent_records_lost']})")
    share = spike_records / sample_records if sample_records else 0.0
    print(f"# sample: records {sample_records}, taxa with spikes {spike_taxa}, spike records {spike_records} "
          f"({share:.4%}), of the ten taxa with the most records {sum(biggest)} with spikes")


if __name__ == "__main__":
    sys.exit(main())
