#!/usr/bin/env python3
"""Species called in several samples of a run on a few reads, always beside the same more abundant relative: the
mark of a reference artefact (a contaminated or transferred gene copy that every present organism with that gene puts
a perfect read on), not of the species being there. At GTDB r226 such calls recurred per reference: 288 of 764
paired-end false positives were 116 taxa called falsely in two samples or more, and a fifth of all false positives
were reads of a present species of another genus at identity 0.99 (docs/claude/2026-10-03-false-positive-anatomy).
protal --build flags such copies in advance (suspect_copies.tsv); this looks at a finished run's outputs for what the
database did not flag.

Reads every <sample>.profile.log of a run (protal -o: every taxon the sample's reads hit, with Predicted, Probability,
Lineage, Name, TaxID and the Summary's Total Hits) and reports, per called species with at most --max-hits hits in a
sample, the samples it was called in, and the species present in all of them with at least --ratio times its hits
that shares the deepest rank with it (its likely source), with that rank. A species called in --min-samples samples
or more, all with the same companion, is marked suspect. Also written: how many of its calls were on few hits.

Usage: recurrent_calls.py OUTPUT_DIR [more dirs or .profile.log files] [-o recurrent_calls.tsv] [--max-hits 3]
       [--min-samples 2] [--ratio 10]
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
PREFIXES = {"d": "domain", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}


def lineage_of(text):
    """{rank: name} of a GTDB lineage string (d__...;p__...;...), or of ';'/'|'-separated names with prefixes."""
    out = {}
    for token in text.replace("|", ";").split(";"):
        token = token.strip()
        if len(token) > 3 and token[1:3] == "__" and token[0] in PREFIXES:
            out[PREFIXES[token[0]]] = token
    return out


def shared_rank(a, b):
    """The deepest rank two lineages share, or None."""
    deepest = None
    for rank in RANKS:
        if rank in a and a.get(rank) == b.get(rank):
            deepest = rank
        else:
            break
    return deepest


def total_hits(summary):
    """The Total Hits of a .profile.log Summary field ("{ Genes: 3, Total Hits: 7, VCOV: ... }"), or None."""
    key = "Total Hits:"
    at = summary.find(key)
    if at < 0:
        return None
    rest = summary[at + len(key):].strip()
    digits = ""
    for c in rest:
        if c.isdigit():
            digits += c
        else:
            break
    return int(digits) if digits else None


def read_profile_log(path):
    """[(taxid, name, lineage dict, hits, called, probability)] of a .profile.log."""
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                called = str(row.get("Predicted", "")).strip().lower() in ("1", "true", "yes")
                probability = float(row.get("Probability", "nan"))
                taxid = row.get("TaxID", "").strip()
                name = row.get("Name", "").strip()
                lineage = lineage_of(row.get("Lineage", ""))
                hits = total_hits(row.get("Summary", "") or "")
            except (TypeError, ValueError):
                continue
            if taxid:
                rows.append((taxid, name, lineage, hits, called, probability))
    return rows


def profile_logs(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            files += sorted(glob.glob(os.path.join(p, "**", "*.profile.log"), recursive=True))
        elif p.endswith(".profile.log"):
            files.append(p)
        else:
            sys.exit(f"{p}: not a folder nor a .profile.log")
    if not files:
        sys.exit("no .profile.log files found")
    return files


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="protal output folders (searched for *.profile.log) or .profile.log files")
    ap.add_argument("-o", "--output", default="recurrent_calls.tsv")
    ap.add_argument("--max-hits", type=int, default=3, help="a call on at most this many hits (reads) is a thin call (default 3)")
    ap.add_argument("--min-samples", type=int, default=2, help="thin calls in this many samples with one companion mark a species suspect (default 2)")
    ap.add_argument("--ratio", type=float, default=10, help="a companion has at least this many times the species' hits in the sample (default 10)")
    args = ap.parse_args()

    files = profile_logs(args.paths)
    # Per species: the samples it was called thinly in, and per such sample the companions (present species with
    # --ratio times its hits), each with the rank shared.
    thin = defaultdict(dict)       # taxid -> {sample: (hits, probability)}
    calls = defaultdict(int)       # taxid -> samples called in at all
    names, lineages = {}, {}
    companions = defaultdict(lambda: defaultdict(dict))  # taxid -> sample -> {companion taxid: (rank, companion hits)}
    for path in files:
        sample = os.path.basename(path)[: -len(".profile.log")]
        rows = read_profile_log(path)
        called = [r for r in rows if r[4]]
        for taxid, name, lineage, hits, _, probability in called:
            names[taxid], lineages[taxid] = name, lineage
            calls[taxid] += 1
            if hits is None or hits > args.max_hits:
                continue
            thin[taxid][sample] = (hits, probability)
            for other, _, other_lineage, other_hits, _, _ in called:
                if other == taxid or other_hits is None or other_hits < args.ratio * max(hits, 1):
                    continue
                rank = shared_rank(lineage, other_lineage)
                if rank is None or rank == "species":
                    continue
                companions[taxid][sample][other] = (rank, other_hits)
    out_rows = []
    for taxid, samples in sorted(thin.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        # Companions present in every thin sample of the species.
        shared = None
        for sample in samples:
            here = set(companions[taxid].get(sample, {}))
            shared = here if shared is None else shared & here
        shared = shared or set()
        best, best_rank, best_hits = "", "", 0
        for c in sorted(shared):
            rank = companions[taxid][next(iter(samples))][c][0]
            hits_sum = sum(companions[taxid][s][c][1] for s in samples)
            # The deepest shared rank first, then the most hits.
            key = (RANKS.index(rank), hits_sum)
            if not best or key > (RANKS.index(best_rank), best_hits):
                best, best_rank, best_hits = c, rank, hits_sum
        suspect = len(samples) >= args.min_samples and bool(best)
        out_rows.append({
            "taxid": taxid, "species": names.get(taxid, ""), "samples_called": calls[taxid], "thin_calls": len(samples),
            "thin_samples": ",".join(sorted(samples)), "hits": ",".join(str(samples[s][0]) for s in sorted(samples)),
            "probabilities": ",".join(f"{samples[s][1]:.3g}" for s in sorted(samples)),
            "companion_taxid": best, "companion": names.get(best, ""), "shared_rank": best_rank,
            "companion_hits": best_hits, "suspect": int(suspect)})
    with open(args.output, "w", newline="") as fh:
        fields = ["taxid", "species", "samples_called", "thin_calls", "thin_samples", "hits", "probabilities", "companion_taxid",
                  "companion", "shared_rank", "companion_hits", "suspect"]
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)
    suspects = [r for r in out_rows if r["suspect"]]
    print(f"{len(files)} samples; {len(out_rows)} species called on at most {args.max_hits} hits in some sample, "
          f"{len(suspects)} of them in {args.min_samples} or more samples always beside the same more abundant relative "
          f"(suspect): {args.output}")
    for r in suspects[:20]:
        print(f"  {r['species']} ({r['taxid']}): thin in {r['thin_calls']} of {r['samples_called']} samples, beside "
              f"{r['companion']} ({r['shared_rank']}, {r['companion_hits']} hits)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
