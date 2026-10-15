#!/usr/bin/env python3
"""What protal saw of the species it missed (the false negatives of a scored run, score_mocks.py's species.tsv): per
missed species and sample, from the sample's SAM
  records, reads, genes   best records on the species' own genes (taxid_gene), their reads and genes
  kept_records, kept_reads  of those, the ones the profiler keeps (MAPQ >= 4)
  alt_reads, alt_tie_reads  reads whose best record is another taxon's but whose ZA lists the species within 5 edits,
                          and within 0 edits (it fits as well)
and from the forest profile (.profile.log of reprofile.sh): its probability and depth (VCov), and, for a truth of
abundances, the depth its share predicts (the share times the sample's median of depth / share over its true positives).
Classes: no_alignment (no record, maybe an alternative), filtered (records, none kept), rejected (kept records, the model
says no). Standard library only.

usage: fn_alignments.py SPECIES_TSV SAM_RUN_FOLDER PROFILE_FOLDER GENOME2TIID OUT_TSV [MIN_MAPQ]
"""
import collections
import csv
import glob
import os
import statistics
import subprocess
import sys

KNOB = 0.5


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def profile_rows(path):
    """{species: (probability, vcov, predicted)} of a .profile.log."""
    rows = {}
    if os.path.exists(path):
        for r in read_tsv(path):
            rows[r["Lineage"].split(";")[-1]] = (float(r["Probability"]), float(r["VCov"]), r["Predicted"])
    return rows


def sam_counts(path, wanted, min_mapq):
    """Per wanted taxid: best records, their reads and genes, kept records and reads; reads listing it in ZA."""
    c = {t: {"records": 0, "reads": set(), "genes": set(), "kept_records": 0, "kept_reads": set(), "alt_reads": set(),
             "alt_tie_reads": set()} for t in wanted}
    text = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True).stdout
    for line in text:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if int(f[1]) & (4 | 256):
            continue
        taxid, gene = f[2].split("_", 1)
        if taxid in c:
            e = c[taxid]
            e["records"] += 1
            e["reads"].add(f[0])
            e["genes"].add(gene)
            if int(f[4]) >= min_mapq:
                e["kept_records"] += 1
                e["kept_reads"].add(f[0])
        for tag in f[11:]:
            if tag.startswith("ZA:Z:") and tag != "ZA:Z:*":
                for item in tag[5:].split(","):
                    t, edits = item.split(":")
                    if t in c and t != taxid:
                        c[t]["alt_reads"].add(f[0])
                        if int(edits) <= 0:
                            c[t]["alt_tie_reads"].add(f[0])
    return c


def main():
    species_tsv, run, profiles, genome2tiid, out = sys.argv[1:6]
    min_mapq = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    taxid_of = {}
    with open(genome2tiid) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            taxid_of.setdefault(f[3].split(";")[-1], f[1])
    rows = read_tsv(species_tsv)
    by_sample = collections.defaultdict(list)
    for r in rows:
        by_sample[r["sample"]].append(r)
    out_rows = []
    for sample, species in sorted(by_sample.items()):
        missed = [r for r in species if r["status"] == "FN"]
        if not missed:
            continue
        sams = glob.glob(os.path.join(run, "**", sample + ".sam.zst"), recursive=True)
        if not sams:
            print(f"{sample}: no SAM", file=sys.stderr)
            continue
        wanted = {taxid_of[r["species"]] for r in missed if r["species"] in taxid_of}
        counts = sam_counts(sams[0], wanted, min_mapq)
        profile = profile_rows(os.path.join(profiles, sample + ".profile.log"))
        # depth per unit of share, from the true positives of a truth of abundances
        ratios = [profile[r["species"]][1] / float(r["truth_share"]) for r in species if r["status"] == "TP"
                  and r["species"] in profile and r["truth_share"] not in ("member", "-") and float(r["truth_share"]) > 0]
        per_share = statistics.median(ratios) if ratios else None
        for r in missed:
            t = taxid_of.get(r["species"])
            e = counts.get(t) if t else None
            p, vcov, _ = profile.get(r["species"], (None, None, None))
            share = r["truth_share"]
            expected = (f"{per_share * float(share):.3g}" if per_share is not None and share not in ("member", "-")
                        else "NA")
            if e is None:
                kind = "not_in_database"
            elif e["records"] == 0:
                kind = "no_alignment"
            elif e["kept_records"] == 0:
                kind = "filtered"
            else:
                kind = "rejected"
            out_rows.append({"sample": sample, "species": r["species"], "truth_share": share, "class": kind,
                             "records": e["records"] if e else 0, "reads": len(e["reads"]) if e else 0,
                             "genes": len(e["genes"]) if e else 0, "kept_records": e["kept_records"] if e else 0,
                             "kept_reads": len(e["kept_reads"]) if e else 0,
                             "alt_reads": len(e["alt_reads"]) if e else 0,
                             "alt_tie_reads": len(e["alt_tie_reads"]) if e else 0,
                             "probability": f"{p:.3f}" if p is not None else "NA",
                             "vcov": f"{vcov:.3g}" if vcov is not None else "NA", "expected_vcov": expected})
    columns = ["sample", "species", "truth_share", "class", "records", "reads", "genes", "kept_records", "kept_reads",
               "alt_reads", "alt_tie_reads", "probability", "vcov", "expected_vcov"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    tally = collections.Counter((r["sample"], r["class"]) for r in out_rows)
    for (sample, kind), n in sorted(tally.items()):
        print(f"{sample}\t{kind}\t{n}")


if __name__ == "__main__":
    main()
