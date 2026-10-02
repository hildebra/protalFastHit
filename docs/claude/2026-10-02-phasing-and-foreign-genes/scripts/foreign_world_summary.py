#!/usr/bin/env python3
"""foreign_world_summary.py OUT - what foreign_world.sh's runs made of the donor's gene A.

Per read type and triple (hgt<i>, its control ctl<i> without the recipient): the links of T's gene A in hgt<i> and
whether it is foreign; its depth there against the control's (VCov, .profile.genes.log); T's depth (TaxVCOV) in the
control and in hgt<i> with the gene kept and left out; the strain MSA's multi-allelic positions of gene A in hgt<i>'s
row (.meta.tsv; none when it is left out); and the other foreign genes of all samples (none expected).
"""
import collections
import csv
import glob
import os
import statistics
import sys

out = sys.argv[1]
triples = list(csv.DictReader(open(os.path.join(out, "world", "triples.tsv")), delimiter="\t"))


def genes_log(run, sample):
    path = os.path.join(out, "runs", run, f"{sample}.profile.genes.log")
    rows = {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path), delimiter="\t"):
            rows[(r["TaxID"], r["GeneID"])] = r
    return rows


def meta_multi(run, species, sample, gene):
    """The multi-allelic positions of the gene in the sample's MSA row(s) (strain rows: <sample>_hap<k>, comma-
    separated, with the number of rows the sample has), '-' if it has none there."""
    for path in glob.glob(os.path.join(out, "runs", run, "**", species.replace(" ", "_") + ".meta.tsv"), recursive=True):
        found, rows = [], set()
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["sample"] == sample or r["sample"].startswith(sample + "_hap"):
                rows.add(r["sample"])
                if r["gene_id"] == gene:
                    found.append(r["multi_allelic"])
        if found:
            return ",".join(found) + (f" ({len(rows)} rows)" if len(rows) > 1 else "")
        if rows:
            return "-" + (f" ({len(rows)} rows)" if len(rows) > 1 else "")
    return "-"


print("| reads | triple | T's gene A: links judged / unlikely | foreign | gene A depth: hgt / ctl | T's depth: ctl | hgt kept | "
      "hgt dropped | gene A multi-allelic in T's MSA row: kept / dropped |")
print("|" + "---|" * 9)
other = collections.Counter()
errors = collections.defaultdict(list)
for kind in ("pe", "pb", "ont"):
    for tr in triples:
        i = tr["sample"].removeprefix("hgt")
        t, gene = tr["t_taxid"], tr["gene"]
        keep, drop = genes_log(f"{kind}.keep", f"hgt{i}"), genes_log(f"{kind}.drop", f"hgt{i}")
        ctl = genes_log(f"{kind}.keep", f"ctl{i}")
        row = keep.get((t, gene))
        if row is None:
            print(f"| {kind} | {i} | no reads | | | | | | |")
            continue
        tax = lambda rows: next((float(r["TaxVCOV"]) for (tx, _), r in rows.items() if tx == t), float("nan"))
        gene_ctl = float(ctl[(t, gene)]["VCov"]) if (t, gene) in ctl else 0.0
        d_ctl, d_keep, d_drop = tax(ctl), tax(keep), tax(drop)
        errors[(kind, "kept")].append(d_keep / d_ctl - 1)
        errors[(kind, "dropped")].append(d_drop / d_ctl - 1)
        species = tr["t_species"].removeprefix("s__")
        multi = (meta_multi(f"{kind}.keep", "s__" + species, f"hgt{i}", gene),
                 meta_multi(f"{kind}.drop", "s__" + species, f"hgt{i}", gene))
        print(f"| {kind} | {i} | {row['LinksJudged']} / {row['LinksUnlikely']} | {'yes' if row['Foreign'] == '1' else 'no'} | "
              f"{float(row['VCov']):.1f} / {gene_ctl:.1f} | {d_ctl:.2f} | {d_keep:.2f} | {d_drop:.2f} | {multi[0]} / {multi[1]} |")
        for run in (f"{kind}.keep",):
            for sample in (f"hgt{i}", f"ctl{i}"):
                for (tx, g), r in genes_log(run, sample).items():
                    if r["Foreign"] == "1" and not (sample == f"hgt{i}" and tx == t and g == gene):
                        other[kind] += 1
print()
for (kind, arm), v in sorted(errors.items()):
    print(f"{kind}, gene A {arm}: T's depth against the control, median {statistics.median(v):+.3f}, "
          f"max {max(v, key=abs):+.3f}")
print("other foreign genes (all samples, kept runs):", dict(other))
