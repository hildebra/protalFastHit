#!/usr/bin/env python3
"""How many GTDB species have a 16S and a 23S rRNA gene, from a release's metadata and its ssu_all FASTA.

    python3 rrna_coverage.py ~/GTDB/r226 > coverage_output.txt      # ~2 minutes, ~3 GB of memory

Reads RELEASE/{bac120,ar53}_metadata_rNNN.tsv.gz (GTDB's per-genome ssu_* and lsu_23s_* columns: the longest hit's
length, the number of hits, the contig it lies on, its best SILVA hit) and RELEASE/genomic_files_all/ssu_all_rNNN.fna.gz
(every 16S hit of >= 200 bp of every genome). GTDB distributes no 23S sequences: the 23S numbers are the metadata's.

A gene counts as full length at MIMAG's high-quality lengths: 16S >= 1200 bp (bacteria) or >= 900 bp (archaea), 23S >=
1900 bp. Genome categories (isolate, MAG, SAG) as scripts/download_gtdb.py tells them apart.
"""

import argparse
import collections
import glob
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
from download_gtdb import genome_category  # noqa: E402

COLS = ["accession", "gtdb_representative", "gtdb_genome_representative", "gtdb_taxonomy", "ncbi_genome_category",
        "ncbi_organism_name", "ncbi_assembly_level", "genome_size", "contig_count", "checkm2_completeness",
        "checkm2_contamination", "ssu_count", "ssu_length", "ssu_contig_len", "ssu_query_id", "ssu_silva_taxonomy",
        "ssu_silva_blast_perc_identity", "lsu_23s_count", "lsu_23s_length", "lsu_23s_contig_len", "lsu_23s_query_id",
        "lsu_silva_23s_taxonomy", "lsu_silva_23s_blast_perc_identity", "lsu_5s_count", "lsu_5s_length"]
FULL16 = {"Bacteria": 1200, "Archaea": 900}
FULL23 = 1900
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def read_metadata(release):
    frames = []
    for path in sorted(glob.glob(os.path.join(release, "*_metadata_r*.tsv.gz"))):
        df = pd.read_csv(path, sep="\t", usecols=COLS, na_values=["none"], keep_default_na=False, low_memory=False)
        frames.append(df)
        print(f"# {os.path.basename(path)}: {len(df):,} genomes", file=sys.stderr)
    df = pd.concat(frames, ignore_index=True)
    ranks = df["gtdb_taxonomy"].str.split(";", expand=True)
    for i, rank in enumerate(RANKS):
        df[rank] = ranks[i].str[3:]
    # "none" was read as missing: give genome_category the text it expects
    df["category"] = [genome_category(c, o) for c, o in zip(df["ncbi_genome_category"].fillna("none").astype(str),
                                                             df["ncbi_organism_name"].fillna("").astype(str))]
    df["is_isolate"] = df["category"] == "isolate"
    df["is_rep"] = df["gtdb_representative"].astype(str).str.lower().isin(["t", "true"])
    for col in ("ssu_count", "lsu_23s_count", "lsu_5s_count"):
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    for col in ("ssu_length", "lsu_23s_length", "lsu_5s_length", "ssu_contig_len", "lsu_23s_contig_len", "genome_size",
                "ssu_silva_blast_perc_identity", "lsu_silva_23s_blast_perc_identity"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    full16 = df["domain"].map(FULL16)
    df["has16"] = df["ssu_length"].notna()
    df["full16"] = df["ssu_length"] >= full16
    df["has23"] = df["lsu_23s_length"].notna()
    df["full23"] = df["lsu_23s_length"] >= FULL23
    df["both_full"] = df["full16"] & df["full23"]
    df["same_contig"] = df["both_full"] & (df["ssu_query_id"] == df["lsu_23s_query_id"])
    return df


def pct(n, d):
    return f"{100.0 * n / d:5.1f}%" if d else "   NA"


def table(rows, header):
    widths = [max(len(str(r[i])) for r in rows + [header]) for i in range(len(header))]
    out = ["  ".join(str(h).ljust(w) for h, w in zip(header, widths))]
    out += ["  ".join(str(c).rjust(w) if i else str(c).ljust(w) for i, (c, w) in enumerate(zip(r, widths))) for r in rows]
    return "\n".join(out)


def ssu_fasta_lengths(path):
    """{accession: (hits, longest ssu_len)} from ssu_all's headers, and the number of records."""
    per = {}
    n = 0
    rx = re.compile(r"\[ssu_len=(\d+)\]")
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] != ">":
                continue
            n += 1
            acc = line[1:line.index("~")]
            m = rx.search(line)
            length = int(m.group(1)) if m else 0
            hits, best = per.get(acc, (0, 0))
            per[acc] = (hits + 1, max(best, length))
    return per, n


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("release", help="GTDB release folder (e.g. ~/GTDB/r226)")
    args = p.parse_args()
    release = os.path.expanduser(args.release)
    df = read_metadata(release)
    reps = df[df["is_rep"]].set_index("gtdb_genome_representative")
    if not reps.index.is_unique:
        sys.exit("representatives are not unique")
    sp = df.groupby("gtdb_genome_representative")
    species = pd.DataFrame({
        "n_genomes": sp.size(),
        "any16": sp["has16"].any(), "anyfull16": sp["full16"].any(),
        "any23": sp["has23"].any(), "anyfull23": sp["full23"].any(),
        "anyboth": sp["both_full"].any(), "anysame": sp["same_contig"].any(),
        "best16": sp["ssu_length"].max(), "best23": sp["lsu_23s_length"].max(),
        "isolates": sp["is_isolate"].sum(),
    })
    species = species.join(reps[["domain", "phylum", "class", "order", "family", "genus", "species", "category",
                                 "has16", "full16", "has23", "full23", "both_full", "same_contig", "ssu_length",
                                 "lsu_23s_length", "ssu_contig_len", "lsu_23s_contig_len", "ssu_silva_taxonomy",
                                 "lsu_silva_23s_taxonomy", "lsu_5s_count", "ssu_count", "lsu_23s_count"]])

    print(f"GTDB release folder: {release}")
    print(f"genomes {len(df):,}, species {len(species):,}")
    print()
    print("## 1. Genomes and species")
    rows = []
    for dom in ("Bacteria", "Archaea"):
        g = df[df["domain"] == dom]
        s = species[species["domain"] == dom]
        cats = s["category"].value_counts()
        rows.append([dom, f"{len(g):,}", f"{len(s):,}", f"{(s['n_genomes'] == 1).sum():,}",
                     f"{cats.get('isolate', 0):,}", f"{cats.get('mag', 0):,}", f"{cats.get('sag', 0):,}"])
    print(table(rows, ["domain", "genomes", "species", "singleton species", "rep isolate", "rep MAG", "rep SAG"]))
    print()

    print("## 2. Species with a 16S / 23S: on the representative, and on any genome of the species")
    print("(full = MIMAG length: 16S >= 1200 bp bacteria / 900 bp archaea, 23S >= 1900 bp; same contig = the longest 16S")
    print(" and the longest 23S on one contig, an assembled operon)")
    rows = []
    for dom in ("Bacteria", "Archaea", "all"):
        s = species if dom == "all" else species[species["domain"] == dom]
        n = len(s)
        rows.append([dom, "representative", pct(s["has16"].sum(), n), pct(s["full16"].sum(), n), pct(s["has23"].sum(), n),
                     pct(s["full23"].sum(), n), pct(s["both_full"].sum(), n), pct(s["same_contig"].sum(), n),
                     pct((s["lsu_5s_count"] > 0).sum(), n)])
        rows.append([dom, "any genome", pct(s["any16"].sum(), n), pct(s["anyfull16"].sum(), n), pct(s["any23"].sum(), n),
                     pct(s["anyfull23"].sum(), n), pct(s["anyboth"].sum(), n), pct(s["anysame"].sum(), n), ""])
    print(table(rows, ["domain", "from", "16S", "16S full", "23S", "23S full", "both full", "both, same contig", "5S"]))
    n_all = len(species)
    print(f"species with no 16S hit at all on any genome: {(~species['any16']).sum():,} ({pct((~species['any16']).sum(), n_all)})")
    print(f"species with no full-length 16S on any genome: {(~species['anyfull16']).sum():,}")
    print(f"species with no 23S hit on any genome: {(~species['any23']).sum():,}; no full-length 23S: {(~species['anyfull23']).sum():,}")
    rescued16 = (~species["full16"] & species["anyfull16"]).sum()
    rescued23 = (~species["full23"] & species["anyfull23"]).sum()
    print(f"representative lacks a full 16S but another genome has one: {rescued16:,}; the same for 23S: {rescued23:,}")
    print()

    print("## 2b. Clades with a full-length rRNA gene (any genome of any of their species)")
    print("('held-out species' = species whose clade at that rank has a full gene in ANOTHER species: were the species")
    print(" missing from the database, its reads could still meet its clade's gene)")
    rows = []
    for rank in ("genus", "family", "order", "class", "phylum"):
        g = species.groupby(rank)
        n16 = g["anyfull16"].transform("sum")
        n23 = g["anyfull23"].transform("sum")
        nboth = g["anyboth"].transform("sum")
        clades = g.size()
        rows.append([rank, f"{len(clades):,}", pct(g["anyfull16"].any().sum(), len(clades)),
                     pct(g["anyfull23"].any().sum(), len(clades)), pct(g["anyboth"].any().sum(), len(clades)),
                     pct(((n16 - species["anyfull16"]) > 0).sum(), n_all), pct(((n23 - species["anyfull23"]) > 0).sum(), n_all),
                     pct(((nboth - species["anyboth"]) > 0).sum(), n_all)])
    print(table(rows, ["rank", "clades", "16S full", "23S full", "both", "held-out species: 16S", "23S", "both"]))
    print()

    print("## 3. By the representative's category (representative / any genome)")
    rows = []
    for cat in ("isolate", "mag", "sag"):
        s = species[species["category"] == cat]
        n = len(s)
        rows.append([cat, f"{n:,}", f"{pct(s['full16'].sum(), n)} / {pct(s['anyfull16'].sum(), n)}",
                     f"{pct(s['full23'].sum(), n)} / {pct(s['anyfull23'].sum(), n)}",
                     f"{pct(s['both_full'].sum(), n)} / {pct(s['anyboth'].sum(), n)}",
                     pct((s["isolates"] > 0).sum(), n)])
    print(table(rows, ["rep category", "species", "16S full", "23S full", "both full", "species with an isolate genome"]))
    print()
    print("By species size (any genome):")
    rows = []
    for label, mask in (("1 genome", species["n_genomes"] == 1), ("2-9 genomes", species["n_genomes"].between(2, 9)),
                        ("10+ genomes", species["n_genomes"] >= 10)):
        s = species[mask]
        n = len(s)
        rows.append([label, f"{n:,}", pct(s["anyfull16"].sum(), n), pct(s["anyfull23"].sum(), n), pct(s["anyboth"].sum(), n)])
    print(table(rows, ["species size", "species", "16S full", "23S full", "both full"]))
    print()

    print("## 4. Lengths of the representative's longest hit (where there is one)")
    rows = []
    for dom in ("Bacteria", "Archaea"):
        s = species[species["domain"] == dom]
        for col, label in (("ssu_length", "16S"), ("lsu_23s_length", "23S")):
            v = s[col].dropna()
            q = np.percentile(v, [5, 25, 50, 75, 95]).round().astype(int)
            rows.append([dom, label, f"{len(v):,}"] + [str(x) for x in q])
    print(table(rows, ["domain", "gene", "n", "p5", "p25", "median", "p75", "p95"]))
    bins16 = [0, 500, 900, 1200, 1400, 10_000]
    bins23 = [0, 1000, 1900, 2500, 10_000]
    for col, bins, label in (("ssu_length", bins16, "16S"), ("lsu_23s_length", bins23, "23S")):
        v = species[col]
        counts = pd.cut(v.dropna(), bins, right=False).value_counts().sort_index()
        print(f"{label}: none {v.isna().sum():,}; " + "; ".join(f"{iv.left}-{iv.right - 1}: {c:,}" for iv, c in counts.items()))
    print()

    print("## 5. Phyla (>= 300 species) by representative full 16S, lowest first")
    g = species.groupby("phylum")
    ph = pd.DataFrame({"species": g.size(), "full16": g["full16"].mean(), "anyfull16": g["anyfull16"].mean(),
                       "full23": g["full23"].mean(), "anyfull23": g["anyfull23"].mean(),
                       "mag": 1 - g["category"].agg(lambda c: (c == "isolate").mean())})
    ph = ph[ph["species"] >= 300].sort_values("full16")
    rows = [[name, f"{int(r.species):,}", f"{100 * r.full16:.1f}% / {100 * r.anyfull16:.1f}%",
             f"{100 * r.full23:.1f}% / {100 * r.anyfull23:.1f}%", f"{100 * r.mag:.0f}%"] for name, r in ph.iterrows()]
    print(table(rows, ["phylum", "species", "16S full rep / any", "23S full rep / any", "rep MAG or SAG"]))
    print()

    print("## 6. Copy numbers and read yield, complete isolate genomes")
    comp = df[(df["ncbi_assembly_level"] == "Complete Genome") & (df["category"] == "isolate")]
    print(f"complete isolate genomes: {len(comp):,} of {comp['gtdb_genome_representative'].nunique():,} species")
    rows = []
    for dom in ("Bacteria", "Archaea"):
        c = comp[comp["domain"] == dom]
        # one genome per species (the first), so large species do not dominate
        c1 = c.drop_duplicates("gtdb_genome_representative")
        q16 = np.percentile(c1["ssu_count"], [10, 50, 90]).astype(int)
        q23 = np.percentile(c1["lsu_23s_count"], [10, 50, 90]).astype(int)
        operon = c1["ssu_count"] * 1540 + c1["lsu_23s_count"] * 2900 + c1["lsu_5s_count"] * 115
        share = (operon / c1["genome_size"]) * 100
        rows.append([dom, f"{len(c1):,}", f"{q16[0]}/{q16[1]}/{q16[2]}", f"{q23[0]}/{q23[1]}/{q23[2]}",
                     f"{np.median(operon) / 1000:.1f} kb", f"{np.median(share):.2f}%",
                     f"{np.percentile(share, 10):.2f}-{np.percentile(share, 90):.2f}%"])
    print(table(rows, ["domain", "species", "16S copies p10/50/90", "23S copies p10/50/90", "rRNA bp median",
                       "of genome, median", "p10-p90"]))
    markers = pd.read_csv(os.path.join(REPO, "scripts", "mini_db", "markers_r226.tsv"), sep="\t", comment="#")
    for mset, gsize in (("bac120", "Bacteria"), ("ar53", "Archaea")):
        m = markers[markers["set"] == mset]
        bp = m["length_aa"].sum() * 3
        med = comp[comp["domain"] == gsize].drop_duplicates("gtdb_genome_representative")["genome_size"].median()
        print(f"{mset}: {len(m)} markers, ~{bp / 1000:.0f} kb by HMM length = {100 * bp / med:.2f}% of a median "
              f"{med / 1e6:.2f} Mb genome")
    for rank in ("genus", "family", "order"):
        with_complete = set(comp[rank])
        clades = species[rank].unique()
        print(f"{rank}: {pct(sum(c in with_complete for c in clades), len(clades))} of {len(clades):,} clades have a complete "
              f"isolate genome (a copy number); {pct(species[rank].isin(with_complete).sum(), n_all)} of species are in one")
    fam = comp.drop_duplicates("gtdb_genome_representative").groupby("family")["ssu_count"]
    spread = fam.agg(lambda c: np.percentile(c, 90) / max(np.percentile(c, 10), 1))[fam.size() >= 5]
    print(f"16S copies within a family (families with >= 5 complete species): p90/p10 ratio median {spread.median():.1f}, "
          f"p90 {np.percentile(spread, 90):.1f}")
    mag = df[df["category"] != "isolate"]
    print(f"MAGs/SAGs: 16S hits per genome where any: median {mag.loc[mag['ssu_count'] > 0, 'ssu_count'].median():.0f}, "
          f"share with >1: {pct((mag['ssu_count'] > 1).sum(), (mag['ssu_count'] > 0).sum())} (operons collapse in assemblies)")
    print()

    print("## 7. Can the representative's rRNA be trusted?")
    rows = []
    for cat in ("isolate", "mag", "sag"):
        s = species[species["category"] == cat]
        h16 = s[s["has16"]]
        h23 = s[s["has23"]]
        dom16 = h16["ssu_silva_taxonomy"].astype(str).str.split(";").str[0]
        dom23 = h23["lsu_silva_23s_taxonomy"].astype(str).str.split(";").str[0]
        organelle16 = h16["ssu_silva_taxonomy"].astype(str).str.contains("Mitochondria|Chloroplast", regex=True)
        print(f"{cat}: SILVA 16S domains where they differ: {dict(dom16[dom16 != h16['domain']].value_counts().head(4))}")
        rows.append([cat, f"{len(h16):,}", pct((dom16 != h16["domain"]).sum(), len(h16)), pct(organelle16.sum(), len(h16)),
                     pct((h16["ssu_contig_len"] < 2000).sum(), len(h16)), pct((h16["ssu_contig_len"] < 5000).sum(), len(h16)),
                     f"{len(h23):,}", pct((dom23 != h23["domain"]).sum(), len(h23)),
                     pct((h23["lsu_23s_contig_len"] < 4000).sum(), len(h23))])
    print(table(rows, ["rep category", "with 16S", "SILVA domain differs", "organelle hit", "16S contig < 2 kb",
                       "< 5 kb", "with 23S", "SILVA domain differs", "23S contig < 4 kb"]))
    print("(a 16S on a contig barely longer than the gene was binned without genomic context: the usual way a")
    print(" contaminating 16S enters a MAG)")
    print()

    print("## 8. The 16S sequences GTDB distributes (genomic_files_all/ssu_all)")
    path = glob.glob(os.path.join(release, "genomic_files_all", "ssu_all_r*.fna.gz"))
    if path:
        per, n = ssu_fasta_lengths(path[0])
        acc = df["accession"]
        in_fasta = acc.isin(per.keys())
        print(f"records {n:,} from {len(per):,} genomes; metadata genomes with a 16S: {df['has16'].sum():,}, of them in the "
              f"FASTA: {(in_fasta & df['has16']).sum():,}; FASTA genomes not in the metadata: "
              f"{len(set(per) - set(acc)):,}")
        longest = acc.map(lambda a: per.get(a, (0, 0))[1])
        both = df["has16"] & in_fasta
        agree = (longest[both] == df.loc[both, "ssu_length"]).mean()
        print(f"the FASTA's longest hit equals the metadata's ssu_length for {100 * agree:.1f}% of those genomes")
        hits = acc.map(lambda a: per.get(a, (0, 0))[0])
        print(f"hits per genome (with any): median {hits[hits > 0].median():.0f}, mean {hits[hits > 0].mean():.2f}; "
              f"metadata ssu_count equals the FASTA's hits for {100 * (hits[both] == df.loc[both, 'ssu_count']).mean():.1f}%")
        rep_seq = species.index.isin(per.keys())
        print(f"species whose representative has a 16S sequence in the FASTA: {rep_seq.sum():,} ({pct(rep_seq.sum(), n_all)})")
    else:
        print("no ssu_all FASTA in the release folder")


if __name__ == "__main__":
    main()
