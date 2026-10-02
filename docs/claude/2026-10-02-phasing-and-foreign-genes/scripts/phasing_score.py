#!/usr/bin/env python3
"""phasing_score.py [BENCH_DIR [RUNS]] - scores the phasing benchmark's runs (run_phasing.sh; RUNS: their folder under
BENCH_DIR/phasing, default runs) against the known strains.

Truth: make_mixtures.py's truth.tsv (each sample's strains of each species and their shares) and the strains' genomes
(make_strains.py: substitutions only, so a strain's genes have the representative's coordinates). From each run's raw
MSA of each of the 12 species (protal's own, before qcmsa; its columns map to the genes' positions):

- a row's sample and strain row: <sample> (one row) or <sample>_hap<k>;
- per mixture sample (2 or 3 strains), its discriminating positions: the gene positions where its strains do not all
  have one base. Each of its rows is matched to the strain whose bases it shows at most of them; a strain is resolved
  at the share of the discriminating positions where a row matched to it shows its base (the best such row), and that
  row's wrong calls are those of its A/C/G/T calls there that are another strain's base. A strain is recovered when
  it is resolved at half the positions or more with at most 10% wrong calls;
- SNPs over all gene positions: FP a row's A/C/G/T base that differs from the reference and from the strain the row
  is matched to (every row); TP and FN every strain's SNPs against the reference, called or not by the strain's best
  row (a strain without a row misses all of them), so that recall counts a mixture's minor strain with or without
  phasing;
- placement: whether a mixture strain's row is nearest (p-distance, columns where both have A/C/G/T) to the pure sample
  of that strain;
- pure samples given strain rows (falsely phased), and the strain rows qcmsa keeps in the filtered MSA.
Writes BENCH_DIR/phasing/results/rows.tsv (one line per run, species and mixture sample) and summary.md.
"""
import collections
import csv
import glob
import os
import re
import statistics
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
RUNS_NAME = sys.argv[2] if len(sys.argv) > 2 else "runs"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = [sys.argv[0], B]
sys.path.insert(0, os.path.join(HERE, "..", "..", "2026-10-02-v072-benchmark", "scripts"))
import strain_score as ss  # noqa: E402

P = os.path.join(B, "phasing")
RUNS = os.path.join(P, RUNS_NAME)
OUT = os.path.join(P, "results" if RUNS_NAME == "runs" else "results_" + RUNS_NAME)
ACGT = set("ACGT")


def main():
    truth = collections.defaultdict(lambda: collections.defaultdict(list))  # species -> sample -> [(strain, share)]
    coverage = collections.defaultdict(dict)  # species -> sample -> total coverage
    genome_of = {}
    for r in csv.DictReader(open(os.path.join(P, "truth.tsv")), delimiter="\t"):
        truth[r["species"]][r["sample"]].append((r["strain"], float(r["share"])))
        coverage[r["species"]][r["sample"]] = coverage[r["species"]].get(r["sample"], 0) + float(r["coverage"])
        genome_of[r["strain"]] = r["genome"]
    rep_of = {}
    for r in csv.DictReader(open(os.path.join(ss.S, "strains.tsv")), delimiter="\t"):
        if r["role"] == "strain":
            rep_of[r["species"]] = r["representative"]
    rep_path = {r["accession"]: r["fasta_path"] for r in csv.DictReader(open(os.path.join(ss.W, "genomes.tsv")), delimiter="\t")}
    markers = collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(ss.W, "marker_positions.tsv")), delimiter="\t"):
        markers[r["accession"]][r["marker"]] = r
    marker_of = {}
    with open(os.path.join(B, "V072", "protal_db", "gene2geneid.tsv")) as fh:
        for line in fh:
            name, gid = line.rstrip("\n").split("\t")[:2]
            marker_of[gid] = name
    seq_cache = {}

    def strain_gene(strain, rep, m):
        if strain not in seq_cache:
            order = ss.contig_order(rep_path[rep])
            seqs = ss.read_fasta(genome_of[strain])
            seq_cache[strain] = (seqs, {name: f"{strain}_contig{i + 1}" for i, name in enumerate(order)})
        seqs, index = seq_cache[strain]
        g = seqs[index[m["contig"]]][int(m["start"]) - 1:int(m["end"])]
        return g.translate(ss.COMPLEMENT)[::-1] if m["strand"] == "-" else g

    lines, falsely = [], collections.Counter()
    for run_dir in sorted(d for d in glob.glob(os.path.join(RUNS, "*")) if os.path.isdir(d)):
        run = os.path.basename(run_dir)
        if not os.path.exists(run_dir + ".done"):
            continue
        for raw in sorted(glob.glob(os.path.join(run_dir, "**", "*.raw.msa.fna"), recursive=True)):
            base = raw[:-len(".raw.msa.fna")]
            species = os.path.basename(base).removeprefix("s__").replace("_", " ")
            if species not in truth:
                continue
            msa = ss.read_fasta(raw)
            ref_name = next((n for n in msa if n.endswith("_reference")), None)
            if ref_name is None:
                continue
            ref = msa[ref_name]
            rep = rep_of[species]
            # The gene positions (MSA columns) and every strain's bases there.
            columns, strain_bases = [], collections.defaultdict(list)
            strains_here = {s for sample in truth[species].values() for s, _ in sample}
            for name, lo, hi in ss.parse_partition(base + ".raw.partition.txt"):
                marker = marker_of.get(re.sub(r"\D", "", name))
                m = markers[rep].get(marker) if marker else None
                if m is None:
                    continue
                cols = [i for i in range(lo - 1, hi) if ref[i] != "-"]
                genes = {s: strain_gene(s, rep, m) for s in strains_here}
                if any(len(g) != len(cols) for g in genes.values()):
                    continue
                columns.extend(cols)
                for s, g in genes.items():
                    strain_bases[s].extend(g)
            ref_bases = [ref[i] for i in columns]
            rows = collections.defaultdict(dict)  # sample -> row name -> its bases at the columns
            for name, seq in msa.items():
                if name == ref_name:
                    continue
                sample = re.sub(r"_hap\d+$", "", name)
                rows[sample][name] = [seq[i] for i in columns]
            for sample in rows:
                if sample.startswith("pure") and len(rows[sample]) > 1:
                    falsely[run] += 1
            final_path = base + ".msa.fna"
            final = ss.read_fasta(final_path) if os.path.exists(final_path) else {}
            pure_rows = {s: r[s] for s, r in rows.items() if s.startswith("pure") and s in r}

            def pdist(a, b):
                both = diff = 0
                for x, y in zip(a, b):
                    if x in ACGT and y in ACGT:
                        both += 1
                        diff += x != y
                return diff / both if both else 1.0

            for sample, mix in sorted(truth[species].items()):
                if len(mix) < 2 or sample not in rows:
                    continue
                strains = [s for s, _ in mix]
                disc = [k for k in range(len(columns)) if len({strain_bases[s][k] for s in strains}) > 1]
                best = {s: (0.0, 0.0, None) for s in strains}  # strain -> (resolved, wrong share, row)
                tp = fp = fn = 0
                for name, bases in rows[sample].items():
                    agree = {s: sum(1 for k in disc if bases[k] == strain_bases[s][k]) for s in strains}
                    s = max(strains, key=lambda x: agree[x])
                    called = [k for k in disc if bases[k] in ACGT]
                    wrong = sum(1 for k in called if bases[k] != strain_bases[s][k] and
                                any(bases[k] == strain_bases[o][k] for o in strains if o != s))
                    resolved = agree[s] / len(disc) if disc else 0.0
                    if resolved > best[s][0] or best[s][2] is None:
                        best[s] = (resolved, wrong / len(called) if called else 0.0, name)
                    # Miscalled SNPs: every row's, against the strain it is matched to.
                    for k in range(len(columns)):
                        c, r, t = bases[k], ref_bases[k], strain_bases[s][k]
                        fp += c in ACGT and c != r and c != t
                # Every strain's SNPs, from its best row (none: all missed), so that a mixture's minor strain counts
                # whether it has a row or not.
                for s in strains:
                    name = best[s][2]
                    for k in range(len(columns)):
                        r, t = ref_bases[k], strain_bases[s][k]
                        if t == r:
                            continue
                        hit = name is not None and rows[sample][name][k] == t
                        tp += hit
                        fn += not hit
                for s, share in mix:
                    resolved, wrong, name = best[s]
                    pure = f"pure{int(s.rsplit('_', 1)[1])}"
                    placed = ""
                    if name is not None and pure_rows:
                        nearest = min(pure_rows, key=lambda p: pdist(rows[sample][name], pure_rows[p]))
                        placed = int(nearest == pure)
                    lines.append(dict(run=run, species=species, sample=sample, strains=len(mix), strain=s, share=share,
                                      coverage=round(coverage[species][sample], 2), rows=len(rows[sample]),
                                      discriminating=len(disc), resolved=round(resolved, 4), wrong=round(wrong, 4),
                                      recovered=int(resolved >= 0.5 and wrong <= 0.1), placed=placed,
                                      kept_by_qcmsa=int(name in final) if name else 0, TP=tp, FP=fp, FN=fn))
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "rows.tsv"), "w") as fh:
        cols = list(lines[0])
        fh.write("\t".join(cols) + "\n")
        fh.writelines("\t".join(str(l[c]) for c in cols) + "\n" for l in lines)

    def frac(v):
        v = [x for x in v if x != ""]
        return f"{sum(v) / len(v):.2f}" if v else "-"

    out = ["# Phasing benchmark (phasing_score.py)", "",
           "Per run: the mixture samples' strains (12 species x 8 samples: 4 at 70:30, 2 at 50:50, one 50:30:20, one "
           "85:15). Shares of the strains; SNPs summed over the strains' matched rows.", "",
           "| run | strains | recovered | resolved (mean) | wrong calls (mean) | nearest pure sample is the strain's | "
           "kept by qcmsa | SNP precision | SNP recall | pure samples given strain rows |", "|" + "---|" * 10]
    by_run = collections.defaultdict(list)
    for l in lines:
        by_run[l["run"]].append(l)
    for run, ls in sorted(by_run.items()):
        # SNPs once per (species, sample): the counts are per sample.
        seen, tp, fp, fn = set(), 0, 0, 0
        for l in ls:
            if (l["species"], l["sample"]) not in seen:
                seen.add((l["species"], l["sample"]))
                tp, fp, fn = tp + l["TP"], fp + l["FP"], fn + l["FN"]
        out.append(f"| {run} | {len(ls)} | {frac([l['recovered'] for l in ls])} | "
                   f"{statistics.mean(l['resolved'] for l in ls):.3f} | {statistics.mean(l['wrong'] for l in ls):.3f} | "
                   f"{frac([l['placed'] for l in ls])} | {frac([l['kept_by_qcmsa'] for l in ls])} | "
                   f"{tp / (tp + fp) if tp + fp else float('nan'):.4f} | {tp / (tp + fn) if tp + fn else float('nan'):.4f} | "
                   f"{falsely[run]} |")
    out += ["", "By design and coverage (strains recovered):", "",
            "| run | design | strain share | coverage < 8x | 8-20x | > 20x |", "|---|---|---|---|---|---|"]
    for run, ls in sorted(by_run.items()):
        for design in sorted({(l["strains"], l["share"]) for l in ls}, key=lambda d: (d[0], -d[1])):
            sel = [l for l in ls if (l["strains"], l["share"]) == design]
            bins = [[l["recovered"] for l in sel if l["coverage"] < 8], [l["recovered"] for l in sel if 8 <= l["coverage"] <= 20],
                    [l["recovered"] for l in sel if l["coverage"] > 20]]
            out.append(f"| {run} | {design[0]} strains | {design[1]} | " + " | ".join(
                f"{frac(b)} ({len(b)})" for b in bins) + " |")
    text = "\n".join(out) + "\n"
    with open(os.path.join(OUT, "summary.md"), "w") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
