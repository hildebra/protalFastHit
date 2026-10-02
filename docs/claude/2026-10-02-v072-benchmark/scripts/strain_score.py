#!/usr/bin/env python3
"""strain_score.py [BENCH_DIR] - scores the strain benchmark's runs (strain_runs.sh) against the known strains.

Truth (make_strains.py): each of the 12 species has 8 strains, strain j in sample j, grown from the representative
(the databases' reference) by substitutions along a known tree, so a strain's genes are the representative's
coordinates (world marker_positions.tsv) in its own genome.

Per run and species with an MSA (<species>.raw.msa.fna, protal's own, and <species>.msa.fna, after its qcmsa):
- MSA length: columns of the raw MSA and of the filtered one, their sample rows and genes;
- SNP calling, from the raw MSA (qcmsa drops columns, so its columns no longer map to the genes' positions): at each
  gene position (a column where the reference row has a base), a sample's base is a called SNP when it is A/C/G/T
  and differs from the reference; TP when it is the strain's base there, FP (a miscalled SNP) when the strain has the
  reference's base there or another; FN a strain's SNP not called so (the reference's base, N, an IUPAC code or a gap).
  Precision, recall, F1; miscalled SNPs per 100 kb of called sequence (positions with A/C/G/T);
- phylogeny: pairwise p-distances of the sample rows (columns where both have A/C/G/T), a neighbour-joining tree, its
  splits against the true tree's (pruned to the samples the MSA has): the Robinson-Foulds distance (normalised to
  0-1), whether the topology is the true one, and the Spearman correlation of the p-distances with the true tree's
  patristic distances; for the filtered and for the raw MSA.
Writes BENCH_DIR/results_strains/species.tsv (one row per run and species), summary.md.
"""
import collections
import csv
import glob
import gzip
import math
import os
import re
import statistics
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
S = os.path.join(B, "strains")
RUNS = os.path.join(B, "strain_runs")
OUT = os.path.join(B, "results_strains")
W = os.path.join(B, "world", "full", "simulation")
VERSIONS = {"v060": "0.6.0a", "v070": "0.7.0", "v071": "0.7.1", "v072": "0.7.2"}
ACGT = set("ACGT")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    out, name, seq = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(seq)
                name, seq = line[1:].strip().split()[0], []
            else:
                seq.append(line.strip().upper())
    if name is not None:
        out[name] = "".join(seq)
    return out


def contig_order(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]


def parse_partition(path):
    """The blocks (name, first, last column, 1-based inclusive). 0.6.0a writes 0-based columns (its first block starts
    at 0), 0.7 1-based ones."""
    blocks = []
    with open(path) as fh:
        for line in fh:
            m = re.match(r"\s*\w+\s*,\s*(\S+)\s*=\s*(\d+)\s*-\s*(\d+)", line)
            if m:
                blocks.append((m.group(1), int(m.group(2)), int(m.group(3))))
    if blocks and min(lo for _, lo, _ in blocks) == 0:
        blocks = [(name, lo + 1, hi + 1) for name, lo, hi in blocks]
    return blocks


def parse_newick(text):
    """{leaf: path of (node id, branch length) to the root} and the splits as frozensets of leaves."""
    pos = 0
    counter = [0]
    children, length, leaves_of = {}, {}, {}

    def node():
        nonlocal pos
        kids = []
        if text[pos] == "(":
            pos += 1
            while True:
                kids.append(node())
                if text[pos] == ",":
                    pos += 1
                    continue
                pos += 1  # ')'
                break
            counter[0] += 1
            me = f"n{counter[0]}"
        else:
            m = re.match(r"[^:,();]+", text[pos:])
            me = m.group(0)
            pos += len(me)
        if pos < len(text) and text[pos] == ":":
            m = re.match(r":([0-9.eE+-]+)", text[pos:])
            length[me] = float(m.group(1))
            pos += len(m.group(0))
        children[me] = kids
        return me

    root = node()
    parent = {c: p for p, cs in children.items() for c in cs}

    def leaves(k):
        if k not in leaves_of:
            leaves_of[k] = frozenset([k]) if not children[k] else frozenset().union(*(leaves(c) for c in children[k]))
        return leaves_of[k]

    leaves(root)
    tips = [k for k in children if not children[k]]
    return tips, parent, length, leaves_of, root


def patristic(tips, parent, length):
    def path(k):
        out = {}
        d = 0.0
        while k in parent:
            out[k] = d
            d += length.get(k, 0.0)
            k = parent[k]
        out[k] = d
        return out
    paths = {t: path(t) for t in tips}
    dist = {}
    for a in tips:
        for b in tips:
            if a < b:
                common = min((paths[a][k] + paths[b][k]) for k in paths[a] if k in paths[b])
                dist[(a, b)] = dist[(b, a)] = common
    return dist


def splits_of(leaf_sets, leaves):
    """Non-trivial unrooted splits of the leaf set `leaves`, each as the side without the smallest leaf."""
    leaves = frozenset(leaves)
    first = min(leaves)
    out = set()
    for s in leaf_sets:
        s = s & leaves
        if 1 < len(s) < len(leaves) - 1:
            out.add(s if first not in s else leaves - s)
    return out


def nj(names, d):
    """Neighbour-joining tree of names with distances d[(a, b)]; its clusters (leaf sets of the joined nodes)."""
    nodes = list(names)
    leafset = {n: frozenset([n]) for n in nodes}
    dist = {(a, b): d[(a, b)] for a in nodes for b in nodes if a != b}
    clusters = []
    k = 0
    while len(nodes) > 3:
        n = len(nodes)
        r = {a: sum(dist[(a, b)] for b in nodes if b != a) for a in nodes}
        best, pair = None, None
        for i, a in enumerate(nodes):
            for b in nodes[i + 1:]:
                q = (n - 2) * dist[(a, b)] - r[a] - r[b]
                if best is None or q < best:
                    best, pair = q, (a, b)
        a, b = pair
        k += 1
        u = f"_nj{k}"
        leafset[u] = leafset[a] | leafset[b]
        clusters.append(leafset[u])
        for c in nodes:
            if c not in (a, b):
                dist[(u, c)] = dist[(c, u)] = (dist[(a, c)] + dist[(b, c)] - dist[(a, b)]) / 2
        nodes = [c for c in nodes if c not in (a, b)] + [u]
    return clusters


def spearman(x, y):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for t in range(i, j + 1):
                r[order[t]] = (i + j) / 2
            i = j + 1
        return r
    if len(x) < 3:
        return float("nan")
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def phylogeny(rows, true_tree, sample_strain):
    """RF (normalised), exact topology, Spearman of p-distances with the true patristic distances; None if fewer than
    4 sample rows."""
    names = [n for n in rows if n in sample_strain]
    if len(names) < 4:
        return None
    tips, parent, length, leaves_of, root = true_tree
    d = {}
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            both = diff = 0
            for x, y in zip(rows[a], rows[b]):
                if x in ACGT and y in ACGT:
                    both += 1
                    diff += x != y
            d[(a, b)] = d[(b, a)] = diff / both if both else 0.75
    est = splits_of(nj(names, d), names)
    strain_sample = {sample_strain[n]: n for n in names}
    true_sets = [frozenset(strain_sample[t] for t in s if t in strain_sample) for s in leaves_of.values()]
    true = splits_of(true_sets, names)
    rf = len(est ^ true) / (2 * (len(names) - 3)) if len(names) > 3 else 0.0
    pat = patristic(tips, parent, length)
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
    rho = spearman([d[p] for p in pairs], [pat[(sample_strain[p[0]], sample_strain[p[1]])] for p in pairs])
    return dict(rf=rf, exact=int(rf == 0), rho=rho, leaves=len(names))


def main():
    strains = collections.defaultdict(dict)  # species -> sample -> strain row
    rep_path = {}
    with open(os.path.join(W, "genomes.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rep_path[row["accession"]] = row["fasta_path"]
    with open(os.path.join(S, "strains.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["role"] == "strain":
                strains[row["species"]][row["sample"]] = row
    trees = {}
    with open(os.path.join(S, "trees.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            trees[row["species"]] = parse_newick(row["newick"])
    markers = collections.defaultdict(dict)  # rep accession -> marker -> row
    with open(os.path.join(W, "marker_positions.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            markers[row["accession"]][row["marker"]] = row
    marker_of = {}
    with open(os.path.join(B, "V071", "protal_db", "gene2geneid.tsv")) as fh:
        for line in fh:
            name, gid = line.rstrip("\n").split("\t")[:2]
            marker_of[gid] = name
    genomes = {}

    def strain_genes(species, sample):
        """{geneid: the strain's gene, oriented as the reference}; and the representative's, for the check."""
        row = strains[species][sample]
        rep = row["representative"]
        if (species, sample) not in genomes:
            order = contig_order(rep_path[rep])
            seqs = read_fasta(row["genome"])
            index = {name: f"{row['strain']}_contig{i + 1}" for i, name in enumerate(order)}
            genomes[(species, sample)] = (seqs, index)
        seqs, index = genomes[(species, sample)]
        return seqs, index, markers[rep]

    rep_genomes = {}

    def gene_of(seqs, index, m, contig_name=None):
        seq = seqs[index[m["contig"]] if contig_name is None else contig_name]
        g = seq[int(m["start"]) - 1:int(m["end"])]
        return g.translate(COMPLEMENT)[::-1] if m["strand"] == "-" else g

    rows_out = []
    for run_dir in sorted(d for d in glob.glob(os.path.join(RUNS, "*")) if os.path.isdir(d)):
        run = os.path.basename(run_dir)
        if not os.path.exists(run_dir + ".done"):
            continue
        variant, reads = run.split(".")
        for raw in sorted(glob.glob(os.path.join(run_dir, "**", "*.raw.msa.fna"), recursive=True)):
            base = raw[:-len(".raw.msa.fna")]
            species = os.path.basename(base).removeprefix("s__").replace("_", " ")
            if species not in strains:
                continue  # a background species
            msa = read_fasta(raw)
            ref_name = next((n for n in msa if n.endswith("_reference")), None)
            if ref_name is None:
                continue
            ref = msa[ref_name]
            sample_rows = {n: s for n, s in msa.items() if n in strains[species]}
            sample_strain = {n: strains[species][n]["strain"] for n in sample_rows}
            blocks = parse_partition(base + ".raw.partition.txt")
            tp = fp = fn = called = positions = genes_checked = genes_mismatch = 0
            rep_acc = next(iter(strains[species].values()))["representative"]
            if rep_acc not in rep_genomes:
                rep_genomes[rep_acc] = read_fasta(rep_path[rep_acc])
            for name, lo, hi in blocks:
                gid = re.sub(r"\D", "", name)
                marker = marker_of.get(gid)
                m = markers[rep_acc].get(marker) if marker else None
                if m is None:
                    genes_mismatch += 1
                    continue
                cols = [i for i in range(lo - 1, hi) if ref[i] != "-"]
                ref_gene = "".join(ref[i] for i in cols)
                if gene_of(rep_genomes[rep_acc], None, m, m["contig"]) != ref_gene:
                    genes_mismatch += 1
                    continue
                genes_checked += 1
                for sample, row in sample_rows.items():
                    seqs, index, mk = strain_genes(species, sample)
                    truth = gene_of(seqs, index, mk[marker])
                    for k, i in enumerate(cols):
                        c, r, t = row[i], ref_gene[k], truth[k]
                        positions += 1
                        if c in ACGT:
                            called += 1
                        snp = c in ACGT and c != r
                        if snp and c == t:
                            tp += 1
                        elif snp:
                            fp += 1
                        if t != r and not (snp and c == t):
                            fn += 1
            final = base + ".msa.fna"
            final_msa = read_fasta(final) if os.path.exists(final) else {}
            final_rows = {n: s for n, s in final_msa.items() if n in strains[species]}
            final_blocks = parse_partition(base + ".partition.txt") if os.path.exists(base + ".partition.txt") else []
            ph_final = phylogeny(final_rows, trees[species], sample_strain) if final_rows else None
            ph_raw = phylogeny(sample_rows, trees[species], sample_strain)
            prec = tp / (tp + fp) if tp + fp else float("nan")
            rec = tp / (tp + fn) if tp + fn else float("nan")
            rows_out.append(dict(
                run=run, variant=variant, reads=reads, species=species, raw_columns=len(ref), raw_rows=len(sample_rows),
                raw_genes=len(blocks), genes_checked=genes_checked, genes_not_checked=genes_mismatch,
                final_columns=len(next(iter(final_rows.values()))) if final_rows else 0, final_rows=len(final_rows),
                final_genes=len(final_blocks), positions=positions, called=called, TP=tp, FP=fp, FN=fn,
                precision=prec, recall=rec, F1=2 * prec * rec / (prec + rec) if tp else 0.0,
                miscalled_per_100kb=1e5 * fp / called if called else float("nan"),
                rf_final=ph_final["rf"] if ph_final else float("nan"),
                exact_final=ph_final["exact"] if ph_final else float("nan"),
                rho_final=ph_final["rho"] if ph_final else float("nan"),
                leaves_final=ph_final["leaves"] if ph_final else 0,
                rf_raw=ph_raw["rf"] if ph_raw else float("nan"), exact_raw=ph_raw["exact"] if ph_raw else float("nan"),
                rho_raw=ph_raw["rho"] if ph_raw else float("nan")))
            print(f"{run} {species}: TP {tp} FP {fp} FN {fn}, final RF {rows_out[-1]['rf_final']}", flush=True)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "species.tsv"), "w") as fh:
        cols = list(rows_out[0])
        fh.write("\t".join(cols) + "\n")
        for r in rows_out:
            fh.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    lines = ["# The strain benchmark (strain_score.py)", "",
             "12 species, 8 strains each (one per sample), against each version's finished database. Means over the "
             "species with an MSA; SNPs summed over them.", "",
             "| reads | version | species with an MSA | raw MSA columns | filtered MSA columns | filtered rows | "
             "SNP precision | SNP recall | SNP F1 | miscalled SNPs per 100 kb called | true topology (filtered MSA) | "
             "RF (filtered) | Spearman (filtered) | true topology (raw MSA) | RF (raw) | wall (s) | peak RSS (GB) |",
             "|" + "---|" * 17]
    def mean(v):
        v = [x for x in v if x == x]
        return statistics.mean(v) if v else float("nan")
    for reads in ("pe", "pb", "ont"):
        for variant in VERSIONS:
            rs = [r for r in rows_out if r["variant"] == variant and r["reads"] == reads]
            run = f"{variant}.{reads}"
            if not os.path.exists(os.path.join(RUNS, run + ".done")):
                continue
            tp, fp, fn = sum(r["TP"] for r in rs), sum(r["FP"] for r in rs), sum(r["FN"] for r in rs)
            called = sum(r["called"] for r in rs)
            prec = tp / (tp + fp) if tp + fp else float("nan")
            rec = tp / (tp + fn) if tp + fn else float("nan")
            f1 = 2 * prec * rec / (prec + rec) if tp else 0.0
            wall = rss = float("nan")
            with open(os.path.join(RUNS, run + ".time")) as fh:
                for line in fh:
                    key, _, value = line.strip().rpartition(": ")
                    if key.startswith("Elapsed (wall clock)"):
                        wall = sum(float(p) * 60 ** i for i, p in enumerate(reversed(value.split(":"))))
                    elif key == "Maximum resident set size (kbytes)":
                        rss = int(value) / 1024 ** 2
            fmt = lambda x, d=3: "-" if x != x else f"{x:.{d}f}"
            lines.append(f"| {reads} | {VERSIONS[variant]} | {len(rs)} of 12 | {fmt(mean([r['raw_columns'] for r in rs]), 0)} | "
                         f"{fmt(mean([r['final_columns'] for r in rs]), 0)} | {fmt(mean([r['final_rows'] for r in rs]), 1)} | "
                         f"{fmt(prec, 4)} | {fmt(rec, 4)} | {fmt(f1, 4)} | {fmt(1e5 * fp / called if called else float('nan'), 1)} | "
                         f"{fmt(mean([r['exact_final'] for r in rs]), 2)} | {fmt(mean([r['rf_final'] for r in rs]), 3)} | "
                         f"{fmt(mean([r['rho_final'] for r in rs]), 3)} | {fmt(mean([r['exact_raw'] for r in rs]), 2)} | "
                         f"{fmt(mean([r['rf_raw'] for r in rs]), 3)} | {fmt(wall, 0)} | {fmt(rss, 2)} |")
    text = "\n".join(lines) + "\n"
    with open(os.path.join(OUT, "summary.md"), "w") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
