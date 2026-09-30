#!/usr/bin/env python3
"""Run the justfile's IQ-TREE command on MSAs and compare branch lengths with the simulated star tree.

usage: treeeval.py <manifest.tsv> <outdir> <label>=<msa> [...]
Truth: each distinct genome is a leaf at distance strain_divergence from the species ancestor; a leaf
whose genome is unique among the MSA rows should get a terminal branch of that length."""
import sys, os, re, subprocess, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from truth import load_truth, genome_of
from consistency import read_fasta

def leaf_lengths(newick):
    return {m.group(1): float(m.group(2)) for m in re.finditer(r"[(,]([^(),:;]+):([0-9.eE-]+)", newick)}

def run(msa, pref):
    os.makedirs(os.path.dirname(pref), exist_ok=True)
    cmd = ["iqtree2", "-s", msa, "-m", "GTR+G", "-B", "1000", "-T", "2", "--seqtype", "DNA", "--prefix", pref, "-redo"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    open(pref + ".stdout.log", "w").write(p.stdout + p.stderr)
    return p.returncode, p.stdout + p.stderr

def main(manifest, outdir, *items):
    div, rep, samp = load_truth(manifest)
    for item in items:
        label, msa = item.split("=", 1)
        sp = os.path.basename(msa).split(".")[0]
        names, seqs = read_fasta(msa)
        if len(names) < 4:
            print(f"{label:16} {sp:22} {len(names)} seqs: skipped (<4, as the justfile does)")
            continue
        pref = os.path.join(outdir, label, sp)
        rc, log = run(msa, pref)
        if rc != 0 or not os.path.exists(pref + ".treefile"):
            err = [l for l in log.splitlines() if "ERROR" in l][:2]
            print(f"{label:16} {sp:22} IQ-TREE rc={rc} {err}")
            continue
        tree = open(pref + ".treefile").read()
        tl = float(re.search(r"Total tree length \(sum of branch lengths\): ([0-9.]+)", open(pref + ".iqtree").read()).group(1))
        m = re.search(r"(\d+) parsimony-informative, (\d+) singleton sites, (\d+) constant sites", log)
        genomes = {n: genome_of(n, sp, rep, samp) for n in names}
        distinct = {tuple(g) for g in genomes.values()}
        true_tl = sum(div[g[0]] for g in distinct if len(g) == 1 and g[0] in div)
        ll = leaf_lengths(tree)
        uniq = [n for n in names if list(genomes.values()).count(genomes[n]) == 1 and len(genomes[n]) == 1 and genomes[n][0] in div]
        ratios = [ll[n] / div[genomes[n][0]] for n in uniq if n in ll]
        rs = f"{min(ratios):.2f}-{max(ratios):.2f} (mean {sum(ratios)/len(ratios):.2f})" if ratios else "-"
        warns = sorted({l.strip()[:90] for l in log.splitlines() if l.startswith("WARNING") or "identical" in l})
        print(f"{label:16} {sp:22} n={len(names)} sites PI/singleton/const={m.groups() if m else '?'} tree_len={tl:.4f} "
              f"true={true_tl:.4f} ({tl/true_tl if true_tl else float('nan'):.2f}x) unique-leaf terminal/true={rs}")
        for w in warns:
            print(f"{'':16}   {w}")

if __name__ == "__main__":
    main(*sys.argv[1:])
