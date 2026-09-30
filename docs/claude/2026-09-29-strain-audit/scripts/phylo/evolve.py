#!/usr/bin/env python3
"""Evolve a species' representative genome along a random tree to make a known strain phylogeny.

Model: K80 substitutions (kappa) with Gamma(shape) site-rate heterogeneity, plus a few small
indels (rate = indel_ratio x substitution rate per site, length 1-6, geometric). The root is the
representative genome (= the protal database reference). Tree: Yule topology with branch lengths
jittered by a lognormal rate factor (non-clocklike), scaled so the MEAN root-to-tip distance is
--height substitutions/site.

Outputs (in --outdir):
  <tip>.fna.gz          tip genomes
  true.nwk              true tree, tips named by tip (= sample) name, lengths in subst/site
  true_with_ref.nwk     same, plus the reference as a zero-length tip at the root
  truth_markers.tsv     tip, gene_id, marker, aligned tip sequence in reference gene coordinates
                        (gene orientation; '-' = deleted base), n_insertions
  genomes.tsv           genome table rows (accession, taxonomy, fasta_path, genome_length, rep)
  evolve_log.tsv        realised substitutions/indels per branch
"""
import argparse
import gzip
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phylo_utils as pu  # noqa: E402

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
CODE = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate(b"ACGT"):
    CODE[b] = i
    CODE[b + 32] = i
# transitions: A<->G (0<->2), C<->T (1<->3)
TRANSITION = np.array([2, 3, 0, 1], dtype=np.uint8)
TRANSV = np.array([[1, 3], [0, 2], [1, 3], [0, 2]], dtype=np.uint8)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "-": "-", "N": "N"}


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    names, seqs, cur = [], [], []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line[0] == ">":
                if names:
                    seqs.append("".join(cur))
                names.append(line[1:].split()[0])
                cur = []
            else:
                cur.append(line)
    if names:
        seqs.append("".join(cur))
    return names, seqs


def revcomp(s):
    return "".join(COMP.get(c, "N") for c in reversed(s))


def yule_tree(ntips, rng, tipnames):
    """Pure-birth tree: returns root Node with branch lengths in time units."""
    # simulate forward: lineages split at rate 1 each
    root = pu.Node()
    active = [root]
    t_start = {id(root): 0.0}
    t = 0.0
    while len(active) < ntips:
        t += rng.exponential(1.0 / len(active))
        k = rng.integers(len(active))
        parent = active.pop(k)
        parent.length = t - t_start[id(parent)] if parent is not root else None
        for _ in range(2):
            c = pu.Node()
            parent.add(c)
            t_start[id(c)] = t
            active.append(c)
    t += rng.exponential(1.0 / len(active))
    for i, leaf in enumerate(active):
        leaf.length = t - t_start[id(leaf)]
    order = rng.permutation(len(active))
    for leaf, j in zip(active, order):
        leaf.name = tipnames[j]
    return root


def root_to_tip(root):
    out = {}

    def rec(nd, d):
        if nd.is_leaf():
            out[nd.name] = d
        for c in nd.children:
            rec(c, d + c.length)
    rec(root, 0.0)
    return out


def mutate(seq, origin, rates_root, bl, kappa, indel_ratio, rng, shape):
    """Evolve one contig along a branch of length bl. seq: uint8 codes 0-3, origin: int64 root coords
    (-1 = inserted). Returns (seq, origin, n_subs, n_ins, n_del)."""
    L = len(seq)
    r = np.where(origin >= 0, rates_root[np.clip(origin, 0, None)], 1.0)
    lam = bl * r
    nev = rng.poisson(lam)
    idx = np.nonzero(nev)[0]
    n_subs = int(nev.sum())
    seq = seq.copy()
    p_ts = kappa / (kappa + 2.0)
    for _round in range(int(nev.max()) if len(idx) else 0):
        sel = idx[nev[idx] > _round]
        u = rng.random(len(sel))
        cur = seq[sel]
        ts = u < p_ts
        which = rng.integers(0, 2, len(sel))
        new = np.where(ts, TRANSITION[cur], TRANSV[cur, which])
        seq[sel] = new
    # indels
    n_indel = rng.poisson(indel_ratio * bl * L)
    n_ins = n_del = 0
    if n_indel:
        seq_l = seq
        org_l = origin
        for _ in range(n_indel):
            ln = min(6, int(rng.geometric(0.6)))
            pos = int(rng.integers(0, len(seq_l) - ln))
            if rng.random() < 0.5:
                ins = rng.integers(0, 4, ln).astype(np.uint8)
                seq_l = np.concatenate([seq_l[:pos], ins, seq_l[pos:]])
                org_l = np.concatenate([org_l[:pos], np.full(ln, -1, dtype=np.int64), org_l[pos:]])
                n_ins += 1
            else:
                seq_l = np.concatenate([seq_l[:pos], seq_l[pos + ln:]])
                org_l = np.concatenate([org_l[:pos], org_l[pos + ln:]])
                n_del += 1
        seq, origin = seq_l, org_l
    return seq, origin, n_subs, n_ins, n_del


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--accession", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--markers", required=True, help="marker_positions.tsv")
    ap.add_argument("--gene2id", required=True, help="protal_db gene2geneid.tsv")
    ap.add_argument("--dbref", required=True, help="protal_db full_reference.fna")
    ap.add_argument("--tiid", required=True, help="protal taxon id of the species in full_reference.fna")
    ap.add_argument("--ntips", type=int, default=12)
    ap.add_argument("--height", type=float, default=0.01)
    ap.add_argument("--kappa", type=float, default=3.0)
    ap.add_argument("--shape", type=float, default=0.8)
    ap.add_argument("--rate_sd", type=float, default=0.4, help="lognormal sd of per-branch rate factor")
    ap.add_argument("--indel_ratio", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tip_prefix", default="t")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--refname", required=True, help="MSA reference row name, e.g. s__Mockella_alpha_reference")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    names, seqs = read_fasta(args.genome)
    root_seqs = []
    for s in seqs:
        a = CODE[np.frombuffer(s.encode(), dtype=np.uint8)]
        if (a == 255).any():
            bad = (a == 255)
            a[bad] = rng.integers(0, 4, bad.sum())
        root_seqs.append(a.astype(np.uint8))
    rates = [rng.gamma(args.shape, 1.0 / args.shape, len(s)) for s in root_seqs]

    # --- markers of this genome, checked against the DB genes ---
    gene2id = {}
    with open(args.gene2id) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2:
                gene2id[f[0]] = f[1]
    dbseq = {}
    cur = None
    with open(args.dbref) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                nm = line[1:].split()[0]
                cur = nm if nm not in dbseq else None
                if cur is not None:
                    dbseq[cur] = []
            elif cur is not None:
                dbseq[cur].append(line)
    dbseq = {k: "".join(v) for k, v in dbseq.items()}
    contig_idx = {n: i for i, n in enumerate(names)}
    markers = []
    with open(args.markers) as fh:
        header = fh.readline()
        for line in fh:
            acc, marker, contig, start, end, strand = line.rstrip("\n").split("\t")
            if acc != args.accession:
                continue
            gid = gene2id.get(marker)
            if gid is None:
                continue
            markers.append((gid, marker, contig_idx[contig], int(start) - 1, int(end), strand))
    ok = mism = missing = 0
    rootstr = ["".join(chr(BASES[x]) for x in s) for s in root_seqs]
    for gid, marker, ci, s, e, strand in markers:
        g = rootstr[ci][s:e]
        if strand == "-":
            g = revcomp(g)
        ref = dbseq.get(f"{args.tiid}_{gid}")
        if ref is None:
            missing += 1
        elif ref.upper() == g:
            ok += 1
        else:
            mism += 1
    sys.stderr.write(f"markers: {len(markers)}; match DB gene copy1: {ok}, mismatch: {mism}, not in DB: {missing}\n")

    # --- tree ---
    tipnames = [f"{args.tip_prefix}{i+1:02d}" for i in range(args.ntips)]
    root = yule_tree(args.ntips, rng, tipnames)
    for nd in root.nodes():
        if nd is not root:
            nd.length *= float(np.exp(rng.normal(0, args.rate_sd)))
    rtt = root_to_tip(root)
    scale = args.height / (sum(rtt.values()) / len(rtt))
    for nd in root.nodes():
        if nd is not root:
            nd.length *= scale
    with open(os.path.join(args.outdir, "true.nwk"), "w") as fh:
        fh.write(pu.to_newick(root) + "\n")
    # reference = root sequence: attach it as a zero-length tip at the root
    import copy
    r2 = copy.deepcopy(root)
    pu.fix_parents(r2)
    r2.add(pu.Node(args.refname, 0.0))
    with open(os.path.join(args.outdir, "true_with_ref.nwk"), "w") as fh:
        fh.write(pu.to_newick(r2) + "\n")

    # --- evolve ---
    log = open(os.path.join(args.outdir, "evolve_log.tsv"), "w")
    log.write("node\tbranch_length\tsubs\tins\tdel\tgenome_len\n")
    tips = {}
    counter = [0]

    def rec(nd, state):
        if nd is not root:
            new = []
            ts = ti = td = 0
            for ci, (sq, og) in enumerate(state):
                sq2, og2, a, b, c = mutate(sq, og, rates[ci], nd.length, args.kappa, args.indel_ratio, rng, args.shape)
                new.append((sq2, og2))
                ts += a
                ti += b
                td += c
            state = new
            counter[0] += 1
            lab = nd.name if nd.is_leaf() else f"internal{counter[0]}"
            log.write(f"{lab}\t{nd.length:.6g}\t{ts}\t{ti}\t{td}\t{sum(len(s) for s,_ in state)}\n")
        if nd.is_leaf():
            tips[nd.name] = state
        for c in nd.children:
            rec(c, state)

    rec(root, [(s, np.arange(len(s), dtype=np.int64)) for s in root_seqs])
    log.close()

    gt = open(os.path.join(args.outdir, "genomes.tsv"), "w")
    tm = open(os.path.join(args.outdir, "truth_markers.tsv"), "w")
    tm.write("tip\tgene_id\tmarker\tseq\tn_ins\n")
    for tip in tipnames:
        state = tips[tip]
        acc = f"{args.accession.replace('GCF_', 'SIM_')}_{tip}"
        path = os.path.join(os.path.abspath(args.outdir), f"{acc}.fna.gz")
        glen = 0
        with gzip.open(path, "wt", compresslevel=3) as fh:
            for ci, (sq, og) in enumerate(state):
                s = BASES[sq].tobytes().decode()
                glen += len(s)
                fh.write(f">{acc}_contig{ci+1}\n")
                for i in range(0, len(s), 80):
                    fh.write(s[i:i + 80] + "\n")
        gt.write(f"{acc}\t{args.taxonomy}\t{path}\t{glen}\tf\n")
        # truth per marker in reference coordinates
        for gid, marker, ci, s, e, strand in markers:
            sq, og = state[ci]
            inv = np.full(e - s, -1, dtype=np.int64)
            m = (og >= s) & (og < e)
            pos = np.nonzero(m)[0]
            inv[og[pos] - s] = pos
            chars = np.where(inv >= 0, BASES[sq[np.clip(inv, 0, None)]], ord("-")).astype(np.uint8).tobytes().decode()
            # insertions inside the gene: inserted bases between first and last mapped positions
            if len(pos):
                lo, hi = pos.min(), pos.max()
                n_ins = int((og[lo:hi + 1] < 0).sum())
            else:
                n_ins = 0
            if strand == "-":
                chars = revcomp(chars)
            tm.write(f"{tip}\t{gid}\t{marker}\t{chars}\t{n_ins}\n")
    gt.close()
    tm.close()
    sys.stderr.write(f"wrote {len(tipnames)} tips to {args.outdir}; mean root-to-tip {args.height}\n")


if __name__ == "__main__":
    main()
