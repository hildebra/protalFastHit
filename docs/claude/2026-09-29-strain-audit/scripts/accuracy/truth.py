#!/usr/bin/env python3
"""truth.py - true marker-gene sequences of every genome, aligned to the representative's genes.

For every genome of the synthetic world and every marker gene it carries, extract the gene from
the genome FASTA (simulation/marker_positions.tsv), check it against the world's marker FASTA
(genomic_files_all/*/fna), and globally align it (affine-gap Needleman-Wunsch, numpy row sweep)
to the species representative's copy of that gene (the DB reference).

Output: truth.pkl with
  genes[(species, geneid)] = reference gene sequence (the rep's)
  truth[(accession, geneid)] = dict(tb=list of true base per ref position ('-' = deleted),
                                    ins=dict ref_pos -> inserted string placed BEFORE ref_pos
                                        (protal's convention: an insertion recorded at rpos),
                                    seq=true gene sequence, n_snp, n_ins, n_del)
  meta: genome -> species, species -> rep accession, domain, geneid <-> marker
Usage: truth.py [out.pkl]
"""
import gzip, os, pickle, sys
import numpy as np

W = os.path.expanduser("~/audit5/world/gtdb_r226")
DB = os.path.expanduser("~/audit5/world/protal_db")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/audit5/accuracy/truth.pkl")

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, name, cur = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    seqs[name] = "".join(cur)
                name = line[1:].split()[0]
                cur = []
            else:
                cur.append(line.upper())
    if name is not None:
        seqs[name] = "".join(cur)
    return seqs


# --- scoring for the global alignment
MATCH, MISMATCH, GO, GE = 2, -4, 10, 1  # gap of length L costs GO + (L-1)*GE
NEG = -10 ** 9


def align(a, b):
    """Wrapper storing the G-source table needed during traceback."""
    n, m = len(a), len(b)
    A = np.frombuffer(a.encode(), dtype=np.uint8)
    B = np.frombuffer(b.encode(), dtype=np.uint8)
    jj = np.arange(m + 1, dtype=np.int64)
    Hptr = np.zeros((n + 1, m + 1), dtype=np.int8)
    Gsrc_t = np.zeros((n + 1, m + 1), dtype=np.int8)
    Fext = np.zeros((n + 1, m + 1), dtype=bool)
    Eext = np.zeros((n + 1, m + 1), dtype=bool)
    H = np.empty(m + 1, dtype=np.int64)
    F = np.full(m + 1, NEG, dtype=np.int64)
    H[0] = 0
    H[1:] = -(GO + GE * (jj[1:] - 1))
    Hptr[0, 1:] = 2
    Eext[0, 2:] = True
    for i in range(1, n + 1):
        sc = np.where(B == A[i - 1], MATCH, MISMATCH).astype(np.int64)
        diag = np.full(m + 1, NEG, dtype=np.int64)
        diag[1:] = H[:-1] + sc
        Fo = H - GO
        Fe = F - GE
        Fnew = np.maximum(Fo, Fe)
        Fext[i] = Fe > Fo
        G = np.maximum(diag, Fnew)
        gs = (Fnew > diag).astype(np.int8)
        Gsrc_t[i] = gs
        Acc = G + GE * jj
        C = np.maximum.accumulate(Acc)
        E = np.full(m + 1, NEG, dtype=np.int64)
        E[1:] = C[:-1] - GO - GE * (jj[1:] - 1)
        Eo = np.full(m + 1, NEG, dtype=np.int64)
        Eo[1:] = G[:-1] - GO
        Eprev = np.full(m + 1, NEG, dtype=np.int64)
        Eprev[1:] = E[:-1]
        Eext[i] = (Eprev - GE) > Eo
        Hptr[i] = np.where(E > G, 2, gs)
        H, F = np.maximum(G, E), Fnew
    score = int(H[m])
    i, j, state = n, m, int(Hptr[n, m])
    out = []
    while i > 0 or j > 0:
        if i == 0:
            state = 2
        elif j == 0:
            state = 1
        if state == 0:
            out.append((a[i - 1], b[j - 1]))
            i -= 1; j -= 1
            state = int(Hptr[i, j])
        elif state == 1:
            ext = bool(Fext[i, j])
            out.append((a[i - 1], "-"))
            i -= 1
            state = 1 if ext else int(Hptr[i, j])
        else:
            ext = bool(Eext[i, j])
            out.append(("-", b[j - 1]))
            j -= 1
            state = 2 if ext else int(Gsrc_t[i, j])
    out.reverse()
    return out, score


def truth_from_alignment(pairs, ref_len):
    tb = []
    ins = {}
    pending = []
    for r, q in pairs:
        if r == "-":
            pending.append(q)
        else:
            if pending:
                ins[len(tb)] = "".join(pending)
                pending = []
            tb.append(q)
    if pending:
        ins[len(tb)] = "".join(pending)
    assert len(tb) == ref_len
    return tb, ins


def left_normalise(ref, tb, ins):
    """Shift each deletion / insertion as far left as the sequence allows (both VCF and most
    aligners' convention), so that indel placement is canonical."""
    # deletions: runs of '-' in tb
    tb = list(tb)
    changed = True
    while changed:
        changed = False
        p = 0
        L = len(tb)
        while p < L:
            if tb[p] == "-":
                q = p
                while q < L and tb[q] == "-":
                    q += 1
                # run p..q-1 deleted; can shift left if ref[p-1] == ref[q-1] and tb[p-1]==ref[p-1]
                if p > 0 and tb[p - 1] == ref[p - 1] and ref[p - 1] == ref[q - 1] and (p) not in ins and q not in ins:
                    tb[p - 1] = "-"
                    tb[q - 1] = ref[q - 1]
                    changed = True
                    break
                p = q
            else:
                p += 1
    new_ins = {}
    for pos, s in sorted(ins.items()):
        # inserted s before ref pos; shift left while the base before equals last inserted base
        while pos > 0 and tb[pos - 1] != "-" and tb[pos - 1] == s[-1] and (pos - 1) not in new_ins and (pos - 1) not in ins:
            s = tb[pos - 1] + s[:-1]
            pos -= 1
        new_ins[pos] = s
    return tb, new_ins


def main():
    gene2id = {}
    with open(os.path.join(DB, "gene2geneid.tsv")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].isdigit():
                gene2id[f[0]] = int(f[1])
    genomes = {}
    with open(os.path.join(W, "simulation/genomes.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            tax = f["gtdb_taxonomy"]
            sp = tax.split(";s__")[-1]
            genomes[f["accession"]] = dict(species=sp, taxonomy=tax, fasta=f["fasta_path"],
                                           length=int(f["genome_length"]), rep=f["gtdb_representative"] == "t",
                                           domain="ar53" if "d__Archaea" in tax else "bac120")
    rep_of = {g["species"]: a for a, g in genomes.items() if g["rep"]}
    pos = {}
    with open(os.path.join(W, "simulation/marker_positions.tsv")) as fh:
        fh.readline()
        for line in fh:
            acc, marker, contig, s, e, strand = line.rstrip("\n").split("\t")
            pos[(acc, marker)] = (contig, int(s), int(e), strand)
    # marker FASTAs of all genomes
    mfa = {}
    for dom in ("bac120", "ar53"):
        d = os.path.join(W, f"genomic_files_all/{dom}_marker_genes_all_r226/fna")
        for fn in os.listdir(d):
            marker = fn[len(dom) + 1:-4]
            for name, seq in read_fasta(os.path.join(d, fn)).items():
                mfa[(dom, marker, name[3:])] = seq
    genes, truth = {}, {}
    mism_marker_fasta = 0
    checked = 0
    for acc, g in sorted(genomes.items()):
        gseq = read_fasta(g["fasta"])
        for (a2, marker), (contig, s, e, strand) in pos.items():
            if a2 != acc:
                continue
            seq = gseq[contig][s - 1:e]
            if strand == "-":
                seq = revcomp(seq)
            ref_fa = mfa.get((g["domain"], marker, acc))
            checked += 1
            if ref_fa != seq:
                mism_marker_fasta += 1
            gid = gene2id[marker]
            if g["rep"]:
                genes[(g["species"], gid)] = seq
            g.setdefault("genes", {})[gid] = seq
    print(f"marker genes checked against marker FASTAs: {checked}, mismatching: {mism_marker_fasta}", file=sys.stderr)
    stats = dict(n=0, ungapped=0, nw=0)
    for acc, g in sorted(genomes.items()):
        sp = g["species"]
        for gid, seq in sorted(g["genes"].items()):
            ref = genes.get((sp, gid))
            if ref is None:
                continue
            stats["n"] += 1
            if g["rep"]:
                assert seq == ref
                truth[(acc, gid)] = dict(tb=list(ref), ins={}, seq=seq, n_snp=0, n_ins=0, n_del=0)
                continue
            use_nw = len(seq) != len(ref)
            if not use_nw:
                mm = np.frombuffer(seq.encode(), np.uint8) != np.frombuffer(ref.encode(), np.uint8)
                # a local frame shift shows as a mismatch cluster
                win = np.convolve(mm.astype(int), np.ones(15, int), mode="valid")
                use_nw = mm.mean() > 0.05 or (win.size and win.max() >= 4)
            if use_nw:
                pairs, _ = align(ref, seq)
                tb, ins = truth_from_alignment(pairs, len(ref))
                tb, ins = left_normalise(ref, tb, ins)
                stats["nw"] += 1
            else:
                tb, ins = list(seq), {}
                stats["ungapped"] += 1
            # consistency: rebuild the sequence
            rebuilt = []
            for p in range(len(ref)):
                if p in ins:
                    rebuilt.append(ins[p])
                if tb[p] != "-":
                    rebuilt.append(tb[p])
            if len(ref) in ins:
                rebuilt.append(ins[len(ref)])
            assert "".join(rebuilt) == seq, (acc, gid)
            n_snp = sum(1 for p in range(len(ref)) if tb[p] != "-" and tb[p] != ref[p])
            n_del = sum(1 for p in range(len(ref)) if tb[p] == "-" and (p == 0 or tb[p - 1] != "-"))
            truth[(acc, gid)] = dict(tb=tb, ins=ins, seq=seq, n_snp=n_snp, n_ins=len(ins), n_del=n_del)
    print(f"alignments: {stats}", file=sys.stderr)
    meta = dict(genomes={a: {k: v for k, v in g.items() if k != "genes"} for a, g in genomes.items()},
                rep_of=rep_of, gene2id=gene2id)
    with open(OUT, "wb") as fh:
        pickle.dump(dict(genes=genes, truth=truth, meta=meta), fh)
    # summary per genome
    for acc, g in sorted(genomes.items()):
        rows = [t for (a, _), t in truth.items() if a == acc]
        L = sum(len(t["tb"]) for t in rows)
        s = sum(t["n_snp"] for t in rows)
        i = sum(t["n_ins"] for t in rows)
        d = sum(t["n_del"] for t in rows)
        print(f"{acc}\t{g['species']}\tgenes={len(rows)}\tref_bp={L}\tsnps={s}\t({100.0 * s / max(L, 1):.2f}%)\tins={i}\tdel={d}")


if __name__ == "__main__":
    main()
