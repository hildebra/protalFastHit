"""Are the column weights' alignments right? Their copies' mappings onto the family columns (column_weights.tsv)
against GTDB's own alignments of the same proteins (the profile-HMM alignments of GTDB-Tk,
<set>_msa_marker_genes_reps_r226.tar.gz: each representative's protein, its insertions trimmed, on the marker's
columns): for every pair of representative copies of one family, the residue pairs GTDB aligns, and whether the
table aligns them too.

Per pair: the residue pairs GTDB aligns (both residues on one of its columns), how many of those the table maps
(both residues' codons on family columns: coverage), how many of those it puts on the same family column
(agreement), and the share of GTDB's pairs whose amino acids differ (the pair's protein divergence). Pairs of one
genus are a species against its genus reference or two species through it; pairs of two genera go through the family
reference.

usage: validate_alignments.py <unpacked database dir> <msa tar.gz> <out prefix>
Writes <out>.pairs.tsv and <out>.summary.txt (also on stdout).
"""
import collections
import difflib
import itertools
import os
import sys
import tarfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "scripts"))
import ancestry_oracle as oracle  # noqa: E402

DB, MSA, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
AMINO = oracle.AMINO


def translate(nt):
    return "".join(oracle.amino_acid(oracle.CODE[nt[i]], oracle.CODE[nt[i + 1]], oracle.CODE[nt[i + 2]])
                   for i in range(0, len(nt) - 2, 3))


def plain(acc):
    return acc[3:] if acc[:3] in ("GB_", "RS_") else acc


taxid_of, genus_of, family_of = {}, {}, {}
with open(f"{DB}/genome2tiid.tsv") as f:
    for line in f:
        acc, taxid, rep, lineage = line.rstrip("\n").split("\t")[:4]
        if acc == rep:
            taxid_of[acc] = int(taxid)
            genus_of[int(taxid)] = lineage.split(";")[5]
            family_of[int(taxid)] = lineage.split(";")[4]
gene_of = {}
with open(f"{DB}/gene2geneid.tsv") as f:
    for line in f:
        marker, geneid = line.split()[:2]
        gene_of[marker] = int(geneid)
copies_nt = {}
with open(f"{DB}/reference.fna") as f:
    name = None
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            t, g = line[1:].split("_")[:2]
            name = (int(t), int(g))
        elif name:
            copies_nt[name] = copies_nt.get(name, "") + line
families, copies = oracle.read_column_weights(f"{DB}/column_weights.tsv")


def our_columns(taxid, gene):
    """Per residue of the copy's protein the family's codon column (first base on a column divisible by 3 and the
    codon's three bases on consecutive columns), else None; and the family row."""
    if (taxid, gene) not in copies:
        return None, None
    family, runs = copies[(taxid, gene)]
    nt_col = {}
    for begin, column, length in runs:
        for k in range(length):
            nt_col[begin + k] = column + k
    n = len(copies_nt[(taxid, gene)]) // 3
    out = [None] * n
    for r in range(n):
        c = nt_col.get(3 * r)
        if c is not None and c % 3 == 0 and nt_col.get(3 * r + 1) == c + 1 and nt_col.get(3 * r + 2) == c + 2:
            out[r] = c // 3
    return out, family


def gtdb_columns(row, protein):
    """Per GTDB column the protein's residue index (or None): the row's residues found in the protein through
    their matching blocks (the row is the protein less its insertions)."""
    residues = [(col, a) for col, a in enumerate(row) if a not in "-."]
    seq = "".join(a for _, a in residues)
    at = [None] * len(row)
    for blk in difflib.SequenceMatcher(None, seq, protein, autojunk=False).get_matching_blocks():
        for k in range(blk.size):
            at[residues[blk.a + k][0]] = blk.b + k
    return at


rows_out = []
with tarfile.open(MSA, "r|gz") as tar:
    for member in tar:
        if not member.isfile() or not member.name.endswith(".faa"):
            continue
        marker = os.path.basename(member.name).split("_reps_")[-1][:-4]
        gene = gene_of.get(marker)
        if gene is None:
            continue
        msa, name = {}, None
        for line in tar.extractfile(member).read().decode().splitlines():
            if line.startswith(">"):
                name = plain(line[1:].split()[0])
                msa[name] = []
            elif name in msa:
                msa[name].append(line.strip())
        per_copy = {}
        for acc, parts in msa.items():
            taxid = taxid_of.get(acc)
            if taxid is None or (taxid, gene) not in copies_nt:
                continue
            protein = translate(copies_nt[(taxid, gene)])
            ours, family = our_columns(taxid, gene)
            per_copy[taxid] = (gtdb_columns("".join(parts), protein), ours, family, protein)
        for a, b in itertools.combinations(sorted(per_copy), 2):
            if family_of[a] != family_of[b]:
                continue  # the families map onto columns of their own: no table aligns two families
            ga, oa, fa, pa = per_copy[a]
            gb, ob, fb, pb = per_copy[b]
            pairs = [(ra, rb) for ra, rb in zip(ga, gb) if ra is not None and rb is not None]
            if not pairs:
                continue
            differ = sum(pa[ra] != pb[rb] for ra, rb in pairs)
            if oa is None or ob is None or fa != fb:
                mapped = agree = 0
            else:
                both = [(ra, rb) for ra, rb in pairs if ra < len(oa) and rb < len(ob) and oa[ra] is not None and ob[rb] is not None]
                mapped = len(both)
                agree = sum(oa[ra] == ob[rb] for ra, rb in both)
            kind = "same genus" if genus_of[a] == genus_of[b] else "two genera"
            rows_out.append((marker, a, b, kind, len(pairs), mapped, agree, differ / len(pairs)))

with open(f"{OUT}.pairs.tsv", "w") as o:
    o.write("marker\ttaxid_a\ttaxid_b\tkind\tgtdb_pairs\tmapped\tagree\tprotein_divergence\n")
    for r in rows_out:
        o.write("\t".join(str(round(x, 4)) if isinstance(x, float) else str(x) for x in r) + "\n")

lines = []
say = lines.append
say(f"{len(rows_out)} copy pairs of {len({r[0] for r in rows_out})} markers; GTDB's aligned residue pairs: "
    f"coverage = the table maps both residues, agreement = on the same family column (of those mapped)")
for kind in ("same genus", "two genera", None):
    rs = [r for r in rows_out if kind is None or r[3] == kind]
    n, m, a = sum(r[4] for r in rs), sum(r[5] for r in rs), sum(r[6] for r in rs)
    say(f"  {kind or 'all':10s}: {len(rs)} pairs, {n} residue pairs, coverage {m / max(n, 1):.3f}, agreement {a / max(m, 1):.4f}")
say("  by the pair's protein divergence (GTDB's columns): pairs, coverage, agreement")
for lo, hi in ((0, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.6), (0.6, 1.01)):
    rs = [r for r in rows_out if lo <= r[7] < hi]
    n, m, a = sum(r[4] for r in rs), sum(r[5] for r in rs), sum(r[6] for r in rs)
    say(f"    [{lo:.2f}, {hi:.2f}): {len(rs):5d} pairs, coverage {m / max(n, 1):.3f}, agreement {a / max(m, 1):.4f}")
text = "\n".join(lines)
print(text)
with open(f"{OUT}.summary.txt", "w") as o:
    o.write(text + "\n")
