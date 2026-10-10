"""Which reference would keep the column weights' alignments short? On GTDB's own alignments of the representatives'
proteins (<set>_msa_marker_genes_reps_r226.tar.gz) and GTDB's tree (<set>_r226.tree), per chosen genus and marker:
the protein divergence (p-distance over the columns both carry) of each species to
  - a random member (what --build does now: a member chosen by hash),
  - the medoid (the member with the least summed divergence to the others),
  - the genus's ancestral sequence, by Fitch parsimony on the tree pruned to the genus (per column; a tie between
    states resolved by the states' counts among the members, then alphabetically; a gap is a state),
and per family and marker the same for the genus references (the medoids) against a random genus reference, their
medoid and the family's ancestral sequence. The share beyond the amino-acid limits (0.3 within, 0.5 among) is what
the alignment loses.

usage: reference_choice.py <taxa.tsv> <msa tar.gz> <tree> <out prefix>
"""
import collections
import os
import random
import statistics
import sys
import tarfile

TAXA, MSA, TREE, OUT = sys.argv[1:5]
WITHIN, AMONG = 0.3, 0.5
rng = random.Random(1)

rep_of, genus_of, family_of = {}, {}, {}
with open(TAXA) as f:
    next(f)
    for line in f:
        acc, role, sp, gen, fam = line.rstrip("\n").split("\t")
        acc = acc[3:] if acc[:3] in ("GB_", "RS_") else acc
        if role == "representative":
            rep_of[sp] = acc
            genus_of[acc], family_of[acc] = gen, fam
members = collections.defaultdict(list)
for acc, gen in genus_of.items():
    members[gen].append(acc)
genera_of = collections.defaultdict(list)
for gen in members:
    genera_of[family_of[members[gen][0]]].append(gen)


def parse_newick(text):
    """The tree as nested lists: a leaf is its name, an inner node a list of children."""
    pos = 0

    def node():
        nonlocal pos
        if text[pos] == "(":
            pos += 1
            kids = [node()]
            while text[pos] == ",":
                pos += 1
                kids.append(node())
            pos += 1  # ")"
            label(); length()
            return kids
        name = label()
        length()
        return name

    def label():
        nonlocal pos
        if pos < len(text) and text[pos] == "'":
            end = text.index("'", pos + 1)
            s = text[pos + 1:end]
            pos = end + 1
            return s
        start = pos
        while pos < len(text) and text[pos] not in ",():;":
            pos += 1
        return text[start:pos]

    def length():
        nonlocal pos
        if pos < len(text) and text[pos] == ":":
            pos += 1
            while pos < len(text) and text[pos] not in ",();":
                pos += 1

    return node()


def prune(tree, keep):
    """The tree on the leaves in `keep`, unary nodes removed (None if none is left)."""
    if isinstance(tree, str):
        name = tree[3:] if tree[:3] in ("GB_", "RS_") else tree
        return name if name in keep else None
    kids = [k for k in (prune(c, keep) for c in tree) if k is not None]
    if not kids:
        return None
    return kids[0] if len(kids) == 1 else kids


def fitch(tree, column):
    """The Fitch state set at the root (column: leaf -> state)."""
    if isinstance(tree, str):
        return {column[tree]}
    sets = [fitch(c, column) for c in tree]
    common = set.intersection(*sets)
    return common if common else set.union(*sets)


def ancestor(tree, rows, leaves):
    out = []
    for i in range(len(next(iter(rows.values())))):
        column = {l: rows[l][i] for l in leaves}
        states = fitch(tree, column)
        counts = collections.Counter(column.values())
        out.append(min(states, key=lambda s: (-counts[s], s)))
    return "".join(out)


def distance(a, b):
    both = [(x, y) for x, y in zip(a, b) if x not in "-." and y not in "-."]
    return sum(x != y for x, y in both) / len(both) if len(both) >= 30 else None


def medoid(names, rows):
    best, best_sum = None, None
    for a in names:
        ds = [distance(rows[a], rows[b]) for b in names if b != a]
        s = sum(d for d in ds if d is not None) if ds else 0
        if best_sum is None or s < best_sum:
            best, best_sum = a, s
    return best


with open(TREE) as f:
    tree = parse_newick(f.read().strip())
within = collections.defaultdict(list)   # choice -> divergences of species to their genus reference
among = collections.defaultdict(list)    # choice -> divergences of genus references to the family reference
with tarfile.open(MSA, "r|gz") as tar:
    for member in tar:
        if not member.isfile() or not member.name.endswith(".faa"):
            continue
        rows, name = {}, None
        for line in tar.extractfile(member).read().decode().splitlines():
            if line.startswith(">"):
                acc = line[1:].split()[0]
                name = acc[3:] if acc[:3] in ("GB_", "RS_") else acc
                name = name if name in genus_of else None
                if name:
                    rows[name] = []
            elif name:
                rows[name].append(line.strip())
        rows = {k: "".join(v) for k, v in rows.items()}
        if len(rows) < 3:
            continue
        for fam, gens in genera_of.items():
            gen_refs = {}
            for gen in gens:
                names = [a for a in members[gen] if a in rows]
                if len(names) < 2:
                    if names:
                        gen_refs[gen] = names[0]
                    continue
                sub = prune(tree, set(names))
                anc = ancestor(sub, rows, names)
                med = medoid(names, rows)
                rnd = rng.choice(names)
                gen_refs[gen] = med
                for a in names:
                    for choice, ref in (("random member", rows[rnd]), ("medoid", rows[med]), ("ancestor", anc)):
                        if choice != "ancestor" and a == (rnd if choice == "random member" else med):
                            continue
                        d = distance(rows[a], ref)
                        if d is not None:
                            within[choice].append(d)
            refs = list(gen_refs.values())
            if len(refs) < 2:
                continue
            sub = prune(tree, set(refs))
            anc = ancestor(sub, rows, refs)
            med = medoid(refs, rows)
            rnd = rng.choice(refs)
            for a in refs:
                for choice, ref in (("random member", rows[rnd]), ("medoid", rows[med]), ("ancestor", anc)):
                    if choice != "ancestor" and a == (rnd if choice == "random member" else med):
                        continue
                    d = distance(rows[a], ref)
                    if d is not None:
                        among[choice].append(d)

lines = []
for title, data, limit in (("species against their genus reference", within, WITHIN),
                           ("genus references (medoids) against the family reference", among, AMONG)):
    lines.append(f"{title} (protein p-distance on GTDB's columns; share beyond {limit}):")
    for choice in ("random member", "medoid", "ancestor"):
        xs = sorted(data[choice])
        if not xs:
            continue
        q = lambda p: xs[min(len(xs) - 1, int(p * len(xs)))]
        lines.append(f"  {choice:14s} n={len(xs):6d}  quartiles {q(0.25):.3f} {q(0.5):.3f} {q(0.75):.3f}  "
                     f"90% {q(0.9):.3f}  beyond {sum(x > limit for x in xs) / len(xs):.3f}")
text = "\n".join(lines)
print(text)
with open(f"{OUT}.txt", "w") as o:
    o.write(text + "\n")
