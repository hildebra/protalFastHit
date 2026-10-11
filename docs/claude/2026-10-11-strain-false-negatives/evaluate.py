"""Read-level and call-level results of the strain world's runs (run_samples.sh), by how far each test strain is from
its species' representative and from the nearest genome of the allele half.

Divergence: per genome pair, the 21-mer containment of the strain's marker genes in the other genome's copies of the
same genes, pooled over the genes both have, as ANI = containment^(1/21) (cached in ~/strainfn/divergence.tsv).
Classes by ANI to the representative: easy >= 99.5%, medium 98.5-99.5%, hard < 98.5%.

Read level (deep samples, 60 pairs per strain): a pair's primary records; aligned (either mate), on the right species
(the mate with the higher MAPQ), counted (right species with MAPQ >= 4 and over 50 aligned bases, the profiler's
filter), counted on another species, and the identity of the right-species alignments (CIGAR M and = over M, =, X,
I, D). Calls (shallow samples): the strain's species predicted in the profile, by the strain's pairs; false positives
are predicted species absent from the sample.

usage: evaluate.py <strain dir> <out prefix> <run>...
"""
import collections
import glob
import os
import re
import subprocess
import sys

S, OUTP, RUNS = sys.argv[1], sys.argv[2], sys.argv[3:]
K = 21
CIG = re.compile(r"(\d+)([MIDNSHP=X])")

name_of, taxid_of = {}, {}
with open(f"{S}/converted/internal_taxonomy.dmp") as f:
    next(f)
    for line in f:
        r = line.rstrip("\n").split("\t")
        if r[4] == "species":
            name_of[r[0]] = r[3]
            taxid_of[r[3]] = r[0]

truth = [line.rstrip("\n").split("\t") for line in open(f"{S}/reads/truth.tsv")][1:]
taxa = [line.rstrip("\n").split("\t") for line in open(f"{S}/world/taxa.tsv")][1:]
halves = [line.rstrip("\n").split("\t") for line in open(f"{S}/world/halves.tsv")][1:]
rep_of = {r[2]: r[0] for r in taxa if r[1] == "representative"}
stratum_of = {r[2]: r[5] for r in taxa}
alleles_of = collections.defaultdict(list)
for acc, half, sp, st, rank in halves:
    if half == "allele":
        alleles_of[sp].append((int(rank), acc))

# ---- divergence ----
DIV = f"{S}/divergence.tsv"
strains = sorted({(r[2], r[3]) for r in truth if r[4] == "strain"})
div = {}
if os.path.exists(DIV):
    for line in open(DIV):
        acc, ani_rep, ani5, ani20 = line.split()
        div[acc] = (float(ani_rep), float(ani5), float(ani20))
if any(acc not in div for acc, _ in strains):
    wanted = {acc for acc, _ in strains} | {rep_of[sp] for _, sp in strains} | \
             {a for _, sp in strains for _, a in alleles_of[sp]}
    seqs = collections.defaultdict(dict)  # accession -> gene -> sequence
    for path in glob.glob(f"{S}/world/genomic_files_all/**/fna/*.fna", recursive=True):
        gene = os.path.basename(path)
        acc, buf = None, []
        for line in open(path):
            if line.startswith(">"):
                if acc in wanted:
                    seqs[acc][gene] = "".join(buf).upper()
                acc, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if acc in wanted:
            seqs[acc][gene] = "".join(buf).upper()
    kmers = {acc: {g: {s[i:i + K] for i in range(len(s) - K + 1)} for g, s in gs.items()} for acc, gs in seqs.items()}

    def ani(a, b):
        hit = tot = 0
        for g, ka in kmers.get(a, {}).items():
            kb = kmers.get(b, {}).get(g)
            if kb:
                hit += len(ka & kb)
                tot += len(ka)
        return (hit / tot) ** (1 / K) if tot and hit else 0.0

    with open(DIV, "w") as o:
        for acc, sp in strains:
            r = ani(acc, rep_of[sp])
            near = sorted(((rank, ani(acc, a)) for rank, a in alleles_of[sp]), key=lambda x: x[0])
            a5 = max([x for rank, x in near if rank <= 5] or [0.0])
            a20 = max([x for _, x in near] or [0.0])
            div[acc] = (r, a5, a20)
            o.write(f"{acc}\t{r:.5f}\t{a5:.5f}\t{a20:.5f}\n")


def klass(acc):
    r = div[acc][0]
    return "easy" if r >= 0.995 else "medium" if r >= 0.985 else "hard"


def near_allele(acc):
    r, _, a20 = div[acc]
    return "allele nearer" if a20 >= r + 0.003 else "rep as near"


# ---- read level ----
def sam_pairs(path):
    """qname -> {mate: (species taxid, mapq, aligned bases, identity)} of the primary records."""
    out = collections.defaultdict(dict)
    proc = subprocess.Popen(["zstdcat", path], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        if line[0] == "@":
            continue
        f = line.split("\t", 6)
        flag = int(f[1])
        if flag & 0x904:
            continue
        ops = collections.Counter()
        for n, op in CIG.findall(f[5]):
            ops[op] += int(n)
        match = ops["M"] + ops["="]
        aligned = match + ops["X"] + ops["I"]
        ident = match / max(1, match + ops["X"] + ops["I"] + ops["D"])
        out[f[0]][2 if flag & 0x80 else 1] = (f[2].split("_")[0], int(f[4]), aligned, ident)
    proc.wait()
    return out


def read_level(run, kind="deep"):
    """(class, sample, accession) -> Counter of the strain's pairs; class -> identities of right-species pairs."""
    stats = collections.defaultdict(lambda: collections.Counter())
    ident = collections.defaultdict(list)
    for sample in sorted({r[0] for r in truth if r[1] == kind}):
        members = {r[2]: r for r in truth if r[0] == sample}
        for qname, mates in sam_pairs(f"{S}/run_{run}/alignments/{sample}.sam.zst").items():
            acc = qname.rsplit("__", 1)[0]
            m = members.get(acc)
            if not m or m[4] != "strain":
                continue
            best = max(mates.values(), key=lambda x: (x[1], x[2]))
            key = (klass(acc), sample, acc)
            stats[key]["aligned"] += 1
            if best[0] == taxid_of.get(m[3]):
                stats[key]["right"] += 1
                ident[klass(acc)].append(best[3])
                if best[1] >= 4 and best[2] > 50:
                    stats[key]["counted"] += 1
                else:
                    stats[key]["right_low_mapq"] += 1
            elif best[1] >= 4 and best[2] > 50:
                stats[key]["counted_wrong"] += 1
        for acc, m in members.items():
            if m[4] == "strain":
                stats[(klass(acc), sample, acc)]["pairs"] += int(m[5])
    return stats, ident


# ---- calls ----
def profile_rows(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(head)}
        for line in f:
            r = line.rstrip("\n").split("\t")
            if len(r) > ci["Name"]:
                yield r[ci["Name"]], r[ci["Predicted"]] == "1", float(r[ci["Probability"]] or 0)


def calls(run, kind, knob=None):
    """[(class, near, pairs, called, p, stratum, sample, accession)], false positives, background calls; with knob,
    called = p >= knob instead of the profile's prediction."""
    res = []
    fp = tp_bg = 0
    for sample in sorted({r[0] for r in truth if r[1] == kind}):
        members = [r for r in truth if r[0] == sample]
        present = {r[3] for r in members}
        rows = {name: (pred if knob is None else p >= knob, p)
                for name, pred, p in profile_rows(f"{S}/run_{run}/profiles/{sample}.profile.log")}
        for r in members:
            pred, p = rows.get(r[3], (False, 0.0))
            if r[4] == "strain":
                res.append((klass(r[2]), near_allele(r[2]), int(r[5]), pred, p, stratum_of[r[3]], sample, r[2]))
            else:
                tp_bg += pred
        fp += sum(1 for name, (pred, _) in rows.items() if pred and name not in present)
    return res, fp, tp_bg


classes = ("easy", "medium", "hard")
n_class = collections.Counter(klass(acc) for acc, _ in strains)
print("test strains by class (ANI to rep): " + ", ".join(f"{c} {n_class[c]}" for c in classes))
print("strains whose nearest allele-half genome is >= 0.3% nearer than the rep: "
      f"{sum(near_allele(a) == 'allele nearer' for a, _ in strains)} of {len(strains)}")
for c in classes:
    xs = [div[a] for a, _ in strains if klass(a) == c]
    if xs:
        print(f"  {c}: median ANI to rep {sorted(x[0] for x in xs)[len(xs) // 2]:.4f}, to nearest of 5 allele genomes "
              f"{sorted(x[1] for x in xs)[len(xs) // 2]:.4f}, of 20 {sorted(x[2] for x in xs)[len(xs) // 2]:.4f}")

with open(f"{OUTP}_reads.tsv", "w") as o:
    o.write("run\tclass\tstrains\tpairs\taligned\tright\tcounted\tright_low_mapq\tcounted_wrong\tidentity\n")
    print("\n## read level (deep samples): share of the strains' pairs")
    print("| run | class | pairs | aligned | right species | counted (MAPQ>=4) | right, MAPQ<4 | counted on another | identity |")
    print("|---|---|---|---|---|---|---|---|---|")
    for run in RUNS:
        stats, ident = read_level(run)
        for c in classes:
            t = collections.Counter()
            n = 0
            for (k, _, _), v in stats.items():
                if k == c:
                    t.update(v)
                    n += 1
            if not t["pairs"]:
                continue
            P = t["pairs"]
            mi = sum(ident[c]) / max(1, len(ident[c]))
            o.write(f"{run}\t{c}\t{n}\t{P}\t{t['aligned']}\t{t['right']}\t{t['counted']}\t{t['right_low_mapq']}\t"
                    f"{t['counted_wrong']}\t{mi:.4f}\n")
            print(f"| {run} | {c} | {P} | {t['aligned'] / P:.3f} | {t['right'] / P:.3f} | {t['counted'] / P:.3f} | "
                  f"{t['right_low_mapq'] / P:.3f} | {t['counted_wrong'] / P:.4f} | {mi:.4f} |")

bins = ((1, 2), (3, 6), (8, 16), (24, 32))
with open(f"{OUTP}_calls.tsv", "w") as o:
    o.write("run\tclass\tpairs\tn\tcalled\tmean_p\n")
    print("\n## calls (shallow samples): share of test strains whose species is called, by pairs of the strain")
    print("| run | class | " + " | ".join(f"{a}-{b} pairs" for a, b in bins) + " | all | FP (20 samples) | deep: called | deep FP |")
    print("|---|---|" + "---|" * (len(bins) + 4))
    for run in RUNS:
        res, fp, _ = calls(run, "shallow")
        dres, dfp, _ = calls(run, "deep")
        for c in classes + ("all",):
            cells = []
            for a, b in bins:
                xs = [x for x in res if (c == "all" or x[0] == c) and a <= x[2] <= b]
                cells.append(f"{sum(x[3] for x in xs) / len(xs):.3f} ({len(xs)})" if xs else "-")
                if xs:
                    o.write(f"{run}\t{c}\t{a}-{b}\t{len(xs)}\t{sum(x[3] for x in xs)}\t{sum(x[4] for x in xs) / len(xs):.4f}\n")
            xs = [x for x in res if c == "all" or x[0] == c]
            dx = [x for x in dres if c == "all" or x[0] == c]
            print(f"| {run} | {c} | " + " | ".join(cells) + f" | {sum(x[3] for x in xs) / max(1, len(xs)):.3f} | "
                  f"{fp if c == 'all' else ''} | {sum(x[3] for x in dx) / max(1, len(dx)):.3f} | {dfp if c == 'all' else ''} |")
    print("\n## knobs (shallow samples): called share of test strains (all classes, hard) and false positives")
    print("| run | " + " | ".join(f"knob {k}: all / hard / FP" for k in (0.5, 0.35, 0.2)) + " |")
    print("|---|---|---|---|")
    for run in RUNS:
        cells = []
        for k in (0.5, 0.35, 0.2):
            res, fp, _ = calls(run, "shallow", k)
            hard = [x for x in res if x[0] == "hard"]
            cells.append(f"{sum(x[3] for x in res) / len(res):.3f} / {sum(x[3] for x in hard) / max(1, len(hard)):.3f} / {fp}")
        print(f"| {run} | " + " | ".join(cells) + " |")
    print("\n## why a test strain is missed (shallow samples): its pairs' fate")
    print("| run | missed | no pair counted: none aligned | aligned elsewhere only | right species, MAPQ < 4 only | "
          "1 pair counted | 2+ pairs counted (model said no) |")
    print("|---|---|---|---|---|---|---|")
    for run in RUNS:
        res, _, _ = calls(run, "shallow")
        stats, _ = read_level(run, "shallow")
        why = collections.Counter()
        for c, _, pairs, pred, p, _, sample, acc in res:
            if pred:
                continue
            t = stats.get((c, sample, acc), collections.Counter())
            why["missed"] += 1
            if t["counted"] >= 2:
                why["model"] += 1
            elif t["counted"] == 1:
                why["one"] += 1
            elif t["right"]:
                why["mapq"] += 1
            elif t["aligned"]:
                why["elsewhere"] += 1
            else:
                why["none"] += 1
        print(f"| {run} | {why['missed']} | {why['none']} | {why['elsewhere']} | {why['mapq']} | {why['one']} | "
              f"{why['model']} |")
    print("\n## evidence pairs per species with and without the MAPQ filter (shallow samples, best record of each pair)")
    print("| run | missed strains | with 2+ right pairs at MAPQ>=4 | with 2+ right pairs at any MAPQ | absent species with "
          "2+ pairs at MAPQ>=4 | absent species with 2+ pairs at any MAPQ |")
    print("|---|---|---|---|---|---|")
    for run in RUNS:
        res, _, _ = calls(run, "shallow")
        missed = {(x[6], x[7]) for x in res if not x[3]}
        m4 = many = a4 = aany = 0
        for sample in sorted({r[0] for r in truth if r[1] == "shallow"}):
            members = {r[2]: r for r in truth if r[0] == sample}
            present = {taxid_of.get(r[3]) for r in members.values()}
            per = collections.defaultdict(collections.Counter)  # (sample-level) taxid -> counts; strain -> counts
            own = collections.defaultdict(collections.Counter)
            for qname, mates in sam_pairs(f"{S}/run_{run}/alignments/{sample}.sam.zst").items():
                acc = qname.rsplit("__", 1)[0]
                best = max(mates.values(), key=lambda x: (x[1], x[2]))
                ok = best[1] >= 4 and best[2] > 50
                per[best[0]]["any"] += 1
                per[best[0]]["q4"] += ok
                if acc in members and best[0] == taxid_of.get(members[acc][3]):
                    own[acc]["any"] += 1
                    own[acc]["q4"] += ok
            for s_, acc in missed:
                if s_ == sample:
                    m4 += own[acc]["q4"] >= 2
                    many += own[acc]["any"] >= 2
            for t, c in per.items():
                if t not in present:
                    a4 += c["q4"] >= 2
                    aany += c["any"] >= 2
        print(f"| {run} | {len(missed)} | {m4} | {many} | {a4} | {aany} |")
    print("\n## calls by whether an allele-half genome is nearer than the rep (shallow, all depths)")
    print("| run | " + " | ".join(f"{c}: {n}" for c in classes for n in ("allele nearer", "rep as near")) + " |")
    print("|---|" + "---|" * 6)
    for run in RUNS:
        res, _, _ = calls(run, "shallow")
        cells = []
        for c in classes:
            for n in ("allele nearer", "rep as near"):
                xs = [x for x in res if x[0] == c and x[1] == n]
                cells.append(f"{sum(x[3] for x in xs) / len(xs):.3f} ({len(xs)})" if xs else "-")
        print(f"| {run} | " + " | ".join(cells) + " |")
