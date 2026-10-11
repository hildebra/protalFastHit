"""Per focal species of the strain world: its stratum, the test strains' median ANI to the representative (evaluate.py's
divergence.tsv), the share of its strains' deep-sample pairs on the right species but under MAPQ 4 (first run given:
the congeners' pull), and the share of its strains called in the shallow samples, per run.

usage: species_table.py <strain dir> <run>...
"""
import collections
import subprocess
import sys

S, RUNS = sys.argv[1], sys.argv[2:]
taxid_of = {}
with open(f"{S}/converted/internal_taxonomy.dmp") as f:
    next(f)
    for line in f:
        r = line.rstrip("\n").split("\t")
        if r[4] == "species":
            taxid_of[r[3]] = r[0]
truth = [line.rstrip("\n").split("\t") for line in open(f"{S}/reads/truth.tsv")][1:]
stratum = {r[2]: r[5] for r in (line.rstrip("\n").split("\t") for line in open(f"{S}/world/taxa.tsv")) if len(r) > 5}
div = {line.split()[0]: float(line.split()[1]) for line in open(f"{S}/divergence.tsv")}


def low_mapq_share(run):
    """species -> (right-species pairs under MAPQ 4, right-species pairs) of its strains in the deep samples."""
    out = collections.defaultdict(lambda: [0, 0])
    for sample in sorted({r[0] for r in truth if r[1] == "deep"}):
        members = {r[2]: r[3] for r in truth if r[0] == sample and r[4] == "strain"}
        best = {}
        proc = subprocess.Popen(["zstdcat", f"{S}/run_{run}/alignments/{sample}.sam.zst"], stdout=subprocess.PIPE, text=True)
        for line in proc.stdout:
            if line[0] == "@":
                continue
            f = line.split("\t", 6)
            if int(f[1]) & 0x904:
                continue
            q = (int(f[4]), f[2].split("_")[0])
            if f[0] not in best or q > best[f[0]]:
                best[f[0]] = q
        proc.wait()
        for qname, (mapq, taxid) in best.items():
            sp = members.get(qname.rsplit("__", 1)[0])
            if sp and taxid == taxid_of[sp]:
                out[sp][1] += 1
                out[sp][0] += mapq < 4
    return out


def called(run):
    out = collections.defaultdict(lambda: [0, 0])
    for sample in sorted({r[0] for r in truth if r[1] == "shallow"}):
        with open(f"{S}/run_{run}/profiles/{sample}.profile.log") as f:
            head = f.readline().rstrip("\n").split("\t")
            ci = {c: i for i, c in enumerate(head)}
            pred = {r[ci["Name"]]: r[ci["Predicted"]] == "1" for r in (l.rstrip("\n").split("\t") for l in f)
                    if len(r) > ci["Name"]}
        for r in truth:
            if r[0] == sample and r[4] == "strain":
                out[r[3]][0] += pred.get(r[3], False)
                out[r[3]][1] += 1
    return out


low = low_mapq_share(RUNS[0])
calls = {run: called(run) for run in RUNS}
species = sorted({r[3] for r in truth if r[4] == "strain"}, key=lambda s: ("tight mid wide".split().index(stratum[s]), s))
print("| species | stratum | median ANI to rep | right pairs under MAPQ 4 (" + RUNS[0] + ") | " +
      " | ".join(f"called: {r}" for r in RUNS) + " |")
print("|---|---|---|---|" + "---|" * len(RUNS))
for sp in species:
    anis = sorted(div[r[2]] for r in truth if r[3] == sp and r[4] == "strain" and r[1] == "shallow")
    lo, n = low[sp]
    print(f"| {sp[3:]} | {stratum[sp]} | {anis[len(anis) // 2]:.4f} | {lo / max(1, n):.3f} | " +
          " | ".join(f"{calls[r][sp][0]}/{calls[r][sp][1]}" for r in RUNS) + " |")
