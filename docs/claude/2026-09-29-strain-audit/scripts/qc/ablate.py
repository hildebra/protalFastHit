#!/usr/bin/env python3
"""What each qcmsa default removes: run qcmsa on every species of a strains dir with the defaults
(plus --reapply-hcov 1000, as protal passes) and with one filter switched off at a time.

Per config and species: sample rows in/out, genes in/out, columns in/out, variable (ACGT) sites
in/out, and how much of the pairwise p-distance between the kept sequences survives."""
import sys, os, subprocess, itertools, re, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from consistency import read_fasta, parse_part

QCMSA = os.path.expanduser("~/audit5/src/scripts/qcmsa.py")
CONFIGS = [
    ("default", []),
    ("gene-min-samples 0", ["--gene-min-samples", "0"]),
    ("gene-min-samples 1", ["--gene-min-samples", "1"]),
    ("gene-min-hcov 0", ["--gene-min-hcov", "0"]),
    ("gene-min-mean-depth 0", ["--gene-min-mean-depth", "0"]),
    ("reapply-hcov 0", ["--reapply-hcov", "0"]),
    ("mrate2 off", ["--min-bad", "1000000", "--no-mask-cell-outliers"]),
    ("min-parsimony 0", ["--min-parsimony-samples", "0"]),
    ("all off", ["--gene-min-samples", "0", "--gene-min-hcov", "0", "--gene-min-mean-depth", "0",
                 "--reapply-hcov", "0", "--min-bad", "1000000", "--no-mask-cell-outliers",
                 "--min-parsimony-samples", "0"]),
]

def pd(a, b):
    n = d = 0
    for x, y in zip(a, b):
        if x in "ACGT" and y in "ACGT":
            n += 1
            d += x != y
    return d / n if n else None

def varsites(seqs):
    v = 0
    for col in zip(*seqs):
        s = {c for c in col if c in "ACGT"}
        v += len(s) > 1
    return v

def main(strains, out):
    os.makedirs(out, exist_ok=True)
    results = []
    for msa in sorted(glob.glob(os.path.join(strains, "*.raw.msa.fna"))):
        sp = os.path.basename(msa)[:-len(".raw.msa.fna")]
        rn, rs = read_fasta(msa)
        raw = dict(zip(rn, rs))
        rpart = parse_part(os.path.join(strains, sp + ".raw.partition.txt"))
        rvar = varsites(rs)
        rd = {}
        for a, b in itertools.combinations(rn, 2):
            rd[(a, b)] = pd(raw[a], raw[b])
        for cname, cargs in CONFIGS:
            pref = os.path.join(out, cname.replace(" ", "_"), sp)
            os.makedirs(os.path.dirname(pref), exist_ok=True)
            for ext in (".msa.fna", ".partition.txt", ".qcmsa_summary.tsv"):
                if os.path.exists(pref + ext):
                    os.remove(pref + ext)
            args = [sys.executable, QCMSA, msa, os.path.join(strains, sp + ".raw.partition.txt"),
                    os.path.join(strains, sp + ".meta.tsv"), "--prefix", pref, "--reapply-hcov", "1000"] + cargs
            p = subprocess.run(args, capture_output=True, text=True)
            r = dict(species=sp, config=cname, rc=p.returncode, samples_in=len(rn) - 1, genes_in=len(rpart),
                     cols_in=len(rs[0]), var_in=rvar)
            if os.path.exists(pref + ".msa.fna"):
                n, s = read_fasta(pref + ".msa.fna")
                part = parse_part(pref + ".partition.txt")
                r.update(samples_out=sum(1 for x in n if not x.endswith("_reference")), genes_out=len(part),
                         cols_out=len(s[0]), var_out=varsites(s))
                # partition consistency
                pos = 1; ok = True
                for g, a, b in part:
                    ok &= a == pos; pos = b + 1
                r["part_ok"] = ok and pos - 1 == len(s[0])
                ratios = []
                d = dict(zip(n, s))
                for a, b in itertools.combinations(n, 2):
                    key = (a, b) if (a, b) in rd else (b, a)
                    if rd.get(key):
                        q = pd(d[a], d[b])
                        if q is not None:
                            ratios.append(q / rd[key])
                r["dist_kept"] = sum(ratios) / len(ratios) if ratios else None
            else:
                r.update(samples_out=0, genes_out=0, cols_out=0, var_out=0, dist_kept=None, part_ok=None)
                r["msg"] = [l for l in p.stderr.splitlines() if "qcmsa.py:" in l][-1:]
            m = re.search(r"reapply-hcov=\d+: dropped (\d+)", p.stderr)
            r["reapply_dropped"] = int(m.group(1)) if m else 0
            results.append(r)
    json.dump(results, open(os.path.join(out, "ablate.json"), "w"), indent=1)
    for r in results:
        dk = f"{r['dist_kept']:.3f}" if r.get("dist_kept") is not None else "-"
        print(f"{r['species']:22} {r['config']:22} samples {r['samples_in']}->{r['samples_out']} genes {r['genes_in']}->{r['genes_out']} "
              f"cols {r['cols_in']}->{r['cols_out']} var {r['var_in']}->{r['var_out']} dist_kept {dk} "
              f"part_ok {r.get('part_ok')} reapply_drop {r['reapply_dropped']} {r.get('msg', '')}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
