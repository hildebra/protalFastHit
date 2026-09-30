#!/usr/bin/env python3
"""Build MSA variants for a run, infer IQ-TREE trees and score them against the truth.

usage: analyze.py <run_dir> [--species Malpha,Tone,Cferv] [--variants raw,filt,...] [--mix mix=t03]
       [--protal_subdir protal] [--tag TAG]

Variants (qcmsa = the traced copy of protal's qcmsa.py, byte-identical output):
  raw            protal .raw.msa.fna, no qcmsa                              GTR+G
  filt           qcmsa defaults (as protal runs it: --reapply-hcov 1000)   GTR+G   <- justfile default
  filt_keepsing  filt + --min-parsimony-samples 0                           GTR+G
  filt_nomr2     filt, MRate2 gene/sample/cell filter off                   GTR+G
  filt_nocov     filt, coverage gates off                                   GTR+G
  filt_none      all qcmsa filters off (only --reapply-hcov)                GTR+G
  filt_varonly   filt + --discard-constant                                  GTR+G (no ASC)
  filt_varasc    filt + --discard-constant                                  GTR+G+ASC
  filt_part      filt, partitioned (-p qcmsa partition, GTR+G per gene)
  raw_part       raw, partitioned (-p raw partition)
  filt_snpaf     (needs a protal run with other SNP settings: use --protal_subdir)
Results: <run_dir>/an[_TAG]/<sp>/<variant>/ and one JSON line per tree in <run_dir>/results[_TAG].jsonl
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import treecmp  # noqa: E402

PHYLO = os.path.dirname(HERE)
SPECIES = {"Malpha": "s__Mockella_alpha", "Cferv": "s__Calidella_fervens", "Tone": "s__Testella_one"}
QCMSA = os.path.join(HERE, "qcmsa_trace.py")
MR2_OFF = ["--min-bad", "1000000", "--no-mask-cell-outliers"]  # --iqr-mult cannot disable it: fence = Q3 when IQR = 0
COV_OFF = ["--gene-min-hcov", "0", "--gene-min-mean-depth", "0", "--gene-min-samples", "0"]
VARIANTS = {
    "filt": ([], "GTR+G", False),
    "filt_keepsing": (["--min-parsimony-samples", "0"], "GTR+G", False),
    "filt_nomr2": (MR2_OFF, "GTR+G", False),
    "filt_nocov": (COV_OFF, "GTR+G", False),
    "filt_none": (["--min-parsimony-samples", "0"] + MR2_OFF + COV_OFF, "GTR+G", False),
    "filt_varonly": (["--discard-constant"], "GTR+G", False),
    "filt_varasc": (["--discard-constant"], "GTR+G+ASC", False),
    "filt_part": ([], "GTR+G", True),
    "filt_strict": (["--preset", "strict"], "GTR+G", False),
    "filt_sabs2": (["--sample-abs-min-bad", "2"], "GTR+G", False),
    "filt_zeros": (["--mrate2-include-zeros"], "GTR+G", False),
}
DEFAULT_VARIANTS = ["raw", "filt", "filt_keepsing", "filt_nomr2", "filt_nocov", "filt_none",
                    "filt_varonly", "filt_varasc", "filt_part", "raw_part"]


def nseq(path):
    return sum(1 for line in open(path) if line.startswith(">"))


def run_iqtree(msa, prefix, model, partition=None, threads=2):
    if os.path.exists(prefix + ".treefile"):
        return True
    cmd = ["iqtree2", "-s", msa, "-m", model, "-B", "1000", "-T", str(threads), "--seqtype", "DNA",
           "--prefix", prefix, "-redo", "-seed", "12345"]
    if partition:
        cmd[3:3] = []
        cmd += ["-p", partition]
    with open(prefix + ".stdout.log", "w") as fh:
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if rc != 0 and model.endswith("+ASC") and os.path.exists(prefix + ".varsites.phy"):
        # IQ-TREE refuses +ASC when columns are invariant once ambiguity is resolved; it writes the
        # truly variable sites, which is what one reruns on.
        with open(prefix + ".asc_retry.txt", "w") as fh:
            fh.write("first +ASC run failed; rerun on .varsites.phy\n")
        vs = prefix + ".varsites.phy"
        shutil.copy(vs, prefix + ".input_varsites.phy")
        cmd = ["iqtree2", "-s", prefix + ".input_varsites.phy", "-m", model, "-B", "1000", "-T", str(threads),
               "--seqtype", "DNA", "--prefix", prefix, "-redo", "-seed", "12345"]
        with open(prefix + ".stdout.log", "a") as fh:
            rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
    return rc == 0 and os.path.exists(prefix + ".treefile")


def iqtree_log_info(prefix):
    info = {"warnings": [], "notes_identical": 0}
    p = prefix + ".log"
    if not os.path.exists(p):
        return info
    txt = open(p, errors="replace").read()
    m = re.findall(r"Alignment has (\d+) sequences with (\d+) columns, (\d+) distinct patterns", txt)
    if m:
        info["iq_seqs"], info["iq_cols"], info["iq_patterns"] = map(int, m[-1])
    m = re.findall(r"(\d+) parsimony-informative, (\d+) singleton sites, (\d+) constant sites", txt)
    if m:
        info["iq_pinf"], info["iq_single"], info["iq_const"] = map(int, m[-1])
    info["warnings"] = sorted(set(re.findall(r"WARNING: ([^\n]+)", txt)))[:20]
    info["notes_identical"] = len(re.findall(r"is identical to", txt))
    info["asc_retry"] = os.path.exists(prefix + ".asc_retry.txt")
    return info


def truth_scores(raw, part, truth, msa=None, colmap=None, mix=None):
    cmd = [sys.executable, os.path.join(HERE, "msa_truth.py"), raw, part, truth, "--json"]
    if msa:
        cmd += ["--msa", msa, "--colmap", colmap]
    if mix:
        cmd += ["--mix", mix]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        return {"error": out.stderr[-500:]}
    res = json.loads(out.stdout)
    agg = {}
    for row in res["rows"].values():
        for k, v in row.items():
            agg[k] = agg.get(k, 0) + v
    return {"cells": agg, "columns": res["columns"], "dropped_columns": res.get("dropped_columns"),
            "rows_detail": res["rows"]}


def qcmsa_summary(prefix):
    p = prefix + ".qcmsa_summary.tsv"
    out = {}
    if not os.path.exists(p):
        return out
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if f[0] == "count":
            out[f[1]] = int(f[2])
        elif f[0] == "sample_filtered":
            out.setdefault("samples_removed", []).append(f[1])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--species", default="Malpha,Tone,Cferv")
    ap.add_argument("--variants", default=",".join(DEFAULT_VARIANTS))
    ap.add_argument("--mix", default=None)
    ap.add_argument("--protal_subdir", default="protal")
    ap.add_argument("--tag", default="")
    ap.add_argument("--reapply_hcov", default="1000")
    a = ap.parse_args()
    run = os.path.abspath(a.run_dir)
    tag = f"_{a.tag}" if a.tag else ""
    an = os.path.join(run, "an" + tag)
    res_path = os.path.join(run, f"results{tag}.jsonl")
    done = set()
    if os.path.exists(res_path):
        for line in open(res_path):
            r = json.loads(line)
            done.add((r["species"], r["variant"]))
    resf = open(res_path, "a")
    for sp in a.species.split(","):
        name = SPECIES[sp]
        sdir = os.path.join(run, a.protal_subdir, "strains")
        raw = os.path.join(sdir, name + ".raw.msa.fna")
        rpart = os.path.join(sdir, name + ".raw.partition.txt")
        meta = os.path.join(sdir, name + ".meta.tsv")
        truth = os.path.join(PHYLO, "strains", sp, "truth_markers.tsv")
        true_tree = os.path.join(PHYLO, "strains", sp, "true.nwk")
        ref = name + "_reference"
        if not os.path.exists(raw):
            resf.write(json.dumps({"run": os.path.basename(run), "tag": a.tag, "species": sp, "variant": "raw",
                                   "status": "no raw MSA"}) + "\n")
            resf.flush()
            continue
        for v in a.variants.split(","):
            if (sp, v) in done:
                continue
            d = os.path.join(an, sp, v)
            os.makedirs(d, exist_ok=True)
            rec = {"run": os.path.basename(run), "tag": a.tag, "species": sp, "variant": v}
            partition = None
            colmap = None
            if v in ("raw", "raw_part"):
                msa = raw
                model = "GTR+G"
                if v == "raw_part":
                    partition = rpart
                rec["n_seqs"] = nseq(raw)
            elif v in ("raw_noref", "filt_noref"):
                # the same MSA without the database reference row
                src = raw if v == "raw_noref" else os.path.join(an, sp, "filt", "q.msa.fna")
                if not os.path.exists(src):
                    rec["status"] = "no source MSA"
                    resf.write(json.dumps(rec) + "\n")
                    resf.flush()
                    continue
                msa = os.path.join(d, "noref.msa.fna")
                with open(src) as fi, open(msa, "w") as fo:
                    keep = True
                    for line in fi:
                        if line.startswith(">"):
                            keep = not line[1:].strip().endswith("_reference")
                        if keep:
                            fo.write(line)
                model = "GTR+G"
                rec["n_seqs"] = nseq(msa)
            else:
                qargs, model, part_flag = VARIANTS[v]
                prefix = os.path.join(d, "q")
                if not os.path.exists(prefix + ".colmap.tsv"):
                    cmd = [sys.executable, QCMSA, raw, rpart, meta, "--prefix", prefix,
                           "--reapply-hcov", a.reapply_hcov] + qargs
                    with open(prefix + ".qcmsa.log", "w") as fh:
                        subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
                msa = prefix + ".msa.fna"
                colmap = prefix + ".colmap.tsv"
                rec["qcmsa"] = qcmsa_summary(prefix)
                if not os.path.exists(msa):
                    rec["status"] = "qcmsa produced no MSA"
                    rec["qcmsa_log_tail"] = open(prefix + ".qcmsa.log").read()[-300:]
                    resf.write(json.dumps(rec) + "\n")
                    resf.flush()
                    continue
                if part_flag:
                    partition = prefix + ".partition.txt"
                rec["n_seqs"] = nseq(msa)
                # default filtered MSA must equal protal's own
                if v == "filt" and a.protal_subdir == "protal":
                    own = os.path.join(sdir, name + ".msa.fna")
                    rec["same_as_protal_filtered"] = (os.path.exists(own) and
                                                      open(own).read() == open(msa).read())
            rec["model"] = model + (" partitioned" if partition else "")
            if rec["n_seqs"] < 4:
                rec["status"] = f"only {rec['n_seqs']} sequences (<4): justfile skips"
                resf.write(json.dumps(rec) + "\n")
                resf.flush()
                continue
            iqp = os.path.join(d, "iq")
            ok = run_iqtree(msa, iqp, model, partition)
            rec.update(iqtree_log_info(iqp))
            if not ok:
                rec["status"] = "iqtree failed"
                rec["iq_tail"] = open(iqp + ".stdout.log").read()[-600:]
            else:
                rec["status"] = "ok"
                rec.update(treecmp.compare(iqp + ".treefile", true_tree, ref, a.mix))
            if v in ("raw", "filt", "filt_keepsing", "filt_none", "filt_nomr2", "filt_nocov"):
                ts = truth_scores(raw, rpart, truth, None if v == "raw" else msa, colmap, a.mix)
                rec["truth"] = ts
            resf.write(json.dumps(rec) + "\n")
            resf.flush()
            print(f"[analyze] {os.path.basename(run)}{tag} {sp} {v}: {rec.get('status')} "
                  f"rf={rec.get('rf')} nrf={rec.get('nrf')}", flush=True)
    resf.close()


if __name__ == "__main__":
    main()
