#!/usr/bin/env python3
"""Bytes on disk per read pair or base of each read type the collector writes, as it writes them (simulate_metagenomes:
BGZF level 6; hifi_reads.py and the pbsim3 renaming: gzip level 1), and the same reads recompressed (gzip -6, zstd -3,
zstd -19) for comparison; plus pbsim3's temporary files per Mb.

usage: disk_costs.py --scripts PROTAL/scripts --genomes genomes.tsv --simulator simulate_metagenomes --out DIR
"""

import argparse
import glob
import gzip
import json
import os
import random
import shutil
import subprocess
import sys


def sizes(path, work):
    """{as written, gzip -6, zstd -3, zstd -19} bytes of a gzipped FASTQ."""
    plain = os.path.join(work, "plain.fq")
    with gzip.open(path, "rb") as fin, open(plain, "wb") as fout:
        shutil.copyfileobj(fin, fout, 16 << 20)
    out = {"written": os.path.getsize(path), "plain": os.path.getsize(plain)}
    out["gzip6"] = len(subprocess.run(["gzip", "-6", "-c", plain], capture_output=True, check=True).stdout)
    for level in (3, 19):
        out[f"zstd{level}"] = len(subprocess.run(["zstd", f"-{level}", "-T4", "-c", "-q", plain], capture_output=True,
                                                 check=True).stdout)
    os.remove(plain)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scripts", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--simulator", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pbsim", default="pbsim")
    ap.add_argument("--pbsim_models")
    ap.add_argument("--pairs", type=int, default=200_000)
    ap.add_argument("--bases", type=int, default=30_000_000)
    a = ap.parse_args()
    sys.path.insert(0, a.scripts)
    import collect_training_data as ctd
    import hifi_reads
    shutil.rmtree(a.out, ignore_errors=True)
    os.makedirs(a.out)
    results = {}

    # Paired-end reads: the design's three setups, and the scenarios' (150 bp HSXt at Q35: -qs -5 -qs2 -3).
    setups = [("pe100_HS20", "100", "HS20", "300", "40", ""), ("pe150_HSXt", "150", "HSXt", "350", "50", ""),
              ("pe250_MSv3", "250", "MSv3", "550", "50", ""), ("pe150_HSXt_Q35", "150", "HSXt", "350", "50", "-qs -5 -qs2 -3")]
    for name, length, profile, fmean, fsd, extra in setups:
        sim = os.path.join(a.out, name)
        command = [a.simulator, "--genome_table", a.genomes, "-o", sim, "-n", "1", "--sample_prefix", name + "_s",
                   "--total_read_pairs", str(a.pairs), "--species_per_sample", "20", "--read_length", length,
                   "--sequencer", profile, "--fragment_mean", fmean, "--fragment_stdev", fsd, "--seed", "3", "-t", "4",
                   "--protal_metafile", os.path.join(sim, "protal")] + (["--extra_art_args", extra] if extra else [])
        with open(os.path.join(a.out, name + ".log"), "w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        files = sorted(glob.glob(os.path.join(sim, "reads", "*_R[12].fq.gz")))
        total = {}
        for f in files:
            for k, v in sizes(f, a.out).items():
                total[k] = total.get(k, 0) + v
        results[name] = {"unit": "pair", "count": a.pairs, **{k: v / a.pairs for k, v in total.items()}}
        shutil.rmtree(sim)
        print(name, json.dumps(results[name]), flush=True)

    # Drawn reads: Ultima, HiFi (hifi_reads.py), Nanopore (pbsim3 + renaming), per Mb.
    rows = [line.rstrip("\n").split("\t") for line in open(a.genomes)]
    rng = random.Random(3)
    genomes = [{"genome": n, "fasta": p, "weight": rng.lognormvariate(0, 1.5) * int(l)} for n, _, p, l in rows[:60]]

    class Opts:
        pbsim_models = a.pbsim_models
        pbsim = a.pbsim
    for name, text in (("ultima", "ultima:300:40:25:2"), ("hifi", "hifi:15000:3000:3"),
                       ("ont", "qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36")):
        setup = ctd.parse_long_setup(text)
        tmp = os.path.join(a.out, name + "_tmp")
        task = {"sample": name, "out": os.path.join(a.out, name + ".fq.gz"), "bases": a.bases, "setup": setup,
                "model": ctd.pbsim_model(Opts, setup["model"]) if setup["method"] in ctd.PBSIM_METHODS else None,
                "pbsim": a.pbsim, "seed": 9, "tmp": tmp, "genomes": genomes}
        temporary = {}
        if setup["method"] in ctd.PBSIM_METHODS:  # pbsim3's own files, measured before the collector removes them
            os.makedirs(tmp, exist_ok=True)
            names, err = ctd.long_read_templates({**task, "templates": os.path.join(tmp, "templates.fa")})
            assert not err
            command = [a.pbsim, "--strategy", "templ", "--method", setup["method"], f"--{setup['method']}",
                       task["model"], "--template", os.path.join(tmp, "templates.fa"), "--accuracy-mean",
                       str(setup["accuracy"]), "--seed", "9", "--prefix", os.path.join(tmp, "r"), "--id-prefix", "r",
                       "--difference-ratio", setup["ratio"]]
            subprocess.run(command, capture_output=True, check=True)
            temporary = {os.path.basename(f): os.path.getsize(f) / (a.bases / 1e6) for f in glob.glob(os.path.join(tmp, "*"))}
            shutil.rmtree(tmp)
        error = ctd.long_read_sample(task)
        assert not error, error
        s = sizes(task["out"], a.out)
        results[name] = {"unit": "Mb", "count": a.bases / 1e6, **{k: v / (a.bases / 1e6) for k, v in s.items()},
                         "temporary_per_mb": temporary}
        print(name, json.dumps(results[name]), flush=True)
    with open(os.path.join(a.out, "disk_costs.json"), "w") as fh:
        json.dump(results, fh, indent=1)


if __name__ == "__main__":
    main()
