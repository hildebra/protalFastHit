#!/usr/bin/env python3
"""What a long-read sample of collect_training_data.py costs, step by step: the collector's long_read_sample
redone with a clock around each step (drawing the templates in Python, then for ONT pbsim3 --strategy templ and
the Python renaming and gzip of its reads, for PacBio HiFi hifi_reads.simulate), for a community of a paired-end
sample simulated by sim_cost.sh (its manifest: genome, fasta_path, relative_abundance, genome_length).

usage: long_cost.py MANIFEST OUT --bases 30000000 [--pbsim PATH --models DIR] [--types ont,pb]
Prints one tab-separated line per step: type, bases, step, wall seconds, CPU seconds of this process, of children.
"""
import argparse
import csv
import glob
import gzip
import os
import resource
import subprocess
import sys
import time

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import collect_training_data as ctd  # noqa: E402


def clock():
    own, kids = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.time(), own.ru_utime + own.ru_stime, kids.ru_utime + kids.ru_stime


def report(kind, bases, step, before):
    after = clock()
    print(f"{kind}\t{bases}\t{step}\t{after[0] - before[0]:.2f}\t{after[1] - before[1]:.2f}\t{after[2] - before[2]:.2f}",
          flush=True)
    return clock()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest")
    p.add_argument("out")
    p.add_argument("--bases", type=int, default=30_000_000)
    p.add_argument("--types", default="ont,pb")
    p.add_argument("--pbsim", default=os.path.expanduser("~/micromamba/envs/protal-db-build/bin/pbsim"))
    p.add_argument("--models", help="folder of QSHMM-ONT-HQ.model (default: found next to pbsim's env)")
    opts = p.parse_args()
    os.makedirs(opts.out, exist_ok=True)
    with open(opts.manifest) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")]
    sample = rows[0]["sample"]
    rows = [r for r in rows if r["sample"] == sample]
    genomes = [{"genome": r["genome"], "fasta": r["fasta_path"],
                "weight": float(r["relative_abundance"]) * float(r["genome_length"])} for r in rows]
    models = opts.models or os.path.dirname(next(iter(glob.glob(os.path.join(
        os.path.dirname(os.path.dirname(opts.pbsim)), "**", "QSHMM-ONT-HQ.model"), recursive=True))))
    print("type\tbases\tstep\twall_s\tcpu_s\tchild_cpu_s", flush=True)
    print(f"#\tcommunity {sample}: {len(genomes)} genomes", flush=True)
    for kind in opts.types.split(","):
        setup = ctd.parse_long_setup("hifi:15000:3000:3" if kind == "pb" else "qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36")
        tmp = os.path.join(opts.out, kind)
        os.makedirs(tmp, exist_ok=True)
        task = {"seed": 11, "genomes": genomes, "setup": setup, "bases": opts.bases,
                "templates": os.path.join(tmp, "templates.fa")}
        t = clock()
        names, error = ctd.long_read_templates(task)
        if error:
            sys.exit(error)
        t = report(kind, opts.bases, f"templates ({len(names)} reads, Python)", t)
        out = os.path.join(tmp, "reads.fq.gz")
        if kind == "pb":
            import hifi_reads
            hifi_reads.simulate(task["templates"], out, setup["q_sd"], task["seed"])
            report(kind, opts.bases, "hifi_reads.simulate (numpy, gzip)", t)
            continue
        prefix = os.path.join(tmp, "r")
        command = [opts.pbsim, "--strategy", "templ", "--method", "qshmm", "--qshmm",
                   os.path.join(models, "QSHMM-ONT-HQ.model"), "--template", task["templates"],
                   "--accuracy-mean", str(setup["accuracy"]), "--seed", "11", "--prefix", prefix, "--id-prefix", "r"]
        if setup.get("ratio"):
            command += ["--difference-ratio", setup["ratio"]]
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        t = report(kind, opts.bases, "pbsim3 --strategy templ", t)
        fastq = sorted(glob.glob(prefix + ".fq*"))[0]
        reads, lines = 0, 0
        with (gzip.open(fastq, "rb") if fastq.endswith(".gz") else open(fastq, "rb")) as fin, \
                gzip.open(out, "wb", compresslevel=1) as fout:
            for line in fin:
                if lines % 4 == 0:
                    reads += 1
                    line = b"@" + names[reads - 1].encode() + b"\n"
                fout.write(line)
                lines += 1
        report(kind, opts.bases, "rename + gzip in Python", t)


if __name__ == "__main__":
    main()
