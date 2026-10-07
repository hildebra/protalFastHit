#!/usr/bin/env python3
"""The read models of simulate_metagenomes --long_templates against the programs they port, on the same templates:
hifi_reads.py (HiFi and its flow model for Ultima reads) and pbsim3 --strategy templ --method qshmm (Nanopore). Per
model and program: reads, read length over template length, base qualities (mean, per-read quantiles, the share of
reads at Q93 throughout, the total variation distance of the quality histograms), the error rate the qualities expect,
the errors made (as each program counts them) and the time on one core.

    compare_models.py --scripts SCRIPTS --simulator BIN --pbsim PBSIM --model QSHMM-ONT-HQ.model --genomes DIR --out DIR

SCRIPTS: protal's scripts folder (hifi_reads.py, compressed.py); DIR: *.fna.gz genomes to draw templates from.
"""
import argparse
import collections
import glob
import gzip
import os
import random
import re
import resource
import subprocess
import sys
import time


def read_genomes(folder, n):
    """The contigs of 100 bases or more of the first n genomes."""
    genomes = []
    for path in sorted(glob.glob(os.path.join(folder, "*.fna.gz")))[:n]:
        data = gzip.decompress(open(path, "rb").read())
        contigs = [b"".join(r.partition(b"\n")[2].split()).upper() for r in (b"\n" + data).split(b"\n>")[1:]]
        genomes.append([c for c in contigs if len(c) >= 100])
    return genomes


def write_templates(genomes, path, count, mean, sd, seed):
    """`count` templates as the collector draws them: a genome, a gamma length in [100, 1e6], a uniform start on a
    contig (cut at its end), either strand."""
    rng = random.Random(seed)
    complement = bytes.maketrans(b"ACGTN", b"TGCAN")
    shape, scale = (mean / sd) ** 2, sd * sd / mean
    with open(path, "wb") as out:
        for i in range(count):
            contigs = genomes[rng.randrange(len(genomes))]
            while True:
                length = int(round(rng.gammavariate(shape, scale)))
                if 100 <= length <= 1000000:
                    break
            weights = [len(c) - 99 for c in contigs]
            contig = rng.choices(contigs, weights)[0]
            start = rng.randrange(len(contig) - 99)
            seq = contig[start:start + length]
            if rng.random() < 0.5:
                seq = seq.translate(complement)[::-1]
            out.write(b">t%d\n%s\n" % (i + 1, seq))


def read_fasta_lengths(path):
    lengths, name = {}, None
    for line in open(path, "rb"):
        if line.startswith(b">"):
            name = line[1:].split()[0].decode()
            lengths[name] = 0
        else:
            lengths[name] += len(line.strip())
    return lengths


def fastq_records(path, scripts):
    sys.path.insert(0, scripts)
    import compressed
    with compressed.open_read(path) as fh:
        while True:
            name = fh.readline()
            if not name:
                return
            seq, _, qual = fh.readline().rstrip(), fh.readline(), fh.readline().rstrip()
            yield name[1:].split()[0].decode(), seq, qual


def timed(command, **kwargs):
    """Runs a command; -> (wall seconds, CPU seconds of it, its output)."""
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    began = time.perf_counter()
    done = subprocess.run(command, capture_output=True, text=True, **kwargs)
    wall = time.perf_counter() - began
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    if done.returncode != 0:
        sys.exit(f"{' '.join(command)} failed: {done.stderr[-2000:]}")
    return wall, after.ru_utime - before.ru_utime + after.ru_stime - before.ru_stime, done.stdout + done.stderr


def summary(path, templates, scripts):
    """Statistics of a FASTQ of one read per template (named after it)."""
    ratios, read_q, hist = [], [], collections.Counter()
    bases = expected = q93 = 0
    for name, seq, qual in fastq_records(path, scripts):
        q = [c - 33 for c in qual]
        ratios.append(len(seq) / templates[name])
        read_q.append(sum(q) / len(q) if q else 0)
        hist.update(q)
        bases += len(seq)
        expected += sum(10 ** (-x / 10) for x in q)
        q93 += all(x == 93 for x in q)
    n = len(ratios)
    mean = sum(ratios) / n
    sd = (sum((r - mean) ** 2 for r in ratios) / n) ** 0.5
    qs = sorted(read_q)
    return {"reads": n, "ratio": mean, "ratio_sd": sd, "base_q": sum(k * v for k, v in hist.items()) / bases,
            "read_q": [qs[int(p * (n - 1))] for p in (0.05, 0.25, 0.5, 0.75, 0.95)], "q93": q93 / n,
            "expected": expected / bases, "hist": hist, "bases": bases}


def tvd(a, b):
    na, nb = sum(a.values()), sum(b.values())
    return 0.5 * sum(abs(a[k] / na - b[k] / nb) for k in set(a) | set(b))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--scripts", required=True)
    p.add_argument("--simulator", required=True)
    p.add_argument("--pbsim", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--genomes", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--n_genomes", type=int, default=40)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sys.path.insert(0, a.scripts)
    genomes = read_genomes(a.genomes, a.n_genomes)
    setups = [("hifi", "hifi:15000:3000:3", 2000, 15000, 3000),
              ("ultima", "ultima:300:40:25:2", 60000, 300, 40),
              ("ont", "qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36", 4000, 8000, 6000)]
    rows = []
    for kind, setup, count, mean, sd in setups:
        templates_path = os.path.join(a.out, f"{kind}_templates.fa")
        write_templates(genomes, templates_path, count, mean, sd, seed=11)
        templates = read_fasta_lengths(templates_path)
        mb = sum(templates.values()) / 1e6
        cpp_out = os.path.join(a.out, f"{kind}_cpp.fq.zst")
        stats = os.path.join(a.out, f"{kind}_cpp.stats")
        command = [a.simulator, "--long_templates", templates_path, "--long_setup", setup, "--long_out", cpp_out,
                   "--long_stats", stats, "--seed", "5"] + (["--long_model", a.model] if kind == "ont" else [])
        wall, cpu, _ = timed(command)
        cpp = summary(cpp_out, templates, a.scripts)
        made = dict(zip(*[line.rstrip("\n").split("\t") for line in open(stats)]))
        cpp["errors"] = int(made["errors"]) / int(made["read_bases"])
        cpp["cpu_per_mb"] = cpu / mb
        if kind == "ont":
            prefix = os.path.join(a.out, "pbsim")
            wall, cpu, log = timed([a.pbsim, "--strategy", "templ", "--method", "qshmm", "--qshmm", a.model, "--template",
                                    templates_path, "--accuracy-mean", "0.97", "--difference-ratio", "39:24:36",
                                    "--seed", "5", "--prefix", prefix, "--id-prefix", "r"])
            fq = [f for f in glob.glob(prefix + ".fq*") + glob.glob(prefix + ".fastq*")][0]
            renamed = os.path.join(a.out, "ont_pbsim_named.fq.gz")  # pbsim names reads r_<n>, n the template's place
            names = list(templates)
            with gzip.open(renamed, "wt") as out:
                for name, seq, qual in fastq_records(fq, a.scripts):
                    out.write(f"@{names[int(name.split('_')[1]) - 1]}\n{seq.decode()}\n+\n{qual.decode()}\n")
            other = summary(renamed, templates, a.scripts)
            rates = {k: float(v) for k, v in re.findall(r"(substitution|insertion|deletion) rate\. : ([0-9.]+)", log)}
            other["errors"] = sum(rates.values())
            other["cpu_per_mb"] = cpu / mb
            name = "pbsim3 3.0.5"
        else:
            import hifi_reads
            py_out = os.path.join(a.out, f"{kind}_py.fq.gz")
            q_mean = 25.0 if kind == "ultima" else None
            q_sd = 2.0 if kind == "ultima" else 3.0
            began = time.process_time()
            hifi_reads.simulate(templates_path, py_out, q_sd, 5, q_mean)
            cpu = time.process_time() - began
            other = summary(py_out, templates, a.scripts)
            import numpy as np  # hifi_reads.mutate's own error count, on the same templates (another draw)
            seqs = [s for _, s in hifi_reads.read_fasta(templates_path)]
            events = 0
            rng = np.random.default_rng(6)
            for i in range(0, len(seqs), 500):
                _, st = hifi_reads.mutate(seqs[i:i + 500], rng, q_sd, q_mean)
                events += int(st["events"].sum())
            other["errors"] = events / other["bases"]
            other["cpu_per_mb"] = cpu / mb
            name = "hifi_reads.py"
        for program, s in (("simulate_metagenomes", cpp), (name, other)):
            rows.append((kind, program, s))
        rows.append((kind, "TVD of the quality histograms", {"tvd": tvd(cpp["hist"], other["hist"])}))
    print("model\tprogram\treads\tlength ratio (SD)\tmean base Q\tper-read Q 5/25/50/75/95%\tQ93 reads\t"
          "errors expected by Q\terrors made\tCPU s per Mb")
    for kind, program, s in rows:
        if "tvd" in s:
            print(f"{kind}\t{program}\t{s['tvd']:.4f}")
            continue
        print(f"{kind}\t{program}\t{s['reads']}\t{s['ratio']:.4f} ({s['ratio_sd']:.4f})\t{s['base_q']:.2f}\t"
              f"{'/'.join(f'{q:.1f}' for q in s['read_q'])}\t{s['q93']:.4f}\t{s['expected']:.5f}\t{s['errors']:.5f}\t"
              f"{s['cpu_per_mb']:.4f}")


if __name__ == "__main__":
    main()
