#!/usr/bin/env python3
"""Unit costs of the collector's long-read code on one core: reading a genome (zlib, then the Python parse, and a
one-pass translate parse), and per read the template drawing, hifi_reads.py (HiFi, Ultima flow model) and pbsim3's
renaming, each as wall and CPU of this process.

usage: unit_costs.py --scripts PROTAL/scripts --genomes genomes.tsv [--pbsim pbsim --pbsim_models DIR] [--n 40]
"""

import argparse
import gzip
import os
import random
import resource
import sys
import tempfile
import time


def cpu():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def children_cpu():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime


UPPER = bytes.maketrans(b"abcdefghijklmnopqrstuvwxyz", b"ABCDEFGHIJKLMNOPQRSTUVWXYZ")
WHITESPACE = b" \t\n\r\x0b\x0c"


def parse_translate(data):
    """read_contigs' parse in one C pass per record (whitespace deleted and upper case at once)."""
    return [record.partition(b"\n")[2].translate(UPPER, WHITESPACE) for record in (b"\n" + data).split(b"\n>")[1:]]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scripts", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--pbsim", default="pbsim")
    ap.add_argument("--pbsim_models")
    ap.add_argument("--n", type=int, default=40, help="genomes for the parse timings")
    ap.add_argument("--community", type=int, default=200, help="genomes of the community the reads are drawn from")
    a = ap.parse_args()
    sys.path.insert(0, a.scripts)
    import collect_training_data as ctd
    import hifi_reads

    rows = [line.rstrip("\n").split("\t") for line in open(a.genomes)]
    print(f"python {sys.version.split()[0]}")

    # 1. A genome read as the collector reads it: file, zlib, parse.
    t = {"read": 0.0, "zlib": 0.0, "parse": 0.0, "translate": 0.0}
    mb = 0
    for name, _, path, _ in rows[:a.n]:
        c0 = cpu()
        with open(path, "rb") as fh:
            data = fh.read()
        c1 = cpu()
        plain = gzip.decompress(data)
        c2 = cpu()
        contigs = [b"".join(record.partition(b"\n")[2].split()).upper() for record in (b"\n" + plain).split(b"\n>")[1:]]
        c3 = cpu()
        fast = parse_translate(plain)
        c4 = cpu()
        assert fast == contigs, name
        assert contigs == ctd.read_contigs(path)
        t["read"] += c1 - c0
        t["zlib"] += c2 - c1
        t["parse"] += c3 - c2
        t["translate"] += c4 - c3
        mb += len(plain) / 1e6
    n = min(a.n, len(rows))
    print(f"genome read ({n} genomes, {mb / n:.2f} MB of FASTA each): file {1e3 * t['read'] / n:.1f} ms, "
          f"zlib {1e3 * t['zlib'] / n:.1f} ms (GIL released), parse {1e3 * t['parse'] / n:.1f} ms (GIL held), "
          f"translate parse {1e3 * t['translate'] / n:.1f} ms (GIL held; same contigs)")

    # 2. Per read: one long-read sample's steps, each timed, on a community of `community` genomes.
    rng = random.Random(7)
    community = [{"genome": name, "fasta": path, "weight": rng.lognormvariate(0, 1.5) * int(length)}
                 for name, _, path, length in rows[:a.community]]
    setups = [("ultima", "ultima:300:40:25:2", 30_000_000), ("pb", "hifi:15000:3000:3", 60_000_000),
              ("ont", "qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36", 30_000_000)]
    reads_of = {}

    class Opts:
        pbsim_models = a.pbsim_models
        pbsim = a.pbsim

    # read_contigs counted and timed (it is called once per genome and drawing round)
    calls = {"n": 0, "cpu": 0.0}
    original = ctd.read_contigs

    def counted(path):
        c = cpu()
        out = original(path)
        calls["n"] += 1
        calls["cpu"] += cpu() - c
        return out
    ctd.read_contigs = counted
    for kind, text, bases in setups:
        setup = ctd.parse_long_setup(text)
        with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(a.genomes))) as tmp:
            task = {"sample": kind, "out": os.path.join(tmp, "reads.fq.gz"), "bases": bases, "setup": setup,
                    "model": ctd.pbsim_model(Opts, setup["model"]) if setup["method"] in ctd.PBSIM_METHODS else None,
                    "pbsim": a.pbsim, "seed": 11, "tmp": os.path.join(tmp, "t"), "genomes": community,
                    "templates": os.path.join(tmp, "templates.fa")}
            os.makedirs(task["tmp"])
            calls.update(n=0, cpu=0.0)
            w0, c0, k0 = time.time(), cpu(), children_cpu()
            names, error = ctd.long_read_templates(task)
            w1, c1, k1 = time.time(), cpu(), children_cpu()
            assert not error, error
            reads = len(names)
            parse_n, parse_cpu = calls["n"], calls["cpu"]
            if setup["method"] in ("hifi", "ultima"):
                made = hifi_reads.simulate(task["templates"], task["out"], setup["q_sd"], task["seed"], setup.get("q_mean"))
                assert made == reads
                w2, c2, k2 = time.time(), cpu(), children_cpu()
                print(f"{kind}: {bases / 1e6:.0f} Mb, {reads} reads, {parse_n} genome reads ({parse_n / len(community):.2f} "
                      f"per genome): templates {w1 - w0:.1f} s wall, {c1 - c0:.1f} s CPU (genome reads {parse_cpu:.1f} s, "
                      f"the rest {1e6 * (c1 - c0 - parse_cpu) / reads:.1f} us per read); hifi_reads {w2 - w1:.1f} s wall, "
                      f"{c2 - c1:.1f} s CPU ({1e6 * (c2 - c1) / reads:.1f} us per read, {1e3 * (c2 - c1) / (bases / 1e6):.0f} "
                      f"ms per Mb)")
            else:
                prefix = os.path.join(task["tmp"], "r")
                command = [a.pbsim, "--strategy", "templ", "--method", setup["method"], f"--{setup['method']}",
                           task["model"], "--template", task["templates"], "--accuracy-mean", str(setup["accuracy"]),
                           "--seed", "11", "--prefix", prefix, "--id-prefix", "r", "--difference-ratio", setup["ratio"]]
                import subprocess
                with open(os.path.join(tmp, "pbsim.log"), "w") as log:
                    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
                w2, c2, k2 = time.time(), cpu(), children_cpu()
                import glob
                fq = sorted(glob.glob(prefix + ".fq*"))[0]
                count = lines = 0
                with (gzip.open(fq, "rb") if fq.endswith(".gz") else open(fq, "rb")) as fin, \
                        gzip.open(task["out"], "wb", compresslevel=1) as fout:
                    for line in fin:
                        if lines % 4 == 0:
                            count += 1
                            line = b"@" + names[count - 1].encode() + b"\n"
                        fout.write(line)
                        lines += 1
                w3, c3 = time.time(), cpu()
                print(f"{kind}: {bases / 1e6:.0f} Mb, {reads} reads, {parse_n} genome reads ({parse_n / len(community):.2f} "
                      f"per genome): templates {w1 - w0:.1f} s wall, {c1 - c0:.1f} s CPU (genome reads {parse_cpu:.1f} s, "
                      f"the rest {1e6 * (c1 - c0 - parse_cpu) / reads:.1f} us per read); pbsim {w2 - w1:.1f} s wall, "
                      f"{k2 - k1:.1f} s CPU ({1e3 * (k2 - k1) / (bases / 1e6):.0f} ms per Mb); renaming {w3 - w2:.1f} s "
                      f"wall, {c3 - c2:.1f} s CPU ({1e6 * (c3 - c2) / reads:.1f} us per read, "
                      f"{1e3 * (c3 - c2) / (bases / 1e6):.0f} ms per Mb)")
            reads_of[kind] = reads

    # 3. The host scenario's paired-end fragments (scenarios.host_pe_chunk's loop, without ART): per fragment.
    import scenarios
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(a.genomes))) as tmp:
        seq_path = os.path.join(tmp, "host.fa")
        with open(seq_path, "wb") as fh:
            for name, _, path, _ in rows[:20]:
                fh.write(gzip.open(path).read())
        folder = scenarios.prepare_host(seq_path, os.path.join(tmp, "host"))
        host = scenarios.Host.of(folder)
        r = random.Random(3)
        n = 200_000
        c0 = cpu()
        with open(os.path.join(tmp, "frag.fa"), "wb") as fh:
            for i in range(n):
                size = max(151, int(round(r.gauss(350, 50))))
                seq = host.draw(r, size, min_length=151)
                if r.random() < 0.5:
                    seq = seq.translate(scenarios.COMPLEMENT)[::-1]
                fh.write(b">h%d_%d\n%s\n" % (1, i + 1, seq))
        print(f"host fragments: {1e6 * (cpu() - c0) / n:.2f} us per fragment (GIL held), {n} fragments")


if __name__ == "__main__":
    main()
