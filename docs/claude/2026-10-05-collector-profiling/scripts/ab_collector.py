#!/usr/bin/env python3
"""A/B of the collector's long reads: run its simulate_long (the jobs of long_unit_jobs on its Scheduler) from one
scripts folder on long-read units that replay a fake paired-end point's communities (Ultima, HiFi, Nanopore; samples
in one piece and in chunks; one unit with a host share), timed; --compare checks that two runs made the same reads
(decompressed). With --processes the run starts the collector's worker processes (Workers.started), as its main does.

usage: ab_collector.py --scripts DIR --genomes genomes.tsv --out DIR/label [--slots 6] [--processes]
       ab_collector.py --compare DIR/a DIR/b
"""

import argparse
import glob
import gzip
import json
import os
import random
import resource
import sys
import time


def cpu(who):
    r = resource.getrusage(who)
    return r.ru_utime + r.ru_stime


def run(a):
    sys.path.insert(0, a.scripts)
    import collect_training_data as collect
    import scenarios
    root = os.path.abspath(a.out)
    point = os.path.join(root, "points", "rl150_p1000")
    os.makedirs(os.path.join(point, "sim"), exist_ok=True)
    rows = [line.rstrip("\n").split("\t") for line in open(a.genomes)]
    rng = random.Random(11)
    samples = [f"rl150_p1000_s_{s + 1}" for s in range(a.samples)]
    with open(os.path.join(point, "sim", "manifest.tsv"), "w") as fh:
        fh.write("sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\tvertical_coverage\trelative_abundance"
                 "\tfastq_r1\tfastq_r2\tfasta_path\tart_seed\n")
        for sample in samples:
            for name, lineage, path, length in rng.sample(rows, a.community):
                fh.write(f"{sample}\t{name}\t{name}\t{lineage}\t{length}\t10\t1\t{rng.lognormvariate(0, 1.5):.6f}"
                         f"\tr1\tr2\t{path}\t1\n")
    with open(os.path.join(point, "sim", "protal.meta"), "w") as fh:
        fh.write(f"#OUTPUT_DIR\t{point}/protal\n#INPUT_DIR\t{point}/sim/reads\n"
                 "#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\n")
        fh.writelines(f"{s}\ta\tb\tc\td\te\t/truth/{s}\n" for s in samples)
    host_fa = os.path.join(root, "host.fa")
    with open(host_fa, "wb") as out:
        for _, _, path, _ in rows[-8:]:
            out.write(gzip.open(path).read())
    host_folder = scenarios.prepare_host(host_fa, os.path.join(root, "host"))
    opts = argparse.Namespace(out=root, seed=1, pbsim="pbsim", pbsim_models=a.pbsim_models, samples=a.samples,
                              long_read_chunk=a.chunk, host_folder=host_folder)
    units = []
    for kind, text in (("ultima", "ultima:300:40:25:2"), ("hifi", "hifi:15000:3000:3"),
                       ("ont", "qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36")):
        setup = collect.parse_long_setup(text)
        for bases in (a.chunk // 2, a.chunk * 4):  # in one piece; in 4 chunks
            name = f"{kind}_b{bases}"
            units.append({"type": "ont" if kind == "ont" else "pb", "name": name, "bases": bases, "setup": setup,
                          "samples": a.samples, "communities": [{"name": "rl150_p1000"}],
                          "point": {"name": name, "read_length": str(setup["length_mean"]), "read_pairs": str(bases)}})
    hosted = dict(units[3], name="hifi_host",  # HiFi in chunks, half of its bases from the host host_share=0.5, seed_index=1_000_123,
                  point={"name": "hifi_host", "read_length": "15000", "read_pairs": str(a.chunk * 4)})
    units.append(hosted)
    w0, s0, c0 = time.time(), cpu(resource.RUSAGE_SELF), cpu(resource.RUSAGE_CHILDREN)
    if a.processes:
        with collect.Workers.started(a.slots):
            failures = collect.simulate_long(list(enumerate(units)), opts, a.slots)
    else:
        failures = collect.simulate_long(list(enumerate(units)), opts, a.slots)
    wall = time.time() - w0
    result = {"scripts": a.scripts, "processes": a.processes, "slots": a.slots, "wall": round(wall, 1),
              "self_cpu": round(cpu(resource.RUSAGE_SELF) - s0, 1),
              "children_cpu": round(cpu(resource.RUSAGE_CHILDREN) - c0, 1), "failures": failures}
    print(json.dumps(result), flush=True)
    if failures:
        sys.exit(1)


def compare(a, b):
    def reads(root):
        return {os.path.relpath(f, root): f for f in glob.glob(os.path.join(root, "points", "*", "sim", "reads", "*.fq.gz"))}
    first, second = reads(a), reads(b)
    if sorted(first) != sorted(second):
        print(f"different samples: {sorted(set(first) ^ set(second))}")
        return False
    same = True
    for name in sorted(first):
        x, y = gzip.open(first[name]).read(), gzip.open(second[name]).read()
        if x != y:
            print(f"DIFFERENT: {name} ({len(x)} and {len(y)} bytes)")
            same = False
    left = glob.glob(os.path.join(b, "points", "*", "sim", "tmp"))
    print(f"{len(first)} samples compared: {'identical' if same else 'different'}"
          f"{'; temporary folders left: ' + ', '.join(left) if left else ''}")
    return same and not left


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scripts")
    ap.add_argument("--genomes")
    ap.add_argument("--out")
    ap.add_argument("--pbsim_models")
    ap.add_argument("--slots", type=int, default=6)
    ap.add_argument("--samples", type=int, default=3)
    ap.add_argument("--community", type=int, default=150)
    ap.add_argument("--chunk", type=int, default=4_000_000)
    ap.add_argument("--processes", action="store_true")
    ap.add_argument("--compare", nargs=2)
    a = ap.parse_args()
    if a.compare:
        sys.exit(0 if compare(*a.compare) else 1)
    run(a)


if __name__ == "__main__":
    main()
