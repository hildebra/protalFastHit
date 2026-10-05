#!/usr/bin/env python3
"""How the collector's long-read jobs scale with the Scheduler's slots: the current code (jobs on threads of one
process, each chunk of a sample drawing its own templates), and three prototypes:

- processes: the same jobs, each run in a worker process (ProcessPoolExecutor) by the Scheduler's thread;
- onepass:   a sample's templates for all its chunks drawn in one pass over its genomes (each genome read once per
             drawing round instead of once per chunk and round; the same reads, see --check), then the chunks' reads
             (hifi_reads.py, pbsim3) as jobs of their own;
- onepass+processes: both.

Workload: --samples samples of --bases each of one read setup, in chunks of --chunk bases, from a community of
--community genomes (lognormal abundances, sigma 1.5). Prints per run: wall, the main process's CPU, its children's
(pbsim3; worker processes count as children once they end), and their sum over wall (cores used).

usage: gil_scaling.py --scripts PROTAL/scripts --genomes genomes.tsv --out DIR --mode threads --slots 6 [--setup ultima]
"""

import argparse
import bisect
import collections
import concurrent.futures
import gzip
import itertools
import json
import os
import random
import resource
import shutil
import sys
import time

SETUPS = {"ultima": "ultima:300:40:25:2", "pb": "hifi:15000:3000:3", "ont": "qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36"}


def rusage(who):
    r = resource.getrusage(who)
    return r.ru_utime + r.ru_stime


def setup_path(scripts):
    if scripts not in sys.path:
        sys.path.insert(0, scripts)


# ---- the prototype: every chunk's templates in one pass over the sample's genomes ---------------------------

def templates_once(chunks, ctd):
    """long_read_templates for every chunk of one sample (the same genomes, their own seeds, bases, name steps and
    template files), each genome read once per drawing round for all chunks. Each chunk's random stream is consumed
    in the same order as long_read_templates consumes it (its round's plan, then its draws genome by genome in
    sorted order), so the templates are the same. -> ([names of each chunk], None) or (None, why)."""
    genomes = chunks[0]["genomes"]
    cumulative, total = [], 0.0
    for genome in genomes:
        total += genome["weight"]
        cumulative.append(total)
    if total <= 0:
        return None, "no genome with reads to simulate (relative abundances and lengths are 0)"
    states = []
    for task in chunks:
        states.append({"task": task, "rng": random.Random(task["seed"]), "names": [], "bases": 0,
                       "step": task.get("name_step", 1), "offset": task.get("name_offset", 0),
                       "out": open(task["templates"], "wb")})
    try:
        while True:
            plans = []
            for s in states:
                task, rng = s["task"], s["rng"]
                if s["bases"] >= task["bases"]:
                    plans.append(None)
                    continue
                planned, need = collections.defaultdict(list), task["bases"] - s["bases"]
                setup = task["setup"]
                while need > 0:
                    g = bisect.bisect_right(cumulative, rng.random() * total)
                    length = ctd.read_length(rng, setup["length_mean"], setup["length_sd"])
                    planned[g].append(length)
                    need -= length
                plans.append(planned)
            if all(p is None for p in plans):
                break
            for g in sorted(set().union(*(p for p in plans if p))):
                contigs = starts = None
                for s, planned in zip(states, plans):
                    if not planned or g not in planned:
                        continue
                    rng, names, out = s["rng"], s["names"], s["out"]
                    if contigs is None:
                        contigs = [c for c in ctd.read_contigs(genomes[g]["fasta"]) if len(c) >= ctd.PBSIM_MIN_LENGTH]
                        if not contigs:
                            return None, f"{genomes[g]['genome']}: no sequence of {ctd.PBSIM_MIN_LENGTH} bases or more"
                        starts = list(itertools.accumulate(len(c) - ctd.PBSIM_MIN_LENGTH + 1 for c in contigs))
                    for length in planned[g]:
                        at = rng.randrange(starts[-1])
                        k = bisect.bisect_right(starts, at)
                        start = at - (starts[k - 1] if k else 0)
                        seq = contigs[k][start:start + length]
                        if rng.random() < 0.5:
                            seq = seq.translate(ctd.COMPLEMENT)[::-1]
                        names.append(f"g{g}x_{len(names) * s['step'] + s['offset'] + 1}")
                        out.write(b">" + names[-1].encode() + b"\n" + seq + b"\n")
                        s["bases"] += len(seq)
    finally:
        for s in states:
            s["out"].close()
    return [s["names"] for s in states], None


def reads_of_templates(task, names, ctd):
    """long_read_sample after its templates: the reads of task["templates"] into task["out"]."""
    import subprocess
    import glob
    setup = task["setup"]
    tmp = task["tmp"]
    if setup["method"] in ("hifi", "ultima"):
        import hifi_reads
        reads = hifi_reads.simulate(task["templates"], task["out"] + ".partial", setup["q_sd"], task["seed"],
                                    setup.get("q_mean"))
        if reads != len(names):
            return f"{task['sample']}: hifi_reads.py made {reads} reads of {len(names)} templates"
        os.replace(task["out"] + ".partial", task["out"])
        shutil.rmtree(tmp, ignore_errors=True)
        return None
    prefix = os.path.join(tmp, "r")
    command = [task["pbsim"], "--strategy", "templ", "--method", setup["method"], f"--{setup['method']}", task["model"],
               "--template", task["templates"], "--accuracy-mean", str(setup["accuracy"]), "--seed", str(task["seed"]),
               "--prefix", prefix, "--id-prefix", "r"]
    if setup.get("ratio"):
        command += ["--difference-ratio", setup["ratio"]]
    with open(os.path.join(tmp, "pbsim.log"), "w") as log:
        rc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc:
        return f"pbsim failed {rc}"
    fq = sorted(set(glob.glob(prefix + ".fq*") + glob.glob(prefix + ".fastq*")))[0]
    reads = lines = 0
    with (gzip.open(fq, "rb") if fq.endswith(".gz") else open(fq, "rb")) as fin, \
            gzip.open(task["out"] + ".partial", "wb", compresslevel=1) as fout:
        for line in fin:
            if lines % 4 == 0:
                reads += 1
                line = b"@" + names[reads - 1].encode() + b"\n"
            fout.write(line)
            lines += 1
    if reads != len(names):
        return f"pbsim made {reads} reads of {len(names)}"
    os.replace(task["out"] + ".partial", task["out"])
    shutil.rmtree(tmp, ignore_errors=True)
    return None


# ---- worker-process entry points (module level: picklable) ----------------------------------------------------

def worker_sample(scripts, task):
    setup_path(scripts)
    import collect_training_data as ctd
    return ctd.long_read_sample(task), rusage(resource.RUSAGE_SELF)


def worker_templates(scripts, chunks):
    setup_path(scripts)
    import collect_training_data as ctd
    for c in chunks:
        os.makedirs(c["tmp"], exist_ok=True)
    names, error = templates_once(chunks, ctd)
    return names, error


def worker_reads(scripts, task, names):
    setup_path(scripts)
    import collect_training_data as ctd
    return reads_of_templates(task, names, ctd)


# ---- the runs -----------------------------------------------------------------------------------------------

def sample_tasks(a, ctd, out):
    rows = [line.rstrip("\n").split("\t") for line in open(a.genomes)]
    setup = ctd.parse_long_setup(SETUPS[a.setup])

    class Opts:
        pbsim_models = a.pbsim_models
        pbsim = a.pbsim
    model = ctd.pbsim_model(Opts, setup["model"]) if setup["method"] in ctd.PBSIM_METHODS else None
    tasks = []
    for s in range(a.samples):
        rng = random.Random(100 + s)
        community = [{"genome": n, "fasta": p, "weight": rng.lognormvariate(0, 1.5) * int(length)}
                     for n, _, p, length in rng.sample(rows, a.community)]
        sample = f"{a.setup}_s{s + 1}"
        tasks.append({"sample": sample, "out": os.path.join(out, "reads", sample + ".fq.gz"), "bases": a.bases,
                      "setup": setup, "model": model, "pbsim": a.pbsim, "seed": (1 * 1000003 + 7 * 1009 + s) * 101,
                      "tmp": os.path.join(out, "tmp", sample), "genomes": community})
    return tasks


def run(a):
    setup_path(a.scripts)
    import collect_training_data as ctd
    out = os.path.join(a.out, f"{a.setup}_{a.mode}_{a.slots}")
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(os.path.join(out, "reads"))
    tasks = sample_tasks(a, ctd, out)
    need = 2 if tasks[0]["setup"]["method"] in ctd.PBSIM_METHODS else 1
    scheduler = ctd.Scheduler(a.slots)
    pool = concurrent.futures.ProcessPoolExecutor(a.slots) if "processes" in a.mode else None
    parses = collections.Counter()
    w0, s0, c0 = time.time(), rusage(resource.RUSAGE_SELF), rusage(resource.RUSAGE_CHILDREN)
    for task in tasks:
        os.makedirs(task["tmp"], exist_ok=True)
        chunks = ctd.long_read_chunks(task, a.chunk)
        names = [f"long:{task['sample']}:{c + 1}" for c in range(len(chunks))]
        if a.mode.startswith("onepass"):
            for c in chunks:
                c["templates"] = os.path.join(c["tmp"], "templates.fa")
            box = {}

            def draw(chunks=chunks, box=box):
                if pool:
                    box["names"], error = pool.submit(worker_templates, a.scripts, chunks).result()
                else:
                    for c in chunks:
                        os.makedirs(c["tmp"], exist_ok=True)
                    box["names"], error = templates_once(chunks, ctd)
                return error, []
            scheduler.add(f"draw:{task['sample']}", draw, priority=3e9)
            for c, (name, part) in enumerate(zip(names, chunks)):
                def reads(c=c, part=part, box=box):
                    if pool:
                        return pool.submit(worker_reads, a.scripts, part, box["names"][c]).result(), []
                    return reads_of_templates(part, box["names"][c], ctd), []
                scheduler.add(name, reads, need=need, after=[f"draw:{task['sample']}"],
                              priority=ctd.long_read_seconds(part))
        else:
            for name, part in zip(names, chunks):
                if pool:
                    job = lambda part=part: (pool.submit(worker_sample, a.scripts, part).result()[0], [])
                else:
                    job = lambda part=part: (ctd.long_read_sample(part), [])
                scheduler.add(name, job, need=need, priority=ctd.long_read_seconds(part))
        scheduler.add(f"join:{task['sample']}", lambda task=task, chunks=chunks: (ctd.join_chunks(task, chunks), []),
                      after=names, priority=2e9)
    failures = scheduler.run()
    if pool:
        pool.shutdown()
    wall = time.time() - w0
    self_cpu = rusage(resource.RUSAGE_SELF) - s0
    child_cpu = rusage(resource.RUSAGE_CHILDREN) - c0
    if failures:
        sys.exit(f"failures: {failures}")
    result = {"setup": a.setup, "mode": a.mode, "slots": a.slots, "samples": a.samples, "bases": a.bases,
              "chunk": a.chunk, "community": a.community, "wall": round(wall, 2), "self_cpu": round(self_cpu, 2),
              "children_cpu": round(child_cpu, 2), "cores_used": round((self_cpu + child_cpu) / wall, 2),
              "python": sys.version.split()[0], "load": os.getloadavg()[0]}
    print(json.dumps(result), flush=True)
    with open(os.path.join(a.out, "results.jsonl"), "a") as fh:
        fh.write(json.dumps(result) + "\n")
    return out


def check(a):
    """The reads of mode threads and mode onepass (same setup, slots) are the same, decompressed."""
    def reads(mode):
        folder = os.path.join(a.out, f"{a.setup}_{mode}_{a.slots}", "reads")
        return {f: gzip.open(os.path.join(folder, f)).read() for f in sorted(os.listdir(folder))}
    first, second = reads("threads"), reads("onepass")
    same = first == second
    print(f"check {a.setup}: threads and onepass reads {'identical' if same else 'DIFFER'} "
          f"({len(first)} samples, {sum(len(v) for v in first.values()) / 1e6:.1f} MB)")
    return same


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scripts", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", default="threads", choices=["threads", "processes", "onepass", "onepass+processes"])
    ap.add_argument("--slots", type=int, default=6)
    ap.add_argument("--setup", default="ultima", choices=sorted(SETUPS))
    ap.add_argument("--samples", type=int, default=2)
    ap.add_argument("--bases", type=int, default=60_000_000)
    ap.add_argument("--chunk", type=int, default=10_000_000)
    ap.add_argument("--community", type=int, default=300)
    ap.add_argument("--pbsim", default="pbsim")
    ap.add_argument("--pbsim_models")
    ap.add_argument("--check", action="store_true", help="compare the threads and onepass reads of a finished pair")
    ap.add_argument("--keep", action="store_true", help="keep the reads (for --check)")
    a = ap.parse_args()
    if a.check:
        sys.exit(0 if check(a) else 1)
    out = run(a)
    if not a.keep:
        shutil.rmtree(out, ignore_errors=True)


if __name__ == "__main__":
    main()
