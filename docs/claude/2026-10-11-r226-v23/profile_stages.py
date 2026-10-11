"""Where the r226 build's protal runs spent their time per sample: the aligning (logs/protal_cpu_<collection>.tsv, wall
and CPU on 84 threads), the profiling stages ("Profiling sample ... took", wall) and the SAM parsing threads' time
("SAM read: ... parsing", CPU summed over the threads; the per-record feature evidence runs there), from
logs/protal_runs_*.log, summed by read type and by kind of sample (design or scenario).

usage: profile_stages.py <build dir>
"""
import re
import sys
from collections import defaultdict

D = sys.argv[1]
PROF = re.compile(r"^Profiling sample (\S+) took ([\d.]+)s: (.*)$")
PARSE = re.compile(r"^Sample (\S+) SAM read: .*\(parsing ([\d.]+) s, adding to the taxa ([\d.]+) s")
STAGE = re.compile(r"([a-zA-Z][a-zA-Z ,-]+?) ([\d.]+)s")


def read_type(name):
    if name.startswith("pb_") or "_pb_" in name:
        return "pb"
    if name.startswith("ont_") or "_ont_" in name:
        return "ont"
    if name.endswith("_se") or "_se_" in name:
        return "se"
    return "pe"


def kind(name):
    return "scenario" if name.startswith("sc_") else "design"


align = defaultdict(lambda: [0.0, 0.0, 0])
prof = defaultdict(lambda: defaultdict(float))
for coll in ("training", "test"):
    with open(f"{D}/logs/protal_cpu_{coll}.tsv") as f:
        next(f)
        for line in f:
            run, stage, wall, cpu, *_ = line.rstrip("\n").split("\t")
            if stage.startswith("aligning "):
                s = stage[len("aligning "):]
                a = align[(read_type(s), kind(s))]
                a[0] += float(wall)
                a[1] += float(cpu)
                a[2] += 1
    with open(f"{D}/logs/protal_runs_{coll}.log") as f:
        for line in f:
            q = PARSE.match(line)
            if q:
                p = prof[(read_type(q.group(1)), kind(q.group(1)))]
                p["parse CPU"] += float(q.group(2))
                p["add CPU"] += float(q.group(3))
                continue
            m = PROF.match(line)
            if not m:
                continue
            s, total, rest = m.group(1), float(m.group(2)), m.group(3)
            p = prof[(read_type(s), kind(s))]
            p["total"] += total
            p["samples"] += 1
            for st, sec in STAGE.findall(rest):
                p[st.strip().lstrip(", ")] += float(sec)

stages = ["reading the SAM", "record evidence and sample context", "read EM", "congener distances",
          "SNP post-processing", "scoring", "writing the profile"]
print("per read type and kind of sample, seconds summed over the samples (training + test collections)")
print("aligning: wall and CPU on 84 threads; profiling: wall of the stage (its parts run in parallel, so they may sum "
      "above the total)")
hdr = ["type", "kind", "samples", "align wall", "align CPU", "profile wall", "parse CPU", "parse/align CPU"] + stages
print("\t".join(hdr))
for key in sorted(set(align) | set(prof)):
    a, p = align[key], prof[key]
    row = [key[0], key[1], str(a[2]), f"{a[0]:.0f}", f"{a[1]:.0f}", f"{p['total']:.0f}", f"{p['parse CPU']:.0f}",
           f"{p['parse CPU'] / max(a[1], 1):.2f}"] + [f"{p[s]:.0f}" for s in stages]
    print("\t".join(row))
