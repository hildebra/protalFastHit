"""Copy a training table (collect_training_data.py) keeping only a random share of one class of present taxa, the other
rows unchanged and byte for byte (the trainer's features are protal's, per row, so the rows left keep their values).

usage: thin_table.py <table> <out> <real|insilico> <share kept> [seed]
  real: present taxa simulated from another real genome of their species (meta_rep_genome 0, meta_insilico_strain 0)
  insilico: present taxa simulated from an in-silico strain (meta_insilico_strain 1)
"""
import random
import sys

src, out, which, keep = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
rng = random.Random(int(sys.argv[5]) if len(sys.argv) > 5 else 1)

kept = dropped = 0
with open(src) as f, open(out, "w") as o:
    header = f.readline()
    o.write(header)
    col = {c: i for i, c in enumerate(header.rstrip("\n").split("\t"))}
    t, rep, ins = col["truth"], col["meta_rep_genome"], col["meta_insilico_strain"]
    for line in f:
        v = line.split("\t", max(t, rep, ins) + 1)
        if v[t] == "1":
            is_ins = v[ins] == "1"
            is_real = not is_ins and v[rep] == "0"
            if (is_real if which == "real" else is_ins):
                if rng.random() >= keep:
                    dropped += 1
                    continue
                kept += 1
        o.write(line)
print(f"{which}: kept {kept}, dropped {dropped} present rows ({kept / max(kept + dropped, 1):.3f})")
