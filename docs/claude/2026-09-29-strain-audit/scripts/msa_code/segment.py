#!/usr/bin/env python3
# print ref/truth/row of one gene segment: segment.py <out> <manifest> <species> <sample> <gene> <from> <to>
import sys, os, importlib.util
out, man, species, sample, gid, a0, b0 = sys.argv[1:8]
gid, a0, b0 = int(gid), int(a0), int(b0)
sys.argv = [sys.argv[0], out, man]
spec = importlib.util.spec_from_file_location("tc", os.path.join(os.path.dirname(os.path.abspath(__file__)), "tc_lib.py"))
tc = importlib.util.module_from_spec(spec); spec.loader.exec_module(tc)
msa = tc.read_fasta(os.path.join(out, "strains", species + ".raw.msa.fna"))
parts = {g: (a, b) for g, a, b in tc.read_parts(os.path.join(out, "strains", species + ".raw.partition.txt"))}
a, b = parts[gid]
ref = msa[species + "_reference"][a:b]; row = msa[sample][a:b]
keep = [i for i, c in enumerate(ref) if c != "-"]
ref = "".join(ref[i] for i in keep); row = "".join(row[i] for i in keep)
tg = tc.true_gene(tc.truth[(sample, species)][0], gid)
print("ref  ", ref[a0:b0]); print("truth", tg[a0:b0]); print("row  ", row[a0:b0])
