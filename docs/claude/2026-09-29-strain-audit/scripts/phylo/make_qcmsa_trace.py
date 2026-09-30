#!/usr/bin/env python3
"""Make qcmsa_trace.py: a copy of protal's scripts/qcmsa.py (unchanged filtering) that also writes
<prefix>.colmap.tsv, the raw-MSA column behind every output column, for the audit's truth
comparisons. The repo's qcmsa.py is not modified."""
import os
import sys

src = os.path.expanduser("~/audit5/src/scripts/qcmsa.py")
dst = os.path.join(os.path.dirname(os.path.abspath(__file__)), "qcmsa_trace.py")
s = open(src).read()
anchor = '    final_seqs = ["".join(row[j] for j in surviving) for row in msa_rows]\n'
assert s.count(anchor) == 1
add = anchor + (
    "    # --- audit tracer (not in protal): raw-MSA column of every output column ---\n"
    "    with open(prefix + '.colmap.tsv', 'w') as _fh:\n"
    "        _fh.write('out_col\\traw_col\\tgene\\n')\n"
    "        for _k, _j in enumerate(surviving):\n"
    "            _fh.write(f'{_k}\\t{col_orig[_j]}\\t{col_gene[_j]}\\n')\n"
)
open(dst, "w").write(s.replace(anchor, add))
print("wrote", dst)
