#!/usr/bin/env python3
"""Write qcmsa_colmap.py: a verbatim copy of ~/audit5/src/scripts/qcmsa.py that additionally writes
<prefix>.colmap.txt (kept row names, then the 0-based raw-MSA column of every output column), so
the filtered MSA can be mapped back to reference coordinates. The filtering logic is untouched;
evaluate.py checks that the copy's .msa.fna equals the one protal's qcmsa run wrote."""
import os

SRC = os.path.expanduser("~/audit5/src/scripts/qcmsa.py")
DST = os.path.expanduser("~/audit5/accuracy/qcmsa_colmap.py")
src = open(SRC).read()
anchor = '    final_seqs = ["".join(row[j] for j in surviving) for row in msa_rows]\n'
assert src.count(anchor) == 1
add = anchor + (
    '    # [accuracy audit] kept rows and the original (0-based) raw-MSA column of every output column\n'
    '    with open(prefix + ".colmap.txt", "w") as _fh:\n'
    '        _fh.write("#names\\t" + "\\t".join(kept_names) + "\\n")\n'
    '        _fh.write("\\n".join(str(col_orig[j]) for j in surviving) + "\\n")\n'
)
open(DST, "w").write(src.replace(anchor, add))
