#!/usr/bin/env python3
"""make_colmap.py QCMSA_PY OUT_PY - an instrumented copy of a qcmsa.py that also writes
<prefix>.colmap.txt (kept rows, then the raw-MSA column of every output column), as the strain
audit's make_qcmsa_colmap.py does for one fixed qcmsa; the filtering is untouched, and evaluate.py
checks that the copy's .msa.fna equals the run's."""
import sys

src, dst = sys.argv[1], sys.argv[2]
text = open(src).read()
anchor = '    final_seqs = ["".join(row[j] for j in surviving) for row in msa_rows]\n'
assert text.count(anchor) == 1, f"{src}: the anchor is not there once"
add = anchor + (
    '    # [accuracy audit] kept rows and the original (0-based) raw-MSA column of every output column\n'
    '    with open(prefix + ".colmap.txt", "w") as _fh:\n'
    '        _fh.write("#names\\t" + "\\t".join(kept_names) + "\\n")\n'
    '        _fh.write("\\n".join(str(col_orig[j]) for j in surviving) + "\\n")\n'
)
open(dst, "w").write(text.replace(anchor, add))
