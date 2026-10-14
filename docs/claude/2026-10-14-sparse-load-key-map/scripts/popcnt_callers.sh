#!/bin/bash
# The functions of protal (HEAD build) that call libgcc's __popcountdi2, demangled.
objdump -d --no-show-raw-insn ~/sparseprof/bin/protal_base | awk '/^[0-9a-f]+ <.*>:$/ {f = $0} /call.*<__popcountdi2/ {print f}' | sort | uniq -c | sort -rn | c++filt | cut -c1-220 | head -30
