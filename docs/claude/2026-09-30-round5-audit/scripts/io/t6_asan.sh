#!/bin/bash
# T6: the ASan/UBSan build on reader edge cases, the three SAM formats, and --profile_only; the
# I/O unit tests under the sanitizers.
set -u
source $(dirname "$0")/lib.sh
P=$W/src/build-asan/protal
export ASAN_OPTIONS=detect_leaks=0:abort_on_error=0 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=0
T=$W/t6; rm -rf $T; mkdir -p $T/in; cd $T
I=$T/in; A1=$R/sa_R1.fq; A2=$R/sa_R2.fq
# FASTA with blank lines after long lines (StripString on an empty line), CRLF FASTA, empty record
python3 - <<'EOF'
src = open("/home/falk/audit6/io/reads/sa_R1.fq").read().split("\n")
src2 = open("/home/falk/audit6/io/reads/sa_R2.fq").read().split("\n")
def fa(lines, blank):
    out = []
    for i in range(0, 800, 4):
        out.append(">" + lines[i][1:])
        out.append(lines[i + 1])
        if blank and i % 8 == 0:
            out.append("")          # blank line inside the sequence part
    return "\n".join(out) + "\n\n\n"
open("/home/falk/audit6/io/t6/in/blank_R1.fa", "w").write(fa(src, True))
open("/home/falk/audit6/io/t6/in/blank_R2.fa", "w").write(fa(src2, True))
open("/home/falk/audit6/io/t6/in/crlf_R1.fa", "w").write(fa(src, False).replace("\n", "\r\n"))
open("/home/falk/audit6/io/t6/in/crlf_R2.fa", "w").write(fa(src2, False).replace("\n", "\r\n"))
# FASTQ: record 3 has a quality shorter than its sequence; record 5 a header "@" alone
def fq(lines, qshort, bare):
    out = []
    for n, i in enumerate(range(0, 400, 4)):
        h, s, p, q = lines[i:i + 4]
        if n == 3 and qshort: q = q[:50]
        if n == 5 and bare: h = "@"
        out += [h, s, p, q]
    return "\n".join(out) + "\n"
open("/home/falk/audit6/io/t6/in/qshort_R1.fq", "w").write(fq(src, True, False))
open("/home/falk/audit6/io/t6/in/qshort_R2.fq", "w").write(fq(src2, False, False))
open("/home/falk/audit6/io/t6/in/bare_R1.fq", "w").write(fq(src, False, True))
open("/home/falk/audit6/io/t6/in/bare_R2.fq", "w").write(fq(src2, False, False))
# last record cut to 2 lines (plain file, no compression to detect it)
open("/home/falk/audit6/io/t6/in/cut_R1.fq", "w").write("\n".join(src[:398]) + "\n")
EOF
head -n 400 $A2 > $I/cut_R2.fq
head -n 400 $A1 > $I/full_R1.fq
asan_check() { grep -a -E "ERROR: AddressSanitizer|runtime error:" $1 | grep -v -c "IndexCodec.h:203"; }
run() { # run NAME ARGS...
  local name=$1; shift
  timeout 900 $P --db $DB --no_qcmsa --no_strains -t 2 -o out --prefix $name "$@" > $name.log 2>&1
  local rc=$? s=$(ls out/$name.sam* 2>/dev/null | grep -v -e err -e partial | head -1)
  printf '%-10s rc=%s sanitizer_reports=%s records=%s | %s\n' $name $rc "$(asan_check $name.log)" \
    "$([ -n "$s" ] && samtext $s | grep -vc '^@' || echo none)" \
    "$(grep -a -E 'ERROR: AddressSanitizer|runtime error|\[ERROR\]|malformed|unrecog|differ|skipped' $name.log | grep -a -v IndexCodec | head -3 | tr '\n' '|' | cut -c1-300)"
}
run fablank -1 $I/blank_R1.fa -2 $I/blank_R2.fa --sam_format sam
run facrlf -1 $I/crlf_R1.fa -2 $I/crlf_R2.fa --sam_format sam
run fablse -1 $I/blank_R1.fa --model_se $DB/model_pe.xml --sam_format sam
run qshort -1 $I/qshort_R1.fq -2 $I/qshort_R2.fq --sam_format sam
run bare -1 $I/bare_R1.fq -2 $I/bare_R2.fq --sam_format sam
run cutse -1 $I/cut_R1.fq --model_se $DB/model_pe.xml --sam_format sam
run fullse -1 $I/full_R1.fq --model_se $DB/model_pe.xml --sam_format sam
run cutpe -1 $I/cut_R1.fq -2 $I/cut_R2.fq --sam_format sam
echo "cutse vs fullse qnames: $(samtext out/cutse.sam | grep -v '^@' | cut -f1 | sort -u | wc -l) vs $(samtext out/fullse.sam | grep -v '^@' | cut -f1 | sort -u | wc -l)"
grep -a -v '^@' out/qshort.sam | awk -F'\t' 'length($10)!=length($11){c++} END{print "qshort: records with QUAL length != SEQ length:", c+0}'
for f in sam gz zst; do run pe_$f -1 $A1 -2 $A2 --sam_format $f; done
for f in sam gz zst; do
  s=$(ls out/pe_$f.sam* | grep -v err)
  timeout 900 $P --db $DB --no_qcmsa -t 2 -o po --prefix po_$f --profile_only $s > po_$f.log 2>&1
  echo "po_$f rc=$? sanitizer_reports=$(asan_check po_$f.log) profile $(cmp -s po/po_$f.profile $W/t1/sam_t1/sa.profile && echo SAME || echo diff)"
done
echo "--- unit tests under the sanitizers (I/O ones)"
cd $W/src/build-asan/tests
ASAN_OPTIONS=detect_leaks=1 timeout 1200 ./protal_tests --gtest_filter='ThreadedGzStream*:SamFile*:SeqReader*:GeneTables*:ReferenceMap*:UniqueKmers*:FromSam*:SamReader*:SamRoundTrip*' > $T/unit.log 2>&1
echo "unit rc=$? $(grep -a -E '^\[  (PASSED|FAILED)' $T/unit.log | tr '\n' ' ') sanitizer_reports=$(asan_check $T/unit.log)"
grep -a -E 'FAILED|ERROR: AddressSanitizer|runtime error' $T/unit.log | head -5
