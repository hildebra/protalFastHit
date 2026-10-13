#!/bin/sh
# Fetch each run of runs_in.tsv from the ENA portal API and write runs.tsv: the curated columns
# plus ENA's study, instrument, layout, read count, bases (Gb) and sample title. A run ENA does not
# return, or returns with 0 bases (raw signal only, or a deposit ENA never processed), is marked in
# the "check" column. Needs curl and awk; ENA answers in about a minute.
#   sh fetch_runs.sh [runs_in.tsv] [runs.tsv]
set -eu
in=${1:-runs_in.tsv}
out=${2:-runs.tsv}
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT

ids=$(awk -F'\t' 'NR > 1 { print $5 }' "$in" | paste -sd, -)
curl -sS --max-time 300 -G "https://www.ebi.ac.uk/ena/portal/api/search" \
  --data-urlencode "result=read_run" \
  --data-urlencode "includeAccessions=$ids" \
  --data-urlencode "fields=run_accession,study_accession,instrument_model,library_layout,library_strategy,read_count,base_count,sample_title" \
  --data-urlencode "format=tsv" > "$tmp"

awk -F'\t' -v OFS='\t' '
  FNR == NR { if (FNR > 1) ena[$1] = $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS sprintf("%.2f", $7 / 1e9) OFS $8; next }
  FNR == 1 { print $0, "study", "instrument", "layout", "strategy", "reads", "gbases", "sample_title", "check"; next }
  {
    if (!($5 in ena)) { print $0, "-", "-", "-", "-", "-", "-", "-", "not in ENA"; next }
    split(ena[$5], f, OFS)
    print $0, ena[$5], (f[6] + 0 > 0 ? "ok" : "0 bases in ENA")
  }' "$tmp" "$in" > "$out"

awk -F'\t' 'NR > 1 && $NF != "ok" { print "check: " $5 " " $NF }' "$out" >&2
echo "$(($(wc -l < "$out") - 1)) runs written to $out" >&2
