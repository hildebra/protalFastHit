#!/bin/bash
# Real-sized copies of a synthetic world's genomes, for timing the simulators.
#
# The synthetic worlds (simulate_gtdb_release.py) have genomes of ~270 kb in 3 contigs; GTDB's are ~3.4 Mb
# (median) in 1 to hundreds of contigs. ART, pbsim3, the genome copies and the length scan cost per base and
# per contig, so they are timed on copies padded to GTDB-like sizes: each genome keeps its own contigs (and so
# its marker genes) and gains random contigs (pseudo-random bytes from AES-CTR with the accession as key, so the
# copies are the same on every run) up to a size drawn per genome, lognormal with median 3.4 Mb and sigma 0.35
# (1.5-8 Mb), in a number of contigs drawn lognormal with median 40 and sigma 1.2 (3-600). The padding is
# outside every marker gene, as most of a real genome is: its reads do not align.
#
# usage: pad_genomes.sh GENOME_TABLE OUT_DIR [JOBS]   writes OUT_DIR/genomes/*.fna.gz and OUT_DIR/genomes.tsv
set -euo pipefail
table=$1 out=$2 jobs=${3:-4}
mkdir -p "$out/genomes"
acgt=$(printf 'ACGT%.0s' $(seq 64))
export out acgt

# The plan: accession, path, target bases, contigs (awk's srand: the same plan on every run).
awk -F'\t' -v seed=7 'BEGIN { srand(seed) }
  function gauss() { return sqrt(-2 * log(1 - rand())) * cos(6.283185307 * rand()) }
  NF >= 3 && $1 !~ /^#/ {
    size = int(3.4e6 * exp(0.35 * gauss())); if (size < 1.5e6) size = 1.5e6; if (size > 8e6) size = 8e6
    contigs = int(40 * exp(1.2 * gauss()) + 0.5); if (contigs < 3) contigs = 3; if (contigs > 600) contigs = 600
    print $1 "\t" $3 "\t" size "\t" contigs }' "$table" > "$out/plan.tsv"

pad_one() {
  local acc=$1 path=$2 size=$3 contigs=$4
  local dest="$out/genomes/$(basename "$path")"
  [ -s "$dest" ] && return 0
  local own
  own=$(zcat "$path" | grep -v '^>' | tr -d '\n\r' | wc -c)
  local pad=$(( size > own ? size - own : 0 ))
  local extra=$(( contigs > 3 ? contigs - 3 : 1 ))
  local per=$(( pad / extra + 1 ))
  local lines=$(( (per + 79) / 80 ))
  {
    zcat "$path"
    openssl enc -aes-128-ctr -nosalt -pass "pass:$acc" -pbkdf2 < /dev/zero 2>/dev/null | head -c "$pad" | tr '\000-\377' "$acgt" |
      fold -w 80 | awk -v lines="$lines" -v acc="$acc" '(NR - 1) % lines == 0 { print ">" acc "_pad" int((NR - 1) / lines) + 1 } { print }'
  } | gzip -4 > "$dest.partial"
  mv "$dest.partial" "$dest"
}
export -f pad_one
tr '\t' '\n' < "$out/plan.tsv" | xargs -d '\n' -n 4 -P "$jobs" bash -c 'pad_one "$@"' _

# The genome table of the copies: the same rows, the FASTA path replaced.
awk -F'\t' -v OFS='\t' -v dir="$out/genomes" '$1 !~ /^#/ && NF >= 3 { n = split($3, p, "/"); $3 = dir "/" p[n] } { print }' \
  "$table" > "$out/genomes.tsv"
echo "padded $(wc -l < "$out/plan.tsv") genomes into $out/genomes ($(du -sh "$out/genomes" | cut -f1))"
