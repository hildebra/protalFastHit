#!/usr/bin/env bash
# Each arm's records (ss_compare.sh's WORK/<arm>.pairs: read, reference) scored against the read's true source
# (docs/claude/2026-10-13-sparse-seeding): simulate_metagenomes names a long read g<k>x_<n>, k the genome's place among
# the sample's rows of the simulation's genomes table (sf_hifi_sim.sh's OUT/tmp/genomes.tsv); the genome's lineage comes
# from the paired-end manifest the sample replays, a record's from its reference's taxid (<taxid>_<gene id>) through the
# database's genome2tiid.tsv. Per arm: the records (not secondary), those of the read's genome, of its species and of
# its genus; the reads with a record, and those with one of their species.
#   bash ss_read_truth.sh WORK GENOMES_TSV SAMPLE MANIFEST GENOME2TIID ARMS...
set -uo pipefail
work=$1 genomes=$2 sample=$3 manifest=$4 g2t=$5; shift 5
printf 'arm\trecords\tgenome\tspecies\tgenus\treads\treads_species\n'
for arm in "$@"; do
  awk -F'\t' -v s="$sample" -v a="$arm" '
    FILENAME == ARGV[1] { if (FNR > 1 && $1 == s) place[n++] = $2; next }
    FILENAME == ARGV[2] { if (FNR == 1) { for (i = 1; i <= NF; i++) col[$i] = i; next } truth[$col["genome"]] = $col["taxonomy"]; next }
    FILENAME == ARGV[3] { lineage[$2] = $4; tiid_of[$1] = $2; next }
    function genus(l) { sub(/;s__.*/, "", l); return l }
    { k = $1; sub(/^g/, "", k); sub(/x_.*/, "", k); g = place[k + 0]; t = $2; sub(/_.*/, "", t)
      records++; reads[$1]
      if ((g in tiid_of) && tiid_of[g] == t) same_genome++
      if (lineage[t] == truth[g]) { same_species++; reads_species[$1] }
      if (genus(lineage[t]) == genus(truth[g])) same_genus++ }
    END { nr = 0; ns = 0; for (r in reads) nr++; for (r in reads_species) ns++
          printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", a, records, same_genome, same_species, same_genus, nr, ns }' \
      "$genomes" "$manifest" "$g2t" "$work/$arm.pairs"
done
