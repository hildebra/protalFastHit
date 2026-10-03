#!/usr/bin/env bash
# Experiment tree perf4/tim: HEAD plus wall-clock prints around the profiling passes (stderr "[tim] name seconds").
set -uo pipefail
W=$HOME/mt-work/perf4; T=$W/tim
rm -rf $T; mkdir -p $T; (cd $W/ref && tar -c --exclude=./build .) | tar -x -C $T
cd $T
P=src/Profiling/Profiler.h; R=src/RunProtal.h
perl -0pi -e 's|#include "Haplotypes.h"|#include "Haplotypes.h"\n#include <chrono>\n#define TIM_NOW() std::chrono::steady_clock::now()\n#define TIM_PRINT(name, t0) do { std::cerr << "[tim] " << name << " " << std::chrono::duration<double>(TIM_NOW() - (t0)).count() << std::endl; } while (0)\n|' $P
perl -0pi -e 's|(void ApplyRecordEvidence\(size_t threads = 1\) \{\n\s*)FinishLink\(\);|${1}auto const tim_fl = TIM_NOW(); FinishLink(); TIM_PRINT("finish_link", tim_fl);|' $P
perl -0pi -e 's|(\n\s*)ApplySampleContext\(threads\);\n|${1}auto const tim_sc = TIM_NOW(); ApplySampleContext(threads); TIM_PRINT("sample_context_total", tim_sc);\n|' $P
perl -0pi -e 's|(\n\s*)auto const shares = context::AbundanceWeightedShares\(m_evidence.Ambiguity\(\), m_evidence.RecordCountsByTaxon\(\)\);|${1}auto const tim_em = TIM_NOW(); auto const shares = context::AbundanceWeightedShares(m_evidence.Ambiguity(), m_evidence.RecordCountsByTaxon()); TIM_PRINT("em_shares", tim_em); std::cerr << "[tim] em_classes " << m_evidence.Ambiguity().size() << std::endl;|' $P
perl -0pi -e 's|(\n\s*)sam_chunks::ParallelFor\(pairs.size\(\), std::max<size_t>\(threads, 1\), \[&\]\(size_t p\) \{\n(\s*)pair_distance\[p\] = m_distances->Between\(pairs\[p\].first, pairs\[p\].second\);\n(\s*)\}\);|${1}auto const tim_d = TIM_NOW(); sam_chunks::ParallelFor(pairs.size(), std::max<size_t>(threads, 1), [&](size_t p) {\n${2}pair_distance[p] = m_distances->Between(pairs[p].first, pairs[p].second);\n${3}}); TIM_PRINT("congener_distances", tim_d); std::cerr << "[tim] congener_pairs " << pairs.size() << std::endl;|' $P
perl -0pi -e 's|(\n\s*)profile.ApplyRecordEvidence\(threads\);\n(\s*)m_post_process_bm.Start\(\);\n(\s*)profile.PostProcessSNPs\(snp_min_cov, snp_min_obs_fwdrev, snp_min_af, snp_min_mean_qual, snp_min_phred_sum, snp_require_strand, threads\);\n(\s*)m_post_process_bm.Stop\(\);|${1}auto const tim_ev = TIM_NOW(); profile.ApplyRecordEvidence(threads); TIM_PRINT("record_evidence_total", tim_ev);\n${2}m_post_process_bm.Start(); auto const tim_pp = TIM_NOW();\n${3}profile.PostProcessSNPs(snp_min_cov, snp_min_obs_fwdrev, snp_min_af, snp_min_mean_qual, snp_min_phred_sum, snp_require_strand, threads);\n${4}m_post_process_bm.Stop(); TIM_PRINT("post_process_snps", tim_pp);|' $P
perl -0pi -e 's|(\n\s*)auto const error = ProfileSamParallel\(file_path, profile, rejected, threads, serial\);|${1}auto const tim_ps = TIM_NOW(); auto const error = ProfileSamParallel(file_path, profile, rejected, threads, serial); TIM_PRINT("parse_and_add_parallel", tim_ps);|' $P
perl -0pi -e 's|(\n\s*)m_reads = read_id;\n(\s*)profile.ApplyRecordEvidence\(threads\);|${1}m_reads = read_id; TIM_PRINT("parse_and_add_serial", tim_ser);\n${2}auto const tim_ev2 = TIM_NOW(); profile.ApplyRecordEvidence(threads); TIM_PRINT("record_evidence_total", tim_ev2);|' $P
perl -0pi -e 's|(\n\s*)size_t read_id = 0;\n(\s*)// One link per read|${1}size_t read_id = 0; auto const tim_ser = TIM_NOW();\n${2}// One link per read|' $P
perl -0pi -e 's|(\n\s*)if \(threads_per_sample > 1\) profile.ScoreTaxa\(filter, threads_per_sample\);|${1}{ auto const tim_st = std::chrono::steady_clock::now(); if (threads_per_sample > 1) profile.ScoreTaxa(filter, threads_per_sample); std::cerr << "[tim] score_taxa " << std::chrono::duration<double>(std::chrono::steady_clock::now() - tim_st).count() << std::endl; }|' $R
perl -0pi -e 's|(\n\s*)profile.WriteSparseProfile\(taxonomy, filter, os, &os_total, &os_dismissed, threads_per_sample\);\n(\s*)profile.WriteGeneProfile\(taxonomy, filter, &os_genes, threads_per_sample\);|${1}{ auto const tim_w = std::chrono::steady_clock::now(); profile.WriteSparseProfile(taxonomy, filter, os, &os_total, &os_dismissed, threads_per_sample); std::cerr << "[tim] write_sparse_profile " << std::chrono::duration<double>(std::chrono::steady_clock::now() - tim_w).count() << std::endl; }\n${2}{ auto const tim_g = std::chrono::steady_clock::now(); profile.WriteGeneProfile(taxonomy, filter, &os_genes, threads_per_sample); std::cerr << "[tim] write_gene_profile " << std::chrono::duration<double>(std::chrono::steady_clock::now() - tim_g).count() << std::endl; }|' $R
grep -c 'TIM_PRINT\|\[tim\]' $P $R
grep -q '#include <chrono>' $R || perl -0pi -e 's|#include "RunStatus.h"|#include "RunStatus.h"\n#include <chrono>|' $R
nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=OFF > $W/tim.configure.log 2>&1 && nice cmake --build build --target protal -j 5 > $W/tim.build.log 2>&1 && echo "OK tim build" || { echo "FAIL tim build"; grep -E 'error' -A3 $W/tim.build.log | head -40; exit 1; }
B=$T/build/protal; DB=$HOME/bench071/V073/protal_db
while ! grep -q "Run protal took" $W/o.map4.log 2>/dev/null; do sleep 5; done
for spec in "5M 6 $W/o.pe5M_t6/s.sam.zst" "5M 1 $W/o.pe5M_t6/s.sam.zst" "500k 1 $W/o.pe500k_t1/s.sam.zst" "500k 6 $W/o.pe500k_t1/s.sam.zst" "ont90M 6 $W/o.ont90M_t6/s.sam.zst"; do
  set -- $spec; n=$1; t=$2; sam=$3; rm -rf $W/o.tim_${n}_t$t
  $B --db $DB --profile_only $sam --prefix s -o $W/o.tim_${n}_t$t -t $t --no_qcmsa --verbose > $W/o.tim_${n}_t$t.log 2> $W/o.tim_${n}_t$t.err
  echo "== tim $n t=$t"; grep -E '^\[tim\]' $W/o.tim_${n}_t$t.err; grep -E "^(Thread 0 Profile sample|Profiling|Run protal) took" $W/o.tim_${n}_t$t.log
done
