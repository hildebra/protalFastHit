#!/usr/bin/env bash
# One row per run: wall, alignment and profiling stage times, and the MTAUDIT counters.
cd ${MT_AUDIT:-$HOME/mt-audit}/runs
printf "run\twall_s\tcpu_pct\talign_s\tprofile_s\treader_s_sum\toutput_s_sum"
for k in reader_wait reader_hold underflow_wait inflate_busy inflate_idle sam_pack sam_wait sam_hold prof_total prof_group prof_addsam prof_post; do printf "\t$k"; done
printf "\treader_n\tsam_n\tinflated_MB\tsam_MB\tload\n"
for t in *.time; do
  n=${t%.time}
  wall=$(grep Elapsed $t | awk '{print $NF}'); cpu=$(grep 'Percent of CPU' $t | awk '{print $NF}')
  al=$(grep 'Aligning reads took' $n.log | sed 's/.*took //'); pr=$(grep 'Profiling took' $n.log | sed 's/.*took //')
  rd=$(awk -F'\t' '$1=="Sequence reader"{print $2}' $n/misc/*_runtime.tsv 2>/dev/null)
  ou=$(awk -F'\t' '$1=="Output handler"{print $2}' $n/misc/*_runtime.tsv 2>/dev/null)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s" $n "$wall" "$cpu" "${al// /}" "${pr// /}" "$rd" "$ou"
  for k in reader_wait reader_hold underflow_wait inflate_busy inflate_idle sam_pack sam_wait sam_hold prof_total prof_group prof_addsam prof_post; do
    v=$(awk -F'\t' -v k=$k '$1=="MTAUDIT" && $2==k{print $3}' $n.err); printf "\t%s" "${v:--}"
  done
  rn=$(awk -F'\t' '$1=="MTAUDIT" && $2=="reader_hold"{print $4}' $n.err); sn=$(awk -F'\t' '$1=="MTAUDIT" && $2=="sam_hold"{print $4}' $n.err)
  ib=$(awk -F'\t' '$1=="MTAUDIT" && $2=="inflated_bytes"{printf "%.0f", $3/1e6}' $n.err); sb=$(awk -F'\t' '$1=="MTAUDIT" && $2=="sam_bytes"{printf "%.0f", $3/1e6}' $n.err)
  printf "\t%s\t%s\t%s\t%s\t%s\n" "${rn:--}" "${sn:--}" "${ib:--}" "${sb:--}" "$(cut -d' ' -f1 $n.load)"
done
