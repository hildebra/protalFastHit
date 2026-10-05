#!/usr/bin/env bash
# Runs build <name> (~/det-order/<name>) N times (default 4) on the v0.7.3 benchmark world at 6 threads: 500k simulated
# read pairs (pe) and 90 Mb of simulated PacBio HiFi (pb). Outputs in ~/det-order/runs/<name>/{pe,pb}<i>.
#   runs.sh <name> [N] [extra protal options]
set -uo pipefail
NAME=$1; N=${2:-4}; shift; shift || true
W=$HOME/det-order; BIN=$W/$NAME/tree/build/protal; O=$W/runs/$NAME; mkdir -p $O
DB=$HOME/bench071/V073/protal_db/database.protal
P=$HOME/bench071/samples/points/rl150_p500000/sim/reads
R1=$P/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000_s_1_R2.fq.gz
PB=$HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz
for i in $(seq 1 $N); do
  rm -rf $O/pe$i $O/pb$i
  $BIN --db $DB -1 $R1 -2 $R2 --read_type pe --prefix s -o $O/pe$i -t 6 --no_qcmsa "$@" > $O/pe$i.log 2>&1 || echo "pe$i failed"
  $BIN --db $DB -1 $PB --read_type pb --prefix s -o $O/pb$i -t 6 --no_qcmsa "$@" > $O/pb$i.log 2>&1 || echo "pb$i failed"
done
echo "done $NAME"
