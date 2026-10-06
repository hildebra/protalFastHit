#!/bin/bash
# Alternate ISA-L and zlib-ng runs so both see the same machine load.
cd ~/isal/bench
GZ=~/zstdbench/reads.gzip6.fq.gz
BGZF=~/zstdbench/reads.bgzf.fq.gz
out=~/isal/bench/inflate_results.tsv
: > $out
uptime
for r in 1 2 3 4 5; do
  ./bench_inflate isal $GZ 1 >> $out
  ./bench_inflate zng $GZ 1 >> $out
done
for r in 1 2 3 4 5; do
  ./bench_inflate isal $BGZF 1 | sed 's/^isal/isal_bgzf/' >> $out
  ./bench_inflate zng $BGZF 1 | sed 's/^zng/zng_bgzf/' >> $out
done
uptime
# min and median of wall and cpu MB/s per method
perl -e '
  while (<>) { my @f = split /\t/; my ($m) = $f[0]; my %h = map { split /=/, $_, 2 } @f[1..$#f];
    push @{$w{$m}}, $h{wall_MBps}; push @{$c{$m}}, $h{cpu_MBps}; push @order, $m unless $seen{$m}++; }
  printf "%-10s %5s %12s %12s %12s %12s\n", "method", "runs", "wall_max", "wall_median", "cpu_max", "cpu_median";
  for my $m (@order) { my @a = sort { $a <=> $b } @{$w{$m}}; my @b = sort { $a <=> $b } @{$c{$m}};
    printf "%-10s %5d %12.1f %12.1f %12.1f %12.1f\n", $m, scalar @a, $a[-1], $a[$#a/2], $b[-1], $b[$#b/2]; }
' $out
