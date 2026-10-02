#!/usr/bin/env perl
# Species-level F1 as lr_score.pl computes it, on the v0.7.2 benchmark's own long-read runs (0.7.2 binary, V072).
use strict; use warnings;
my $B = "$ENV{HOME}/bench071";
for my $scen (qw(ont_b3000000 pb_b3000000 ont_b90000000 pb_b90000000)) {
  my %truth_of; open my $t, '<', "$B/samples_lr072/points/$scen/sim/samples.tsv" or die; <$t>;
  while (<$t>) { chomp; my @f = split /\t/; $truth_of{$f[0]} = $f[2] }
  my @f1;
  for my $n (1 .. 4) { my $s = "${scen}_s_$n"; my $t = $scen =~ /^pb/ ? "pb" : "ont"; my $p = "$B/runs_lr072/v072.full.$t.$s/$s.profile"; open my $f, '<', $p or next;
    my %sp; while (<$f>) { chomp; my @c = split /\t/; $sp{$c[1]} = 1 if @c >= 3 }
    my %tr; open my $tf, '<', $truth_of{$s} or die; while (<$tf>) { chomp; $tr{$_} = 1 if length }
    my $tp = grep { $tr{$_} } keys %sp; my $fp = keys(%sp) - $tp; my $fn = keys(%tr) - $tp; push @f1, 2 * $tp / (2 * $tp + $fp + $fn) }
  my $m = 0; $m += $_ for @f1; printf "%-14s 0.7.2's runs, F1 as lr_score.pl: %.4f (%d samples)\n", $scen, @f1 ? $m / @f1 : 0, scalar @f1;
}
