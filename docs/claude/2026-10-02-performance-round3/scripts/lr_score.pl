#!/usr/bin/env perl
# lr_score.pl: per scenario of lr_bench.sh, species-level F1 of each build against the truth, how far the two builds'
# profiles are from each other (Bray-Curtis over genomes; identical profiles counted), and the alignment times.
use strict; use warnings;
my $B = "$ENV{HOME}/bench071"; my $W = "$ENV{HOME}/mt-work/perf3/lrbench";
my %truth_of;
for my $scen (qw(ont_b3000000 pb_b3000000 ont_b90000000 pb_b90000000)) {
  open my $t, '<', "$B/samples_lr072/points/$scen/sim/samples.tsv" or die; <$t>;
  while (<$t>) { chomp; my @f = split /\t/; $truth_of{$f[0]} = $f[2] } }
sub profile { my $p = shift; my %g; open my $f, '<', $p or return undef; while (<$f>) { chomp; my @c = split /\t/; next if @c < 3; $g{$c[0]} = [$c[1], $c[2]] } return \%g }
sub species { my $g = shift; my %s; $s{$_->[0]} = 1 for values %$g; return \%s }
my %times; open my $tt, '<', "$W/times.tsv" or die; <$tt>; while (<$tt>) { chomp; my @f = split /\t/; push @{$times{"$f[0]|$f[1]"}}, [@f[2..4]] }
printf "%-14s %8s %8s %8s %6s %10s %10s %10s %10s\n", 'scenario', 'F1 ref', 'F1 lr', 'BC r-l', 'same', 'align ref', 'align lr', 'user ref', 'user lr';
for my $scen (qw(ont_b3000000 pb_b3000000 ont_b90000000 pb_b90000000)) {
  my (@f1r, @f1l, @bc, $same, @ar, @al, @ur, @ul); $same = 0;
  for my $n (1 .. 4) {
    my $s = "${scen}_s_$n"; my $r = profile("$W/ref.$s/s.profile"); my $l = profile("$W/lr.$s/s.profile"); next unless $r && $l;
    my %t; open my $tf, '<', $truth_of{$s} or die "$truth_of{$s}"; while (<$tf>) { chomp; $t{$_} = 1 if length }
    for my $pair ([$r, \@f1r], [$l, \@f1l]) { my $sp = species($pair->[0]); my $tp = grep { $t{$_} } keys %$sp; my $fp = keys(%$sp) - $tp; my $fn = keys(%t) - $tp;
      push @{$pair->[1]}, 2 * $tp / (2 * $tp + $fp + $fn) }
    my %all = (%$r, %$l); my ($num, $den) = (0, 0);
    for my $g (keys %all) { my $a = $r->{$g} ? $r->{$g}[1] : 0; my $b = $l->{$g} ? $l->{$g}[1] : 0; $num += abs($a - $b); $den += $a + $b }
    push @bc, $den ? $num / $den : 0; $same++ if $num == 0;
    push @ar, $times{"ref|$s"}[0][2]; push @al, $times{"lr|$s"}[0][2]; push @ur, $times{"ref|$s"}[0][1]; push @ul, $times{"lr|$s"}[0][1];
  }
  my $mean = sub { my @x = grep { defined } @_; return @x ? (eval { my $s = 0; $s += $_ for @x; $s }) / @x : 0 };
  printf "%-14s %8.4f %8.4f %8.4f %4d/%d %10.1f %10.1f %10.1f %10.1f\n", $scen, $mean->(@f1r), $mean->(@f1l), $mean->(@bc), $same, scalar(@bc), $mean->(@ar), $mean->(@al), $mean->(@ur), $mean->(@ul);
}
