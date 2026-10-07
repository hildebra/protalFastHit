#!/usr/bin/env perl
# Estimate the error SAM bytes per kind (FP, FN, unseen) under a cap of N fragments per reason.
# usage: zstd -dc S.sam.zst | caps.pl S.taxa.tsv N1,N2,...
use strict; use warnings;
my ($taxa, $caps) = @ARGV;
open my $t, '<', $taxa or die; my $h = <$t>; my %kind;
while (<$t>) { my @f = split /\t/; $kind{$f[3]} = $f[2]; }
my (%bytes, %reasons, @order);
while (my $line = <STDIN>) {
  next if $line =~ /^@/;
  my ($q) = $line =~ /^(\S+)/;
  my ($xe) = $line =~ /\txe:Z:(\S+)/;
  push @order, $q unless exists $bytes{$q};
  $bytes{$q} += length $line;
  $reasons{$q}{$_} = 1 for split /,/, $xe;
}
my %by;
for my $q (@order) { push @{$by{$_}}, $q for keys %{$reasons{$q}}; }
sub kind_of { my ($r) = @_; my ($w, $t) = split /:/, $r; return $w eq 'FP' ? 'FP' : $w eq 'FN' ? 'FN' : ($kind{$t} // 'unseen') eq 'FN' ? 'FN' : 'unseen'; }
for my $n (0, split /,/, $caps) {
  my %kept;
  for my $r (keys %by) { my @q = @{$by{$r}}; @q = @q[0 .. $n - 1] if $n && @q > $n; $kept{kind_of($r)}{$_} = 1 for @q; }
  printf "cap %4s:", $n || 'all';
  for my $k (qw(FP FN unseen)) { my $b = 0; $b += $bytes{$_} for keys %{$kept{$k} // {}}; printf "  %s %6d fragments %7.1f MB", $k, scalar(keys %{$kept{$k} // {}}), $b / 1e6; }
  print "\n";
}
