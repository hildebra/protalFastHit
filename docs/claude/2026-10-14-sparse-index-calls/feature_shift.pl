#!/usr/bin/env perl
# How far each feature of the training dumps moves between two runs of the same reads (another read length, another
# seeding), against how far it separates present from absent species. Rows joined by sample and taxon.
# usage: feature_shift.pl <profiles A> <profiles B> [min |shift|]
# shift = (mean in B - mean in A) / (mean of present - mean of absent, in A), over the present (true) and the absent
# species with a row in both; ratio = median B/A of the present ones (rows with A != 0).
use strict; use warnings;
use List::Util qw(sum max);
my ($da, $db, $min) = @ARGV; $min //= 0.1;
my %skip = map { $_ => 1 } qw(truth truth_raw prediction probability taxon taxon_name dataset);
my @cols;
sub load { my ($dir) = @_; my %r;
    for my $f (glob("$dir/*.truth_annotated")) { (my $s = $f) =~ s{.*/}{}; $s =~ s/\.profile\.truth_annotated$//;
        open(my $fh, "<", $f) or die; chomp(my $h = <$fh>); my @h = split /\t/, $h; @cols = grep { !$skip{$_} } @h;
        while (<$fh>) { chomp; my @v = split /\t/; my %x; @x{@h} = @v; $r{$s}{$x{taxon}} = \%x; } }
    return \%r; }
my ($A, $B) = (load($da), load($db));
my (@t, @f);
for my $s (keys %$A) { for my $k (keys %{$A->{$s}}) { my $y = $B->{$s}{$k} or next; my $x = $A->{$s}{$k};
    push @{$x->{truth} ? \@t : \@f}, [$x, $y]; } }
sub mean { my ($rows, $i, $c) = @_; return @$rows ? sum(map { $_->[$i]{$c} } @$rows) / @$rows : 0; }
sub median { my @s = sort { $a <=> $b } @_; return @s ? $s[$#s / 2] : 0; }
printf "rows in both: %d present, %d absent\n", scalar @t, scalar @f;
my @out;
for my $c (@cols) {
    next if grep { !defined $_->[0]{$c} || $_->[0]{$c} !~ /^-?[\d.eE+-]+$/ } @t[0 .. 0];
    my ($at, $bt, $af, $bf) = (mean(\@t, 0, $c), mean(\@t, 1, $c), mean(\@f, 0, $c), mean(\@f, 1, $c));
    my $sep = $at - $af;
    my ($st, $sf) = abs($sep) > 1e-12 ? (($bt - $at) / abs($sep), ($bf - $af) / abs($sep)) : (0, 0);
    my $ratio = median(map { $_->[1]{$c} / $_->[0]{$c} } grep { $_->[0]{$c} != 0 } @t);
    push @out, [$c, $at, $bt, $af, $bf, $st, $sf, $ratio];
}
printf "%-36s %10s %10s %10s %10s %8s %8s %8s\n", "feature", "present A", "present B", "absent A", "absent B", "shift p", "shift a", "ratio p";
for my $o (sort { max(abs($b->[5]), abs($b->[6])) <=> max(abs($a->[5]), abs($a->[6])) } @out) {
    my $ratio_moved = $o->[7] > 0 && abs(log($o->[7])) >= 0.1;
    next if max(abs($o->[5]), abs($o->[6])) < $min && !$ratio_moved;
    printf "%-36s %10.4g %10.4g %10.4g %10.4g %8.2f %8.2f %8.3f\n", @$o;
}
