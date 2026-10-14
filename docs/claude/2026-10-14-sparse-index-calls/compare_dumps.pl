#!/usr/bin/env perl
# Join two runs' training dumps (<profile>.truth_annotated) by sample and taxon: calls of each, the true species the
# second run loses or gains, and per model feature how the second run's values compare with the first's (true species
# with records in both; and the lost ones). usage: compare_dumps.pl <profiles A> <profiles B> <model.xml>
use strict; use warnings;
use List::Util qw(sum);
my ($da, $db, $model) = @ARGV;
open(my $mh, "<", $model) or die; my %feat;
while (<$mh>) { $feat{$1} = 1 while /MiningField name="([^"]+)"/g; } delete $feat{truth};
sub load { my ($dir) = @_; my %r;
    for my $f (glob("$dir/*.truth_annotated")) { (my $s = $f) =~ s{.*/}{}; $s =~ s/\.profile\.truth_annotated$//;
        open(my $fh, "<", $f) or die; chomp(my $h = <$fh>); my @h = split /\t/, $h;
        while (<$fh>) { chomp; my @v = split /\t/; my %x; @x{@h} = @v; $r{$s}{$x{taxon}} = \%x; } }
    return \%r; }
my ($A, $B) = (load($da), load($db));
my (%c, @true_both, @lost, @gained);
for my $s (sort keys %$A) {
    for my $t (sort keys %{$A->{$s}}) { my $x = $A->{$s}{$t}; my $y = $B->{$s}{$t};
        $c{"A_tp"}++ if $x->{truth} && $x->{prediction}; $c{"A_fp"}++ if !$x->{truth} && $x->{prediction};
        next unless $y;
        push @true_both, [$x, $y] if $x->{truth};
        push @lost, [$x, $y, $s] if $x->{truth} && $x->{prediction} && !$y->{prediction};
        push @gained, [$x, $y, $s] if $x->{truth} && !$x->{prediction} && $y->{prediction}; }
    for my $t (sort keys %{$B->{$s}}) { my $y = $B->{$s}{$t};
        $c{"B_tp"}++ if $y->{truth} && $y->{prediction}; $c{"B_fp"}++ if !$y->{truth} && $y->{prediction};
        $c{"B_only_true"}++ if $y->{truth} && !$A->{$s}{$t}; }
}
print "calls: ", join(" ", map { "$_=" . ($c{$_} // 0) } sort keys %c), "; true in both ", scalar(@true_both),
      ", lost ", scalar(@lost), ", gained ", scalar(@gained), "\n";
sub median { my @s = sort { $a <=> $b } @_; return @s ? $s[$#s / 2] : 0; }
my @rows;
for my $f (sort keys %feat) {
    my @pairs = grep { defined $_->[0]{$f} } @true_both; next unless @pairs;
    my $ma = sum(map { $_->[0]{$f} } @pairs) / @pairs; my $mb = sum(map { $_->[1]{$f} } @pairs) / @pairs;
    my @ratio = map { $_->[1]{$f} / $_->[0]{$f} } grep { $_->[0]{$f} != 0 } @pairs;
    my $changed = grep { abs($_->[1]{$f} - $_->[0]{$f}) > 1e-9 * (1 + abs($_->[0]{$f})) } @pairs;
    my ($la, $lb) = @lost ? (sum(map { $_->[0]{$f} } @lost) / @lost, sum(map { $_->[1]{$f} } @lost) / @lost) : (0, 0);
    push @rows, [$f, $changed / @pairs, $ma, $mb, median(@ratio), $la, $lb];
}
printf "%-34s %8s %10s %10s %8s | %10s %10s\n", "feature (true species)", "changed", "mean A", "mean B", "med B/A", "lost: A", "lost: B";
printf "%-34s %8.3f %10.4g %10.4g %8.3f | %10.4g %10.4g\n", @$_ for sort { $b->[1] <=> $a->[1] } grep { $_->[1] > 0.01 } @rows;
printf "probability (true, both): mean A %.4f, B %.4f; lost: A %.4f, B %.4f\n",
    sum(map { $_->[0]{probability} } @true_both) / @true_both, sum(map { $_->[1]{probability} } @true_both) / @true_both,
    @lost ? (sum(map { $_->[0]{probability} } @lost) / @lost, sum(map { $_->[1]{probability} } @lost) / @lost) : (0, 0);
