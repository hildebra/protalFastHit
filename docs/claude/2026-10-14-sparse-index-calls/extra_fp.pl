#!/usr/bin/env perl
# The false positives of run B that run A did not call (same sample and taxon), with some of their features in both runs.
# usage: extra_fp.pl <profiles A> <profiles B> [feature,...]
use strict; use warnings;
my ($da, $db, $list) = @ARGV;
my @show = split /,/, ($list // "fragments,identity,top_identity,mean_mapq,uniqueness,failed_candidate_rate,congener_fit_share,lu_per_kb,split_fragment_share");
sub load { my ($dir) = @_; my %r;
    for my $f (glob("$dir/*.truth_annotated")) { (my $s = $f) =~ s{.*/}{}; $s =~ s/\.profile\.truth_annotated$//;
        open(my $fh, "<", $f) or die; chomp(my $h = <$fh>); my @h = split /\t/, $h;
        while (<$fh>) { chomp; my @v = split /\t/; my %x; @x{@h} = @v; $r{$s}{$x{taxon}} = \%x; } }
    return \%r; }
my ($A, $B) = (load($da), load($db));
printf "%-10s %-34s %s\n", "sample", "taxon", join(" ", map { sprintf "%14s", $_ } "probability", @show);
for my $s (sort keys %$B) { for my $t (sort keys %{$B->{$s}}) { my $y = $B->{$s}{$t};
    next unless !$y->{truth} && $y->{prediction}; my $x = $A->{$s}{$t};
    next if $x && $x->{prediction};
    printf "%-10s %-34s %s\n", $s, substr($y->{taxon_name}, 0, 34),
        join(" ", map { my $c = $_; sprintf "%14s", ($x ? sprintf("%.3g", $x->{$c}) : "-") . "/" . sprintf("%.3g", $y->{$c}) } "probability", @show);
} }
