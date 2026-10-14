#!/usr/bin/env perl
# The strain world's runs compared at the same false positives: the share of test strains (evaluate.py's: truth.tsv rows
# of kind "strain") whose species a run calls when its knob lets through at most K false positives over the samples of
# a depth (shallow or deep), for several K; by the profiles' probabilities (*.profile.log), without the singleton rule.
# usage: matched_fp.pl <strain world dir> <run> ...   (run: the name after run_, e.g. v_l150_mix_new)
use strict;
use warnings;

my ($world, @runs) = @ARGV;
my %budget = (shallow => [2, 5, 10, 20], deep => [5, 10, 20, 40]);
my (%members, %kind_of);
open(my $tf, "<", "$world/reads/truth.tsv") or die "$world/reads/truth.tsv: $!\n";
<$tf>;
while (<$tf>) {
    chomp;
    my @r = split /\t/;    # sample, kind, accession, species, role, pairs
    push @{$members{$r[0]}}, \@r;
    $kind_of{$r[0]} = $r[1];
}
close $tf;

printf "%-24s %s\n", "run", join(" | ", map { my $k = $_; join(" ", map { sprintf "%s<=%-3d", substr($k, 0, 1), $_ } @{$budget{$k}}) } qw(shallow deep));
for my $run (@runs) {
    my @cells;
    for my $kind (qw(shallow deep)) {
        my (@strain_p, @fp_p);
        for my $sample (sort grep { $kind_of{$_} eq $kind } keys %members) {
            my %p;
            open(my $fh, "<", "$world/run_$run/profiles/$sample.profile.log") or die "$run $sample: $!\n";
            chomp(my $h = <$fh>);
            my @h = split /\t/, $h;
            my %ci = map { $h[$_] => $_ } 0 .. $#h;
            while (<$fh>) {
                chomp;
                my @r = split /\t/;
                next unless defined $r[$ci{Name}];
                $p{$r[$ci{Name}]} = $r[$ci{Probability}] || 0;
            }
            close $fh;
            my %present = map { $_->[3] => 1 } @{$members{$sample}};
            push @strain_p, map { $p{$_->[3]} // 0 } grep { $_->[4] eq "strain" } @{$members{$sample}};
            push @fp_p, map { $p{$_} } grep { !$present{$_} } keys %p;
        }
        my @fp = sort { $b <=> $a } @fp_p;
        for my $k (@{$budget{$kind}}) {
            # The knob just above the (k+1)-th highest false positive: at most k false positives called.
            my $knob = $k < @fp ? $fp[$k] : -1;
            my $called = grep { $_ > $knob } @strain_p;
            push @cells, sprintf "%.3f", $called / @strain_p;
        }
        push @cells, "|";
    }
    pop @cells;
    printf "%-24s %s\n", $run, join(" ", @cells);
}
