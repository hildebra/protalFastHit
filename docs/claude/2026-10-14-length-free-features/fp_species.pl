#!/usr/bin/env perl
# The false positives of strain-world runs by species: per run, how many shallow and deep samples call each absent
# species. usage: fp_species.pl <run dir> ...   (each with profiles/*.profile.truth_annotated)
use strict;
use warnings;

my (%count, %runs);
for my $run (@ARGV) {
    (my $name = $run) =~ s{.*/run_}{};
    $runs{$name} = 1;
    for my $f (glob("$run/profiles/*.truth_annotated")) {
        (my $sample = $f) =~ s{.*/}{};
        my $depth = $sample =~ /^deep/ ? "deep" : "shallow";
        open(my $fh, "<", $f) or die "$f: $!\n";
        chomp(my $h = <$fh>);
        my @h = split /\t/, $h;
        while (<$fh>) {
            chomp;
            my %x;
            @x{@h} = split /\t/;
            next unless !$x{truth} && $x{prediction};
            $count{$x{taxon_name}}{$name}{$depth}++;
        }
    }
}
my @names = map { (my $n = $_) =~ s{.*/run_}{}; $n } @ARGV;
printf "%-36s %s\n", "species (shallow/deep samples)", join(" ", map { sprintf "%-14s", substr($_, 2) } @names);
for my $sp (sort { total($b) <=> total($a) || $a cmp $b } keys %count) {
    printf "%-36s %s\n", substr($sp, 0, 36),
        join(" ", map { my $c = $count{$sp}{$_} // {}; sprintf "%-14s", ($c->{shallow} // 0) . "/" . ($c->{deep} // 0) } @names);
}
printf "%-36s %s\n", "total",
    join(" ", map { my $n = $_; my ($s, $d) = (0, 0); for my $sp (keys %count) { $s += $count{$sp}{$n}{shallow} // 0; $d += $count{$sp}{$n}{deep} // 0 } sprintf "%-14s", "$s/$d" } @names);

sub total { my $sp = shift; my $t = 0; for my $n (keys %{$count{$sp}}) { $t += $_ for values %{$count{$sp}{$n}} } $t }
