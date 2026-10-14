#!/usr/bin/env perl
# The seed tags of two runs' SAMs (alignments/*.sam.zst) record by record: for the records of the same alignment (read,
# flag, gene, position, CIGAR) in both, how many have the same ZS, ZU and ZT, and each run's mean ZS, ZU and ZT and
# ZU / ZS over its own records. usage: uc_tags.pl <run dir A> <run dir B>
use strict;
use warnings;

my ($a, $b) = @ARGV;
die "usage: $0 <run A> <run B>\n" unless defined $b;

sub load {
    my ($dir) = @_;
    my (%tags, @sums);
    my $n = 0;
    for my $f (glob("$dir/alignments/*.sam.zst")) {
        (my $sample = $f) =~ s{.*/}{};
        open(my $fh, "-|", "zstdcat", $f) or die "$f: $!\n";
        while (<$fh>) {
            next if /^@/;
            chomp;
            my @r = split /\t/;
            next if $r[2] eq "*";
            my %t = map { /^(Z[USTL]):i:(\d+)$/ ? ($1 => $2) : () } @r[11 .. $#r];
            my @v = map { $t{$_} // 0 } qw(ZS ZU ZT);
            $tags{"$sample\t$r[0]\t$r[1]\t$r[2]\t$r[3]\t$r[5]"} = join(",", @v);
            $sums[$_] += $v[$_] for 0 .. 2;
            $sums[3] += defined $t{ZL};
            $n++;
        }
        close $fh;
    }
    return (\%tags, \@sums, $n);
}
my ($ta, $sa, $na) = load($a);
my ($tb, $sb, $nb) = load($b);
my ($shared, $same) = (0, 0);
my @examples;
for my $k (keys %$ta) {
    next unless exists $tb->{$k};
    $shared++;
    if ($ta->{$k} eq $tb->{$k}) { $same++ } elsif (@examples < 3) { push @examples, "$k: $ta->{$k} / $tb->{$k}" }
}
for my $x ([$a, $sa, $na], [$b, $sb, $nb]) {
    my ($dir, $s, $n) = @$x;
    printf "%s: %d mapped records; mean ZS %.2f, ZU %.2f, ZT %.2f, ZU/ZS %.3f; with ZL %d\n", $dir, $n, map({ $s->[$_] / ($n || 1) } 0 .. 2),
        $s->[1] / ($s->[0] || 1), $s->[3];
}
printf "the same alignment in both: %d records, %d (%.4f) with the same ZS, ZU and ZT\n", $shared, $same, $same / ($shared || 1);
print "  differ: $_\n" for @examples;
