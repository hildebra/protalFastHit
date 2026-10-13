#!/usr/bin/env perl
# mini_reads.pl REFERENCE.fna PAIRS SEED PREFIX: paired 100 bp reads of 220-320 bp fragments of the reference's genes,
# both orientations, 0.5% substitutions, into PREFIX_R1.fq and PREFIX_R2.fq.
use strict;
use warnings;

my ($fasta, $pairs, $seed, $prefix) = @ARGV;
srand($seed);
my (@genes, $seq);
open(my $in, '<', $fasta) or die "$fasta: $!";
while (<$in>) {
    chomp;
    if (/^>/) { push @genes, $seq if defined $seq && length $seq >= 320; $seq = ''; next; }
    $seq .= uc $_;
}
push @genes, $seq if defined $seq && length $seq >= 320;
close $in;
die "no genes" unless @genes;

sub revcomp { my $s = reverse $_[0]; $s =~ tr/ACGT/TGCA/; return $s; }
sub mutate {
    my @b = split //, $_[0];
    for (@b) { $_ = (grep { $_ ne $b[0] } qw(A C G T))[int rand 3] if rand() < 0.005; }
    return join '', @b;
}
open(my $r1, '>', "${prefix}_R1.fq") or die;
open(my $r2, '>', "${prefix}_R2.fq") or die;
for my $n (1 .. $pairs) {
    my $g = $genes[int rand @genes];
    my $len = 220 + int rand 101;
    my $frag = substr($g, int rand(length($g) - $len + 1), $len);
    $frag = revcomp($frag) if rand() < 0.5;
    my $q = 'I' x 100;
    print $r1 "\@p.$n/1\n", mutate(substr($frag, 0, 100)), "\n+\n$q\n";
    print $r2 "\@p.$n/2\n", mutate(revcomp(substr($frag, -100))), "\n+\n$q\n";
}
close $r1;
close $r2;
