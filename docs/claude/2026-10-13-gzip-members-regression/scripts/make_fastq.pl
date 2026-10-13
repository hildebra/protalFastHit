#!/usr/bin/env perl
# make_fastq.pl RECORDS SEED > out.fq
# Illumina-like FASTQ (NovaSeq style: 151 bp, binned qualities F : , #, long headers), as a real run's R1 would look,
# for timing gzip readers. Bases come from random offsets of a random pool, qualities from a pool of binned strings.
use strict;
use warnings;

my ($records, $seed) = @ARGV;
$records //= 1_000_000;
srand($seed // 1);

my @b = qw(A C G T);
my $pool = join '', map { $b[int rand 4] } 1 .. 4_000_000;
my $plen = length $pool;

# Binned qualities: mostly F, some ':' and ',' runs, '#' tails on a tenth of the reads.
my @quals;
for my $i (1 .. 5000) {
    my $q = '';
    while (length $q < 151) {
        my $r = rand;
        my $c = $r < 0.86 ? 'F' : $r < 0.94 ? ':' : $r < 0.99 ? ',' : '#';
        $q .= $c x (1 + int rand($c eq 'F' ? 30 : 3));
    }
    $q = substr $q, 0, 151;
    if (rand() < 0.1) {
        my $t = 1 + int rand 60;
        substr($q, 151 - $t) = '#' x $t;
    }
    push @quals, $q;
}

my @idx = map { join '', map { $b[int rand 4] } 1 .. 8 } 1 .. 4;
my $out = '';
my ($x, $y, $tile) = (1000, 1000, 1101);
for my $n (1 .. $records) {
    $x += 1 + int rand 40;
    if ($x > 32000) { $x = 1000 + int rand 100; $y += 1 + int rand 30; }
    if ($y > 37000) { $y = 1000; $tile++; }
    my $o = int rand($plen - 151);
    $out .= "\@A00261:712:HNKGWDSX5:1:$tile:$x:$y 1:N:0:$idx[$n % 4]+$idx[($n + 1) % 4]\n"
          . substr($pool, $o, 151) . "\n+\n" . $quals[int rand @quals] . "\n";
    if (length $out > 1 << 20) { print $out; $out = ''; }
}
print $out;
