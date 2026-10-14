#!/usr/bin/env perl
# How the ancestry sites of an ancestry_sites.bin spread over copies by density (sites over the span from the first to
# the last site), and what a cap on the density would keep: per cap (a site per D bases at most), the copies above it
# and the share of all sites kept. usage: anc_density.pl ancestry_sites.bin
use strict;
use warnings;
my $file = shift or die "usage: $0 ancestry_sites.bin\n";
open(my $fh, "<:raw", $file) or die "$file: $!\n";
my $c = do { local $/; <$fh> };
close $fh;
my $at = length "protal ancestry sites 1\n";
my ($n) = unpack("Q<", substr($c, $at, 8));
$at += 16;
my ($offsets_at, $payload_at) = ($at + 8 * $n, $at + 16 * $n);
sub varint { my ($p) = @_; my ($v, $s, $st) = (0, 0, $p); while (1) { my $b = ord(substr($c, $p++, 1)); $v |= ($b & 0x7f) << $s; last unless $b & 0x80; $s += 7 } return ($v, $p - $st) }
my (@copies, $all);
for my $i (0 .. $n - 1) {
    my $p = $payload_at + unpack("Q<", substr($c, $offsets_at + 8 * $i, 8));
    (undef, my $b) = varint($p); $p += $b;
    my $count;
    for (1 .. 5) { (my $v, $b) = varint($p); $p += $b; $count = $v }
    my ($prev, $first, $last) = (0, undef, 0);
    for (1 .. $count) { (my $z, $b) = varint($p); $p += $b; $prev += ($z >> 1) ^ -($z & 1); $first //= $prev; $last = $prev }
    next unless $count;
    push @copies, [$count, $last - $first + 1];
    $all += $count;
}
printf "%d copies, %d sites; density (sites a base of the span) median %.3f\n", scalar @copies, $all,
    (sort { $a <=> $b } map { $_->[0] / $_->[1] } @copies)[@copies / 2];
for my $d (4, 8, 12, 16, 24, 32) {
    my ($kept, $above) = (0, 0);
    for my $x (@copies) {
        my $cap = int($x->[1] / $d) + 1;
        $kept += $x->[0] < $cap ? $x->[0] : $cap;
        $above++ if $x->[0] > $cap;
    }
    printf "a site per %2d bases at most: %5.1f%% of the copies thinned, %5.1f%% of the sites kept\n", $d, 100 * $above / @copies, 100 * $kept / $all;
}
