#!/usr/bin/env perl
# What the bytes of an ancestry_sites.bin (AncestrySites.h, "protal ancestry sites 1") hold, and what leaner codings
# would take: the index (8-byte key and offset per copy), each entry's length, its header fields (congener, congeners,
# compared, identity, site count, packed flag, indel count), the positions (zigzag steps), the bases and the indels;
# then the positions as Rice-coded gaps (the best parameter per copy, 3 bits to name it), the index as 4-byte offsets
# by dense gene row, and a lean header (congener varint, congeners 1 byte, identity 1 byte, counts varint).
# usage: anc_breakdown.pl ancestry_sites.bin
use strict;
use warnings;

my $file = shift or die "usage: $0 ancestry_sites.bin\n";
open(my $fh, "<:raw", $file) or die "$file: $!\n";
my $content = do { local $/; <$fh> };
close $fh;
my $magic = "protal ancestry sites 1\n";
substr($content, 0, length $magic) eq $magic or die "not an ancestry sites table\n";
my $at = length $magic;
my ($n, $payload_size) = unpack("Q<Q<", substr($content, $at, 16));
$at += 16;
my $offsets_at = $at + 8 * $n;
my $payload_at = $at + 16 * $n;

my %bytes = (index => 16 * $n, head => length($magic) + 16);
my ($sites, $indels, $copies_with_sites, $rice_bits, $lean_head, $v2) = (0, 0, 0, 0, 0, length("protal ancestry sites 2
") + 24);
my @sites_per_copy;

sub varint {    # (value, bytes) at position p of $content
    my ($p) = @_;
    my ($v, $shift, $start) = (0, 0, $p);
    while (1) {
        my $b = ord(substr($content, $p++, 1));
        $v |= ($b & 0x7f) << $shift;
        last unless $b & 0x80;
        $shift += 7;
    }
    return ($v, $p - $start);
}
sub vlen { my $v = shift; my $l = 1; while ($v >= 0x80) { $v >>= 7; $l++ } return $l }

for my $i (0 .. $n - 1) {
    my $offset = unpack("Q<", substr($content, $offsets_at + 8 * $i, 8));
    my $p = $payload_at + $offset;
    my ($len, $lb) = varint($p);
    $bytes{length} += $lb;
    $p += $lb;
    my $end = $p + $len;
    my @f;
    for my $field (qw(congener congeners compared identity count)) {
        my ($v, $b) = varint($p);
        $p += $b;
        $bytes{"h_$field"} += $b;
        push @f, $v;
    }
    my ($congener, $congeners, $compared, $identity, $count) = @f;
    # the positions: zigzag steps (already counted by reading them)
    my ($prev, @pos) = (0);
    my $pos_start = $p;
    for (1 .. $count) {
        my ($z, $b) = varint($p);
        $p += $b;
        my $step = ($z >> 1) ^ -($z & 1);
        $prev += $step;
        push @pos, $prev;
    }
    my $positions_v1 = $p - $pos_start;
    $bytes{positions} += $positions_v1;
    my $packed = ord(substr($content, $p, 1));
    $bytes{h_flag} += 1;
    $p += 1;
    my $base_bytes = $packed == 1 ? int(($count + 3) / 4) : $count;
    $bytes{bases} += $base_bytes;
    $p += $base_bytes;
    my ($ni, $b) = varint($p);
    $bytes{h_indels} += $b;
    $p += $b;
    my $indel_start = $p;
    for (1 .. $ni) {
        my (undef, $b1) = varint($p);
        $p += $b1;
        my (undef, $b2) = varint($p);
        $p += $b2;
    }
    my $indel_bytes = $p - $indel_start;
    $bytes{indels} += $indel_bytes;
    $p == $end or die "entry $i: decoded to $p, ends at $end\n";
    $sites += $count;
    $indels += $ni;
    $copies_with_sites++ if $count;
    push @sites_per_copy, $count;
    # Rice-coded gaps: the best k per copy (gap - 1 >= 0; the first gap from 0), 3 bits to name k.
    my @gaps;
    my $q = -1;
    for my $x (@pos) { push @gaps, $x - $q - 1; $q = $x; }
    my $best;
    for my $k (0 .. 7) {
        my $bits = 3;
        $bits += ($_ >> $k) + 1 + $k for @gaps;
        $best = $bits if !defined $best || $bits < $best;
    }
    $rice_bits += $best // 3;
    # A lean header: congener varint, congeners 1 byte, identity 1 byte, compared varint, count varint, indels varint.
    $lean_head += vlen($congener) + 1 + 1 + vlen($compared) + vlen($count) + vlen($ni);
    # Version 2 as EncodeCompact writes it: varints congener, congeners, compared; identity 2 bytes; count varint; with
    # sites a flags byte, the Rice bits padded to a byte (ascending positions) and the bases; the indels as before; the
    # entry's length varint; 8 index bytes (4-byte key and offset).
    my $ascending = 1;
    for my $j (1 .. $#pos) { $ascending = 0 if $pos[$j] <= $pos[$j - 1] }
    my $positions_v2 = 0;
    if ($count) {
        if ($ascending) {
            my ($bestbits, @g);
            my $q2 = -1;
            for my $x (@pos) { push @g, $x - $q2 - 1; $q2 = $x }
            for my $k (0 .. 15) { my $bits = 0; $bits += ($_ >> $k) + 1 + $k for @g; $bestbits = $bits if !defined $bestbits || $bits < $bestbits }
            $positions_v2 = int(($bestbits + 7) / 8);
        } else {
            $positions_v2 = $positions_v1;  # out of order: the zigzag steps of version 1
        }
    }
    my $entry_v2 = vlen($congener) + vlen($congeners) + vlen($compared) + 2 + vlen($count) + ($count ? 1 + $positions_v2 + $base_bytes : 0)
        + vlen($ni) + $indel_bytes;
    $v2 += 8 + vlen($entry_v2) + $entry_v2;
}
my $total = 0;
$total += $_ for values %bytes;
printf "file %s: %d bytes, %d copies (%d with base sites), %d base sites (%.1f a copy), %d indels\n", $file,
    length($content), $n, $copies_with_sites, $sites, $sites / ($n || 1), $indels;
printf "%-14s %12s %6s %8s\n", "part", "bytes", "share", "per copy";
for my $k (sort { $bytes{$b} <=> $bytes{$a} } keys %bytes) {
    printf "%-14s %12d %5.1f%% %8.1f\n", $k, $bytes{$k}, 100 * $bytes{$k} / $total, $bytes{$k} / ($n || 1);
}
printf "%-14s %12d\n", "sum", $total;
my $rice = int($rice_bits / 8 + 0.5);
my $lean = 4 * $n + $bytes{length} + $lean_head + $rice + $bytes{bases} + $bytes{indels};
printf "positions as Rice-coded gaps: %d bytes (%.0f%% of the zigzag varints)\n", $rice, 100 * $rice / ($bytes{positions} || 1);
printf "lean layout (4-byte offsets by row, lean header, Rice gaps): %d bytes, %.0f%% of the file, %.1f a copy\n", $lean,
    100 * $lean / length($content), $lean / ($n || 1);
my @s = sort { $a <=> $b } @sites_per_copy;
printf "version 2 (EncodeCompact, 4-byte index): %d bytes, %.1f%% of the file, %.1f a copy
", $v2, 100 * $v2 / length($content), $v2 / ($n || 1);
printf "sites a copy: median %d, 90%% %d, max %d\n", $s[@s / 2], $s[int(0.9 * @s)], $s[-1];
