#!/usr/bin/env perl
# Dumps a raw deflate stream: each block's type and code lengths, then its tokens (L <byte> or M <length> <distance>),
# one per line, with the input position. usage: deflate_dump.pl FILE OFFSET SIZE
use strict;
use warnings;

my ($file, $offset, $size) = @ARGV;
open(my $h, '<:raw', $file) or die "$file: $!";
local $/;
my $all = <$h>;
my @bytes = unpack('C*', substr($all, $offset, $size));
my ($pos, $bit) = (0, 0);    # byte and bit of the reader

sub bits {
    my ($n) = @_;
    my $v = 0;
    for my $i (0 .. $n - 1) {
        die "out of data\n" if $pos >= @bytes;
        $v |= (($bytes[$pos] >> $bit) & 1) << $i;
        if (++$bit == 8) { $bit = 0; $pos++; }
    }
    return $v;
}

# canonical Huffman decoding table from code lengths: {length => {code => symbol}}
sub table {
    my @len = @_;
    my %count;
    $count{$_}++ for grep { $_ } @len;
    my ($code, %next) = (0);
    for my $l (1 .. 15) {
        $code = ($code + ($count{$l - 1} // 0)) << 1;
        $next{$l} = $code;
    }
    my %t;
    for my $s (0 .. $#len) {
        next unless $len[$s];
        $t{$len[$s]}{ $next{ $len[$s] }++ } = $s;
    }
    return \%t;
}

sub decode {
    my ($t) = @_;
    my ($code, $l) = (0, 0);
    while (1) {
        $code = ($code << 1) | bits(1);
        $l++;
        return $t->{$l}{$code} if exists $t->{$l} && exists $t->{$l}{$code};
        die "bad code\n" if $l > 15;
    }
}

my @lbase = (3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258);
my @lext = (0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0);
my @dbase = (1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073,
             4097, 6145, 8193, 12289, 16385, 24577);
my @dext = (0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13);

my $out = 0;    # input position
my $last = 0;
until ($last) {
    $last = bits(1);
    my $type = bits(2);
    print "BLOCK type $type final $last at bit ", $pos * 8 + $bit, "\n";
    if ($type == 0) {
        $bit = 0, $pos++ if $bit;
        my $len = $bytes[$pos] | $bytes[$pos + 1] << 8;
        $pos += 4;
        print "STORED $len\n";
        $pos += $len;
        $out += $len;
        next;
    }
    my ($lt, $dt);
    if ($type == 1) {
        my @l = ((8) x 144, (9) x 112, (7) x 24, (8) x 8);
        ($lt, $dt) = (table(@l), table((5) x 30));
    } else {
        my $hlit = bits(5) + 257;
        my $hdist = bits(5) + 1;
        my $hclen = bits(4) + 4;
        my @order = (16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15);
        my @cl = (0) x 19;
        $cl[ $order[$_] ] = bits(3) for 0 .. $hclen - 1;
        my $ct = table(@cl);
        my @lens;
        while (@lens < $hlit + $hdist) {
            my $s = decode($ct);
            if ($s < 16) { push @lens, $s; }
            elsif ($s == 16) { my $r = 3 + bits(2); push @lens, ($lens[-1]) x $r; }
            elsif ($s == 17) { push @lens, (0) x (3 + bits(3)); }
            else { push @lens, (0) x (11 + bits(7)); }
        }
        my @ll = @lens[0 .. $hlit - 1];
        my @dl = @lens[$hlit .. $#lens];
        print "HLIT $hlit HDIST $hdist HCLEN $hclen\n";
        print "CL @cl\n";
        print "LL ", join(',', map { "$_:$ll[$_]" } grep { $ll[$_] } 0 .. $#ll), "\n";
        print "DL ", join(',', map { "$_:$dl[$_]" } grep { $dl[$_] } 0 .. $#dl), "\n";
        ($lt, $dt) = (table(@ll), table(@dl));
    }
    while (1) {
        my $s = decode($lt);
        if ($s < 256) { print "$out L $s\n"; $out++; next; }
        last if $s == 256;
        my $len = $lbase[ $s - 257 ] + bits($lext[ $s - 257 ]);
        my $d = decode($dt);
        my $dist = $dbase[$d] + bits($dext[$d]);
        print "$out M $len $dist\n";
        $out += $len;
    }
}
print "END $out bytes, ", $pos + ($bit ? 1 : 0), " of $size compressed bytes\n";
