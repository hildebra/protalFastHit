#!/usr/bin/env perl
# Pool bench_inflate result lines: max/median MB/s (wall and cpu) per method.
use strict;
use warnings;
my (%w, %c, @order, %seen);
while (<>) {
    chomp;
    my @f = split /\t/;
    my $m = $f[0];
    my %h = map { split /=/, $_, 2 } @f[1 .. $#f];
    push @{ $w{$m} }, $h{wall_MBps};
    push @{ $c{$m} }, $h{cpu_MBps};
    push @order, $m unless $seen{$m}++;
}
printf "%-10s %5s %12s %12s %12s %12s\n", "method", "runs", "wall_max", "wall_median", "cpu_max", "cpu_median";
for my $m (@order) {
    my @a = sort { $a <=> $b } @{ $w{$m} };
    my @b = sort { $a <=> $b } @{ $c{$m} };
    my $med = sub { my @x = @_; @x % 2 ? $x[$#x / 2] : ($x[@x / 2 - 1] + $x[@x / 2]) / 2 };
    printf "%-10s %5d %12.1f %12.1f %12.1f %12.1f\n", $m, scalar @a, $a[-1], $med->(@a), $b[-1], $med->(@b);
}
