#!/usr/bin/env perl
# Compare two protal SAMs (zstd) record by record (key: read name + flag): which fields and tags differ, and how.
use strict; use warnings;
my ($a, $b) = @ARGV;
sub load {
    my ($f) = @_; my %r;
    open(my $fh, "-|", "zstdcat", $f) or die;
    while (<$fh>) { next if /^@/; chomp; my @f = split /\t/; my %t;
        for my $x (@f[11..$#f]) { my ($k, undef, $v) = split /:/, $x, 3; $t{$k} = $v; }
        $r{"$f[0]\t$f[1]"} = { rname => $f[2], pos => $f[3], mapq => $f[4], cigar => $f[5], tags => \%t }; }
    return \%r;
}
my ($ra, $rb) = (load($a), load($b));
my %n; my %tagdiff; my %mapq; my %tagsum;
for my $k (keys %$ra) {
    unless (exists $rb->{$k}) { $n{only_a}++; next; }
    my ($x, $y) = ($ra->{$k}, $rb->{$k}); $n{both}++;
    $n{rname_diff}++ if $x->{rname} ne $y->{rname};
    $n{pos_diff}++ if $x->{pos} ne $y->{pos};
    $n{cigar_diff}++ if $x->{cigar} ne $y->{cigar};
    if ($x->{mapq} != $y->{mapq}) { $n{mapq_diff}++; $mapq{$x->{mapq} < $y->{mapq} ? "up" : "down"}++;
        $mapq{"cross4"}++ if ($x->{mapq} >= 4) != ($y->{mapq} >= 4); }
    my %tk = map { $_ => 1 } (keys %{$x->{tags}}, keys %{$y->{tags}});
    for my $t (keys %tk) {
        my ($u, $v) = ($x->{tags}{$t} // "-", $y->{tags}{$t} // "-");
        $tagdiff{$t}++ if $u ne $v;
        if ($t =~ /^Z[UTN]$/) { $tagsum{$t}[0] += ($u eq "-" ? 0 : $u); $tagsum{$t}[1] += ($v eq "-" ? 0 : $v); }
    }
}
for my $k (keys %$rb) { $n{only_b}++ unless exists $ra->{$k}; }
print join(" ", map { "$_=$n{$_}" } sort keys %n), "\n";
print "mapq: ", join(" ", map { "mapq_$_=$mapq{$_}" } sort keys %mapq), "\n";
print "tags differing: ", join(" ", map { "differ_$_=$tagdiff{$_}" } sort keys %tagdiff), "\n";
print "tag sums a/b: ", join(" ", map { "$_=$tagsum{$_}[0]/$tagsum{$_}[1]" } sort keys %tagsum), "\n";
