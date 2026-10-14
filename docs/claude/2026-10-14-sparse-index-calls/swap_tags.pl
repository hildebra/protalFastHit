#!/usr/bin/env perl
# Counterfactual SAM: the base run's records with some of their fields taken from a donor run's same record (read name +
# flag), or scaled. usage: swap_tags.pl <base.sam.zst> <donor.sam.zst|-> <spec> <out.sam.zst>
# spec, comma-separated: MAPQ or a tag (ZU ZT ZA ZF ZC ZN ZR) from the donor (absent there: removed); x<N> multiplies ZU
# and ZT by N (capped at 65535); "none" copies the base.
use strict; use warnings;
my ($base, $donor, $spec, $out) = @ARGV;
my @items = grep { $_ ne "none" } split /,/, $spec;
my $same = grep { $_ eq "same" } @items;  # "same": MAPQ only from a donor record on the same gene (RNAME)
@items = grep { $_ ne "same" } @items;
my %take = map { $_ => 1 } grep { !/^x\d+$/ } @items;
my ($scale) = map { /^x(\d+)$/ ? $1 : () } @items;
my @order = qw(ZU ZT ZA ZF ZC ZR ZN);
my %rank; @rank{@order} = (0 .. $#order);

sub split_rec {
    my ($line) = @_;
    my @f = split /\t/, $line;
    my (%t, %type, @seen);
    for my $x (@f[11 .. $#f]) { my ($k, $ty, $v) = split /:/, $x, 3; $t{$k} = $v; $type{$k} = $ty; push @seen, $k; }
    return ([@f[0 .. 10]], \%t, \%type, \@seen);
}
my %d;
if (%take) {
    open(my $dh, "-|", "zstdcat", $donor) or die "donor $donor";
    while (my $l = <$dh>) { next if $l =~ /^@/; chomp $l; my ($f, $t, $ty) = split_rec($l);
        $d{"$f->[0]\t$f->[1]"} = [$f->[4], $t, $ty, $f->[2]]; }
    close $dh;
}
open(my $bh, "-|", "zstdcat", $base) or die "base $base";
open(my $oh, "|-", "zstd -q -f -o $out") or die "out $out";
my %n;
while (my $l = <$bh>) {
    if ($l =~ /^@/) { print $oh $l; next; }
    chomp $l;
    my ($f, $t, $ty, $seen) = split_rec($l);
    my $key = "$f->[0]\t$f->[1]";
    if (%take && exists $d{$key}) {
        my ($mapq, $dt, $dty, $drname) = @{$d{$key}};
        $n{swapped}++;
        if ($take{MAPQ} && (!$same || $drname eq $f->[2])) { $f->[4] = $mapq; $n{mapq_taken}++; }
        for my $tag (grep { $_ ne "MAPQ" } keys %take) {
            if (exists $dt->{$tag}) { $t->{$tag} = $dt->{$tag}; $ty->{$tag} = $dty->{$tag}; push @$seen, $tag; }
            else { delete $t->{$tag}; }
        }
    } elsif (%take) { $n{no_donor}++; }
    if ($scale) { for my $tag (qw(ZU ZT)) { $t->{$tag} = $t->{$tag} * $scale > 65535 ? 65535 : $t->{$tag} * $scale if exists $t->{$tag}; } }
    my %u; my @keys = grep { exists $t->{$_} && !$u{$_}++ } @$seen;
    @keys = sort { ($rank{$a} // 99) <=> ($rank{$b} // 99) } @keys;
    print $oh join("\t", @$f, map { "$_:$ty->{$_}:$t->{$_}" } @keys), "\n";
}
close $bh; close $oh;
print STDERR "$out: ", join(" ", map { "$_=$n{$_}" } sort keys %n), "\n";
