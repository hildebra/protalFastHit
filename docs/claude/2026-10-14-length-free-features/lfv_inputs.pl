#!/usr/bin/env perl
# A model's inputs changed, as a PMML stub for machine_learning_cmdline.py --reference-pmml (which trains on the inputs
# of a model's MiningSchema): the first MiningSchema of MODEL (the forest's own, every feature in order), without the
# features of --drop, with those of --swap renamed (old:new) and those of --add appended, written to OUT with no trees.
# Usage: lfv_inputs.pl MODEL OUT [--drop a,b,...] [--swap old:new,...] [--add a,b,...]
use strict;
use warnings;
use Getopt::Long;

my ($drop, $swap, $add) = ("", "", "");
GetOptions("drop=s" => \$drop, "swap=s" => \$swap, "add=s" => \$add)
    or die "usage: $0 MODEL OUT [--drop a,b] [--swap old:new] [--add a,b]\n";
my ($model, $out) = @ARGV;
die "usage: $0 MODEL OUT [--drop a,b] [--swap old:new]\n" unless defined $out;
my %drop = map { $_ => 1 } grep { length } split /,/, $drop;
my %swap = map { split /:/, $_, 2 } grep { length } split /,/, $swap;

open my $in, "<", $model or die "$model: $!\n";
my $text = do { local $/; <$in> };
close $in;
my ($schema) = $text =~ m{<MiningSchema>(.*?)</MiningSchema>}s or die "$model: no MiningSchema\n";
my @fields;
while ($schema =~ m{<MiningField\s+([^>]*?)/?>}g) {
    my $attrs = $1;
    my ($name) = $attrs =~ /name="([^"]*)"/ or next;
    next if $attrs =~ /usageType="(?!active)/;   # the predicted field
    push @fields, $name;
}
my %have = map { $_ => 1 } @fields;
for my $f (keys %drop, keys %swap) {
    die "$model: no input $f\n" unless $have{$f};
}
my @kept = map { $swap{$_} // $_ } grep { !$drop{$_} } @fields;
push @kept, grep { length } split /,/, $add;
my %seen;
my @dup = grep { $seen{$_}++ } @kept;
die "inputs twice: @dup\n" if @dup;

open my $o, ">", $out or die "$out: $!\n";
print $o qq{<?xml version="1.0" encoding="UTF-8"?>\n<PMML version="4.4">\n <DataDictionary numberOfFields="},
    scalar(@kept) + 1, qq{">\n};
print $o qq{  <DataField name="truth" optype="categorical" dataType="string"><Value value="FALSE"/><Value value="TRUE"/></DataField>\n};
print $o qq{  <DataField name="$_" optype="continuous" dataType="double"/>\n} for @kept;
print $o qq{ </DataDictionary>\n <MiningModel functionName="classification">\n  <MiningSchema><MiningField name="truth" usageType="predicted"/>};
print $o qq{<MiningField name="$_"/>} for @kept;
print $o qq{</MiningSchema>\n  <Segmentation multipleModelMethod="average"/>\n </MiningModel>\n</PMML>\n};
close $o;
printf STDERR "%s: %d inputs (of %d; dropped %d, swapped %d, added %d)\n", $out, scalar @kept, scalar @fields,
    scalar keys %drop, scalar keys %swap, scalar grep { length } split /,/, $add;
