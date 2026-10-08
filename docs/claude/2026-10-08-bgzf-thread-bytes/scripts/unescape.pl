#!/usr/bin/env perl
# Turns a gtest "Which is: "..." "..."" line (an escaped std::string) back into its bytes.
# usage: unescape.pl LOG_EXCERPT LINE_NUMBER OUT
use strict;
use warnings;

my ($file, $want, $out) = @ARGV;
open(my $in, '<:raw', $file) or die "$file: $!";
my $line;
while (<$in>) {
    if ($. == $want) { $line = $_; last; }
}
close $in;
die "no line $want\n" unless defined $line;
$line =~ s/\r?\n\z//;
$line =~ s/^\s*Which is: // or die "not a Which is line\n";

my %simple = ('0' => "\0", "'" => "'", '\\' => '\\', 'a' => "\a", 'b' => "\b", 'f' => "\f",
              'n' => "\n", 'r' => "\r", 't' => "\t", 'v' => "\x0b", '"' => '"');
my $bytes = '';
my $i = 0;
my $n = length $line;
while ($i < $n) {
    my $c = substr($line, $i, 1);
    if ($c eq ' ') { $i++; next; }
    die "expected a quote at $i\n" unless $c eq '"';
    $i++;
    while (1) {
        die "unterminated at $i\n" if $i >= $n;
        $c = substr($line, $i, 1);
        if ($c eq '"') { $i++; last; }
        if ($c eq '\\') {
            my $e = substr($line, $i + 1, 1);
            if ($e eq 'x') {
                my $j = $i + 2;
                $j++ while $j < $n && substr($line, $j, 1) =~ /[0-9A-Fa-f]/;
                $bytes .= chr(hex(substr($line, $i + 2, $j - $i - 2)));
                $i = $j;
            } elsif (exists $simple{$e}) {
                $bytes .= $simple{$e};
                $i += 2;
            } else {
                die "unknown escape \\$e at $i\n";
            }
            next;
        }
        $bytes .= $c;
        $i++;
    }
}
open(my $o, '>:raw', $out) or die "$out: $!";
print $o $bytes;
close $o;
printf "%s: %d bytes\n", $out, length $bytes;
