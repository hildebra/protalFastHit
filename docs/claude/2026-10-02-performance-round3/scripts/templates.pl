#!/usr/bin/env perl
# templates.pl MANIFEST BASES SEED > templates.fa: read templates for long reads from a community manifest: genomes drawn
# by relative abundance x genome length, a contig by its length, a start uniformly, a length ~ N(15000, 3000) (at least
# 1,000, at most the contig), either strand.
use strict; use warnings;
my ($manifest, $bases, $seed) = @ARGV; srand($seed // 1);
open my $m, '<', $manifest or die; my $h = <$m>; chomp $h; my @c = split /\t/, $h; my %i; @i{@c} = 0 .. $#c;
my (@g, @w); my $tw = 0;
while (<$m>) { chomp; my @f = split /\t/; push @g, $f[$i{fasta_path}]; my $w = $f[$i{relative_abundance}] * $f[$i{genome_length}]; push @w, $w; $tw += $w }
my %cache;
sub contigs { my $p = shift; return $cache{$p} if $cache{$p}; my $fh; if ($p =~ /\.gz$/) { open $fh, '-|', 'zcat', $p or die } else { open $fh, '<', $p or die }
  my (@s, $cur); while (<$fh>) { chomp; if (/^>/) { push @s, $cur if defined $cur; $cur = '' } else { $cur .= uc $_ } } push @s, $cur if defined $cur; return $cache{$p} = \@s }
my ($done, $n) = (0, 0);
while ($done < $bases) {
  my $r = rand($tw); my $k = 0; $r -= $w[$k], $k++ while $k < $#w && $r >= $w[$k];
  my $s = contigs($g[$k]); my $tot = 0; $tot += length $_ for @$s; my $x = rand($tot); my $j = 0; $x -= length $s->[$j], $j++ while $j < $#$s && $x >= length $s->[$j];
  my $contig = $s->[$j]; my $len = int(15000 + 3000 * sqrt(-2 * log(rand() || 1e-12)) * cos(6.283185307 * rand()));
  $len = 1000 if $len < 1000; $len = length $contig if $len > length $contig; my $start = int(rand(length($contig) - $len + 1));
  my $t = substr($contig, $start, $len); if (rand() < 0.5) { $t = reverse $t; $t =~ tr/ACGT/TGCA/ }
  printf ">t%d %s\n%s\n", $n++, $g[$k], $t; $done += $len;
}
