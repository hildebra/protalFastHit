use strict; use warnings;
my (%runs, @hdr);
for my $d (glob "perf_*_[0-9]") {
  my ($lvl) = $d =~ /perf_(avx2|avx512)_/;
  open my $f, "<", "$d/runs.tsv" or die; chomp(my $h = <$f>); @hdr = split /\t/, $h;
  while (<$f>) { chomp; my @v = split /\t/; my %r; @r{@hdr} = @v; push @{$runs{$r{sample}}{$lvl}}, \%r; }
}
sub med { my @a = sort { $a <=> $b } grep { defined && $_ ne "NA" && $_ ne "" } @_; return "NA" unless @a; my $n = @a; $n % 2 ? $a[$n/2] : ($a[$n/2-1] + $a[$n/2]) / 2 }
my @cols = qw(wall_s user_s aligning_s instructions cycles ipc);
for my $s (sort keys %runs) {
  for my $l (qw(avx2 avx512)) {
    my @r = @{$runs{$s}{$l} || []};
    printf "%-4s %-6s n=%d", $s, $l, scalar @r;
    for my $c (@cols) { my $m = med(map { $_->{$c} } @r); printf "  %s %s", $c, $m =~ /^[\d.]+$/ && $m > 1e9 ? sprintf("%.3fT", $m/1e12) : $m; }
    print "\n";
  }
  # counts must agree
  my @count_cols = grep { my $i = $_; $i =~ /^(reads|anchored_reads|tried|screened|from_anchors|whole_windows|made|written|kmers|kmers_in_index|blocks_scanned|flex_cells|seeds|paired_seeds|dropped_lookups|anchors)$/ } @hdr;
  my %vals; for my $l (qw(avx2 avx512)) { for my $r (@{$runs{$s}{$l}}) { $vals{$_}{join ",", map { $r->{$_} // "" } $_}++ for @count_cols } }
  my @diff = grep { keys %{$vals{$_}} > 1 } @count_cols;
  print "     counts (", scalar(@count_cols), " columns) ", (@diff ? "DIFFER: @diff" : "equal across all runs"), "\n";
}
