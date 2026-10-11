#!/usr/bin/env perl
# The cost probe: a local patch of src/Profiling/Profiler.h (8d8c89b) that skips feature groups' evidence at run time,
# for measuring what each costs (instructions saved under callgrind). Not for production: the skipped features take
# their empty values. PROTAL_SKIP is a bitmask: 1 ancestry sites (also the fixed-site part of polymorphic), 2 column
# weights, 4 strain alleles (alleles + polymorphic), 8 congener gaps, 16 congener distances (distance group's
# relative_* and close share), 32 untried candidates, 64 consistency's unexpected congener fit and gene overlap.
# usage: cost_probe.pl <protal source dir>   (patches <dir>/src/Profiling/Profiler.h in place; dies if a hunk misses)
use strict;
use warnings;

my $f = "$ARGV[0]/src/Profiling/Profiler.h";
open(my $in, '<', $f) or die "$f: $!";
local $/;
my $s = <$in>;
close $in;

sub patch {
    my ($name, $from, $to) = @_;
    my $n = ($s =~ s/\Q$from\E/$to/);
    die "hunk $name did not apply\n" unless $n;
}

patch('include', "\nnamespace protal {\n", "\n#include <cstdlib>\n\nnamespace protal {\n");
patch('helper', "    namespace profiler {\n", <<'EOF');
    namespace profiler {
        // Cost probe (local measurement patch): PROTAL_SKIP bitmask of feature groups whose evidence is skipped.
        inline unsigned CostSkip() {
            static unsigned const skip = [] {
                char const* e = std::getenv("PROTAL_SKIP");
                return e ? static_cast<unsigned>(std::strtoul(e, nullptr, 0)) : 0u;
            }();
            return skip;
        }
EOF
patch('ancestry',
      "                if (!sites) sites = m_genome_loader->AncestrySitesOf(taxid, geneid);\n" .
      "                auto& shared = m_species_sites[key];\n" .
      "                if (!shared) shared = m_genome_loader->SpeciesAncestrySitesOf(taxid, geneid);\n",
      "                static auto const no_sites = std::make_shared<ancestry::Sites const>();\n" .
      "                if (!sites) sites = (CostSkip() & 1) ? no_sites : m_genome_loader->AncestrySitesOf(taxid, geneid);\n" .
      "                auto& shared = m_species_sites[key];\n" .
      "                if (!shared) shared = (CostSkip() & 1) ? no_sites : m_genome_loader->SpeciesAncestrySitesOf(taxid, geneid);\n");
patch('weights', "                if (!m_genome_loader->GetColumnWeights().Empty()) {\n                    e.weights_known = true;",
      "                if (!(CostSkip() & 2) && !m_genome_loader->GetColumnWeights().Empty()) {\n                    e.weights_known = true;");
patch('alleles', "if (kept) NoteAlleles(e, taxid, geneid, sam, sites.get());",
      "if (kept && !(CostSkip() & 4)) NoteAlleles(e, taxid, geneid, sam, sites.get());");
patch('gaps', "if (kept) NoteCopy(e, taxid, geneid, differences, aligned);",
      "if (kept && !(CostSkip() & 8)) NoteCopy(e, taxid, geneid, differences, aligned);");
patch('distance sketches',
      "                m_bm_distances.Start();\n                sam_chunks::ParallelFor(to_sketch.size()",
      "                m_bm_distances.Start();\n                if (!(CostSkip() & 16)) sam_chunks::ParallelFor(to_sketch.size()");
patch('distance pairs', "                sam_chunks::ParallelFor(pairs.size(), std::max<size_t>(threads, 1), [&](size_t p) {",
      "                if (!(CostSkip() & 16)) sam_chunks::ParallelFor(pairs.size(), std::max<size_t>(threads, 1), [&](size_t p) {");
patch('close share', "if (best >= 0) s.relative_close_share = CloseShare(",
      "if (best >= 0 && !(CostSkip() & 16)) s.relative_close_share = CloseShare(");
patch('untried', "if (!sam->m_untried.empty()) evidence.NoteUntriedCandidates(",
      "if (!sam->m_untried.empty() && !(CostSkip() & 32)) evidence.NoteUntriedCandidates(");
patch('unexpected fit', "if (same_genus && !neighbours.Empty() && alternative != taxid",
      "if (!(CostSkip() & 64) && same_genus && !neighbours.Empty() && alternative != taxid");
patch('gene overlap', "found != partners.end() && m_taxa.contains(found->second.first)) {",
      "found != partners.end() && !(CostSkip() & 64) && m_taxa.contains(found->second.first)) {");

open(my $out, '>', $f) or die "$f: $!";
print $out $s;
close $out;
print "patched $f\n";
