#!/usr/bin/env perl
# Experiment only, never committed: patches a source tree (path as the argument) so that the hash maps keyed by taxid
# or gene id iterate in another order. With std::hash (the identity on integers) a tsl::sparse_map of dense taxids
# iterates in taxid order whatever the order of insertion, so these maps get a scrambling hash instead:
#   - GenomeLoader's genome map (src/SequenceUtils/GenomeLoader.h), whose genomes are also inserted again in descending
#     taxid order into a map with room for four times as many;
#   - the profile's taxon map and, where it is still a hash map, its gene map (src/Profiling/Profiler.h).
# Outputs must not change.
use strict;
use warnings;
my $root = shift or die "usage: $0 SOURCE_TREE\n";
my $hash = 'struct ScrambledHash { size_t operator()(size_t k) const { uint64_t x = k * 0x9E3779B97F4A7C15ull; return x ^ (x >> 29); } };';

sub patch {
    my ($file, @edits) = @_;
    open my $in, '<', $file or die "$file: $!\n";
    local $/;
    my $text = <$in>;
    close $in;
    while (my ($from, $to, $optional) = splice(@edits, 0, 3)) {
        my $n = ($text =~ s/\Q$from\E/$to/);
        die "$file: not found: $from\n" unless $n || $optional;
    }
    open my $out, '>', $file or die "$file: $!\n";
    print $out $text;
    close $out;
}

my $call = 'if (!gene_table || !LoadGeneTable(*gene_table, unique_size, m_threads)) LoadPositionMap(m_map, m_threads);';
my $method = <<'CPP';
        void ReorderForTest() {
            std::vector<GenomeKey> keys;
            for (auto const& [key, _] : m_genomes) keys.push_back(key);
            std::sort(keys.rbegin(), keys.rend());
            GenomeMap other;
            other.reserve(keys.size() * 4);
            for (auto key : keys) other.insert({ key, std::move(m_genomes.at(key)) });
            m_genomes = std::move(other);
            std::cerr << "TEST: genome map with a scrambling hash, made in descending taxid order; it iterates";
            size_t n = 0;
            for (auto const& [key, _] : m_genomes) if (n++ < 8) std::cerr << ' ' << key;
            std::cerr << " ..." << std::endl;
        }

CPP
patch("$root/src/SequenceUtils/GenomeLoader.h",
      'using GenomeMap = tsl::sparse_map<GenomeKey, Genome>;',
      "$hash\n        using GenomeMap = tsl::sparse_map<GenomeKey, Genome, ScrambledHash>;", 0,
      $call, "$call\n            ReorderForTest();", 0,
      '        bool FromGeneTable() const', $method . '        bool FromGeneTable() const', 0);
patch("$root/src/Profiling/Profiler.h",
      '        using GeneMap = ', "        $hash\n        using GeneMap = ", 0,
      'using TaxonMap = tsl::sparse_map<uint32_t, Taxon>;',
      'using TaxonMap = tsl::sparse_map<uint32_t, Taxon, ScrambledHash>;', 0,
      'using GeneMap = tsl::sparse_map<uint32_t, Gene>;',
      'using GeneMap = tsl::sparse_map<uint32_t, Gene, ScrambledHash>;', 1);
