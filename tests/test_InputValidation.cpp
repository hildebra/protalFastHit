// Unit tests for the database and input files protal reads: malformed files must stop with a clear
// message (exit 8) instead of being half-read, and the index must know which reference it was built for.
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <unistd.h>
#include <zstd.h>
#include "Options.h"
#include "Hash/Seedmap.h"
#include "Profiling/Profiler.h"
#include "SequenceUtils/GenomeLoader.h"
#include "SequenceUtils/KmerUtils.h"
#include "Taxonomy/Taxonomy.h"
#include "Utilities/ReferenceFingerprint.h"
#include "SequenceUtils/ReadTypeDetection.h"
#include "gzstream/gzstream.h"
#include "TestUtil.h"

namespace fs = std::filesystem;
using protal::test::ScratchDir;

namespace {
    // Two genera with two and one species, as the mini DB lays them out.
    constexpr char kTaxonomy[] =
            "id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n"
            "1\t6\t0\ts__Mockella alpha\tspecies\t7\tGCF_1\n"
            "2\t6\t0\ts__Mockella beta\tspecies\t7\tGCF_2\n"
            "3\t7\t0\ts__Fakibacter gamma\tspecies\t7\tGCF_3\n"
            "4\t4\t0\troot\tno rank\t0\t\n"
            "5\t4\t0\td__Bacteria\tdomain\t1\t\n"
            "6\t5\t0\tg__Mockella\tgenus\t6\t\n"
            "7\t5\t0\tg__Fakibacter\tgenus\t6\t\n";

    void LoadTaxonomy(std::string const& path) { protal::taxonomy::IntTaxonomy taxonomy(path); }
}

TEST(Taxonomy, LoadsWithOrWithoutHeaderAndWithCRLF) {
    ScratchDir dir;
    protal::taxonomy::IntTaxonomy with_header(dir.Write("a.dmp", kTaxonomy));
    EXPECT_EQ(with_header.Get("s__Mockella beta"), 2u);
    EXPECT_EQ(with_header.root_id, 4u);

    std::string body = kTaxonomy;
    body = body.substr(body.find('\n') + 1);
    protal::taxonomy::IntTaxonomy without_header(dir.Write("b.dmp", body));
    EXPECT_EQ(without_header.map.size(), with_header.map.size());

    std::string crlf;
    for (char c : std::string(kTaxonomy)) crlf += c == '\n' ? std::string("\r\n") : std::string(1, c);
    protal::taxonomy::IntTaxonomy windows(dir.Write("c.dmp", crlf + "\r\n"));
    EXPECT_EQ(windows.map.at(3).rank, "species");
    EXPECT_EQ(windows.map.at(3).rep_genome, "GCF_3");
}

TEST(Taxonomy, RejectsMalformedFiles) {
    ScratchDir dir;
    std::string header = "id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n";
    EXPECT_EXIT(LoadTaxonomy((dir.path / "none.dmp").string()),
                testing::ExitedWithCode(8), "cannot open");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("empty.dmp", header)), testing::ExitedWithCode(8), "no taxa");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("short.dmp", header + "1\t1\t0\troot\tno rank\n")),
                testing::ExitedWithCode(8), "line 2: expected at least 6");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("text.dmp", header + "1\tone\t0\troot\tno rank\t0\n")),
                testing::ExitedWithCode(8), "must be integers");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("twice.dmp", std::string(kTaxonomy) + "2\t6\t0\ts__Other\tspecies\t7\t\n")),
                testing::ExitedWithCode(8), "id 2 is defined twice");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("orphan.dmp", std::string(kTaxonomy) + "8\t9\t0\ts__Lost\tspecies\t7\t\n")),
                testing::ExitedWithCode(8), "parent id 9 is used but never defined");
    // Lineages are walked up to the root, the taxon that is its own parent: a cycle, or a root
    // given a parent, would make the walk go on for ever.
    EXPECT_EXIT(LoadTaxonomy(dir.Write("cycle.dmp", std::string(kTaxonomy) + "8\t9\t0\ts__Round\tspecies\t7\t\n9\t8\t0\tg__Round\tgenus\t6\t\n")),
                testing::ExitedWithCode(8), "never reaches the root .a cycle through id [89]");
    std::string rooted_below = kTaxonomy;
    rooted_below.replace(rooted_below.find("4\t4\t0\troot"), 3, "4\t5");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("rootless.dmp", rooted_below)), testing::ExitedWithCode(8), "never reaches the root");
    // Truth files and --msa_species name taxa, so a name must be one taxon's.
    EXPECT_EXIT(LoadTaxonomy(dir.Write("same.dmp", std::string(kTaxonomy) + "8\t6\t0\ts__Mockella alpha\tspecies\t7\t\n")),
                testing::ExitedWithCode(8), "line 9: the name 's__Mockella alpha' is given to ids 1 and 8; names must be unique");
}

TEST(Truth, ReadsLineagesAndTaxids) {
    ScratchDir dir;
    protal::taxonomy::IntTaxonomy taxonomy(dir.Write("t.dmp", kTaxonomy));

    // simulate_metagenomes writes the lineage in a later column; any field may hold it.
    auto lineage = protal::GetTruth(dir.Write("lineage.tsv",
            "genome\tlineage\treads\n"
            "GCF_1\td__Bacteria;g__Mockella;s__Mockella alpha\t100\r\n"
            "\n"
            "# a comment\n"
            "GCF_3\td__Bacteria;g__Fakibacter;s__Fakibacter gamma\t50\n"), taxonomy);
    EXPECT_EQ(lineage.size(), 2u);
    EXPECT_TRUE(lineage.contains(1));
    EXPECT_TRUE(lineage.contains(3));

    // Single-column species names and internal taxids.
    auto names = protal::GetTruth(dir.Write("names.tsv", "s__Mockella beta\n2\n3\t0.4\n"), taxonomy);
    EXPECT_EQ(names.size(), 2u);
    EXPECT_TRUE(names.contains(2));
    EXPECT_TRUE(names.contains(3));
}

TEST(Truth, SkipsLinesWithoutAKnownSpecies) {
    ScratchDir dir;
    protal::taxonomy::IntTaxonomy taxonomy(dir.Write("t.dmp", kTaxonomy));
    testing::internal::CaptureStderr();
    auto truth = protal::GetTruth(dir.Write("partial.tsv",
            "GCF_9\td__Bacteria;g__Mockella;s__\n"           // unclassified at species level
            "GCF_8\td__Bacteria;s__Unknown species\n"
            "99\n"                                           // taxid not in the taxonomy
            "just text\n"
            "GCF_2\td__Bacteria;g__Mockella;s__Mockella beta\n"), taxonomy);
    auto log = testing::internal::GetCapturedStderr();
    EXPECT_EQ(truth.size(), 1u);
    EXPECT_TRUE(truth.contains(2));
    EXPECT_NE(log.find("No species classification"), std::string::npos);
    EXPECT_NE(log.find("Species unknown to taxonomy: s__Unknown species"), std::string::npos);
    EXPECT_NE(log.find("taxid 99 is not in the taxonomy"), std::string::npos);
    EXPECT_NE(log.find("line 4: no lineage or taxid"), std::string::npos);
}

TEST(Truth, MissingFileStops) {
    ScratchDir dir;
    protal::taxonomy::IntTaxonomy taxonomy(dir.Write("t.dmp", kTaxonomy));
    EXPECT_EXIT(protal::GetTruth((dir.path / "none.tsv").string(), taxonomy), testing::ExitedWithCode(8), "Cannot read truth file");
}

namespace {
    // Two genes of taxon 1 and one of taxon 2; the map holds the byte range of each sequence line.
    struct Reference {
        std::string fna;
        std::string map;
        Reference() {
            std::ostringstream fna_os, map_os;
            std::pair<int, int> genes[] = { { 1, 1 }, { 1, 2 }, { 2, 1 } };
            std::string sequences[] = { "ACGTACGTAA", "ccggttaa", "GGGGCCCCAT" };
            for (int i = 0; i < 3; i++) {
                fna_os << '>' << genes[i].first << '_' << genes[i].second << '\n';
                size_t start = fna_os.str().size();
                fna_os << sequences[i] << '\n';
                map_os << genes[i].first << '\t' << genes[i].second << '\t' << start << '\t' << start + sequences[i].size() << '\n';
            }
            fna = fna_os.str();
            map = map_os.str();
        }
    };

    void LoadMap(std::string const& fna, std::string const& map) { protal::GenomeLoader loader(fna, map); }
}

TEST(ReferenceMap, LoadsGenesAndUppercasesThem) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    protal::GenomeLoader loader(fna, dir.Write("reference.map", ref.map + "\n"));
    EXPECT_EQ(loader.GeneCount(), 3u);
    EXPECT_TRUE(loader.HasGene(1, 2));
    EXPECT_FALSE(loader.HasGene(2, 2));
    EXPECT_FALSE(loader.HasGene(3, 1));
    EXPECT_EQ(loader.GeneLength(1, 2), 8u);

    auto& gene = loader.GetGenome(1).GetGeneOMP(2);
    EXPECT_EQ(gene.Sequence(), "CCGGTTAA");
}

// LoadAllGenomes puts the genes not yet loaded into one arena; a genome loaded gene by gene
// before keeps its sequences.
TEST(ReferenceMap, LoadAllGenomesAfterAGenomeWasLoadedOnItsOwn) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);
    for (int threads : { 1, 3 }) {
        protal::GenomeLoader loader(fna, map);
        EXPECT_EQ(loader.GetGenome(1).GetGeneOMP(2).Sequence(), "CCGGTTAA");
        loader.LoadAllGenomes(threads);
        EXPECT_TRUE(loader.AllGenomesLoaded());
        EXPECT_EQ(loader.GetGenome(1).GetGene(1).Sequence(), "ACGTACGTAA") << threads << " threads";
        EXPECT_EQ(loader.GetGenome(1).GetGene(2).Sequence(), "CCGGTTAA") << threads << " threads";
        EXPECT_EQ(loader.GetGenome(2).GetGene(1).Sequence(), "GGGGCCCCAT") << threads << " threads";

        protal::GenomeLoader fresh(fna, map);
        fresh.LoadAllGenomes(threads);
        EXPECT_EQ(fresh.GetGenome(1).GetGene(2).Sequence(), "CCGGTTAA") << threads << " threads";
        EXPECT_EQ(fresh.GetGenome(2).GetGene(1).Sequence(), "GGGGCCCCAT") << threads << " threads";
    }
}

TEST(ReferenceMap, RejectsMalformedLines) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto size = std::to_string(ref.fna.size());
    EXPECT_EXIT(LoadMap(fna, dir.Write("m1", ref.map + "3\t1\t5\n")), testing::ExitedWithCode(8), "line 4: expected 4 tab-separated columns, found 3");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m2", ref.map + "3\t1\t5\tx\n")), testing::ExitedWithCode(8), "column 4 is not");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m3", ref.map + "0\t1\t5\t8\n")), testing::ExitedWithCode(8), "taxid must be between 1 and");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m4", ref.map + "3\t0\t5\t8\n")), testing::ExitedWithCode(8), "gene id must be between 1 and");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m5", ref.map + "3\t1\t8\t8\n")), testing::ExitedWithCode(8), "end byte must be after start byte");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m6", ref.map + "3\t1\t5\t" + size + "0\n")), testing::ExitedWithCode(8), "is past the end of");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m7", ref.map + "1\t2\t5\t8\n")), testing::ExitedWithCode(8), "gene 1_2 is listed twice");
    EXPECT_EXIT(LoadMap(fna, (dir.path / "none.map").string()), testing::ExitedWithCode(8), "cannot open the file");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m8", "\n")), testing::ExitedWithCode(8), "the file lists no genes");
}

TEST(ReferenceFingerprint, TracksMapContentAndReferenceSize) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);
    auto original = protal::ReferenceFingerprint::Of(map, fna);
    EXPECT_EQ(original.fna_size, ref.fna.size());
    EXPECT_EQ(original, protal::ReferenceFingerprint::Of(map, fna));

    dir.Write("reference.map", ref.map + "3\t1\t5\t8\n");
    EXPECT_NE(original, protal::ReferenceFingerprint::Of(map, fna));
    dir.Write("reference.map", ref.map);
    dir.Write("reference.fna", ref.fna + "A");
    EXPECT_NE(original, protal::ReferenceFingerprint::Of(map, fna));
}

TEST(KmerUtils, LowercaseBasesEncodeLikeUppercase) {
    for (auto [upper, lower] : { std::pair{ 'A', 'a' }, { 'C', 'c' }, { 'G', 'g' }, { 'T', 't' } }) {
        EXPECT_EQ(KmerUtils::BaseToInt(lower), KmerUtils::BaseToInt(upper));
        EXPECT_EQ(KmerUtils::BaseToIntC(lower), KmerUtils::BaseToIntC(upper));
        EXPECT_EQ(KmerUtils::BaseToInt(lower, 0), KmerUtils::BaseToInt(upper, 0));
        EXPECT_EQ(KmerUtils::BaseToIntC(lower, 0), KmerUtils::BaseToIntC(upper, 0));
    }
    EXPECT_EQ(KmerUtils::BaseToInt('n'), 4u);
}

TEST(UniqueKmers, RejectsMalformedLines) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);
    auto load = [&](std::string const& content) {
        protal::GenomeLoader loader(fna, map);
        loader.LoadUniqueKmers(dir.Write("bad.tsv", content));
    };
    EXPECT_EXIT(load("1\t1\t3\n"), testing::ExitedWithCode(8), "line 1: expected 9 tab-separated columns, found 3");
    EXPECT_EXIT(load("1\t1\tx\t0\t0\t0\t0\t0\t10\n"), testing::ExitedWithCode(8), "column 3 is not");
    EXPECT_EXIT(load("7\t1\t3\t0.3\t0\t0\t0\t0\t10\n"), testing::ExitedWithCode(8), "gene 7_1 is not in reference.map");
    EXPECT_EXIT(load(""), testing::ExitedWithCode(8), "the file lists no genes .rebuild the database with --build.");
}

namespace {
    // A reference.map of 600 taxa x 100 genes (~1.7 MB: several chunks of the parallel parser) over a
    // sparse reference.fna, and a unique_kmers.tsv for it. `edit` changes lines (0-based) by index.
    struct BigTables {
        std::string fna, map, unique;
        size_t genes = 60000;
        // interleaved: the lines gene-major (gene 1 of every genome, then gene 2, ...), so that a genome's rows are
        // spread over the whole file and every chunk, as protal's tables never are. taxa, per_genome: the table's size.
        BigTables(ScratchDir const& dir, std::map<size_t, std::string> const& map_edits = {},
                  std::map<size_t, std::string> const& unique_edits = {}, bool interleaved = false,
                  int taxa = 600, int per_genome = 100) {
            genes = static_cast<size_t>(taxa) * static_cast<size_t>(per_genome);
            std::string m, u;
            uint64_t position = 0;
            size_t line = 0;
            m.reserve(genes * 36);
            u.reserve(genes * 40);
            for (size_t index = 0; index < genes; index++, line++) {
                int const taxid = interleaved ? 1 + static_cast<int>(index % static_cast<size_t>(taxa)) : 1 + static_cast<int>(index / static_cast<size_t>(per_genome));
                int const gene = interleaved ? 1 + static_cast<int>(index / static_cast<size_t>(taxa)) : 1 + static_cast<int>(index % static_cast<size_t>(per_genome));
                {
                    uint64_t const length = 900 + (taxid * 7 + gene) % 300;
                    std::string map_line = std::to_string(taxid) + '\t' + std::to_string(gene) + '\t' + std::to_string(position) + '\t' + std::to_string(position + length);
                    std::string unique_line = std::to_string(taxid) + '\t' + std::to_string(gene) + '\t' + std::to_string(gene % 3) + "\t0.1\t" +
                                              std::to_string(gene % 2) + "\t0.1\t1\t0.1\t" + std::to_string(length - 30);
                    if (auto it = map_edits.find(line); it != map_edits.end()) map_line = it->second;
                    if (auto it = unique_edits.find(line); it != unique_edits.end()) unique_line = it->second;
                    m += map_line + '\n';
                    u += unique_line + '\n';
                    position += length + 1;
                }
            }
            fna = dir.Write("big.fna", "");
            std::filesystem::resize_file(fna, position);
            map = dir.Write("big.map", m);
            unique = dir.Write("big_unique.tsv", u);
        }
    };

    protal::GenomeLoader LoadBig(BigTables const& t, int threads) {
        protal::GenomeLoader loader(protal::db::DbFile::OnDisk(t.fna), protal::db::DbFile::OnDisk(t.map), threads);
        loader.LoadUniqueKmers(t.unique, threads);
        return loader;
    }
}

TEST(GeneTables, ParsedInParallelAsLineByLine) {
    ScratchDir dir;
    BigTables tables(dir);
    ASSERT_GT(std::filesystem::file_size(tables.map), size_t{1} << 20);  // several chunks
    auto one = LoadBig(tables, 1);
    auto eight = LoadBig(tables, 8);
    EXPECT_EQ(one.GeneCount(), tables.genes);
    EXPECT_EQ(eight.GeneCount(), tables.genes);
    for (int taxid : { 1, 299, 600 }) {
        for (int gene : { 1, 50, 100 }) {
            EXPECT_EQ(eight.GeneLength(taxid, gene), one.GeneLength(taxid, gene));
            EXPECT_EQ(eight.GetGenome(taxid).GetGene(gene).GetStartByte(), one.GetGenome(taxid).GetGene(gene).GetStartByte());
            EXPECT_EQ(eight.GetGenome(taxid).GetGene(gene).GetUniqueKmerCounts(), one.GetGenome(taxid).GetGene(gene).GetUniqueKmerCounts());
            EXPECT_EQ(eight.GetGenome(taxid).IsGeneHittable(gene), one.GetGenome(taxid).IsGeneHittable(gene));
        }
        EXPECT_EQ(eight.GetGenome(taxid).GetUniqueKmerCounts(), one.GetGenome(taxid).GetUniqueKmerCounts());
    }
}

TEST(GeneTables, TheFirstProblemInTheFileIsReportedWithItsLine) {
    ScratchDir dir;
    // A problem deep in the file (a late chunk) keeps its line number.
    BigTables late(dir, { { 55000, "7\t8\t9" } });
    EXPECT_EXIT(LoadBig(late, 8), testing::ExitedWithCode(8), "line 55001: expected 4 tab-separated columns, found 3");
    // A gene listed twice early wins over a malformed line later, and the other way round.
    BigTables duplicate_first(dir, { { 1000, "1\t1\t0\t10" }, { 50000, "x" } });
    EXPECT_EXIT(LoadBig(duplicate_first, 8), testing::ExitedWithCode(8), "line 1001: gene 1_1 is listed twice");
    BigTables malformed_first(dir, { { 1000, "1\t1\tx\t10" }, { 50000, "1\t1\t0\t10" } });
    EXPECT_EXIT(LoadBig(malformed_first, 8), testing::ExitedWithCode(8), "line 1001: column 3 is not a non-negative integer");
    BigTables huge(dir, { { 3, "1\t4\t0\t99999999999999999999999" } });
    EXPECT_EXIT(LoadBig(huge, 8), testing::ExitedWithCode(8), "line 4: column 4 is too large");
    // unique_kmers.tsv likewise.
    BigTables unknown_gene(dir, {}, { { 40000, "1\t700\t1\t0.1\t0\t0.1\t1\t0.1\t9" } });
    EXPECT_EXIT(LoadBig(unknown_gene, 8), testing::ExitedWithCode(8), "line 40001: gene 1_700 is not in reference.map");
    BigTables short_line(dir, {}, { { 59999, "1\t2" } });
    EXPECT_EXIT(LoadBig(short_line, 8), testing::ExitedWithCode(8), "line 60000: expected 9 tab-separated columns, found 2");
}

// The rows are added by genome, each genome's by one thread (gene_table::AddPieceByGenome). A table whose genomes
// are spread over the file (gene-major order, as gtdb_to_protal_db.py writes reference.map by default) gives every genome a run in every chunk:
// the same genes as on one thread, and a gene listed twice far from its genome's other rows is still caught.
TEST(GeneTables, GenomesSpreadOverTheFileAreAddedOnceEach) {
    ScratchDir dir;
    BigTables tables(dir, {}, {}, true);
    auto one = LoadBig(tables, 1);
    auto eight = LoadBig(tables, 8);
    EXPECT_EQ(one.GeneCount(), tables.genes);
    EXPECT_EQ(eight.GeneCount(), tables.genes);
    for (int taxid = 1; taxid <= 600; taxid += 7) {
        for (int gene = 1; gene <= 100; gene += 9) {
            ASSERT_EQ(eight.GeneLength(taxid, gene), one.GeneLength(taxid, gene));
            ASSERT_EQ(eight.GetGenome(taxid).GetGene(gene).GetStartByte(), one.GetGenome(taxid).GetGene(gene).GetStartByte());
            ASSERT_EQ(eight.GetGenome(taxid).GetGene(gene).GetUniqueKmerCounts(), one.GetGenome(taxid).GetGene(gene).GetUniqueKmerCounts());
            ASSERT_EQ(eight.GetGenome(taxid).IsGeneHittable(gene), one.GetGenome(taxid).IsGeneHittable(gene));
        }
        ASSERT_EQ(eight.GetGenome(taxid).GetUniqueKmerCounts(), one.GetGenome(taxid).GetUniqueKmerCounts());
    }
    // Line 30001 (0-based 30000) is genome 1's gene 51; made a second gene 1_1, far from gene 1_1's own line.
    BigTables duplicate(dir, { { 30000, "1\t1\t0\t10" } }, {}, true);
    EXPECT_EXIT(LoadBig(duplicate, 8), testing::ExitedWithCode(8), "line 30001: gene 1_1 is listed twice");
}

TEST(UniqueKmers, AGenomeWithoutUniqueKmersHasNoHittableGene) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);

    protal::GenomeLoader without(fna, map);  // no unique_kmers.tsv: every gene counts
    EXPECT_TRUE(without.GetGenome(2).IsGeneHittable(1));
    EXPECT_EQ(without.GetGenome(2).GeneNum(), 1u);

    // Genome 2's gene has no unique k-mers, and genome 1 has one with and one without. The first line ends in CRLF and a
    // blank line follows it, as an editor may leave them.
    protal::GenomeLoader loader(fna, map);
    loader.LoadUniqueKmers(dir.Write("u.tsv", "1\t1\t3\t0.3\t0\t0\t0\t0\t10\r\n\n1\t2\t0\t0\t0\t0\t0\t0\t8\n2\t1\t0\t0\t0\t0\t0\t0\t8\n"));
    EXPECT_EQ(loader.GetGenome(1).GetHittableGenes(), std::vector<uint32_t>{ 1 });
    EXPECT_FALSE(loader.GetGenome(2).IsGeneHittable(1));
    EXPECT_TRUE(loader.GetGenome(2).GetHittableGenes().empty());
    EXPECT_EQ(loader.GetGenome(2).GeneNum(), 0u);
}

namespace {
    // n reads of `length` random bases in FASTQ (or FASTA without quality), named name(i), every base of quality q
    // (Phred+33 character) but every tenth of quality low where low is given.
    std::string Reads(size_t n, size_t length, char q, std::function<std::string(size_t)> const& name, bool fasta = false, char low = 0) {
        std::mt19937 rng(5);
        std::string out;
        for (size_t i = 0; i < n; i++) {
            std::string seq(length, 'A'), qual(length, q);
            for (auto& c : seq) c = "ACGT"[rng() % 4];
            if (low) for (size_t k = 0; k < length; k += 10) qual[k] = low;
            out += (fasta ? ">" : "@") + name(i) + "\n" + seq + "\n";
            if (!fasta) out += "+\n" + qual + "\n";
        }
        return out;
    }
    std::string Plain(size_t i) { return "read" + std::to_string(i); }

    std::optional<protal::ReadTypeGuess> Guess(std::string const& path) {
        return protal::GuessReadType(protal::SampleReads(path), protal::MAX_SHORT_READ_LENGTH);
    }
}

// The kind of reads in a single read file, from its first reads (ReadTypeDetection.h).
TEST(ReadTypeDetection, ShortAndLongReadsByLength) {
    ScratchDir dir;
    auto const short_reads = Guess(dir.Write("short.fq", Reads(300, 150, 'I', Plain)));
    ASSERT_TRUE(short_reads);
    EXPECT_EQ(short_reads->type, protal::ReadType::Single);
    // A stray long read among short ones leaves them short reads (the reader then stops at it).
    auto const stray = Guess(dir.Write("stray.fq", Reads(150, 150, 'I', Plain) + Reads(1, 2000, 'I', Plain) + Reads(150, 150, 'I', Plain)));
    ASSERT_TRUE(stray);
    EXPECT_EQ(stray->type, protal::ReadType::Single);
    // Nothing to read: no guess (a missing file, or a pipe, whose reads would be lost to the alignment).
    EXPECT_FALSE(Guess((dir.path / "missing.fq").string()));
    EXPECT_FALSE(Guess(dir.Write("empty.fq", "")));
}

TEST(ReadTypeDetection, LongReadsByTheirQuality) {
    ScratchDir dir;
    auto const ont = Guess(dir.Write("q18.fq", Reads(50, 4000, '3', Plain)));  // Q18
    ASSERT_TRUE(ont);
    EXPECT_EQ(ont->type, protal::ReadType::ONT);
    EXPECT_NE(ont->evidence.find("median read quality Q18.0"), std::string::npos) << ont->evidence;
    EXPECT_NE(ont->evidence.find("the first 50 reads up to 4.0 kb long, median 4.0 kb"), std::string::npos) << ont->evidence;
    auto const hifi = Guess(dir.Write("q35.fq", Reads(50, 4000, 'D', Plain)));  // Q35
    ASSERT_TRUE(hifi);
    EXPECT_EQ(hifi->type, protal::ReadType::PacBio);
    // A read's quality is its bases' mean error probability: a tenth of the bases at Q5 make Q40 reads Q15 (the mean of
    // the Phred values would be Q36.5).
    auto const mixed = Guess(dir.Write("mixed.fq", Reads(50, 4000, 'I', Plain, false, '&')));
    ASSERT_TRUE(mixed);
    EXPECT_EQ(mixed->type, protal::ReadType::ONT);
    EXPECT_NE(mixed->evidence.find("Q15.0"), std::string::npos) << mixed->evidence;
    // Without qualities: PacBio, and the evidence says why.
    auto const fasta = Guess(dir.Write("long.fa", Reads(50, 4000, 'I', Plain, true)));
    ASSERT_TRUE(fasta);
    EXPECT_EQ(fasta->type, protal::ReadType::PacBio);
    EXPECT_NE(fasta->evidence.find("without base qualities"), std::string::npos) << fasta->evidence;
    // Q0 at every base means unknown (pbsim3's PacBio reads): as without qualities.
    auto const q0 = Guess(dir.Write("q0.fq", Reads(50, 4000, '!', Plain)));
    ASSERT_TRUE(q0);
    EXPECT_EQ(q0->type, protal::ReadType::PacBio);
    EXPECT_NE(q0->evidence.find("(FASTA, or Q0 throughout)"), std::string::npos) << q0->evidence;
    // Gzipped reads read as plain ones.
    std::string const gz = (dir.path / "q18.fq.gz").string();
    {
        ogzstream os(gz.c_str());
        os << Reads(50, 4000, '3', Plain);
    }
    auto const zipped = Guess(gz);
    ASSERT_TRUE(zipped);
    EXPECT_EQ(zipped->type, protal::ReadType::ONT);
    // zstd too.
    std::string const text = Reads(50, 4000, 'D', Plain);
    std::string zst(ZSTD_compressBound(text.size()), '\0');
    zst.resize(ZSTD_compress(zst.data(), zst.size(), text.data(), text.size(), 3));
    auto const zstd = Guess(dir.Write("q35.fq.zst", zst));
    ASSERT_TRUE(zstd);
    EXPECT_EQ(zstd->type, protal::ReadType::PacBio);
}

TEST(ReadTypeDetection, NamesBeforeQualities) {
    ScratchDir dir;
    // MinKNOW: a UUID, and runid=; dorado: a UUID. ONT reads, whatever their quality (duplex reads reach Q30).
    auto minknow = [](size_t i) {
        char id[64];
        std::snprintf(id, sizeof id, "%08zx-1c2d-4e5f-8a9b-0123456789ab runid=6f3a ch=12 start_time=2024-01-01T00:00:00Z", i);
        return std::string(id);
    };
    auto const ont = Guess(dir.Write("minknow.fq", Reads(40, 3000, 'D', minknow)));
    ASSERT_TRUE(ont);
    EXPECT_EQ(ont->type, protal::ReadType::ONT);
    EXPECT_NE(ont->evidence.find("named as MinKNOW and dorado name reads"), std::string::npos) << ont->evidence;
    // PacBio: movie/ZMW, also with one underscore in the movie name and more after the ZMW.
    auto const ccs = Guess(dir.Write("ccs.fq", Reads(40, 3000, '0', [](size_t i) { return "m64011_190830_220126/" + std::to_string(i) + "/ccs"; })));
    ASSERT_TRUE(ccs);
    EXPECT_EQ(ccs->type, protal::ReadType::PacBio);
    auto const revio = Guess(dir.Write("revio.fq", Reads(40, 3000, '0', [](size_t i) { return "m84011_220902_175841_s1/" + std::to_string(i) + "/ccs/fwd"; })));
    ASSERT_TRUE(revio);
    EXPECT_EQ(revio->type, protal::ReadType::PacBio);
    auto const one_underscore = Guess(dir.Write("m.fq", Reads(40, 3000, '0', [](size_t i) { return "m64001_000000/" + std::to_string(i) + "/ccs"; })));
    ASSERT_TRUE(one_underscore);
    EXPECT_EQ(one_underscore->type, protal::ReadType::PacBio);
    // Names like neither: by quality ('0' is Q15).
    auto const other = Guess(dir.Write("other.fq", Reads(40, 3000, '0', [](size_t i) { return "m_" + std::to_string(i); })));
    ASSERT_TRUE(other);
    EXPECT_EQ(other->type, protal::ReadType::ONT);
    EXPECT_EQ(other->basis, protal::ReadTypeBasis::Quality);
    EXPECT_EQ(ont->basis, protal::ReadTypeBasis::Names);
    EXPECT_EQ(ccs->basis, protal::ReadTypeBasis::Names);
}

// A UUID anywhere in the header, between characters that are no letters or digits: after an archive's own id (SRA,
// ENA), in a dorado duplex read's "uuid;uuid", before a suffix; not inside a longer word or run of hex digits.
TEST(ReadTypeDetection, AUuidAnywhereInTheHeader) {
    using protal::read_type_detection::HoldsUuid;
    std::string const uuid = "1b2c3d4e-1c2d-4e5f-8a9b-0123456789ab", other = "0f9e8d7c-6b5a-4f3e-9d2c-1b0a98765432";
    EXPECT_TRUE(HoldsUuid(uuid));
    EXPECT_TRUE(HoldsUuid("SRR12345678.1 " + uuid + " length=5000"));
    EXPECT_TRUE(HoldsUuid("ERR12345.7\t" + uuid));
    EXPECT_TRUE(HoldsUuid(uuid + ";" + other));
    EXPECT_TRUE(HoldsUuid(uuid + "_Basecall_1D_template"));
    EXPECT_TRUE(HoldsUuid("read=" + uuid + ","));
    EXPECT_FALSE(HoldsUuid("x" + uuid));
    EXPECT_FALSE(HoldsUuid(uuid + "0"));
    EXPECT_FALSE(HoldsUuid("SRR12345678.1 1 length=5000"));
    EXPECT_FALSE(HoldsUuid(uuid.substr(0, 35)));
    EXPECT_FALSE(HoldsUuid("1b2c3d4e01c2d-4e5f-8a9b-0123456789ab"));
}

// ONT reads from SRA or ENA keep their name, if at all, after the archive's id; dorado names a duplex read uuid;uuid.
// Their names make them ONT reads whatever their quality ('D' is Q35); a PacBio name after an archive's id makes PacBio
// reads.
TEST(ReadTypeDetection, ArchivedReadsByTheNamesTheyKeep) {
    ScratchDir dir;
    auto uuid = [](size_t i) {
        char id[40];
        std::snprintf(id, sizeof id, "%08zx-1c2d-4e5f-8a9b-0123456789ab", i);
        return std::string(id);
    };
    std::vector<std::pair<std::string, std::function<std::string(size_t)>>> const ont_names = {
            { "sra", [&](size_t i) { return "SRR12345678." + std::to_string(i + 1) + " " + uuid(i) + " length=3000"; } },
            { "ena_runid", [&](size_t i) { return "ERR1234567." + std::to_string(i + 1) + " " + uuid(i) + " runid=6f3a ch=12"; } },
            { "runid_only", [](size_t i) { return "SRR1." + std::to_string(i + 1) + " sampleid=x,runid=6f3a"; } },
            { "duplex", [&](size_t i) { return uuid(i) + ";" + uuid(i + 1000); } },
            { "albacore", [&](size_t i) { return uuid(i) + "_Basecall_1D_template"; } },
    };
    for (auto const& [name, read_name] : ont_names) {
        auto const guess = Guess(dir.Write(name + ".fq", Reads(40, 3000, 'D', read_name)));
        ASSERT_TRUE(guess) << name;
        EXPECT_EQ(guess->type, protal::ReadType::ONT) << name;
        EXPECT_EQ(guess->basis, protal::ReadTypeBasis::Names) << name;
    }
    // Lost names: by quality (Q35, PacBio's), which DescribeReadType makes a warning.
    auto const lost = Guess(dir.Write("lost.fq", Reads(40, 3000, 'D', [](size_t i) { return "SRR1." + std::to_string(i + 1) + " length=3000"; })));
    ASSERT_TRUE(lost);
    EXPECT_EQ(lost->type, protal::ReadType::PacBio);
    EXPECT_EQ(lost->basis, protal::ReadTypeBasis::Quality);
    // A PacBio name after an archive's id, Q15 reads ('0'): PacBio reads by the name.
    auto const ccs = Guess(dir.Write("sra_ccs.fq", Reads(40, 3000, '0', [](size_t i) {
        return "SRR9087597." + std::to_string(i + 1) + " m64011_190830_220126/" + std::to_string(i) + "/ccs length=3000";
    })));
    ASSERT_TRUE(ccs);
    EXPECT_EQ(ccs->type, protal::ReadType::PacBio);
    EXPECT_EQ(ccs->basis, protal::ReadTypeBasis::Names);
}

// Read by read, Q0 at every base is no quality, as the alignment's reader takes it: of 26 reads at Q0 and 24 at Q35,
// most have no qualities, so the reads are taken for PacBio's for want of them (the median of all 50 was Q0, ONT).
TEST(ReadTypeDetection, ReadsAtQ0AreReadsWithoutQualities) {
    ScratchDir dir;
    auto const guess = Guess(dir.Write("mixed_q0.fq", Reads(26, 4000, '!', Plain) + Reads(24, 4000, 'D', Plain)));
    ASSERT_TRUE(guess);
    EXPECT_EQ(guess->type, protal::ReadType::PacBio);
    EXPECT_EQ(guess->basis, protal::ReadTypeBasis::NoQualities);
    // With most reads at Q35, their median decides.
    auto const most = Guess(dir.Write("few_q0.fq", Reads(20, 4000, '!', Plain) + Reads(30, 4000, '3', Plain)));
    ASSERT_TRUE(most);
    EXPECT_EQ(most->type, protal::ReadType::ONT);
    EXPECT_NE(most->evidence.find("median read quality Q18.0"), std::string::npos) << most->evidence;
}

// What protal says of each guess (DescribeReadType): a note where names decided, a warning that names the kind, why and
// --read_type where the base quality, or the lack of it, did; for single-end reads nothing, unless they look like long
// reads (a median length beyond Illumina's and a quality below PacBio's), which a note says.
TEST(ReadTypeDetection, AWarningWhereQualitiesAloneDecide) {
    auto describe = [](protal::ReadSample sample) {
        auto const guess = protal::GuessReadType(sample, protal::MAX_SHORT_READ_LENGTH);
        return protal::DescribeReadType(*guess, "s", "s.fq");
    };
    auto sample = [](std::vector<size_t> const& lengths, std::optional<double> quality, size_t ont_names = 0) {
        protal::ReadSample s;
        s.reads = lengths.size();
        s.lengths = lengths;
        s.longest = *std::max_element(lengths.begin(), lengths.end());
        s.median_length = protal::read_type_detection::Median(lengths);
        s.median_quality = quality;
        s.ont_names = ont_names;
        return s;
    };
    std::string const set_type = "--read_type ont or pb (or a map's READ_TYPE) sets the kind of reads";

    auto const named = describe(sample({ 4000, 4000 }, 35.0, 2));
    ASSERT_TRUE(named);
    EXPECT_FALSE(named->warning);
    EXPECT_NE(named->text.find("Sample s: ONT reads (the first 2 reads up to 4.0 kb long, median 4.0 kb, named as MinKNOW and "
                               "dorado name reads in s.fq), aligned and profiled as such; --read_type (or a map's READ_TYPE) "
                               "sets the kind of reads"), std::string::npos) << named->text;

    auto const pacbio = describe(sample({ 4000, 4000 }, 35.0));
    ASSERT_TRUE(pacbio);
    EXPECT_TRUE(pacbio->warning);
    EXPECT_NE(pacbio->text.find("Sample s: PacBio reads (the first 2 reads up to 4.0 kb long, median 4.0 kb, median read quality "
                                "Q35.0 (PacBio from Q25, ONT below) in s.fq), aligned and profiled as such by their base quality "
                                "alone"), std::string::npos) << pacbio->text;
    EXPECT_NE(pacbio->text.find("as after a download from SRA or ENA"), std::string::npos) << pacbio->text;
    EXPECT_NE(pacbio->text.find("ONT reads reach PacBio's Q25"), std::string::npos) << pacbio->text;
    EXPECT_NE(pacbio->text.find(set_type), std::string::npos) << pacbio->text;

    auto const ont = describe(sample({ 4000, 4000 }, 18.0));
    ASSERT_TRUE(ont);
    EXPECT_TRUE(ont->warning);
    EXPECT_NE(ont->text.find("Sample s: ONT reads ("), std::string::npos) << ont->text;
    EXPECT_NE(ont->text.find("by their base quality alone"), std::string::npos) << ont->text;
    EXPECT_EQ(ont->text.find("reach PacBio's Q25"), std::string::npos) << ont->text;
    EXPECT_NE(ont->text.find(set_type), std::string::npos) << ont->text;

    auto const fasta = describe(sample({ 4000, 4000 }, std::nullopt));
    ASSERT_TRUE(fasta);
    EXPECT_TRUE(fasta->warning);
    EXPECT_NE(fasta->text.find("Sample s: PacBio reads ("), std::string::npos) << fasta->text;
    EXPECT_NE(fasta->text.find("aligned and profiled as such for want of base qualities"), std::string::npos) << fasta->text;
    EXPECT_NE(fasta->text.find(set_type), std::string::npos) << fasta->text;

    // Single-end reads: nothing to say of Illumina's; long-read-like ones (median 600 bp, Q14) get a note.
    EXPECT_FALSE(describe(sample({ 150, 150, 150 }, 35.0)));
    EXPECT_FALSE(describe(sample({ 600, 600, 600 }, 35.0)));          // merged Illumina pairs
    EXPECT_FALSE(describe(sample({ 300, 300, 300 }, 14.0)));          // Illumina's length
    EXPECT_FALSE(describe(sample({ 600, 600, 600 }, std::nullopt)));  // no qualities to tell
    auto const long_like = describe(sample({ 600, 600, 600, 2000 }, 14.0));
    ASSERT_TRUE(long_like);
    EXPECT_FALSE(long_like->warning);
    EXPECT_NE(long_like->text.find("Sample s: single-end reads (the first 4 reads up to 2.0 kb long, median 0.6 kb, at most a "
                                   "quarter over 1.0 kb, median read quality Q14.0 in s.fq), aligned and profiled as such, "
                                   "though their lengths and qualities are those of long reads"), std::string::npos)
            << long_like->text;
    EXPECT_NE(long_like->text.find("--read_type ont or pb (or a map's READ_TYPE) for ONT or PacBio reads"), std::string::npos);
}

// At the boundaries: a median read quality of Q25 (kPacBioMinQuality) is PacBio's, a read of MAX_SHORT_READ_LENGTH bases
// is a short read, and a quarter of the reads longer leaves the sample's reads short.
TEST(ReadTypeDetection, TheBoundariesOfQualityAndLength) {
    auto type = [](std::vector<size_t> const& lengths, double quality) {
        protal::ReadSample sample;
        sample.reads = lengths.size();
        sample.lengths = lengths;
        sample.longest = sample.median_length = *std::max_element(lengths.begin(), lengths.end());
        sample.median_quality = quality;
        return protal::GuessReadType(sample, protal::MAX_SHORT_READ_LENGTH)->type;
    };
    double const q25 = protal::kPacBioMinQuality;
    EXPECT_EQ(type({ 4000, 4000 }, q25), protal::ReadType::PacBio);
    EXPECT_EQ(type({ 4000, 4000 }, std::nextafter(q25, 0.0)), protal::ReadType::ONT);
    size_t const longest_short = protal::MAX_SHORT_READ_LENGTH;
    EXPECT_EQ(type({ longest_short, longest_short }, 18), protal::ReadType::Single);
    EXPECT_EQ(type({ longest_short + 1, longest_short + 1 }, 18), protal::ReadType::ONT);
    EXPECT_EQ(type({ 150, 150, 150, 4000 }, 18), protal::ReadType::Single);
    EXPECT_EQ(type({ 150, 150, 4000 }, 18), protal::ReadType::ONT);
}

// The reads looked at: the first 200, or fewer where they reach 5 Mb, as set (the read that reaches the bases is the last).
TEST(ReadTypeDetection, TheFirstReadsUpToTheCaps) {
    ScratchDir dir;
    auto const path = dir.Write("reads.fq", Reads(300, 150, 'I', Plain));
    EXPECT_EQ(protal::SampleReads(path).reads, 200u);
    EXPECT_EQ(protal::SampleReads(path, 10).reads, 10u);
    auto const capped = protal::SampleReads(path, 200, 1000);
    EXPECT_EQ(capped.reads, 7u);  // 1,050 bases
    EXPECT_EQ(capped.lengths, std::vector<size_t>(7, 150));
}

// The preload's sort of its genes by start byte (gene_table::ParallelSort): what std::sort gives, on any number of threads,
// for sizes below, at and above a part per thread (2^16 values), up to 5 parts (merged over three rounds, a run left over
// in two), and for genes listed genome by genome from a reference written gene by gene.
TEST(GeneTables, TheParallelSortIsStdSort) {
    std::mt19937_64 rng(11);
    auto by_first = [](auto const& a, auto const& b) { return a.first < b.first; };
    for (size_t const n : { size_t{0}, size_t{1}, size_t{1000}, size_t{1} << 16, (size_t{1} << 18) + 7, (size_t{5} << 16) + 3 }) {
        std::vector<std::pair<uint64_t, uint32_t>> values(n);
        for (size_t i = 0; i < n; i++) values[i] = { rng(), static_cast<uint32_t>(i) };
        auto expected = values;
        std::sort(expected.begin(), expected.end(), by_first);
        for (int const threads : { 1, 2, 3, 7, 32 }) {
            auto sorted = values;
            protal::gene_table::ParallelSort(sorted, threads, by_first);
            EXPECT_EQ(sorted, expected) << n << " values on " << threads << " threads";
        }
    }
    // 2000 genomes of up to 120 genes, the reference gene by gene: listed genome by genome, sorted back into its order.
    std::vector<std::pair<uint64_t, uint32_t>> listed;
    std::vector<std::vector<uint64_t>> start(2000, std::vector<uint64_t>(120, UINT64_MAX));
    uint64_t offset = 0;
    for (uint32_t gene = 0; gene < 120; gene++) {
        for (uint32_t genome = 0; genome < 2000; genome++) {
            if (rng() % 100 < 84) start[genome][gene] = (offset += 900 + rng() % 400);
        }
    }
    for (uint32_t genome = 0; genome < 2000; genome++) {
        for (uint32_t gene = 0; gene < 120; gene++) {
            if (start[genome][gene] != UINT64_MAX) listed.emplace_back(start[genome][gene], genome * 120 + gene);
        }
    }
    auto expected = listed;
    std::sort(expected.begin(), expected.end(), by_first);
    protal::gene_table::ParallelSort(listed, 6, by_first);
    EXPECT_EQ(listed, expected);
}
