// Unit tests for the binary gene table (GeneTableFile.h; GenomeLoader::WriteGeneTable and LoadGeneTable): the genes and their
// unique k-mer counts load from it as from reference.map and unique_kmers.tsv, with the same reference fingerprint; the same
// tables give the same bytes; a table made from other tables, or a corrupt one, is not used (the text tables are); and a
// single-file database with one loads it.
#include <gtest/gtest.h>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <unistd.h>
#include "SequenceUtils/GenomeLoader.h"
#include "SequenceUtils/GeneTableFile.h"
#include "Utilities/Database.h"
#include "Utilities/ReferenceFingerprint.h"

using namespace protal;
namespace fs = std::filesystem;

namespace {
    struct TempDir {
        fs::path path;
        TempDir() {
            path = fs::temp_directory_path() / ("protal_gene_table_test_" + std::to_string(::getpid()) + "_" +
                                                std::to_string(reinterpret_cast<uintptr_t>(this)));
            fs::create_directories(path);
        }
        ~TempDir() { std::error_code ec; fs::remove_all(path, ec); }
        std::string Write(std::string const& name, std::string const& content) const {
            std::ofstream os(path / name, std::ios::binary);
            os << content;
            return (path / name).string();
        }
    };

    // reference.map of `taxa` genomes with gene ids up to `per_genome` (every third id left out), a sparse reference.fna
    // that covers them, and unique_kmers.tsv for most genes (every seventh has no line, every fifth no unique k-mers).
    struct Tables {
        std::string fna, map, unique;
        size_t genes = 0;
        // scrambled: the genomes listed in a scrambled order (taxid 1 + i * 7 mod taxa), as a converter may list them.
        Tables(TempDir const& dir, int taxa = 300, int per_genome = 30, std::string const& suffix = "", bool scrambled = false) {
            std::string m, u;
            uint64_t position = 0;
            for (int i = 0; i < taxa; i++) {
                int const taxid = scrambled ? 1 + (i * 7) % taxa : 1 + i;
                for (int gene = 1; gene <= per_genome; gene++) {
                    if (gene % 3 == 0) continue;
                    uint64_t const length = 500 + (taxid * 13 + gene * 7) % 900;
                    m += std::to_string(taxid) + '\t' + std::to_string(gene) + '\t' + std::to_string(position) + '\t' +
                         std::to_string(position + length) + '\n';
                    if (gene % 7 != 0) {
                        int const su = gene % 5 == 0 ? 0 : 1 + gene % 4, lu = gene % 5 == 0 ? 0 : gene % 2;
                        u += std::to_string(taxid) + '\t' + std::to_string(gene) + '\t' + std::to_string(su) + "\t0.1\t" +
                             std::to_string(lu) + "\t0.1\t" + std::to_string(gene % 3) + "\t0.1\t" + std::to_string(length - 30) + '\n';
                    }
                    position += length + 1;
                    genes++;
                }
            }
            fna = dir.Write("reference" + suffix + ".fna", "");
            fs::resize_file(fna, position);
            map = dir.Write("reference" + suffix + ".map", m);
            unique = dir.Write("unique_kmers" + suffix + ".tsv", u);
        }
    };

    GenomeLoader TextLoader(Tables const& t, int threads, bool with_unique = true) {
        GenomeLoader loader(db::DbFile::OnDisk(t.fna), db::DbFile::OnDisk(t.map), threads);
        if (with_unique) loader.LoadUniqueKmers(t.unique, threads);
        return loader;
    }

    std::string WriteTable(GenomeLoader const& loader, Tables const& t, std::string const& path, bool with_unique = true) {
        return loader.WriteGeneTable(path, fs::file_size(t.map), with_unique ? fs::file_size(t.unique) : 0, with_unique);
    }

    GenomeLoader TableLoader(Tables const& t, std::string const& table, int threads, bool with_unique = true) {
        return GenomeLoader(db::DbFile::OnDisk(t.fna), db::DbFile::OnDisk(t.map), threads, db::DbFile::OnDisk(table),
                            with_unique ? std::optional<uint64_t>(fs::file_size(t.unique)) : std::nullopt);
    }

    // Every gene of a and b the same: present, start byte, length, unique k-mer counts, hittable; every genome's sums.
    void ExpectSame(GenomeLoader& a, GenomeLoader& b, int taxa = 300, int per_genome = 30) {
        ASSERT_EQ(a.GeneCount(), b.GeneCount());
        for (int taxid = 1; taxid <= taxa; taxid++) {
            auto& ga = a.GetGenome(taxid);
            auto& gb = b.GetGenome(taxid);
            ASSERT_EQ(ga.GetGeneList().size(), gb.GetGeneList().size()) << taxid;
            for (int gene = 1; gene <= per_genome + 1; gene++) {
                ASSERT_EQ(a.HasGene(taxid, gene), b.HasGene(taxid, gene)) << taxid << "_" << gene;
                if (!a.HasGene(taxid, gene)) continue;
                ASSERT_EQ(ga.GetGene(gene).GetStartByte(), gb.GetGene(gene).GetStartByte());
                ASSERT_EQ(ga.GetGene(gene).GetLength(), gb.GetGene(gene).GetLength());
                ASSERT_EQ(ga.GetGene(gene).GetUniqueKmerCounts(), gb.GetGene(gene).GetUniqueKmerCounts());
                ASSERT_EQ(ga.IsGeneHittable(gene), gb.IsGeneHittable(gene));
            }
            ASSERT_EQ(ga.GetUniqueKmerCounts(), gb.GetUniqueKmerCounts()) << taxid;
            ASSERT_EQ(ga.GetHittableGenes(), gb.GetHittableGenes()) << taxid;
            ASSERT_EQ(ga.GeneNum(), gb.GeneNum()) << taxid;
        }
        auto const [ta, na, la] = a.IndexFieldMaxima();
        auto const [tb, nb, lb] = b.IndexFieldMaxima();
        EXPECT_EQ(ta, tb);
        EXPECT_EQ(na, nb);
        EXPECT_EQ(la, lb);
    }

    std::string Slurp(std::string const& path) {
        std::ifstream is(path, std::ios::binary);
        return { std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>() };
    }
}

TEST(GeneTableFile, LoadsAsTheTextTables) {
    TempDir dir;
    Tables t(dir);
    auto text = TextLoader(t, 1);
    ASSERT_TRUE(text.MapFingerprint().has_value());
    // The fingerprint hashed while reference.map was parsed is the one the index is checked against.
    EXPECT_EQ(*text.MapFingerprint(), ReferenceFingerprint::Of(t.map, t.fna));
    EXPECT_EQ(*TextLoader(t, 6).MapFingerprint(), ReferenceFingerprint::Of(t.map, t.fna));
    std::string const table = (dir.path / gene_table_file::kFileName).string();
    ASSERT_EQ(WriteTable(text, t, table), "");
    for (int threads : { 1, 6 }) {
        auto binary = TableLoader(t, table, threads);
        ASSERT_TRUE(binary.FromGeneTable());
        EXPECT_TRUE(binary.UniqueKmersFromGeneTable());
        EXPECT_EQ(*binary.MapFingerprint(), *text.MapFingerprint());
        ExpectSame(text, binary);
    }
}

TEST(GeneTableFile, WithoutUniqueKmers) {
    TempDir dir;
    Tables t(dir);
    auto text = TextLoader(t, 4, false);
    std::string const table = (dir.path / gene_table_file::kFileName).string();
    ASSERT_EQ(WriteTable(text, t, table, false), "");
    auto binary = TableLoader(t, table, 4, false);
    ASSERT_TRUE(binary.FromGeneTable());
    EXPECT_FALSE(binary.UniqueKmersFromGeneTable());
    ExpectSame(text, binary);
    EXPECT_TRUE(binary.GetGenome(1).IsGeneHittable(1));  // not known: every gene counts
    // A table without the unique k-mer counts is not one of a database that has them.
    EXPECT_FALSE(TableLoader(t, table, 4, true).FromGeneTable());
}

TEST(GeneTableFile, TheSameTablesGiveTheSameBytes) {
    TempDir dir;
    Tables t(dir);
    std::string const a = (dir.path / "a.bin").string(), b = (dir.path / "b.bin").string();
    ASSERT_EQ(WriteTable(TextLoader(t, 1), t, a), "");
    ASSERT_EQ(WriteTable(TextLoader(t, 6), t, b), "");
    EXPECT_EQ(Slurp(a), Slurp(b));
    // The size it should have: the header, the genomes and the genes.
    gene_table_file::Header h;
    std::string const bytes = Slurp(a);
    std::memcpy(&h, bytes.data(), sizeof(h));
    EXPECT_EQ(h.genomes, 300u);
    EXPECT_EQ(h.genes, t.genes);
    EXPECT_EQ(gene_table_file::ExpectedSize(h), bytes.size());
}

TEST(GeneTableFile, AnotherDatabasesTableIsNotUsed) {
    TempDir dir;
    Tables t(dir);
    std::string const table = (dir.path / gene_table_file::kFileName).string();
    ASSERT_EQ(WriteTable(TextLoader(t, 4), t, table), "");
    // reference.map with one gene more: the table was made from another map, and the text tables are loaded.
    Tables other(dir, 301, 30, "_other");
    GenomeLoader loader(db::DbFile::OnDisk(other.fna), db::DbFile::OnDisk(other.map), 4, db::DbFile::OnDisk(table),
                        std::optional<uint64_t>(fs::file_size(other.unique)));
    EXPECT_FALSE(loader.FromGeneTable());
    EXPECT_EQ(loader.GeneCount(), other.genes);
    EXPECT_TRUE(loader.HasGene(301, 1));
    EXPECT_EQ(*loader.MapFingerprint(), ReferenceFingerprint::Of(other.map, other.fna));
    // Another unique_kmers.tsv size likewise.
    GenomeLoader unique_changed(db::DbFile::OnDisk(t.fna), db::DbFile::OnDisk(t.map), 4, db::DbFile::OnDisk(table),
                                std::optional<uint64_t>(fs::file_size(t.unique) + 1));
    EXPECT_FALSE(unique_changed.FromGeneTable());
}

TEST(GeneTableFile, ACorruptTableIsNotUsed) {
    TempDir dir;
    Tables t(dir);
    auto text = TextLoader(t, 4);
    std::string const table = (dir.path / gene_table_file::kFileName).string();
    ASSERT_EQ(WriteTable(text, t, table), "");
    std::string const good = Slurp(table);
    auto check = [&](std::string const& bytes, char const* what) {
        std::string const path = dir.Write("bad.bin", bytes);
        GenomeLoader loader(db::DbFile::OnDisk(t.fna), db::DbFile::OnDisk(t.map), 4, db::DbFile::OnDisk(path),
                            std::optional<uint64_t>(fs::file_size(t.unique)));
        EXPECT_FALSE(loader.FromGeneTable()) << what;
        EXPECT_FALSE(loader.UniqueKmersFromGeneTable()) << what;
        loader.LoadUniqueKmers(t.unique, 4);  // as ProtalDB does when the table gave no counts
        ExpectSame(text, loader);  // the text tables, loaded instead
    };
    check(good.substr(0, good.size() - 5), "truncated");
    check(good.substr(0, 40), "shorter than its header");
    std::string magic = good;
    magic[0] = 'X';
    check(magic, "magic");
    std::string version = good;
    version[8] = 2;
    check(version, "version");
    // A gene record's id set to 0 (the first gene: its id is at byte 8 of the record).
    gene_table_file::Header h;
    std::memcpy(&h, good.data(), sizeof(h));
    std::string gene_id = good;
    size_t const first_gene = sizeof(h) + h.genomes * sizeof(gene_table_file::Genome);
    std::memset(gene_id.data() + first_gene + 8, 0, 4);
    check(gene_id, "gene id 0");
    // A taxid twice (the second genome's taxid made the first's).
    std::string order = good;
    std::memcpy(order.data() + sizeof(h) + sizeof(gene_table_file::Genome), good.data() + sizeof(h), 8);
    check(order, "a taxid twice");
}

TEST(GeneTableFile, ASingleFileDatabaseLoadsIt) {
    TempDir dir;
    Tables t(dir);
    auto text = TextLoader(t, 4);
    std::string const table = (dir.path / gene_table_file::kFileName).string();
    ASSERT_EQ(WriteTable(text, t, table), "");
    std::string error;
    std::string const bundle_path = (dir.path / "database.protal").string();
    ASSERT_TRUE(db::Write(bundle_path, { { "reference.fna", t.fna }, { "reference.map", t.map }, { "unique_kmers.tsv", t.unique },
                                         { gene_table_file::kFileName, table } }, zstd::Params{ 3, 0, 2, uint64_t{1} << 16 }, error)) << error;
    auto const bundle = db::Bundle::Open(bundle_path, error);
    ASSERT_TRUE(bundle) << error;
    auto const unique = db::DbFile::InBundle(*bundle, "unique_kmers.tsv");
    GenomeLoader loader(db::DbFile::InBundle(*bundle, "reference.fna"), db::DbFile::InBundle(*bundle, "reference.map"), 4,
                        db::DbFile::InBundle(*bundle, gene_table_file::kFileName), unique.Size());
    ASSERT_TRUE(loader.FromGeneTable());
    EXPECT_TRUE(loader.UniqueKmersFromGeneTable());
    EXPECT_EQ(*loader.MapFingerprint(), ReferenceFingerprint::Of(db::DbFile::InBundle(*bundle, "reference.map"),
                                                                 db::DbFile::InBundle(*bundle, "reference.fna")));
    ExpectSame(text, loader);
    // The text tables from the single file, for comparison: the same, and the same fingerprint.
    GenomeLoader from_text(db::DbFile::InBundle(*bundle, "reference.fna"), db::DbFile::InBundle(*bundle, "reference.map"), 4);
    from_text.LoadUniqueKmers(unique, 4);
    EXPECT_FALSE(from_text.FromGeneTable());
    EXPECT_EQ(*from_text.MapFingerprint(), *loader.MapFingerprint());
    ExpectSame(from_text, loader);
}

// The genome map is built in the order reference.map lists the genomes, from the text tables and from the binary one alike:
// its iteration order is the same, which some outputs depend on (ties), and so are the runs' outputs.
TEST(GeneTableFile, TheGenomeMapIsBuiltInTheSameOrder) {
    TempDir dir;
    Tables t(dir, 300, 30, "", true);  // 300 is not a multiple of 7: every taxid once
    auto text = TextLoader(t, 4);
    std::string const table = (dir.path / gene_table_file::kFileName).string();
    ASSERT_EQ(WriteTable(text, t, table), "");
    auto binary = TableLoader(t, table, 4);
    ASSERT_TRUE(binary.FromGeneTable());
    std::vector<uint64_t> a, b;
    for (auto const& [key, genome] : text.GetGenomeMap()) a.push_back(key);
    for (auto const& [key, genome] : binary.GetGenomeMap()) b.push_back(key);
    EXPECT_EQ(a, b);
    EXPECT_EQ(text.GetGenomeMap().bucket_count(), binary.GetGenomeMap().bucket_count());
    ExpectSame(text, binary);
}
