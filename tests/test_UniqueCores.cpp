// A record's unique k-mers counted from its alignment (UniqueCores.h, since 2026-10-14): which of a gene's cores an
// alignment covers with exact k-mers, the table of every copy's cores' bits (unique_cores.bin) and a thread's counter.
#include <gtest/gtest.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "SequenceUtils/KmerIterator.h"
#include "SequenceUtils/Minimizer.h"
#include "SequenceUtils/UniqueCores.h"
#include "TestUtil.h"

using namespace protal;
namespace uc = protal::unique_cores;

namespace {
    // `n` matches.
    std::string M(size_t n) { return std::string(n, 'M'); }

    uc::Counts Count(std::vector<uint32_t> const& positions, std::vector<uint8_t> const& flags, std::string const& columns, int64_t start) {
        return uc::CountCovered(std::span<uint32_t const>(positions), [&](size_t i) { return flags[i]; }, columns, start);
    }

    // A gene as the counter reads one: its bases (Sequence) and length.
    struct FakeGene {
        std::string bases;
        std::string_view Sequence() const { return bases; }
        size_t GetLength() const { return bases.size(); }
    };
}

// A core counts where its whole 31-mer (8 bases before its first, 23 from it) lies in a run of matches: a mismatch, an
// insertion or a deletion inside it, the alignment's ends or a soft clip leave it out; its bits give the unique ones.
TEST(UniqueCores, AnAlignmentCoversTheCoresOfItsExactStretches) {
    std::vector<uint32_t> const positions{ 20, 40, 60 };  // k-mers [12, 43), [32, 63), [52, 83)
    std::vector<uint8_t> const flags{ uc::kUnique | uc::kUniqueTwo, uc::kUnique, 0 };
    auto all = Count(positions, flags, M(100), 0);
    EXPECT_EQ(all.cores, 3u);
    EXPECT_EQ(all.unique, 2u);
    EXPECT_EQ(all.unique_two, 1u);
    // A mismatch at gene position 45 lies in the second core's k-mer only.
    EXPECT_EQ(Count(positions, flags, M(45) + "X" + M(54), 0).cores, 2u);
    EXPECT_EQ(Count(positions, flags, M(45) + "X" + M(54), 0).unique, 1u);
    // An insertion after gene position 50: the second core's k-mer spans it; the third's starts after it.
    EXPECT_EQ(Count(positions, flags, M(51) + "II" + M(49), 0).cores, 2u);
    // A deletion of gene position 70 lies in the third core's k-mer.
    EXPECT_EQ(Count(positions, flags, M(70) + "D" + M(29), 0).cores, 2u);
    // From gene position 30: the first core's k-mer starts before the alignment.
    EXPECT_EQ(Count(positions, flags, M(70), 30).cores, 2u);
    // A soft clip takes no gene base: from position 10 after four clipped bases the first core is covered.
    EXPECT_EQ(Count(positions, flags, "SSSS" + M(90), 10).cores, 3u);
    // Ending at gene position 82: the third core's k-mer [52, 83) is cut.
    EXPECT_EQ(Count(positions, flags, M(82), 0).cores, 2u);
    EXPECT_EQ(Count(positions, flags, M(83), 0).cores, 3u);
    EXPECT_EQ(Count(positions, flags, "", 0).cores, 0u);
}

// The table: a row per gene id of each taxon (first_row), each row's cores from a byte; rows read back, a gene beyond
// its taxon's rows or a taxid beyond the table has none; the file round-trips and a damaged one is refused.
TEST(UniqueCores, TheTableHoldsEveryRowsBits) {
    std::vector<uint64_t> const first_row{ 0, 0, 2, 3 };  // taxon 0 no genes, taxon 1 genes 1-2, taxon 2 gene 1
    std::vector<uint32_t> const cores{ 5, 0, 9 };
    ASSERT_TRUE(uc::Table::Fits(cores));
    uc::Table table(first_row, cores);
    EXPECT_EQ(table.Rows(), 3u);
    EXPECT_EQ(table.Row(1, 1), std::optional<size_t>(0));
    EXPECT_EQ(table.Row(1, 2), std::optional<size_t>(1));
    EXPECT_EQ(table.Row(2, 1), std::optional<size_t>(2));
    EXPECT_FALSE(table.Row(2, 2).has_value());
    EXPECT_FALSE(table.Row(0, 1).has_value());
    EXPECT_FALSE(table.Row(1, 0).has_value());
    EXPECT_FALSE(table.Row(7, 1).has_value());
    EXPECT_EQ(table.Capacity(0), 8u);
    EXPECT_EQ(table.Capacity(1), 0u);
    EXPECT_EQ(table.Capacity(2), 12u);
    table.Set(0, 4, uc::kUnique);
    table.Set(2, 0, uc::kUnique | uc::kUniqueTwo);
    table.Set(2, 8, uc::kUniqueTwo);
    table.Set(2, 8, uc::kUnique);  // set again: replaced
    EXPECT_EQ(table.Flags(0, 4), uc::kUnique);
    EXPECT_EQ(table.Flags(0, 3), 0);
    EXPECT_EQ(table.Flags(2, 0), uc::kUnique | uc::kUniqueTwo);
    EXPECT_EQ(table.Flags(2, 8), uc::kUnique);
    EXPECT_EQ(table.Totals(), std::make_pair(uint64_t{ 3 }, uint64_t{ 1 }));
    std::ostringstream os;
    table.Write(os);
    uc::Table back;
    ASSERT_EQ(back.Read(os.str()), "");
    EXPECT_EQ(back.Rows(), 3u);
    EXPECT_EQ(back.Row(2, 1), std::optional<size_t>(2));
    EXPECT_EQ(back.Capacity(2), 12u);
    for (size_t row = 0; row < 3; row++) {
        for (size_t core = 0; core < back.Capacity(row); core++) EXPECT_EQ(back.Flags(row, core), table.Flags(row, core));
    }
    uc::Table bad;
    EXPECT_NE(bad.Read("not a table"), "");
    auto const whole = os.str();
    EXPECT_NE(bad.Read(whole.substr(0, whole.size() - 1)), "");
    EXPECT_TRUE(bad.Empty());
}

// A thread's counter finds a gene's cores with the k-mer handler (as the build does) and reads their bits from the
// table; an alignment of the whole gene covers every core. A table whose row does not fit the gene's cores (another
// database's) leaves the record its anchor's counts (nullopt) and is counted; a copy of the counter counts alike.
TEST(UniqueCores, TheCounterCountsAGenesCoresFromTheTable) {
    constexpr size_t kK = 31, kM = 15;
    ClosedSyncmer syncmer{ kM, 7, 2, true };
    SimpleKmerHandler<ClosedSyncmer> handler{ kK, kM, syncmer };
    std::mt19937 rng(5);
    FakeGene gene{ test::RandomSequence(900, rng) };
    KmerList kmers;
    handler(std::string_view(gene.bases), kmers);
    ASSERT_GT(kmers.size(), 20u);
    auto table = std::make_shared<uc::Table>(std::vector<uint64_t>{ 0, 0, 1 }, std::vector<uint32_t>{ static_cast<uint32_t>(kmers.size()) });
    for (size_t i = 0; i < kmers.size(); i += 3) table->Set(0, i, uc::kUnique);
    size_t const unique = (kmers.size() + 2) / 3;
    uc::Counter<SimpleKmerHandler<ClosedSyncmer>> counter(table, handler);
    ASSERT_TRUE(counter.Active());
    auto const whole = counter.Count(1, 1, gene, M(gene.bases.size()), 0);
    ASSERT_TRUE(whole.has_value());
    EXPECT_EQ(whole->cores, kmers.size());
    EXPECT_EQ(whole->unique, unique);
    EXPECT_EQ(whole->unique_two, 0u);
    // The first 150 bases: the cores whose k-mers end there.
    size_t first = 0;
    for (auto const& k : kmers) first += k.second + kK <= 150;
    EXPECT_EQ(counter.Count(1, 1, gene, M(150), 0)->cores, first);
    EXPECT_FALSE(counter.Count(1, 2, gene, M(150), 0).has_value());  // a gene beyond the table
    auto copy = counter;
    EXPECT_EQ(copy.Count(1, 1, gene, M(gene.bases.size()), 0)->unique, unique);
    // A row with room for fewer cores than the gene has.
    auto small = std::make_shared<uc::Table>(std::vector<uint64_t>{ 0, 0, 1 }, std::vector<uint32_t>{ static_cast<uint32_t>(kmers.size() - 5) });
    uc::Counter<SimpleKmerHandler<ClosedSyncmer>> wrong(small, handler);
    EXPECT_FALSE(wrong.Count(1, 1, gene, M(150), 0).has_value());
    EXPECT_FALSE(wrong.Count(1, 1, gene, M(150), 0).has_value());  // from the cache
    EXPECT_EQ(wrong.Mismatched(), 1u);
    uc::Counter<SimpleKmerHandler<ClosedSyncmer>> none;
    EXPECT_FALSE(none.Active());
    EXPECT_FALSE(none.Count(1, 1, gene, M(150), 0).has_value());
}
