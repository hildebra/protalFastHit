// ColumnWeightsBuild.h - how --build makes column_weights.tsv (ColumnWeights.h): per family and marker gene, the
// conservation of every column from the family's own copies, by alignment (WFA2, protal's scores, both ends partly
// free: gene calls start and end apart).
//
// Per family and gene:
//  1. The genus references: one species per genus, chosen by a hash of the genus (the same species for every gene it
//     has; the next by hash where it lacks the gene). Up to Settings::genera genera, chosen by hash, vote; the first
//     of them gives the family reference, whose copy is the columns' coordinate system. Every genus reference is
//     aligned to the family reference (up to kMaxAmongDivergence), for the mapping of its genus's copies and, for the
//     voting genera, the `among` votes: per column the genus references' bases, their consensus and how many agree
//     (CodeOf: (k + 1)/(n + 2)); per codon column the amino acids alike (aa).
//  2. Within each genus: every species' copy aligned to the genus reference (up to kMaxWithinDivergence): per genus
//     column the species' bases (the reference's too), the consensus and the share agreeing, k/n, where two or more
//     were compared (one sequence says nothing of how a column varies); the genus's shares are carried onto the family
//     columns through the reference's alignment, and a family column's `within` is the mean of its genera's shares
//     (one vote per genus, however many species it has) with the pseudocounts over all the species compared
//     (CodeOfPooled: the evidence of the family's genera accumulates). The amino acids alike, into `aa` together with
//     the genus references' (the mean of the two shares, the pseudocounts over both).
//  3. Every copy's mapping: its runs on the family columns, its alignment to the genus reference composed with the
//     reference's to the family reference. A genus not among the voters still maps (its reference is aligned too).
// A family of one genus has that genus's reference as the family reference and no among votes beyond its own; a copy
// without a family in the taxonomy, or whose alignments fail, gets no row (its columns unknown: the features -1 for a
// taxon without any). One alignment per copy and one per genus reference and gene: at r226 (16.6 M copies) a few
// CPU-hours, the gene's copies in memory while it is compared, as the congener gaps do (Build.h WriteCongenerGaps).
#pragma once

#include "ColumnWeights.h"
#include "WFA2Wrapper2.h"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace protal::column_weights {
    inline constexpr double kFreeEnds = 0.15;              // of each copy's length, free at either end
    inline constexpr double kMaxWithinDivergence = 0.2;    // a species' copy farther than this from its genus reference is not used
    inline constexpr double kMaxAmongDivergence = 0.4;     // a genus reference farther than this from the family reference is not used
    inline constexpr double kMinOverlap = 0.5;             // aligned columns needed, of the shorter copy's length
    inline constexpr uint32_t kMinCompared = 2;            // sequences compared at a column for a genus's (or the genus
                                                           // references') share to vote in within (aa)
    inline constexpr size_t kAminoAcids = 22;              // 20 + stop + unknown
    inline constexpr uint32_t kMinConsensusGenera = 3;     // genus references compared at a column for a consensus base
    inline constexpr uint32_t kMinConsensusVotes = 2;      // of them agreeing, and more than half (a 2026-10-10 lesson: with one
                                                           // genus per family the "consensus" was a congener's or the species'
                                                           // own copy, and polarising by it broke the two-congener genera)

    struct Settings {
        size_t genera = 10;    // genera per family that vote in the among estimate (--column_weights; 0: no table)
        size_t species = 24;   // species per genus aligned for the within estimate (every species is aligned for its mapping)
    };

    // Why MapTo found no alignment: an empty sequence or one longer than kMaxLength, the aligner gave up (its score
    // beyond the divergence's bound: far, or gappy), too little overlap, or more divergent than allowed.
    enum class MapFailure : uint8_t { kNone, kSequence, kGaveUp, kOverlap, kDivergence };
    inline constexpr size_t kFailures = 5;
    inline constexpr size_t kAmong = 0, kWithin = 1;  // the alignments' kinds in Stats
    inline constexpr double kFarFactor = 1.5;         // a divergence failure beyond this times the limit counts as far

    struct Stats {
        size_t genes = 0;             // genes with a family of two copies or more
        size_t families = 0;          // (family, gene) rows
        size_t copies = 0;            // copies mapped
        size_t alignments = 0;        // alignments made
        size_t unaligned = 0;         // of them failed
        size_t single_genus = 0;      // families of one genus
        // By kind (kAmong: a genus reference against the family reference; kWithin: a species against its genus
        // reference): the alignments made, those failed by cause (MapFailure), and of the divergence failures those
        // beyond kFarFactor times the limit; and the species' copies never aligned because their genus reference failed.
        std::array<size_t, 2> made{};
        std::array<std::array<size_t, kFailures>, 2> failed{};
        std::array<size_t, 2> far{};
        size_t skipped = 0;

        void Note(size_t kind, MapFailure why, double divergence, double limit) {
            alignments++;
            made[kind]++;
            if (why == MapFailure::kNone) return;
            unaligned++;
            failed[kind][static_cast<size_t>(why)]++;
            far[kind] += why == MapFailure::kDivergence && divergence > kFarFactor * limit;
        }

        void Add(Stats const& o) {
            genes += o.genes;
            families += o.families;
            copies += o.copies;
            alignments += o.alignments;
            unaligned += o.unaligned;
            single_genus += o.single_genus;
            for (size_t k = 0; k < 2; k++) {
                made[k] += o.made[k];
                far[k] += o.far[k];
                for (size_t c = 0; c < kFailures; c++) failed[k][c] += o.failed[k][c];
            }
            skipped += o.skipped;
        }
    };

    namespace detail {
        inline uint64_t Mix(uint64_t x) {
            x += 0x9E3779B97F4A7C15ull;
            x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
            x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
            return x ^ (x >> 31);
        }

        inline uint8_t AminoIndex(char aa) {
            static constexpr char kLetters[] = "ACDEFGHIKLMNPQRSTVWY*";
            for (uint8_t i = 0; kLetters[i]; i++) {
                if (kLetters[i] == aa) return i;
            }
            return static_cast<uint8_t>(kAminoAcids - 1);
        }

        // Of the votes per class, the best class's count (the consensus) and the class.
        template<size_t N>
        inline std::pair<uint32_t, uint8_t> Best(std::array<uint32_t, N> const& votes) {
            uint8_t best = 0;
            for (uint8_t c = 1; c < N; c++) {
                if (votes[c] > votes[best]) best = c;
            }
            return { votes[best], best };
        }
    }

    // The alignment of `copy` to `reference` as a position map: per reference position the copy's position, or -1
    // where the copy has no base there (a gap, or outside the aligned stretch). Empty when the alignment fails
    // (farther than `max_divergence`, too little overlap, or the aligner gives up). `cigar` is scratch. `why` gets the
    // cause of a failure (kNone on success), `divergence` the differences per aligned column where they were counted.
    inline std::vector<int32_t> MapTo(std::string const& reference, std::string const& copy, double max_divergence,
                                      WFA2Wrapper2& aligner, std::string& cigar, MapFailure* why = nullptr,
                                      double* divergence = nullptr) {
        std::vector<int32_t> map;
        if (why) *why = MapFailure::kNone;
        if (divergence) *divergence = 0;
        if (reference.empty() || copy.empty() || reference.size() > kMaxLength || copy.size() > kMaxLength) {
            if (why) *why = MapFailure::kSequence;
            return map;
        }
        if (reference == copy) {
            map.resize(reference.size());
            for (size_t i = 0; i < map.size(); i++) map[i] = static_cast<int32_t>(i);
            return map;
        }
        int const free_ref = static_cast<int>(kFreeEnds * static_cast<double>(reference.size()));
        int const free_copy = static_cast<int>(kFreeEnds * static_cast<double>(copy.size()));
        size_t const longer = std::max(reference.size(), copy.size());
        int const max_score = static_cast<int>(4.5 * max_divergence * static_cast<double>(longer)) + 64;
        aligner.Alignment(copy, reference, free_copy, free_copy, free_ref, free_ref, max_score);
        if (!aligner.Success()) {
            if (why) *why = MapFailure::kGaveUp;
            return map;
        }
        aligner.CigarInto(cigar);
        map.assign(reference.size(), -1);
        size_t ref = 0, query = 0, columns = 0, differences = 0;
        size_t begin = 0, end = cigar.size();
        while (begin < end && (cigar[begin] == 'I' || cigar[begin] == 'D')) begin++;
        while (end > begin && (cigar[end - 1] == 'I' || cigar[end - 1] == 'D')) end--;
        for (size_t i = 0; i < cigar.size(); i++) {
            char const op = cigar[i];
            bool const inside = i >= begin && i < end;
            if (op == 'M' || op == '=' || op == 'X') {
                if (inside && ref < map.size() && query < copy.size()) {
                    map[ref] = static_cast<int32_t>(query);
                    columns++;
                    differences += op == 'X';
                }
                ref++;
                query++;
            } else if (op == 'I') {   // the copy has a base the reference lacks
                query++;
                if (inside) {
                    columns++;
                    differences++;
                }
            } else if (op == 'D') {   // the reference has a base the copy lacks
                ref++;
                if (inside) {
                    columns++;
                    differences++;
                }
            }
        }
        size_t const shorter = std::min(reference.size(), copy.size());
        if (columns == 0 || static_cast<double>(columns) < kMinOverlap * static_cast<double>(shorter)) {
            if (why) *why = MapFailure::kOverlap;
            map.clear();
            return map;
        }
        double const diverged = static_cast<double>(differences) / static_cast<double>(columns);
        if (divergence) *divergence = diverged;
        if (diverged > max_divergence) {
            if (why) *why = MapFailure::kDivergence;
            map.clear();
        }
        return map;
    }

    // The runs of a copy on the family columns from its map (per copy position the column, -1 none).
    inline std::vector<Run> RunsOf(std::vector<int32_t> const& copy_to_column) {
        std::vector<Run> runs;
        size_t i = 0;
        while (i < copy_to_column.size()) {
            if (copy_to_column[i] < 0) {
                i++;
                continue;
            }
            size_t j = i + 1;
            while (j < copy_to_column.size() && copy_to_column[j] == copy_to_column[i] + static_cast<int32_t>(j - i) &&
                   j - i < UINT16_MAX) {
                j++;
            }
            runs.push_back(Run{ static_cast<uint16_t>(i), static_cast<uint16_t>(copy_to_column[i]), static_cast<uint16_t>(j - i) });
            i = j;
        }
        return runs;
    }

    // The votes of one genus's copies on its reference's columns: per column the bases and per codon the amino acids.
    struct GenusVotes {
        std::vector<std::array<uint32_t, 5>> bases;    // per reference position: A, C, G, T, other
        std::vector<std::array<uint32_t, kAminoAcids>> amino;  // per reference codon
        explicit GenusVotes(size_t length) : bases(length, std::array<uint32_t, 5>{}), amino((length + 2) / 3, std::array<uint32_t, kAminoAcids>{}) {}

        // A copy's bases at the reference positions its map gives (copy_at[ref] = copy position or -1).
        void Add(std::string const& copy, std::vector<int32_t> const& copy_at) {
            size_t const n = std::min(copy_at.size(), bases.size());
            for (size_t i = 0; i < n; i++) {
                if (copy_at[i] < 0) continue;
                bases[i][BaseCode(copy[static_cast<size_t>(copy_at[i])])]++;
            }
            for (size_t c = 0; 3 * c + 2 < n; c++) {
                int32_t const a = copy_at[3 * c], b = copy_at[3 * c + 1], d = copy_at[3 * c + 2];
                if (a < 0 || b != a + 1 || d != a + 2) continue;   // a codon split by a gap, or not aligned
                char const aa = AminoAcid(BaseCode(copy[static_cast<size_t>(a)]), BaseCode(copy[static_cast<size_t>(b)]),
                                          BaseCode(copy[static_cast<size_t>(d)]));
                amino[c][detail::AminoIndex(aa)]++;
            }
        }

        // Per position the bases A-T compared and how many carry the consensus; per codon the amino acids alike (the
        // unknown class does not vote).
        void Raw(std::vector<uint32_t>& agree, std::vector<uint32_t>& compared, std::vector<uint32_t>& aa_agree,
                 std::vector<uint32_t>& aa_compared) const {
            agree.assign(bases.size(), 0);
            compared.assign(bases.size(), 0);
            for (size_t i = 0; i < bases.size(); i++) {
                std::array<uint32_t, 4> four{ bases[i][0], bases[i][1], bases[i][2], bases[i][3] };
                compared[i] = four[0] + four[1] + four[2] + four[3];
                agree[i] = detail::Best(four).first;
            }
            aa_agree.assign(amino.size(), 0);
            aa_compared.assign(amino.size(), 0);
            for (size_t c = 0; c < amino.size(); c++) {
                std::array<uint32_t, kAminoAcids - 1> known{};
                for (size_t k = 0; k + 1 < kAminoAcids; k++) {
                    known[k] = amino[c][k];
                    aa_compared[c] += amino[c][k];
                }
                aa_agree[c] = detail::Best(known).first;
            }
        }

        // Per position the share agreeing with the consensus, (k + 1)/(n + 2) over the bases A-T compared (NaN without
        // one), and the consensus base; per codon the same for the amino acids.
        void Shares(std::vector<double>& base_share, std::vector<uint8_t>& consensus, std::vector<double>& aa_share) const {
            base_share.assign(bases.size(), std::nan(""));
            consensus.assign(bases.size(), kNoBase);
            for (size_t i = 0; i < bases.size(); i++) {
                std::array<uint32_t, 4> four{ bases[i][0], bases[i][1], bases[i][2], bases[i][3] };
                uint32_t const compared = four[0] + four[1] + four[2] + four[3];
                if (compared == 0) continue;
                auto const [agree, best] = detail::Best(four);
                base_share[i] = (static_cast<double>(agree) + 1.0) / (static_cast<double>(compared) + 2.0);
                consensus[i] = best;
            }
            aa_share.assign(amino.size(), std::nan(""));
            for (size_t c = 0; c < amino.size(); c++) {
                uint32_t compared = 0;
                for (size_t k = 0; k + 1 < kAminoAcids; k++) compared += amino[c][k];  // the unknown class does not vote
                if (compared == 0) continue;
                std::array<uint32_t, kAminoAcids - 1> known{};
                for (size_t k = 0; k + 1 < kAminoAcids; k++) known[k] = amino[c][k];
                auto const [agree, best] = detail::Best(known);
                aa_share[c] = (static_cast<double>(agree) + 1.0) / (static_cast<double>(compared) + 2.0);
            }
        }
    };

    // One family's copies of a gene: `taxids[i]` holds `seqs[i]`, `genus_of[i]` its genus (0: none; such copies are
    // left out). Adds the family's row and its copies' mappings to `families` and `copies` (the copy rows refer to the
    // family row by its index in `families`), on the calling thread. Returns whether a row was made.
    inline bool BuildFamily(uint32_t family, uint32_t gene, std::vector<uint32_t> const& taxids, std::vector<std::string> const& seqs,
                            std::vector<uint32_t> const& genus_of, Settings const& settings, WFA2Wrapper2& aligner, Stats& stats,
                            std::vector<Family>& families,
                            std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<Run>>>& copies) {
        // The genera and their members, each member by its hash (the genus reference first).
        std::map<uint32_t, std::vector<std::pair<uint64_t, size_t>>> genera;
        for (size_t i = 0; i < taxids.size(); i++) {
            if (genus_of[i] == 0 || seqs[i].empty()) continue;
            genera[genus_of[i]].emplace_back(detail::Mix(taxids[i] ^ detail::Mix(genus_of[i])), i);
        }
        if (genera.empty()) return false;
        std::vector<std::pair<uint64_t, uint32_t>> by_hash;  // the genera by hash
        for (auto& [genus, members] : genera) {
            std::sort(members.begin(), members.end());
            by_hash.emplace_back(detail::Mix(genus ^ detail::Mix(family)), genus);
        }
        std::sort(by_hash.begin(), by_hash.end());
        size_t const voters = std::min(settings.genera, by_hash.size());
        stats.single_genus += genera.size() == 1;
        uint32_t const family_genus = by_hash[0].second;
        size_t const family_ref = genera[family_genus][0].second;
        std::string const& reference = seqs[family_ref];
        if (reference.size() < 3 || reference.size() > kMaxLength) return false;
        size_t const length = reference.size();

        std::string cigar;
        GenusVotes among(length);
        // Per family column (codon) the genera's raw shares summed, the genera voting, and the sequences they compared.
        std::vector<double> within_sum(length, 0.0);
        std::vector<uint32_t> within_n(length, 0), within_compared(length, 0);
        std::vector<double> aa_within_sum((length + 2) / 3, 0.0);
        std::vector<uint32_t> aa_within_n((length + 2) / 3, 0), aa_within_compared((length + 2) / 3, 0);
        std::vector<std::pair<uint32_t, std::vector<Run>>> mapped;  // (copy index, runs)
        uint16_t voted = 0;
        MapFailure why;
        double divergence;
        for (size_t g = 0; g < by_hash.size(); g++) {
            uint32_t const genus = by_hash[g].second;
            auto const& members = genera[genus];
            size_t const ref = members[0].second;
            std::string const& genus_ref = seqs[ref];
            // The genus reference on the family columns: per family column the genus reference's position.
            std::vector<int32_t> ref_at;
            if (ref == family_ref) {
                ref_at.resize(length);
                for (size_t i = 0; i < length; i++) ref_at[i] = static_cast<int32_t>(i);
            } else {
                ref_at = MapTo(reference, genus_ref, kMaxAmongDivergence, aligner, cigar, &why, &divergence);
                stats.Note(kAmong, why, divergence, kMaxAmongDivergence);
                if (ref_at.empty()) {
                    stats.skipped += members.size() - 1;
                    continue;  // a genus too far from the family reference: its copies get no row
                }
            }
            bool const votes = g < voters;
            if (votes) {
                among.Add(genus_ref, ref_at);
                voted++;
            }
            // The genus's species on the genus reference's columns.
            GenusVotes within(genus_ref.size());
            std::vector<int32_t> identity(genus_ref.size());
            for (size_t i = 0; i < identity.size(); i++) identity[i] = static_cast<int32_t>(i);
            within.Add(genus_ref, identity);
            // The genus reference's columns on the family's: per genus reference position the family column.
            std::vector<int32_t> column_of(genus_ref.size(), -1);
            for (size_t c = 0; c < length; c++) {
                if (ref_at[c] >= 0 && static_cast<size_t>(ref_at[c]) < column_of.size()) column_of[static_cast<size_t>(ref_at[c])] = static_cast<int32_t>(c);
            }
            mapped.emplace_back(static_cast<uint32_t>(ref), RunsOf(column_of));
            size_t aligned_species = 1;
            for (size_t m = 1; m < members.size(); m++) {
                size_t const i = members[m].second;
                auto const copy_at = MapTo(genus_ref, seqs[i], kMaxWithinDivergence, aligner, cigar, &why, &divergence);
                stats.Note(kWithin, why, divergence, kMaxWithinDivergence);
                if (copy_at.empty()) continue;
                if (aligned_species < settings.species) {
                    within.Add(seqs[i], copy_at);
                    aligned_species++;
                }
                // The copy's columns: its position at each genus reference position, through the reference's columns.
                std::vector<int32_t> copy_to_column(seqs[i].size(), -1);
                for (size_t r = 0; r < copy_at.size(); r++) {
                    if (copy_at[r] >= 0 && column_of[r] >= 0) copy_to_column[static_cast<size_t>(copy_at[r])] = column_of[r];
                }
                mapped.emplace_back(static_cast<uint32_t>(i), RunsOf(copy_to_column));
            }
            if (!votes) continue;
            std::vector<uint32_t> agree, compared, aa_agree, aa_compared;
            within.Raw(agree, compared, aa_agree, aa_compared);
            for (size_t r = 0; r < agree.size(); r++) {
                if (compared[r] < kMinCompared || column_of[r] < 0) continue;
                size_t const c = static_cast<size_t>(column_of[r]);
                within_sum[c] += static_cast<double>(agree[r]) / static_cast<double>(compared[r]);
                within_n[c]++;
                within_compared[c] += compared[r];
            }
            for (size_t c = 0; c < aa_agree.size(); c++) {
                if (aa_compared[c] < kMinCompared || 3 * c + 2 >= column_of.size()) continue;
                int32_t const a = column_of[3 * c];
                if (a < 0 || a % 3 != 0 || column_of[3 * c + 1] != a + 1 || column_of[3 * c + 2] != a + 2) continue;
                size_t const fc = static_cast<size_t>(a) / 3;
                aa_within_sum[fc] += static_cast<double>(aa_agree[c]) / static_cast<double>(aa_compared[c]);
                aa_within_n[fc]++;
                aa_within_compared[fc] += aa_compared[c];
            }
        }
        if (mapped.empty()) return false;
        Family row;
        row.family = family;
        row.gene = gene;
        row.reference = taxids[family_ref];
        row.genera = voted;
        row.within.assign(length, kNoCode);
        row.among.assign(length, kNoCode);
        row.aa.assign((length + 2) / 3, kNoCode);
        std::vector<double> among_share, among_aa;
        among.Shares(among_share, row.consensus, among_aa);
        std::vector<uint32_t> among_agree, among_compared, among_aa_agree, among_aa_compared;
        among.Raw(among_agree, among_compared, among_aa_agree, among_aa_compared);
        for (size_t c = 0; c < length; c++) {
            if (within_n[c] > 0) row.within[c] = CodeOfPooled(within_sum[c] / within_n[c], within_compared[c]);
            if (!std::isnan(among_share[c])) row.among[c] = CodeOfShare(among_share[c], 1);
            // The consensus base polarises a small genus's ancestry sites (AncestrySites.h), so it must come from
            // outside that genus: a majority of at least kMinConsensusVotes genus references among kMinConsensusGenera
            // or more compared (one of which may be the species' own genus). A family of one genus, or two that
            // differ, has no consensus there.
            auto const& votes = among.bases[c];
            uint32_t const compared = votes[0] + votes[1] + votes[2] + votes[3];
            auto const [agree, best] = detail::Best(std::array<uint32_t, 4>{ votes[0], votes[1], votes[2], votes[3] });
            if (compared < kMinConsensusGenera || agree < kMinConsensusVotes || 2 * agree <= compared) row.consensus[c] = kNoBase;
        }
        // The amino acid: the mean of the within share (over the genera) and the genus references' share, the
        // pseudocounts over the sequences both compared.
        for (size_t c = 0; c < row.aa.size(); c++) {
            double sum = 0;
            size_t n = 0, compared = 0;
            if (aa_within_n[c] > 0) {
                sum += aa_within_sum[c] / aa_within_n[c];
                n++;
                compared += aa_within_compared[c];
            }
            if (c < among_aa_compared.size() && among_aa_compared[c] >= kMinCompared) {
                sum += static_cast<double>(among_aa_agree[c]) / static_cast<double>(among_aa_compared[c]);
                n++;
                compared += among_aa_compared[c];
            }
            if (n > 0) row.aa[c] = CodeOfPooled(sum / static_cast<double>(n), compared);
        }
        uint32_t const index = static_cast<uint32_t>(families.size());
        families.push_back(std::move(row));
        for (auto& [i, runs] : mapped) {
            if (runs.empty()) continue;
            copies.emplace_back(taxids[i], gene, index, std::move(runs));
            stats.copies++;
        }
        stats.families++;
        return true;
    }

    // How a table's codes fall, for the build log: per kind (within, among, aa) a histogram of the columns' (codons')
    // codes, kNoCode included. Line() gives each kind's columns with an estimate, the median and largest code, and the
    // share at kConservedCode or more (the conserved columns the rate features count; among: at `among_full` or more,
    // the full differences of the weighted site shift). Until 2026-10-10 no column reached kConservedCode (r226 v22).
    struct Codes {
        size_t columns = 0;
        std::array<std::array<size_t, kMaxCode + 1>, 3> histogram{};

        std::string Line(uint8_t among_full) const {
            static constexpr char const* kNames[] = { "within", "among", "aa" };
            std::ostringstream s;
            s.setf(std::ios::fixed);
            for (size_t k = 0; k < 3; k++) {
                auto const& h = histogram[k];
                size_t const total = std::accumulate(h.begin(), h.end(), size_t{ 0 });
                size_t const known = total - h[kNoCode];
                uint8_t const threshold = k == 1 ? among_full : kConservedCode;
                size_t high = 0, seen = 0, largest = 0, median = 0;
                for (size_t c = 1; c <= kMaxCode; c++) {
                    if (h[c] == 0) continue;
                    largest = c;
                    if (c >= threshold) high += h[c];
                    if (median == 0 && 2 * (seen + h[c]) >= known) median = c;
                    seen += h[c];
                }
                double const share = known ? 100.0 * static_cast<double>(high) / static_cast<double>(known) : 0.0;
                s << (k ? "; " : "") << kNames[k] << " " << known << " of " << total << " with an estimate, median "
                  << median << ", largest " << largest << ", " << std::setprecision(2) << share << "% at "
                  << static_cast<int>(threshold) << " or more";
            }
            return s.str();
        }
    };

    inline Codes CodeSummary(std::vector<Family> const& families) {
        Codes codes;
        for (auto const& f : families) {
            codes.columns += f.within.size();
            for (uint8_t const c : f.within) codes.histogram[0][std::min<size_t>(c, kMaxCode)]++;
            for (uint8_t const c : f.among) codes.histogram[1][std::min<size_t>(c, kMaxCode)]++;
            for (uint8_t const c : f.aa) codes.histogram[2][std::min<size_t>(c, kMaxCode)]++;
        }
        return codes;
    }

    // One gene's copies (taxids[i] holds seqs[i]; one copy per taxid), each with its genus and family (0: none), on
    // `threads` threads, one family per task. The rows go to `families` and `copies` in family order (deterministic).
    inline void ScanGene(uint32_t gene, std::vector<uint32_t> const& taxids, std::vector<std::string> const& seqs,
                         std::vector<uint32_t> const& genus_of, std::vector<uint32_t> const& family_of, Settings const& settings,
                         int threads, Stats& stats, std::vector<Family>& families,
                         std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<Run>>>& copies) {
        threads = std::max(threads, 1);
        std::map<uint32_t, std::vector<size_t>> by_family;
        for (size_t i = 0; i < taxids.size(); i++) {
            if (family_of[i] != 0 && genus_of[i] != 0) by_family[family_of[i]].push_back(i);
        }
        std::vector<std::pair<uint32_t, std::vector<size_t>>> tasks(by_family.begin(), by_family.end());
        if (tasks.empty()) return;
        stats.genes++;
        std::vector<std::vector<Family>> fam_out(tasks.size());
        std::vector<std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<Run>>>> copy_out(tasks.size());
        std::vector<Stats> stat_out(tasks.size());
        #pragma omp parallel num_threads(threads)
        {
            WFA2Wrapper2 aligner(4, 6, 2, 0);
            #pragma omp for schedule(dynamic, 1)
            for (int64_t t = 0; t < static_cast<int64_t>(tasks.size()); t++) {
                auto const& [family, members] = tasks[static_cast<size_t>(t)];
                std::vector<uint32_t> ids, genus;
                std::vector<std::string> sequences;
                for (size_t const i : members) {
                    ids.push_back(taxids[i]);
                    genus.push_back(genus_of[i]);
                    sequences.push_back(seqs[i]);
                }
                BuildFamily(family, gene, ids, sequences, genus, settings, aligner, stat_out[static_cast<size_t>(t)],
                            fam_out[static_cast<size_t>(t)], copy_out[static_cast<size_t>(t)]);
            }
        }
        for (size_t t = 0; t < tasks.size(); t++) {
            uint32_t const base = static_cast<uint32_t>(families.size());
            for (auto& f : fam_out[t]) families.push_back(std::move(f));
            for (auto& [taxid, g, row, runs] : copy_out[t]) copies.emplace_back(taxid, g, base + row, std::move(runs));
            stats.Add(stat_out[t]);
        }
    }
}
