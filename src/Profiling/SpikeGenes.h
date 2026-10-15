// SpikeGenes.h - genes on which a sample's reads pile far beyond what the rest of their taxon's genes explain: a PCR
// product or a cloned fragment of one gene carried over into the library, or a contaminating fragment. In a real
// D6331 mock sample (Dohlman et al., SRR38158090) every species of the standard had one to three marker genes with
// 48-97% of its records; Akkermansia muciniphila had 2,108 of its 2,177 records on one stretch of secA and 1-8 on each of
// 37 other genes. The coverage features then read every species as a relative whose reads gather on a few genes: the
// model's probabilities fell to 0.38-0.80 and two species of the standard were missed
// (docs/claude/2026-10-15-mock-near-misses). A run counts each gene's records before it profiles a sample
// (Profiler::CountGeneRecords), finds every taxon's spike genes (SpikesOf) and judges whether the sample shows the
// pattern of carried-over products (Judge: spikes on several of its most abundant taxa, a share of its records). Only
// then is the sample in spike mode, and its spike genes' records are left out of the evidence, as a suspect gene copy's
// are (MicrobialProfile::PrepareMAPQ; --spike_genes auto, the default; on: every sample; off: none). Without the check,
// a relative's reads that pile on the few genes of an absent congener nearest them were taken as spikes too, and in two
// long-read mock samples three such congeners then looked present. Only records the profile counts (on genes with
// unique k-mers, MicrobialProfile::CountedRecords) enter the rule and the check.
//
// The rule (SpikesOf): the taxon's background is its records on its genes other than the kMaxGenes with the most records
// per base, over those genes' summed length (every hittable gene of the taxon counts, hit or not), so that up to
// kMaxGenes spikes do not hide each other. Those top genes are taken in order while each holds at least kMinRecords
// records and kMinFold times what the background gives a gene of its length, and while the taxon keeps more than
// kMinRestGenes genes with records. The genes taken are spikes only if together they hold kMinShare of the taxon's
// records or more. Outside that D6331 sample, on nine other real mock samples (five paired-end, four long-read) and 62
// simulated test samples of a reduced r226 database, it took genes of eight taxa, all absent, and none of a present
// species; no such sample had more than one of its ten taxa with the most records spiked, or 2% of its records on spikes.
#pragma once

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace protal::spike_genes {
    inline constexpr uint64_t kMinRecords = 50;   // a spike gene's records at least
    inline constexpr double kMinFold = 50;        // and this many times what the taxon's other genes give it
    inline constexpr double kMinShare = 0.5;      // the spikes hold this share of the taxon's records or more
    inline constexpr size_t kMinRestGenes = 5;    // the taxon keeps more genes with records than this
    inline constexpr size_t kMaxGenes = 3;        // spikes per taxon at most

    inline uint64_t Key(uint32_t taxid, uint32_t geneid) {
        return uint64_t{taxid} << 32 | geneid;
    }

    // Records per gene of a sample: Key(taxid, geneid) -> records.
    using Counts = std::unordered_map<uint64_t, uint64_t>;

    struct Spike {
        uint32_t taxid = 0;
        uint32_t geneid = 0;
        uint64_t records = 0;
        double expected = 0;  // the records the taxon's other genes give it

        double Fold() const {
            return expected > 0 ? static_cast<double>(records) / expected : 0;
        }
    };

    // The spike genes of taxon `taxid`: `records` are its hittable genes with records (gene id, records), `lengths` all
    // its hittable genes (gene id, length); both in any order. None if it has none.
    inline std::vector<Spike> SpikesOf(uint32_t taxid, std::vector<std::pair<uint32_t, uint64_t>> records,
                                       std::vector<std::pair<uint32_t, uint64_t>> const& lengths) {
        std::vector<Spike> spikes;
        if (records.size() <= kMinRestGenes + 1) return spikes;
        std::unordered_map<uint32_t, uint64_t> length_of;
        uint64_t total_length = 0;
        for (auto const& [gene, length] : lengths) {
            length_of[gene] = length;
            total_length += length;
        }
        // Only genes of known length count; in gene order, so that ties break the same way every time.
        std::erase_if(records, [&](auto const& r) { return length_of.find(r.first) == length_of.end() || length_of[r.first] == 0; });
        std::sort(records.begin(), records.end());
        uint64_t total_records = 0;
        for (auto const& [_, n] : records) total_records += n;
        if (total_length == 0 || total_records == 0) return spikes;

        if (records.size() <= kMinRestGenes + 1) return spikes;

        // The genes by records per base, the most first (ties: the lower gene id, as `records` is sorted).
        std::vector<size_t> order(records.size());
        for (size_t i = 0; i < order.size(); i++) order[i] = i;
        auto const rate = [&](size_t i) {
            return static_cast<double>(records[i].second) / static_cast<double>(length_of[records[i].first]);
        };
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return rate(a) > rate(b); });
        size_t const top = std::min(kMaxGenes, records.size() - kMinRestGenes - 1);
        uint64_t background_records = total_records, background_length = total_length;
        for (size_t k = 0; k < top; k++) {
            background_records -= records[order[k]].second;
            background_length -= length_of[records[order[k]].first];
        }
        double const background = background_length == 0 ? 0
            : static_cast<double>(background_records) / static_cast<double>(background_length);

        uint64_t taken_records = 0;
        for (size_t k = 0; k < top; k++) {
            auto const [gene, n] = records[order[k]];
            double const expected = background * static_cast<double>(length_of[gene]);
            if (n < kMinRecords || static_cast<double>(n) < kMinFold * expected) break;
            spikes.push_back({ taxid, gene, n, expected });
            taken_records += n;
        }
        if (spikes.empty() || static_cast<double>(taken_records) < kMinShare * static_cast<double>(total_records)) return {};
        return spikes;
    }

    // Adds the records of SAM text (whole lines; header and blank lines are skipped) to `counts`: per RNAME
    // (taxid_geneid), the records that are neither unmapped nor secondary (FLAG 0x4, 0x100). A long read's supplementary
    // records, its other genes, count. Lines it cannot read are skipped (the profile reports them).
    inline void CountRecords(std::string_view text, Counts& counts) {
        size_t pos = 0;
        while (pos < text.size()) {
            size_t end = text.find('\n', pos);
            if (end == std::string_view::npos) end = text.size();
            std::string_view line = text.substr(pos, end - pos);
            pos = end + 1;
            if (line.empty() || line.front() == '@') continue;
            size_t const tab1 = line.find('\t');
            if (tab1 == std::string_view::npos) continue;
            size_t const tab2 = line.find('\t', tab1 + 1);
            if (tab2 == std::string_view::npos) continue;
            size_t const tab3 = line.find('\t', tab2 + 1);
            if (tab3 == std::string_view::npos) continue;
            unsigned flag = 0;
            auto const* flag_begin = line.data() + tab1 + 1;
            if (std::from_chars(flag_begin, line.data() + tab2, flag).ec != std::errc{}) continue;
            if (flag & 0x104u) continue;
            std::string_view const rname = line.substr(tab2 + 1, tab3 - tab2 - 1);
            size_t const underscore = rname.find('_');
            if (underscore == std::string_view::npos) continue;
            uint32_t taxid = 0, geneid = 0;
            if (std::from_chars(rname.data(), rname.data() + underscore, taxid).ec != std::errc{}) continue;
            if (std::from_chars(rname.data() + underscore + 1, rname.data() + rname.size(), geneid).ec != std::errc{}) continue;
            counts[Key(taxid, geneid)]++;
        }
    }

    // Adds `from`'s counts to `into`.
    inline void Add(Counts& into, Counts const& from) {
        for (auto const& [key, n] : from) into[key] += n;
    }

    // Whether a sample is profiled in spike mode (--spike_genes): its spike genes left out of the evidence. Auto: if the
    // sample's own coverage shows the pattern (Judge); on: every sample; off: none.
    enum class Mode { kAuto, kOn, kOff };

    inline std::optional<Mode> ParseMode(std::string_view text) {
        if (text == "auto") return Mode::kAuto;
        if (text == "on") return Mode::kOn;
        if (text == "off") return Mode::kOff;
        return std::nullopt;
    }

    inline std::string_view ModeName(Mode mode) {
        return mode == Mode::kAuto ? "auto" : mode == Mode::kOn ? "on" : "off";
    }

    // The sample check (Judge): a library with PCR products or fragments carried over puts spikes on many species, the
    // abundant ones among them, while a relative's reads pile on a few genes of an absent congener only, a thin taxon.
    // A sample is in spike mode if at least kSampleMinTopSpiked of its kSampleTopTaxa taxa with the most records have
    // spike genes and the spikes hold kSampleMinShare of its records or more. On the real D6331 sample of Dohlman et al.:
    // 8 of the top 10, 81% of the records; on five other real paired-end and four long-read mock samples and 30
    // simulated test samples: none of the top 10, at most 1% (docs/claude/2026-10-15-mock-near-misses).
    inline constexpr size_t kSampleTopTaxa = 10;
    inline constexpr size_t kSampleMinTopSpiked = 2;
    inline constexpr double kSampleMinShare = 0.05;

    struct SampleCheck {
        uint64_t records = 0;        // the sample's records with a reference (CountRecords)
        uint64_t spike_records = 0;  // on every taxon's spike genes (SpikesOf)
        size_t spiked_taxa = 0;      // taxa with spike genes
        size_t top_taxa = 0;         // its taxa with the most records looked at (kSampleTopTaxa, or fewer)
        size_t top_spiked = 0;       // of those, the ones with spike genes
        bool spiky = false;          // the pattern of carried-over products: spike mode under Mode::kAuto

        double Share() const {
            return records == 0 ? 0 : static_cast<double>(spike_records) / static_cast<double>(records);
        }
    };

    // What a sample's records per gene (`counts`) and its taxa's spike genes (`spikes`, SpikesOf of every taxon) say
    // about the sample. Ties among the taxa with the most records break by the lower taxid.
    inline SampleCheck Judge(Counts const& counts, std::vector<Spike> const& spikes) {
        SampleCheck check;
        std::unordered_map<uint32_t, uint64_t> per_taxon;
        for (auto const& [key, n] : counts) {
            per_taxon[static_cast<uint32_t>(key >> 32)] += n;
            check.records += n;
        }
        std::unordered_map<uint32_t, bool> spiked;
        for (auto const& spike : spikes) {
            check.spike_records += spike.records;
            spiked[spike.taxid] = true;
        }
        check.spiked_taxa = spiked.size();
        std::vector<std::pair<uint32_t, uint64_t>> taxa(per_taxon.begin(), per_taxon.end());
        std::sort(taxa.begin(), taxa.end(), [](auto const& a, auto const& b) {
            return a.second != b.second ? a.second > b.second : a.first < b.first;
        });
        check.top_taxa = std::min(kSampleTopTaxa, taxa.size());
        for (size_t i = 0; i < check.top_taxa; i++) check.top_spiked += spiked.count(taxa[i].first);
        check.spiky = check.top_spiked >= kSampleMinTopSpiked && check.Share() >= kSampleMinShare;
        return check;
    }
}
