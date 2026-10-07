// GeneIncongruence.h - gene copies that do not belong to their species. --build compares every species' copy of each
// marker gene with the other species' copies (bottom sketches of 12-mers, as a run's congener distances) and flags a
// copy that is near-identical to a copy of another genus, family, order, class, phylum or domain while its own
// congeners' copies are farther, or it has none: a contaminating contig in a MAG, or a transferred gene (gene by gene,
// each gene's copies on all threads: ScanGene). Reads on such
// a copy are no evidence of its species, and every present organism with that gene puts a perfect read on it: at GTDB
// r226 a fifth of the false positives were reads of a present species of another genus at identity 0.99, recurring
// per reference (docs/claude/2026-10-03-false-positive-anatomy). A query leaves records on suspect copies out
// (suspect_copies.tsv in the database; --keep_suspect_copies). gene_incongruence.tsv beside the database lists every
// near pair across genera, for the record and for reporting the references.
#pragma once

#include "GeneConservation.h"

#include <omp.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <cstdint>
#include <fstream>
#include <istream>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace protal::gene_incongruence {
    inline const std::string kFileName = "suspect_copies.tsv";          // in the database: the copies a query leaves out
    inline const std::string kReportFileName = "gene_incongruence.tsv";  // beside it: every near pair across genera

    inline constexpr size_t kSketchSize = 128;              // bottom hashes per copy (gene_conservation::BottomSketch)
    inline constexpr double kDefaultSuspectDistance = 0.02;  // a copy within this of another genus's (--suspect_copy_distance)
    inline constexpr double kCongenerMargin = 0.02;          // and whose nearest congener's copy is this much farther; a copy
                                                             // without congeners: nearer to the partner than the partner's
                                                             // congeners by this much
    inline constexpr double kReportDistance = 0.05;          // pairs across genera reported in gene_incongruence.tsv
    inline constexpr size_t kProbeHashes = 8;                // a copy's smallest hashes looked up for candidates
    inline constexpr size_t kMaxBucket = 2048;               // copies sharing a hash beyond which the hash tells nothing
    inline constexpr uint32_t kGeneBits = 20;                // gene ids below 2^20 (as gene_conservation)

    // The deepest rank two species share: congeners (Genus) down to nothing but the domain, or None when their
    // domains differ or are unknown.
    enum class Rank : uint8_t { Genus = 0, Family, Order, Class, Phylum, Domain, None };

    inline constexpr std::array<std::string_view, 7> kRankNames = { "genus", "family", "order", "class", "phylum", "domain", "none" };

    inline std::string_view RankName(Rank rank) {
        return kRankNames[static_cast<size_t>(rank)];
    }

    inline std::optional<Rank> RankFromName(std::string_view name) {
        for (size_t i = 0; i < kRankNames.size(); i++) {
            if (kRankNames[i] == name) return static_cast<Rank>(i);
        }
        return std::nullopt;
    }

    // A species' ancestors by rank (taxids of internal_taxonomy.dmp; 0: unknown), genus first.
    struct Lineage {
        std::array<uint32_t, 6> ancestor{};

        uint32_t At(Rank rank) const {
            return ancestor[static_cast<size_t>(rank)];
        }
    };

    inline Rank SharedRank(Lineage const& a, Lineage const& b) {
        for (size_t r = 0; r < a.ancestor.size(); r++) {
            if (a.ancestor[r] != 0 && a.ancestor[r] == b.ancestor[r]) return static_cast<Rank>(r);
        }
        return Rank::None;
    }

    // Every node's lineage from internal_taxonomy.dmp (columns id, parent_id, external_id, name, rank, ...): the
    // nearest ancestor of each rank, the node itself included. Returns an error text, or "" when the file was read.
    inline std::string LineagesFromTaxonomy(std::istream& is, std::unordered_map<uint32_t, Lineage>& lineages) {
        std::unordered_map<uint32_t, uint32_t> parent;
        std::unordered_map<uint32_t, Rank> rank_of;
        std::string line;
        size_t nodes = 0;
        while (std::getline(is, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            std::vector<std::string_view> fields;
            std::string_view view(line);
            size_t start = 0;
            for (size_t tab; (tab = view.find('\t', start)) != std::string_view::npos; start = tab + 1) fields.push_back(view.substr(start, tab - start));
            fields.push_back(view.substr(start));
            if (fields.size() < 5) continue;
            uint32_t id = 0, up = 0;
            if (std::from_chars(fields[0].data(), fields[0].data() + fields[0].size(), id).ec != std::errc() ||
                std::from_chars(fields[1].data(), fields[1].data() + fields[1].size(), up).ec != std::errc()) continue;
            parent[id] = up;
            if (auto const rank = RankFromName(fields[4]); rank && *rank != Rank::None) rank_of[id] = *rank;
            nodes++;
        }
        if (nodes == 0) return "no taxonomy nodes";
        for (auto const& [id, _] : parent) {
            Lineage lineage;
            uint32_t node = id;
            for (int depth = 0; depth < 64; depth++) {
                if (auto const r = rank_of.find(node); r != rank_of.end() && lineage.ancestor[static_cast<size_t>(r->second)] == 0) {
                    lineage.ancestor[static_cast<size_t>(r->second)] = node;
                }
                auto const it = parent.find(node);
                if (it == parent.end() || it->second == node || it->second == 0) break;
                node = it->second;
            }
            lineages[id] = lineage;
        }
        return "";
    }

    inline std::string LineagesFromTaxonomy(std::string const& path, std::unordered_map<uint32_t, Lineage>& lineages) {
        std::ifstream is(path);
        if (!is) return "cannot read " + path;
        return LineagesFromTaxonomy(is, lineages);
    }

    // A species' copy of a gene, as its bottom sketch.
    struct Copy {
        uint32_t taxid = 0;
        std::vector<uint32_t> sketch;
    };

    // A suspect copy: near-identical to `partner`'s copy, of the shared rank `rank` (never Genus), at `distance`;
    // its nearest congener's copy at `congener_distance` (2 for none within kReportDistance).
    struct Suspect {
        uint32_t taxid = 0;
        uint32_t geneid = 0;
        uint32_t partner = 0;
        Rank rank = Rank::None;
        float distance = 2;
        float congener_distance = 2;
    };

    // A near pair of copies across genera.
    struct Pair {
        uint32_t geneid = 0;
        uint32_t taxid_a = 0;
        uint32_t taxid_b = 0;
        Rank rank = Rank::None;
        float distance = 2;
    };

    struct Result {
        std::vector<Suspect> suspects;  // by gene id, then taxid
        std::vector<Pair> pairs;        // by gene id, then taxids
        size_t copies = 0;              // copies compared
        size_t genes = 0;               // genes with two copies or more
        size_t candidates = 0;          // pairs of copies (each way) that share one of a copy's kProbeHashes smallest hashes
        size_t compared = 0;            // of them, those whose sketches were merged: the others cannot be within
                                        // kReportDistance (detail::SharedBound)
    };

    // One gene's copies, for ScanGene: copy i is taxid taxids[i]'s, its sketch (sorted, as BottomSketch gives it) the
    // length[i] hashes from Sketch(i), all in one array (`stride` hashes a copy) rather than a vector per copy.
    struct GeneCopies {
        std::vector<uint32_t> taxids;
        std::vector<uint32_t> hashes;
        std::vector<uint32_t> length;
        size_t stride = kSketchSize;

        size_t Size() const { return taxids.size(); }
        uint32_t const* Sketch(size_t i) const { return hashes.data() + i * stride; }

        // n copies of up to `stride_` hashes each, to be Set.
        void Resize(size_t n, size_t stride_) {
            stride = stride_;
            taxids.assign(n, 0);
            length.assign(n, 0);
            hashes.assign(n * stride, 0);
        }

        void Set(size_t i, uint32_t taxid, std::vector<uint32_t> const& sketch) {
            taxids[i] = taxid;
            length[i] = static_cast<uint32_t>(std::min(sketch.size(), stride));
            std::copy_n(sketch.begin(), length[i], hashes.begin() + static_cast<std::ptrdiff_t>(i * stride));
        }
    };

    namespace detail {
        // need[s]: the fewest hashes two sketches whose smaller one has s hashes must share among the s smallest of their
        // union for their distance, as Scan keeps it (SketchDistanceOf as a float), to be within kReportDistance; s + 1
        // if no count is. Taken from the distance itself, so a pair skipped for sharing fewer is one Scan never kept.
        inline std::vector<uint32_t> SharedNeeded(size_t max_s) {
            std::vector<uint32_t> need(max_s + 1);
            for (size_t s = 0; s <= max_s; s++) {
                need[s] = static_cast<uint32_t>(s + 1);
                for (size_t c = 1; c <= s; c++) {
                    if (!(static_cast<float>(gene_conservation::SketchDistanceOf(c, s)) > kReportDistance)) {
                        need[s] = static_cast<uint32_t>(c);
                        break;
                    }
                }
            }
            return need;
        }

        // A copy's signature: its hashes as bits of a 1,024-bit set. Two sketches share at most the bits both have
        // plus the hashes of either that fall on a bit another of its own hashes took (`extra`), which bounds the
        // hashes their merge can find shared (SharedBound).
        inline constexpr size_t kSignatureWords = 16;

        inline uint32_t SharedBound(uint64_t const* a, uint64_t const* b, uint32_t extra_a, uint32_t extra_b) {
            uint32_t bits = 0;
            for (size_t w = 0; w < kSignatureWords; w++) bits += static_cast<uint32_t>(__builtin_popcountll(a[w] & b[w]));
            return bits + std::min(extra_a, extra_b);
        }

        // The hashes in both of the s smallest of the union of sorted a and b (SketchDistance's merge, without branches),
        // or nullopt once `need` of them can no longer be reached.
        inline std::optional<size_t> MergeShared(uint32_t const* a, size_t la, uint32_t const* b, size_t lb, size_t s, size_t need) {
            size_t i = 0, j = 0, seen = 0, shared = 0;
            while (seen < s && i < la && j < lb) {
                uint32_t const x = a[i], y = b[j];
                i += x <= y;
                j += y <= x;
                shared += x == y;
                seen++;
                if (shared + (s - seen) < need) return std::nullopt;
            }
            return shared;
        }

        // Every (hash, copy) of a gene's sketches, sorted, in buckets by a mix of the hash: the copies holding a hash are
        // a run in its bucket. Counted, placed and sorted bucket by bucket on all threads.
        class HashIndex {
        public:
            static constexpr int kBucketBits = 12;

            static size_t BucketOf(uint32_t h) { return (h * 0x9E3779B1u) >> (32 - kBucketBits); }

            HashIndex(GeneCopies const& gene, int threads) {
                size_t constexpr buckets = size_t{1} << kBucketBits;
                size_t const n = gene.Size();
                threads = std::max(threads, 1);
                std::vector<uint64_t> counts(static_cast<size_t>(threads) * buckets, 0);  // per thread and bucket
                m_start.assign(buckets + 1, 0);
                uint64_t total = 0;
                for (size_t i = 0; i < n; i++) total += gene.length[i];
                m_entries.resize(total);
                #pragma omp parallel num_threads(threads)
                {
                    size_t const t = static_cast<size_t>(omp_get_thread_num());
                    uint64_t* mine = counts.data() + t * buckets;
                    // The same static schedule twice: each thread places the copies it counted.
                    #pragma omp for schedule(static)
                    for (int64_t i = 0; i < static_cast<int64_t>(n); i++) {
                        uint32_t const* sketch = gene.Sketch(static_cast<size_t>(i));
                        for (uint32_t k = 0; k < gene.length[static_cast<size_t>(i)]; k++) mine[BucketOf(sketch[k])]++;
                    }
                    #pragma omp single
                    {
                        size_t const team = static_cast<size_t>(omp_get_num_threads());
                        uint64_t position = 0;
                        for (size_t b = 0; b < buckets; b++) {
                            m_start[b] = position;
                            for (size_t u = 0; u < team; u++) {
                                uint64_t const c = counts[u * buckets + b];
                                counts[u * buckets + b] = position;
                                position += c;
                            }
                        }
                        m_start[buckets] = position;
                    }
                    #pragma omp for schedule(static)
                    for (int64_t i = 0; i < static_cast<int64_t>(n); i++) {
                        uint32_t const* sketch = gene.Sketch(static_cast<size_t>(i));
                        for (uint32_t k = 0; k < gene.length[static_cast<size_t>(i)]; k++) {
                            m_entries[mine[BucketOf(sketch[k])]++] = { sketch[k], static_cast<uint32_t>(i) };
                        }
                    }
                    #pragma omp for schedule(dynamic, 16)
                    for (int64_t b = 0; b < static_cast<int64_t>(buckets); b++) {
                        std::sort(m_entries.begin() + static_cast<std::ptrdiff_t>(m_start[static_cast<size_t>(b)]),
                                  m_entries.begin() + static_cast<std::ptrdiff_t>(m_start[static_cast<size_t>(b) + 1]));
                    }
                }
            }

            // The entries of hash h: [first, second).
            std::pair<std::pair<uint32_t, uint32_t> const*, std::pair<uint32_t, uint32_t> const*> Find(uint32_t h) const {
                size_t const b = BucketOf(h);
                auto const begin = m_entries.data() + m_start[b], end = m_entries.data() + m_start[b + 1];
                auto const lo = std::lower_bound(begin, end, std::make_pair(h, uint32_t{0}));
                auto const hi = std::upper_bound(lo, end, std::make_pair(h, UINT32_MAX));
                return { lo, hi };
            }

        private:
            std::vector<std::pair<uint32_t, uint32_t>> m_entries;
            std::vector<uint64_t> m_start;
        };
    }

    // Compares the copies of one gene (gene id `geneid`, in `gene`) and adds its suspects, its near pairs across genera
    // and its counts to result (unsorted: Finish sorts them). For each copy its nearest other copy among its congeners
    // and among the species of other genera: the candidates are the copies sharing one of its kProbeHashes smallest
    // hashes (a hash more than kMaxBucket copies hold tells nothing), their distance from the sketches. A candidate whose
    // sketch cannot share the hashes a distance within kReportDistance takes (detail::SharedBound, and the merge stopped
    // once they can no longer be reached) is left out at once: what the comparison keeps is exactly what it kept when
    // every candidate was merged in full, 95% or more of them to no end at GTDB r226 (docs/claude/2026-10-07-database-
    // build-audit.md, S1). The copies are compared on `threads` threads, each copy's nearest ones set by its thread only.
    //
    // A copy is suspect when a copy of another genus (the partner) is within suspect_distance and: the copy's genus has
    // other copies of the gene but its nearest one is kCongenerMargin farther than the partner's (its own genus does not
    // have this copy, the other does); or the copy's genus has no other copy (a singleton genus) and the partner's genus
    // does, with the partner's nearest congener's copy kCongenerMargin farther than this copy (the copy lies inside the
    // other genus's cluster, tighter than the cluster itself: a slow gene shared by a young family does not qualify).
    inline void ScanGene(uint32_t geneid, GeneCopies const& gene, std::unordered_map<uint32_t, Lineage> const& lineages,
                         double suspect_distance, int threads, Result& result) {
        size_t const n = gene.Size();
        result.copies += n;
        if (n < 2) return;
        result.genes++;
        threads = std::max(threads, 1);
        Lineage const no_lineage;
        std::vector<Lineage const*> lineage(n);
        // Copies of the gene per genus: a copy whose genus has others is judged against them.
        std::unordered_map<uint32_t, size_t> per_genus;
        std::vector<uint32_t> genus_of(n, 0);
        for (size_t i = 0; i < n; i++) {
            auto const it = lineages.find(gene.taxids[i]);
            lineage[i] = it == lineages.end() ? &no_lineage : &it->second;
            genus_of[i] = lineage[i]->At(Rank::Genus);
            if (genus_of[i] != 0) per_genus[genus_of[i]]++;
        }
        std::vector<uint64_t> signature(n * detail::kSignatureWords, 0);
        std::vector<uint32_t> extra(n, 0);
        #pragma omp parallel for schedule(static) num_threads(threads)
        for (int64_t i = 0; i < static_cast<int64_t>(n); i++) {
            uint64_t* sig = signature.data() + static_cast<size_t>(i) * detail::kSignatureWords;
            uint32_t const* sketch = gene.Sketch(static_cast<size_t>(i));
            uint32_t const length = gene.length[static_cast<size_t>(i)];
            for (uint32_t k = 0; k < length; k++) sig[(sketch[k] & 1023u) >> 6] |= uint64_t{1} << (sketch[k] & 63u);
            uint32_t bits = 0;
            for (size_t w = 0; w < detail::kSignatureWords; w++) bits += static_cast<uint32_t>(__builtin_popcountll(sig[w]));
            extra[static_cast<size_t>(i)] = length - bits;
        }
        detail::HashIndex const index(gene, threads);
        std::vector<uint32_t> const need = detail::SharedNeeded(gene.stride);

        // Pass 1: each copy's nearest congener's copy and nearest other genus's copy (the partner).
        std::vector<float> congener(n, 2), foreign(n, 2);
        std::vector<size_t> partner(n, n);
        std::vector<Rank> foreign_rank(n, Rank::None);
        std::vector<Pair> pairs;
        size_t candidates = 0, compared = 0;
        #pragma omp parallel num_threads(threads) reduction(+:candidates, compared)
        {
            std::vector<uint32_t> stamp(n, UINT32_MAX);  // stamp[j] == i: j is a candidate of i already
            std::vector<uint32_t> found;
            std::vector<Pair> mine;
            #pragma omp for schedule(dynamic, 64)
            for (int64_t ii = 0; ii < static_cast<int64_t>(n); ii++) {
                auto const i = static_cast<size_t>(ii);
                uint32_t const* a = gene.Sketch(i);
                uint32_t const la = gene.length[i];
                found.clear();
                size_t const probes = std::min<size_t>(kProbeHashes, la);
                for (size_t p = 0; p < probes; p++) {
                    auto const [lo, hi] = index.Find(a[p]);
                    if (static_cast<size_t>(hi - lo) > kMaxBucket) continue;
                    for (auto it = lo; it != hi; ++it) {
                        uint32_t const j = it->second;
                        if (j != i && stamp[j] != i) {
                            stamp[j] = static_cast<uint32_t>(i);
                            found.push_back(j);
                        }
                    }
                }
                candidates += found.size();
                uint32_t const taxid = gene.taxids[i];
                uint64_t const* sig_a = signature.data() + i * detail::kSignatureWords;
                for (uint32_t const j : found) {
                    if (gene.taxids[j] == taxid) continue;
                    uint32_t const lb = gene.length[j];
                    size_t const s = std::min(la, lb);
                    if (detail::SharedBound(sig_a, signature.data() + j * detail::kSignatureWords, extra[i], extra[j]) < need[s]) continue;
                    compared++;
                    auto const shared = detail::MergeShared(a, la, gene.Sketch(j), lb, s, need[s]);
                    if (!shared) continue;
                    auto const d = static_cast<float>(gene_conservation::SketchDistanceOf(*shared, s));
                    if (d > kReportDistance) continue;
                    Rank const rank = SharedRank(*lineage[i], *lineage[j]);
                    if (rank == Rank::Genus) {
                        congener[i] = std::min(congener[i], d);
                        continue;
                    }
                    if (d < foreign[i] || (d == foreign[i] && gene.taxids[j] < gene.taxids[partner[i]])) {
                        foreign[i] = d;
                        partner[i] = j;
                        foreign_rank[i] = rank;
                    }
                    if (taxid < gene.taxids[j]) mine.push_back({ geneid, taxid, gene.taxids[j], rank, d });
                    else mine.push_back({ geneid, gene.taxids[j], taxid, rank, d });
                }
            }
            #pragma omp critical(gene_incongruence_pairs)
            pairs.insert(pairs.end(), mine.begin(), mine.end());
        }
        result.candidates += candidates;
        result.compared += compared;
        // Pass 2: the verdicts.
        for (size_t i = 0; i < n; i++) {
            if (foreign[i] > suspect_distance || partner[i] >= n) continue;
            bool const has_congeners = genus_of[i] != 0 && per_genus[genus_of[i]] > 1;
            bool suspect;
            if (has_congeners) {
                suspect = congener[i] > foreign[i] + kCongenerMargin;  // 2 when none within kReportDistance
            } else {
                size_t const p = partner[i];
                bool const partner_has_congeners = genus_of[p] != 0 && per_genus[genus_of[p]] > 1;
                suspect = partner_has_congeners && congener[p] < 2 && foreign[i] + kCongenerMargin <= congener[p];
            }
            if (suspect) {
                result.suspects.push_back({ gene.taxids[i], geneid, gene.taxids[partner[i]], foreign_rank[i], foreign[i], congener[i] });
            }
        }
        result.pairs.insert(result.pairs.end(), pairs.begin(), pairs.end());
    }

    // Sorts what ScanGene added: the suspects by gene id and taxid, the pairs by gene id and taxids, each once.
    inline void Finish(Result& result) {
        std::sort(result.suspects.begin(), result.suspects.end(), [](Suspect const& a, Suspect const& b) {
            return a.geneid != b.geneid ? a.geneid < b.geneid : a.taxid < b.taxid;
        });
        auto pair_key = [](Pair const& p) { return std::tie(p.geneid, p.taxid_a, p.taxid_b); };
        std::sort(result.pairs.begin(), result.pairs.end(), [&](Pair const& a, Pair const& b) { return pair_key(a) < pair_key(b); });
        result.pairs.erase(std::unique(result.pairs.begin(), result.pairs.end(), [&](Pair const& a, Pair const& b) { return pair_key(a) == pair_key(b); }),
                           result.pairs.end());
    }

    // ScanGene on every gene of copies_by_gene (copies_by_gene[geneid], in any order), then Finish: the same for any
    // thread count.
    inline Result Scan(std::vector<std::vector<Copy>> const& copies_by_gene, std::unordered_map<uint32_t, Lineage> const& lineages,
                       double suspect_distance, int threads) {
        Result result;
        for (size_t geneid = 0; geneid < copies_by_gene.size(); geneid++) {
            auto const& list = copies_by_gene[geneid];
            size_t stride = 0;
            for (auto const& copy : list) stride = std::max(stride, copy.sketch.size());
            GeneCopies gene;
            gene.Resize(list.size(), stride);
            for (size_t i = 0; i < list.size(); i++) gene.Set(i, list[i].taxid, list[i].sketch);
            ScanGene(static_cast<uint32_t>(geneid), gene, lineages, suspect_distance, threads, result);
        }
        Finish(result);
        return result;
    }

    // The suspect copies of a database (suspect_copies.tsv): a run leaves records on them out.
    class Table {
    public:
        static uint64_t Key(uint32_t taxid, uint32_t geneid) {
            return (static_cast<uint64_t>(taxid) << kGeneBits) | geneid;
        }

        void Add(Suspect const& suspect) {
            if (m_keys.insert(Key(suspect.taxid, suspect.geneid)).second) m_rows.push_back(suspect);
        }

        bool Contains(uint32_t taxid, uint32_t geneid) const {
            return !m_keys.empty() && m_keys.count(Key(taxid, geneid)) > 0;
        }

        bool Empty() const { return m_rows.empty(); }
        size_t Size() const { return m_rows.size(); }
        std::vector<Suspect> const& Rows() const { return m_rows; }

        // The species with a suspect copy.
        size_t Species() const {
            std::unordered_set<uint32_t> species;
            for (auto const& row : m_rows) species.insert(row.taxid);
            return species.size();
        }

        // suspect_copies.tsv: a header line starting with "taxid", then per copy taxid, geneid, partner_taxid,
        // shared_rank, distance, congener_distance (tab-separated; lines starting with '#' are skipped). Returns an
        // error text, or "" when every line was read.
        std::string Read(std::istream& is) {
            std::string line;
            size_t number = 0;
            while (std::getline(is, line)) {
                number++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty() || line[0] == '#' || line.rfind("taxid", 0) == 0) continue;
                std::vector<std::string_view> fields;
                std::string_view view(line);
                size_t start = 0;
                for (size_t tab; (tab = view.find('\t', start)) != std::string_view::npos; start = tab + 1) fields.push_back(view.substr(start, tab - start));
                fields.push_back(view.substr(start));
                if (fields.size() < 6) return "line " + std::to_string(number) + " has " + std::to_string(fields.size()) + " fields, not 6";
                Suspect row;
                auto number_of = [&](std::string_view f, auto& value) {
                    auto const [ptr, ec] = std::from_chars(f.data(), f.data() + f.size(), value);
                    return ec == std::errc() && ptr == f.data() + f.size();
                };
                auto const rank = RankFromName(fields[3]);
                bool const no_congener = fields[5] == "none";
                if (no_congener) row.congener_distance = 2;
                if (!number_of(fields[0], row.taxid) || !number_of(fields[1], row.geneid) || !number_of(fields[2], row.partner) ||
                    !rank || !number_of(fields[4], row.distance) || (!no_congener && !number_of(fields[5], row.congener_distance))) {
                    return "line " + std::to_string(number) + " is not taxid, geneid, partner_taxid, shared_rank, distance, congener_distance: " + line;
                }
                if (row.geneid >= (uint32_t{1} << kGeneBits)) return "line " + std::to_string(number) + ": gene id " + std::to_string(row.geneid) + " too large";
                row.rank = *rank;
                Add(row);
            }
            return "";
        }

        void Write(std::ostream& os) const {
            os << "taxid\tgeneid\tpartner_taxid\tshared_rank\tdistance\tcongener_distance\n";
            for (auto const& row : m_rows) {
                os << row.taxid << '\t' << row.geneid << '\t' << row.partner << '\t' << RankName(row.rank) << '\t' << row.distance << '\t';
                if (row.congener_distance >= 2) os << "none";
                else os << row.congener_distance;
                os << '\n';
            }
        }

    private:
        std::unordered_set<uint64_t> m_keys;
        std::vector<Suspect> m_rows;
    };

    // gene_incongruence.tsv: every near pair across genera (distance at most kReportDistance), with which of the two
    // copies are suspect.
    inline void WriteReport(std::ostream& os, Result const& result) {
        std::unordered_set<uint64_t> suspect;
        for (auto const& s : result.suspects) suspect.insert(Table::Key(s.taxid, s.geneid));
        os << "geneid\ttaxid_a\ttaxid_b\tshared_rank\tdistance\tsuspect_a\tsuspect_b\n";
        for (auto const& p : result.pairs) {
            os << p.geneid << '\t' << p.taxid_a << '\t' << p.taxid_b << '\t' << RankName(p.rank) << '\t' << p.distance << '\t'
               << (suspect.count(Table::Key(p.taxid_a, p.geneid)) ? 1 : 0) << '\t' << (suspect.count(Table::Key(p.taxid_b, p.geneid)) ? 1 : 0) << '\n';
        }
    }
}
