// LongReads.h - aligning long reads (PacBio). A long read spans several marker genes, so each gene
// it hits is aligned on its own, over the read's window on that gene: the gene's hits and those of
// its homologs in other taxa at the same place are a segment, with the MAPQ of its best hit against
// the next. A read's segments are written as one primary and supplementary alignments, hard-clipped
// to their aligned bases (ProtalLongReadOutputHandler); the profiler takes each segment as a read.
//
// Seeding and chaining keep read positions in 16 bits (LookupResult, ChainLink), so a read longer
// than kMaxLongReadChunk is seeded in overlapping chunks (ChunkRead). A hit belongs to the chunk
// whose core holds the centre of its gene: with chunks overlapping by at least a gene's length (and
// margins), that chunk holds the whole gene, and no gene is counted twice.
#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>
#include "AlignmentStrategy.h"
#include "AlignmentOutputHandler.h"
#include "GenomeLoader.h"
#include "KmerUtils.h"

namespace protal {
    inline constexpr size_t kMaxLongReadChunk = 65000;

    // An interval [start, end) of a read, in forward read coordinates.
    struct ReadInterval {
        int64_t start = 0;
        int64_t end = 0;

        int64_t Length() const { return end - start; }
        int64_t TwiceCentre() const { return start + end; }
        int64_t Overlap(ReadInterval const& other) const {
            return std::max<int64_t>(0, std::min(end, other.end) - std::max(start, other.start));
        }
        ReadInterval Within(int64_t read_length) const {
            return { std::clamp<int64_t>(start, 0, read_length), std::clamp<int64_t>(end, 0, read_length) };
        }
    };

    struct ReadChunk {
        size_t offset = 0;
        size_t length = 0;
        ReadInterval core;  // this chunk's hits: those whose gene's centre lies here
    };

    // Where the cores of a read's first and last chunk begin and end: a gene the read starts or ends
    // in has its centre beyond the read.
    inline constexpr int64_t kBeyondRead = int64_t{1} << 40;

    // The chunks a read of read_length bases is seeded in: the whole read up to max_chunk bases,
    // else chunks of max_chunk that overlap by at least `overlap` (at most max_chunk / 2). Their
    // cores, split at the middle of each overlap and reaching beyond the read's ends, tile the
    // read, and an interval of up to `overlap` bases lies in the chunk whose core holds its centre.
    inline std::vector<ReadChunk> ChunkRead(size_t read_length, size_t max_chunk, size_t overlap) {
        std::vector<ReadChunk> chunks;
        if (read_length <= max_chunk) {
            chunks.push_back({ 0, read_length, { -kBeyondRead, kBeyondRead } });
            return chunks;
        }
        overlap = std::min(overlap, max_chunk / 2);
        size_t const step = max_chunk - overlap;
        size_t const count = (read_length - max_chunk + step - 1) / step + 1;
        for (size_t i = 0; i < count; i++) {
            chunks.push_back({ std::min(i * step, read_length - max_chunk), max_chunk, {} });
        }
        for (size_t i = 0; i < count; i++) {
            auto& core = chunks[i].core;
            core.start = i == 0 ? -kBeyondRead : static_cast<int64_t>(chunks[i].offset + chunks[i - 1].offset + chunks[i - 1].length) / 2;
            if (i > 0) chunks[i - 1].core.end = core.start;
        }
        chunks.back().core.end = kBeyondRead;
        return chunks;
    }

    // A place where a gene may lie on a read: an anchor of one chunk.
    struct LongReadCandidate {
        Anchor anchor;          // read positions in the anchor's orientation of its chunk
        ReadChunk chunk;
        ReadInterval gene;      // the gene on the read, as the anchor places it (may reach past the read)
        ReadInterval window;    // the bases aligned: the gene with margins, within the read
    };

    // A read's alignment to one gene. seq and qual are its aligned bases in reference orientation;
    // clip_left and clip_right count the read's other bases before and after them.
    struct LongReadHit {
        AlignmentResult alignment;
        std::string seq;
        std::string qual;
        size_t clip_left = 0;
        size_t clip_right = 0;
        ReadInterval aligned;   // the aligned bases, forward read coordinates
    };

    // The hits of one place of a read: a gene and its homologs in other taxa, best first, each
    // alignment once, and the MAPQ of the best against the next; or, if the gene alone could not
    // tell them apart but the read's other genes could (SettleByRead), the best is the hit of the
    // read's taxon and by_read is set.
    struct LongReadSegment {
        std::vector<LongReadHit> hits;
        int mapq = 0;
        bool by_read = false;
    };

    using LongReadSegments = std::vector<LongReadSegment>;

    // The score hits are ranked and MAPQ computed with, as for short reads (Bitscore).
    inline int LongReadBitscore(LongReadHit const& hit) {
        return hit.alignment.GetAlignmentInfo().Score(2, 3, 1, 2);
    }

    // Ranks a segment's hits, drops repeated alignments (the same gene, strand and position from
    // two anchors, which would give MAPQ 0) and sets its MAPQ.
    inline void RankLongReadSegment(LongReadSegment& segment) {
        auto& hits = segment.hits;
        std::stable_sort(hits.begin(), hits.end(), [](LongReadHit const& a, LongReadHit const& b) {
            return LongReadBitscore(a) > LongReadBitscore(b);
        });
        auto same = [](LongReadHit const& a, LongReadHit const& b) {
            auto const& x = a.alignment;
            auto const& y = b.alignment;
            return x.Taxid() == y.Taxid() && x.GeneId() == y.GeneId() && x.Forward() == y.Forward() &&
                   x.GetAlignmentInfo().gene_alignment_start == y.GetAlignmentInfo().gene_alignment_start;
        };
        std::vector<LongReadHit> distinct;
        for (auto& hit : hits) {
            if (std::none_of(distinct.begin(), distinct.end(), [&](LongReadHit const& d) { return same(d, hit); })) {
                distinct.emplace_back(std::move(hit));
            }
        }
        hits = std::move(distinct);
        int const best = hits.empty() ? 0 : LongReadBitscore(hits.front());
        int const second = hits.size() > 1 ? LongReadBitscore(hits[1]) : 0;
        segment.mapq = best > 0 ? MAPQv2(best, second) : 0;
    }

    // A segment is confident with a MAPQ of at least this (the profiler's minimum).
    inline constexpr int kConfidentMapq = 4;
    // The share of its confident segments that one taxon needs to be a read's taxon.
    inline constexpr double kReadTaxonShare = 2.0 / 3.0;

    // The taxon of a read, from its confident segments: each votes for the taxon of its best hit, and
    // a taxon with at least kReadTaxonShare of the votes is the read's. Votes count equally: a gene
    // that the read's species lacks in the database aligns to a relative with a high MAPQ (it has no
    // competitor), and weighting by MAPQ let one such gene outvote several of the species' own.
    // Returns the taxon with the lowest MAPQ of the segments that voted for it (the weakest evidence
    // it rests on); nullopt without one.
    struct ReadTaxon {
        uint32_t taxid = 0;
        int mapq = 0;
    };

    inline std::optional<ReadTaxon> TaxonOfRead(LongReadSegments const& segments) {
        std::vector<std::pair<uint32_t, int>> votes;  // taxon, confident segments
        int total = 0;
        for (auto const& segment : segments) {
            if (segment.mapq < kConfidentMapq || segment.by_read || segment.hits.empty()) continue;
            auto const taxid = segment.hits.front().alignment.Taxid();
            auto it = std::find_if(votes.begin(), votes.end(), [taxid](auto const& v) { return v.first == taxid; });
            if (it == votes.end()) votes.emplace_back(taxid, 1);
            else it->second++;
            total++;
        }
        if (votes.empty()) return std::nullopt;
        auto const best = std::max_element(votes.begin(), votes.end(), [](auto const& a, auto const& b) { return a.second < b.second; });
        if (best->second < kReadTaxonShare * total) return std::nullopt;
        ReadTaxon taxon{ best->first, std::numeric_limits<int>::max() };
        for (auto const& segment : segments) {
            if (segment.mapq >= kConfidentMapq && !segment.by_read && !segment.hits.empty() &&
                segment.hits.front().alignment.Taxid() == taxon.taxid) {
                taxon.mapq = std::min(taxon.mapq, segment.mapq);
            }
        }
        return taxon;
    }

    // Settles an ambiguous segment (MAPQ below kConfidentMapq) by its read's taxon: if the segment
    // has a hit of the taxon that its best cannot be told apart from (the best's MAPQ against it is
    // below kConfidentMapq), that hit becomes the best, with the taxon's MAPQ. Returns whether it did.
    inline bool SettleByRead(LongReadSegment& segment, ReadTaxon const& taxon) {
        if (segment.mapq >= kConfidentMapq || segment.hits.empty()) return false;
        auto& hits = segment.hits;
        auto it = std::find_if(hits.begin(), hits.end(), [&taxon](LongReadHit const& h) { return h.alignment.Taxid() == taxon.taxid; });
        if (it == hits.end()) return false;
        int const best = LongReadBitscore(hits.front()), own = LongReadBitscore(*it);
        if (it != hits.begin() && own < best && MAPQv2(best, own) >= kConfidentMapq) return false;  // clearly another taxon's gene
        std::rotate(hits.begin(), it, it + 1);
        segment.mapq = taxon.mapq;
        segment.by_read = true;
        return true;
    }

    // Groups candidates into segments: the genes of one segment cover at least half of the shorter
    // one's place on the read. Candidates come longest anchor first; segments list their indices.
    inline std::vector<std::vector<size_t>> LongReadSegmentsOf(std::vector<LongReadCandidate> const& candidates, int64_t read_length) {
        std::vector<std::vector<size_t>> members;
        std::vector<ReadInterval> places;
        for (size_t i = 0; i < candidates.size(); i++) {
            auto const place = candidates[i].gene.Within(read_length);
            size_t s = 0;
            while (s < places.size() && 2 * places[s].Overlap(place) < std::min(places[s].Length(), place.Length())) s++;
            if (s == places.size()) {
                places.push_back(place);
                members.emplace_back();
            }
            members[s].push_back(i);
        }
        return members;
    }

    // Seeds, chains and aligns long reads, one read at a time into its segments. Each thread uses
    // its own copy.
    template<typename KmerHandler, typename AnchorFinder>
    class LongReadAligner {
        KmerHandler m_kmer_handler;
        AnchorFinder m_anchor_finder;
        SimpleAlignmentHandler m_alignment_handler;
        GenomeLoader& m_genomes;
        size_t m_align_top;
        double m_min_ani;
        size_t m_max_chunk;
        size_t m_overlap;

        std::string m_chunk;
        std::string m_window;
        std::string m_window_rev;
        std::string m_id;
        KmerList m_kmers;
        SeedList m_seeds;
        AlignmentAnchorList m_anchors;
        std::vector<LongReadCandidate> m_candidates;

        size_t m_last_seeds = 0;
        size_t m_last_anchors = 0;
        size_t m_chunked_reads = 0;
        size_t m_ambiguous = 0;  // segments the gene alone could not assign (MAPQ below kConfidentMapq)
        size_t m_settled = 0;    // of these, those settled by their read's other genes

        std::vector<bool> m_aligned;                        // per candidate of the read
        std::vector<std::vector<size_t>> m_segment_members; // per segment of the read, its candidates

        // Bases aligned beyond a gene's place on the read at each end, for indels between the seeds.
        static int64_t Margin(int64_t gene_length) {
            return 100 + gene_length / 20;
        }

        void AddCandidate(Anchor const& anchor, ReadChunk const& chunk, int64_t read_length) {
            if (anchor.chain.empty()) return;
            int64_t const gene_length = m_genomes.GetGenome(anchor.taxid).GetGeneOMP(anchor.geneid).Sequence().length();
            int64_t const chunk_length = static_cast<int64_t>(chunk.length);
            auto to_read = [&](ReadInterval interval) {  // anchor orientation of the chunk to forward read coordinates
                if (!anchor.forward) interval = { chunk_length - interval.end, chunk_length - interval.start };
                return ReadInterval{ interval.start + static_cast<int64_t>(chunk.offset), interval.end + static_cast<int64_t>(chunk.offset) };
            };
            // The gene starts where the anchor's first seed puts it.
            int64_t const gene_start = static_cast<int64_t>(anchor.chain.front().readpos) - static_cast<int64_t>(anchor.chain.front().genepos);
            auto const gene = to_read({ gene_start, gene_start + gene_length });
            if (gene.TwiceCentre() < 2 * chunk.core.start || gene.TwiceCentre() >= 2 * chunk.core.end) return;
            auto const& last = anchor.chain.back();
            auto const chain = to_read({ anchor.chain.front().readpos, static_cast<int64_t>(last.readpos) + last.length });
            int64_t const margin = Margin(gene_length);
            ReadInterval const window = ReadInterval{ std::min(gene.start - margin, chain.start),
                                                      std::max(gene.end + margin, chain.end) }.Within(read_length);
            m_candidates.push_back({ anchor, chunk, gene, window });
        }

        std::optional<LongReadHit> Align(LongReadCandidate const& candidate, FastxRecord const& record) {
            auto const& read = record.sequence;
            int64_t const read_length = static_cast<int64_t>(read.size());
            auto const& window = candidate.window;
            size_t const window_length = static_cast<size_t>(window.Length());
            m_window.assign(read, window.start, window_length);
            m_window_rev = KmerUtils::ReverseComplement(m_window);

            // The anchor's positions from its chunk to the window, in the anchor's orientation.
            Anchor anchor = candidate.anchor;
            int64_t const shift = anchor.forward ?
                    static_cast<int64_t>(candidate.chunk.offset) - window.start :
                    window.end - static_cast<int64_t>(candidate.chunk.offset + candidate.chunk.length);
            for (auto& link : anchor.chain) {
                int64_t const position = static_cast<int64_t>(link.readpos) + shift;
                // Anchors hold read positions in 16 bits: a window over 64 kb (a gene of ~59 kb) is not aligned.
                if (position < 0 || position + link.length > static_cast<int64_t>(window_length) ||
                    position + link.length > std::numeric_limits<uint16_t>::max()) return std::nullopt;
                link.readpos = static_cast<uint16_t>(position);
            }

            AlignmentResult result;
            m_id = record.id;
            if (!m_alignment_handler.AlignAnchor(anchor, result, m_window, m_window_rev, false, m_id)) return std::nullopt;
            auto& info = result.GetAlignmentInfo();
            if (info.GetProxyANI() < m_min_ani) return std::nullopt;

            // The window's bases outside the alignment become hard clips, as the rest of the read.
            auto const& cigar = info.cigar;
            size_t const left = std::find_if(cigar.begin(), cigar.end(), [](char c) { return c != 'S'; }) - cigar.begin();
            size_t const right = left == cigar.size() ? 0 :
                                 std::find_if(cigar.rbegin(), cigar.rend(), [](char c) { return c != 'S'; }) - cigar.rbegin();
            if (left + right >= window_length) return std::nullopt;
            info.cigar = cigar.substr(left, cigar.size() - left - right);
            info.GetInstructionCountsAndCompress();

            LongReadHit hit;
            hit.alignment = std::move(result);
            auto const& oriented = anchor.forward ? m_window : m_window_rev;
            size_t const aligned_length = window_length - left - right;
            hit.seq = oriented.substr(left, aligned_length);
            if (record.quality.size() == read.size()) {
                std::string quality = record.quality.substr(window.start, window_length);
                if (!anchor.forward) std::reverse(quality.begin(), quality.end());
                hit.qual = quality.substr(left, aligned_length);
            }
            int64_t const oriented_start = anchor.forward ? window.start : read_length - window.end;
            hit.clip_left = static_cast<size_t>(oriented_start) + left;
            hit.clip_right = read.size() - hit.clip_left - aligned_length;
            hit.aligned = anchor.forward ?
                    ReadInterval{ window.start + static_cast<int64_t>(left), window.end - static_cast<int64_t>(right) } :
                    ReadInterval{ window.start + static_cast<int64_t>(right), window.end - static_cast<int64_t>(left) };
            return hit;
        }

    public:
        LongReadAligner(KmerHandler const& kmer_handler, AnchorFinder const& anchor_finder,
                        SimpleAlignmentHandler const& alignment_handler, GenomeLoader& genomes, size_t align_top,
                        double min_ani, size_t max_chunk = kMaxLongReadChunk) :
                m_kmer_handler(kmer_handler),
                m_anchor_finder(anchor_finder),
                m_alignment_handler(alignment_handler),
                m_genomes(genomes),
                m_align_top(align_top),
                m_min_ani(min_ani),
                m_max_chunk(max_chunk) {
            int64_t const longest = static_cast<int64_t>(genomes.MaxGeneLength());
            m_overlap = static_cast<size_t>(longest + 2 * Margin(longest));
        }

        LongReadAligner(LongReadAligner const& other) :
                m_kmer_handler(other.m_kmer_handler),
                m_anchor_finder(other.m_anchor_finder),
                m_alignment_handler(other.m_alignment_handler),
                m_genomes(other.m_genomes),
                m_align_top(other.m_align_top),
                m_min_ani(other.m_min_ani),
                m_max_chunk(other.m_max_chunk),
                m_overlap(other.m_overlap) {}

        // The overlap of the chunks of a read longer than max_chunk: the longest gene and its margins.
        size_t ChunkOverlap() const { return m_overlap; }
        // Whether every gene lies wholly in the chunk that owns it: ChunkRead caps the overlap at
        // half a chunk, which a gene of ~29 kb and its margins exceed (none of GTDB's markers).
        bool ChunksHoldEveryGene() const { return m_overlap <= m_max_chunk / 2; }

        AnchorFinder& GetAnchorFinder() { return m_anchor_finder; }
        KmerHandler& GetKmerHandler() { return m_kmer_handler; }
        SimpleAlignmentHandler& GetAlignmentHandler() { return m_alignment_handler; }

        // Seeds and anchors of the last read, and reads that were seeded in chunks.
        size_t LastSeeds() const { return m_last_seeds; }
        size_t LastAnchors() const { return m_last_anchors; }
        size_t ChunkedReads() const { return m_chunked_reads; }
        // Segments of all reads so far that their gene alone could not assign, and those of them the
        // read's other genes settled.
        size_t AmbiguousSegments() const { return m_ambiguous; }
        size_t SettledSegments() const { return m_settled; }

        // The segments of `record`: for each place of the read where genes were found, the up to
        // align_top longest anchors (and those as long as the last) are aligned.
        void operator()(FastxRecord const& record, LongReadSegments& segments) {
            segments.clear();
            m_candidates.clear();
            m_last_seeds = 0;
            m_last_anchors = 0;
            auto const& read = record.sequence;
            int64_t const read_length = static_cast<int64_t>(read.size());

            auto const chunks = ChunkRead(read.size(), m_max_chunk, m_overlap);
            m_chunked_reads += chunks.size() > 1;
            for (auto const& chunk : chunks) {
                m_chunk.assign(read, chunk.offset, chunk.length);
                m_kmers.clear();
                m_seeds.clear();
                m_anchors.clear();
                m_kmer_handler(std::string_view(m_chunk), m_kmers);
                m_anchor_finder(m_kmers, m_seeds, m_anchors, m_chunk);
                m_last_seeds += m_seeds.size();
                m_last_anchors += m_anchors.size();
                for (auto const& anchor : m_anchors) AddCandidate(anchor, chunk, read_length);
            }
            if (m_candidates.empty()) return;

            std::stable_sort(m_candidates.begin(), m_candidates.end(), [](LongReadCandidate const& a, LongReadCandidate const& b) {
                return a.anchor.total_length > b.anchor.total_length;
            });
            m_aligned.assign(m_candidates.size(), false);
            m_segment_members.clear();
            for (auto& members : LongReadSegmentsOf(m_candidates, read_length)) {
                LongReadSegment segment;
                int take_top = static_cast<int>(m_align_top);
                Anchor const* last = nullptr;
                for (auto i : members) {
                    auto const& candidate = m_candidates[i];
                    if (--take_top < 0 && last && last->total_length != candidate.anchor.total_length) break;
                    last = &candidate.anchor;
                    m_aligned[i] = true;
                    if (auto hit = Align(candidate, record)) segment.hits.push_back(std::move(*hit));
                }
                if (segment.hits.empty()) continue;
                RankLongReadSegment(segment);
                segments.push_back(std::move(segment));
                m_segment_members.push_back(std::move(members));
            }

            // Genes that alone fit several taxa are settled by the read's other genes: their hit of
            // the read's taxon (aligned now if its anchor was not among the longest) becomes the best.
            for (auto const& segment : segments) m_ambiguous += segment.mapq < kConfidentMapq;
            auto const taxon = TaxonOfRead(segments);
            if (!taxon) return;
            for (size_t s = 0; s < segments.size(); s++) {
                auto& segment = segments[s];
                if (segment.mapq >= kConfidentMapq) continue;
                auto const& hits = segment.hits;
                if (std::none_of(hits.begin(), hits.end(), [&taxon](LongReadHit const& h) { return h.alignment.Taxid() == taxon->taxid; })) {
                    for (auto i : m_segment_members[s]) {  // longest anchor first
                        if (m_candidates[i].anchor.taxid != taxon->taxid || m_aligned[i]) continue;
                        m_aligned[i] = true;
                        if (auto hit = Align(m_candidates[i], record)) {
                            segment.hits.push_back(std::move(*hit));
                            RankLongReadSegment(segment);
                        }
                        break;
                    }
                }
                m_settled += SettleByRead(segment, *taxon);
            }
        }
    };

    // The SAM record of a hit, with its CIGAR without the hard clips (see LongReadCigar).
    inline void LongReadHitToSam(SamEntry& sam, LongReadHit const& hit, std::string const& qname) {
        auto const& ar = hit.alignment;
        auto const& info = ar.GetAlignmentInfo();
        sam.m_qname = qname;
        sam.m_flag = 0;
        sam.m_rname = std::to_string(ar.Taxid()) + "_" + std::to_string(ar.GeneId());
        sam.m_pos = info.gene_alignment_start + 1;
        sam.m_mapq = 0;
        sam.m_cigar = info.compressed_cigar;
        sam.m_rnext = "*";
        sam.m_pnext = 0;
        sam.m_tlen = 0;
        sam.m_seq = hit.seq;
        sam.m_qual = hit.qual;
        sam.m_uniques = ar.Uniques();
        sam.m_uniques_two = ar.UniquesTwo();
        Flag::SetReadReverseComplement(sam.m_flag, !ar.Forward());
    }

    // A hit's CIGAR as written: its alignment between the hard clips of the rest of the read.
    inline std::string LongReadCigar(LongReadHit const& hit) {
        return (hit.clip_left ? std::to_string(hit.clip_left) + "H" : std::string()) + hit.alignment.GetAlignmentInfo().compressed_cigar +
               (hit.clip_right ? std::to_string(hit.clip_right) + "H" : std::string());
    }

    /*
     * Output of long reads: per segment its best hit, the others (up to -m) as secondary alignments
     * (0x100, MAPQ 0). The best segment's best hit is the read's primary alignment, those of the
     * other segments are supplementary (0x800), each with its segment's MAPQ; a segment settled by
     * its read's other genes (SettleByRead) is tagged ZR:i:1. Records hold only the aligned bases
     * (hard clips); readers take each primary or supplementary record and the secondary ones after
     * it as one read's candidates.
     */
    class ProtalLongReadOutputHandler {
    private:
        SamSink& m_sink;
        BufferedStringOutput m_sam_output;
        std::vector<uint64_t> m_genes;  // named by the buffered records (SamGeneKey)
        GenomeLoader& m_genomes;
        size_t m_max_out = 1;
        double m_min_cigar_ani = 0;

        // Hands the buffered records and the genes they name to the sink (which serialises writers).
        void Flush() {
            m_sam_output.Drain([this](char const* data, size_t size) { m_sink.Write(data, size, m_genes); });
            if (!m_genes.empty()) m_sink.Write(nullptr, 0, m_genes);
        }

    public:
        size_t alignments = 0;

        ProtalLongReadOutputHandler(SamSink& sink, size_t max_out, size_t sam_buffer_capacity, GenomeLoader& genomes, double min_cigar_ani=0.0) :
                m_sink(sink),
                m_sam_output(sam_buffer_capacity),
                m_genomes(genomes),
                m_max_out(max_out),
                m_min_cigar_ani(min_cigar_ani) {}

        ProtalLongReadOutputHandler(ProtalLongReadOutputHandler const& other) :
                m_sink(other.m_sink),
                m_sam_output(other.m_sam_output.Capacity()),
                m_genomes(other.m_genomes),
                m_max_out(other.m_max_out),
                m_min_cigar_ani(other.m_min_cigar_ani) {}

        ~ProtalLongReadOutputHandler() {
            Flush();
        }

        void operator () (LongReadSegments& segments, FastxRecord& record) {
            if (segments.empty()) return;
            // The best segment first, then the others in read order.
            std::stable_sort(segments.begin(), segments.end(), [](LongReadSegment const& a, LongReadSegment const& b) {
                return a.hits.front().aligned.start < b.hits.front().aligned.start;
            });
            auto best = std::max_element(segments.begin(), segments.end(), [](LongReadSegment const& a, LongReadSegment const& b) {
                return LongReadBitscore(a.hits.front()) < LongReadBitscore(b.hits.front());
            });
            std::rotate(segments.begin(), best, best + 1);

            auto const qname = ReadQName(record.id);
            SNPList snps;
            SamEntry sam;
            bool primary_written = false;
            std::string read_records;
            for (auto const& segment : segments) {
                if (CigarANI(segment.hits.front().alignment.Cigar()) < m_min_cigar_ani) continue;
                bool first = true;
                size_t written = 0;
                for (auto const& hit : segment.hits) {
                    LongReadHitToSam(sam, hit, qname);
                    auto const& ar = hit.alignment;
                    auto const reference = m_genomes.GetGenome(ar.Taxid()).GetGene(ar.GeneId()).Sequence();
                    if (!ExtractSNPs(sam, reference, snps, ar.Taxid(), ar.GeneId(), 0)) {
#pragma omp critical(err_out)
                        {
                            std::cerr << "Inconsistent alignment of " << record.id << ": " << sam.ToString() << std::endl;
                            PrintAlignment(sam, reference, std::cerr);
                        }
                        continue;  // only this candidate; the segment's others are still written
                    }
                    Flag::SetNotPrimaryAlignment(sam.m_flag, !first);
                    Flag::SetSupplementaryAlignment(sam.m_flag, first && primary_written);
                    sam.m_mapq = first ? segment.mapq : 0;
                    sam.m_cigar = LongReadCigar(hit);
                    if (!read_records.empty()) read_records += '\n';
                    read_records += sam.ToString();
                    // ZR:i:1: the best hit is the read's taxon's, which the gene alone could not tell.
                    if (first && segment.by_read) read_records += "\tZR:i:1";
                    m_genes.push_back(SamGeneKey(ar.Taxid(), ar.GeneId()));
                    primary_written = true;
                    first = false;
                    alignments++;
                    if (++written == m_max_out) break;
                }
            }

            // All records of the read go into the buffer as one unit, so a flush cannot let another
            // thread's records land between them.
            if (!read_records.empty() && !m_sam_output.Write(std::move(read_records))) {
                Flush();
            }
        }
    };
}
