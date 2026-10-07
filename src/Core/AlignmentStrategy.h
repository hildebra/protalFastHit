//
// Created by fritsche on 05/09/22.
//

#pragma once

#include "WFA2Wrapper2.h"
#include "AlignmentOutputHandler.h"
#include "GenomeLoader.h"
#include "SeedingStrategy.h"
#include "ChainingStrategy.h"
#include "KmerUtils.h"
#include "FastAlignment.h"
#include "SNPUtils.h"
#include "AlignmentUtils.h"
#include "AnchoredAlignment.h"
#include "AlignmentScreen.h"
#include "TargetClones.h"

namespace protal {
    struct AlignmentOrientation {
        int query_start = 0;
        int query_end = 0;
        int query_len = 0;
        int reference_start = 0;
        int reference_end = 0;
        int reference_len = 0;
        int overlap = 0;
        int query_dove_left = 0;
        int query_dove_right = 0;
        int reference_dove_left = 0;
        int reference_dove_right = 0;

        int MaxDoveLeft() {
            return std::max(query_dove_left, reference_dove_left);
        }

        int MaxDoveRight() {
            return std::max(query_dove_right, reference_dove_right);
        }

        std::string ToString() {
            std::string ret;
            ret += "query_start:           " + std::to_string(query_start) + '\n';
            ret += "query_end:             " + std::to_string(query_end) + '\n';
            ret += "query_len:             " + std::to_string(query_len) + '\n';
            ret += "reference_start:       " + std::to_string(reference_start) + '\n';
            ret += "reference_end:         " + std::to_string(reference_end) + '\n';
            ret += "reference_len:         " + std::to_string(reference_len) + '\n';

            ret += "overlap:               " + std::to_string(overlap) + '\n';
            ret += "query_dove_left:       " + std::to_string(query_dove_left) + '\n';
            ret += "query_dove_right:      " + std::to_string(query_dove_right) + '\n';
            ret += "reference_dove_left:   " + std::to_string(reference_dove_left) + '\n';
            ret += "reference_dove_right:  " + std::to_string(reference_dove_right);
            return ret;
        }

        void Update(int abs_pos, size_t read_length, size_t gene_length, int max_dove) {
            Update(abs_pos, abs_pos, read_length, gene_length, max_dove);
        }

        // As above, with the read's end placed on the gene by a diagonal of its own (abs_pos_right: where read
        // position 0 would lie by the chain's last link). Chains of long reads hold indels, which shift the
        // diagonal along the read: by the first link's diagonal, a read's bases past the gene's end can be more
        // or fewer than there really are, and the alignment either forces the extra bases into the window or
        // is left to the whole window (docs/claude/2026-10-02-performance-round3, item 7). With one diagonal
        // the two give the same window.
        void Update(int abs_pos, int abs_pos_right, size_t read_length, size_t gene_length, int max_dove) {
            // Starts and left dovetails
            query_start = abs_pos < 0 ? -abs_pos : 0;
            query_dove_left = (query_start - max_dove) < 0 ? query_start : max_dove;
            reference_start = std::max(abs_pos, 0);
            reference_dove_left = (reference_start - max_dove) < 0 ? reference_start : max_dove;

            // Ends: where the read's last base lies by the right diagonal, cut at the gene's end.
            long const read_end_on_gene = static_cast<long>(abs_pos_right) + static_cast<long>(read_length);
            long const past_end = std::max(0L, read_end_on_gene - static_cast<long>(gene_length));
            query_end = std::max<int>(query_start, static_cast<int>(static_cast<long>(read_length) - past_end));
            reference_end = std::max<int>(reference_start, static_cast<int>(std::min<long>(read_end_on_gene, static_cast<long>(gene_length))));
            overlap = query_end - query_start;

            // Right dovetails
            query_dove_right = (query_end + max_dove) > read_length ?
                               (read_length - query_end) : max_dove;
            reference_dove_right = (reference_end + max_dove) > gene_length ?
                                   (gene_length - reference_end) : max_dove;

            // Update starts and ends according to dovetails
            query_start -= query_dove_left;
            query_end += query_dove_right;
            reference_start -= reference_dove_left;
            reference_end += reference_dove_right;

            // Get query and reference lengths based on start and end
            query_len = query_end - query_start;
            reference_len = reference_end - reference_start;

//            std::cout << ToString() << std::endl;
        }
    };

    class SimpleAlignmentHandler {
    private:
        WFA2Wrapper2 m_aligner;
        AlignmentResult m_alignment_result;
        GenomeLoader& m_genome_loader;
        size_t m_kmer_size = 15;

        size_t m_align_top = 3;
        double m_max_score_ani = 0;

        bool m_fastalign = false;
        AlignmentOrientation m_alignment_orientation;

        // Reads are aligned from their anchor's exact matches (AnchoredAligner) where the chain
        // allows; off, every read is aligned as a whole into its window (SetAnchoredAlignment).
        bool m_anchored = false;
        AnchoredAligner m_anchored_aligner;
        // Before either method, a candidate whose read and window share too few k-mers for an alignment within
        // the budget to exist is refused (AlignmentScreen; off with SetAlignmentScreen(false)).
        bool m_screen_on = true;
        // Adaptive candidates (SetAdaptiveCandidates, TryCongeners): anchors beyond align_top a divergent read may try on
        // congeners of its best alignment, and every taxon's genus.
        size_t m_adaptive_candidates = 0;
        std::shared_ptr<std::vector<uint32_t> const> m_genera;
        AlignmentScreen m_screen;
        AlignmentScreen::ReadKmers m_read_kmers;  // operator()'s read, packed for its candidates' screens
        std::string m_ops;     // the window's alignment operations, from either method
        std::string m_reverse; // the reverse complement of the read, when the caller has none
        std::string m_window;  // the window's reference bases, for the whole-window alignment


        static bool IsReverse(LookupResult const& first_anchor, LookupResult const& second_anchor) {
            return first_anchor.readpos > second_anchor.readpos;
        }

        inline void ReverseAnchorReadPos(LookupResult& anchor, size_t& read_length) {
            anchor.readpos = read_length - anchor.readpos - m_kmer_size;
        }

        inline void ReverseAnchorReadPos(ChainLink& seed, size_t& read_length) {
            seed.readpos = read_length - seed.readpos - seed.length;
        }

        inline void ReverseAnchorPairReadPos(LookupResult& anchor_a, LookupResult& anchor_b, size_t& read_length) {
            ReverseAnchorReadPos(anchor_a, read_length);
            ReverseAnchorReadPos(anchor_b, read_length);
        }

    public:
        size_t total_alignments = 0;
        // The taxa of the anchors the last call tried to align (the align-top anchors and their ties), in order and
        // with repeats, each with its anchor's gene; with the alignments' taxa this gives the read's failed candidates
        // (FailedCandidates in AlignmentUtils.h).
        std::vector<FailedCandidate> m_attempted;
        std::vector<FailedCandidate> const& Attempted() const { return m_attempted; }
        // The taxa of the last call's anchors at least kCrowdedLength as long as its longest (exact-match bases), at most
        // UINT16_MAX: how many taxa the read's seeds could not tell apart, also those beyond align_top that were never
        // aligned and that ZA cannot list (its ZN tag; 0 without anchors).
        uint16_t m_crowding = 0;
        uint16_t Crowding() const { return m_crowding; }
        std::vector<uint32_t> m_crowded_scratch;  // CrowdedTaxa's
        // The taxa of the last call's crowd (anchors at least kCrowdedLength as long as the longest) that it never aligned
        // against, the strongest anchor first, at most kUntriedListed: the read's untried candidates (ZC,
        // AddUntriedCandidates).
        std::vector<uint32_t> m_untried;
        std::vector<uint32_t> const& Untried() const { return m_untried; }
        // Anchors aligned beyond align_top by the adaptive candidates (TryCongeners), over all calls.
        size_t m_adaptive_alignments = 0;
        // A read whose best alignment's identity (GetProxyANI) is below this is divergent: its crowd's congeners of that
        // alignment's taxon are tried too (TryCongeners).
        static constexpr double kAdaptiveIdentity = 0.99;

        // An anchor at least this share of the longest anchor's exact-match bases counts for the read's crowding.
        static constexpr double kCrowdedLength = 0.8;

        // The taxa of `anchors` (in any order) with an anchor at least kCrowdedLength as long as the longest (total_length),
        // at most UINT16_MAX; 0 without anchors.
        template<typename Anchors>
        static uint16_t CrowdedTaxa(Anchors const& anchors, std::vector<uint32_t>& taxa) {
            size_t longest = 0;
            for (auto const& anchor : anchors) longest = std::max<size_t>(longest, anchor.total_length);
            if (anchors.empty()) return 0;
            taxa.clear();
            for (auto const& anchor : anchors) {
                if (static_cast<double>(anchor.total_length) >= kCrowdedLength * static_cast<double>(longest)) {
                    taxa.push_back(static_cast<uint32_t>(anchor.taxid));
                }
            }
            std::sort(taxa.begin(), taxa.end());
            size_t const distinct = static_cast<size_t>(std::unique(taxa.begin(), taxa.end()) - taxa.begin());
            return static_cast<uint16_t>(std::min<size_t>(distinct, UINT16_MAX));
        }
        size_t total_tail_alignments = 0;
        size_t total_tail_length = 0;
        Benchmark bm_alignment{ "Alignment", 0, Benchmark::kPerRead};
        Benchmark m_bm_alignment {"Raw alignment", 0, Benchmark::kPerRead};
        Benchmark bm_seedext{ "Seed Extension", 0, Benchmark::kPerRead};
        // The k-mer screen of every candidate (part of the alignment handler's time): at GTDB scale it refuses 90-96% of
        // the candidates and is most of the handler (docs/claude/2026-10-06-performance-profiling).
        Benchmark m_bm_screen{ "K-mer screen", 0, Benchmark::kPerRead };
        size_t dummy = 0;
        // Anchors tried (AlignAnchor calls), of them those the k-mer screen refused before WFA2, those aligned
        // from their exact matches, and those aligned as a whole (anchored alignment off, or a chain it does not
        // handle); joined over threads like the benchmarks (JoinCounts).
        size_t m_attempted_alignments = 0;
        size_t m_screened_alignments = 0;
        size_t m_anchored_alignments = 0;
        size_t m_whole_window_alignments = 0;
        size_t m_validity_checks = 0;  // alignments so far: one in 64 is walked base by base (AlignAnchor)

        // Adds a thread's copy's counts, and its screen's time, to this (the global) handler's.
        void JoinCounts(SimpleAlignmentHandler const& other) {
            m_attempted_alignments += other.m_attempted_alignments;
            m_adaptive_alignments += other.m_adaptive_alignments;
            m_screened_alignments += other.m_screened_alignments;
            m_anchored_alignments += other.m_anchored_alignments;
            m_whole_window_alignments += other.m_whole_window_alignments;
            m_bm_screen.Join(other.m_bm_screen);
        }

//        AlignmentInfo m_info;

        void TotalReset() {
            total_alignments = 0;
            total_tail_alignments = 0;
            bm_alignment = Benchmark{ "Alignment", 0, Benchmark::kPerRead };
            bm_seedext = Benchmark{ "Seed Extension", 0, Benchmark::kPerRead };
            dummy = 0;
        }


        SimpleAlignmentHandler(GenomeLoader& genome_loader, WFA2Wrapper2& aligner, size_t kmer_size, size_t align_top, double max_score_ani, bool fastalign) :
                m_genome_loader(genome_loader),
                m_aligner(aligner),
                m_kmer_size(kmer_size),
                m_align_top(align_top),
                m_max_score_ani(max_score_ani),
                m_fastalign(fastalign) {};

        SimpleAlignmentHandler(SimpleAlignmentHandler const& other) :
                m_genome_loader(other.m_genome_loader),
                m_aligner(other.m_aligner),
                m_kmer_size(other.m_kmer_size),
                m_align_top(other.m_align_top),
                m_max_score_ani(other.m_max_score_ani),
                m_fastalign(other.m_fastalign),
                m_anchored(other.m_anchored),
                m_screen_on(other.m_screen_on),
                m_adaptive_candidates(other.m_adaptive_candidates),
                m_genera(other.m_genera) {
            m_anchored_aligner.AllowIndels(other.m_anchored_aligner.IndelsAllowed());
        };

        // Adaptive candidates (--adaptive_candidates): a divergent read tries up to `candidates` more anchors of its
        // crowd, of taxa of its best alignment's genus (`genera`: taxid -> genus, 0 for none). 0 or no genera: off.
        void SetAdaptiveCandidates(size_t candidates, std::shared_ptr<std::vector<uint32_t> const> genera) {
            m_adaptive_candidates = candidates;
            m_genera = std::move(genera);
        }

        void SetAnchoredAlignment(bool anchored) {
            m_anchored = anchored;
        }

        // The k-mer screen before WFA2 (AlignmentScreen): on by default; off, every candidate is aligned.
        void SetAlignmentScreen(bool on) {
            m_screen_on = on;
        }

        // Anchored alignment through chains whose links lie on different diagonals (long reads, AnchoredAligner).
        void SetAnchoredIndels(bool indels) {
            m_anchored_aligner.AllowIndels(indels);
        }

        bool AnchoredAlignment() const {
            return m_anchored;
        }



        inline bool ReverseAnchor(AlignmentAnchor& anchor, size_t read_len) {
            bool reversed = false;
            if (IsReverse(anchor.a, anchor.b)) {
                reversed = true;
                ReverseAnchorPairReadPos(anchor.a, anchor.b, read_len);
            }
            return reversed;
        }

        inline void ReverseAnchor(CAlignmentAnchor& anchor, size_t read_len) {
            for (auto& chainlink : anchor.chain) {
                ReverseAnchorReadPos(chainlink, read_len);
            }
            std::reverse(anchor.chain.begin(), anchor.chain.end());
        }

        inline int SeedIndel(CAlignmentAnchor& anchor) {
            if (anchor.chain.size() == 1) return 0;
            int rpos_diff = static_cast<int>(anchor.Back().readpos) - static_cast<int>(anchor.Front().readpos);
            int gpos_diff = static_cast<int>(anchor.Back().genepos) - static_cast<int>(anchor.Front().genepos);
            return rpos_diff - gpos_diff;
        }



        size_t ExtendSeedLeft(Seed const& s, std::string const& query, std::string_view const gene) {
            size_t extension = 0;
            size_t max_extension_len = std::min(static_cast<uint32_t>(s.readpos), s.genepos);
            for (int qpos = s.readpos, rpos = s.genepos;
                 extension < max_extension_len && query[qpos] == gene[rpos];
                 qpos--, rpos--) {
                extension++;
            }
            return extension;
        }

        size_t ExtendSeedRight(Seed const& s, size_t k, std::string const& query, std::string_view const gene) {
            size_t extension = 0;
            size_t max_extension_len = std::min(
                    static_cast<uint32_t>(query.length() - s.readpos - k),
                    static_cast<uint32_t>(gene.length() - s.genepos - k));
            for (int qpos = s.readpos + k, rpos = s.genepos + k;
                 extension < max_extension_len && query[qpos] == gene[rpos];
                 qpos++, rpos++) {
                extension++;
            }
            return extension;
        }

        static std::pair<std::string_view, std::string_view>
        GetViewsForAlignment(Seed const& a, std::string const& query, std::string const& target) {

            auto min = std::min(static_cast<uint32_t>(a.readpos), a.genepos);
            auto qstart = a.readpos - min;
            auto tstart = a.genepos - min;
            auto overlap_length = std::min(query.length() - qstart, target.length() - tstart);

            if (tstart + overlap_length > target.length()) exit(9);
            if (qstart + overlap_length > query.length()) exit(10);

            return {
                std::string_view(query.c_str() + qstart, overlap_length),
                std::string_view(target.c_str() + tstart, overlap_length)
            };
        }

        static std::pair<size_t, size_t> ExtendSeed(Seed const& s, size_t k, std::string const& query, std::string_view const gene) {
            size_t extension_left = 0;
            size_t extension_right = 0;

            size_t max_extension_len = std::min(static_cast<uint32_t>(s.readpos), s.genepos);
            for (int qpos = s.readpos - 1, rpos = s.genepos - 1;
                 extension_left < max_extension_len && query[qpos] == gene[rpos];
                 qpos--, rpos--) {
                extension_left++;
            }
            max_extension_len = std::min(
                    static_cast<uint32_t>(query.length() - s.readpos - k),
                    static_cast<uint32_t>(gene.length() - s.genepos - k));
            for (int qpos = s.readpos + k, rpos = s.genepos + k;
                 extension_right < max_extension_len && query[qpos] == gene[rpos];
                 qpos++, rpos++) {
                extension_right++;
            }
            return { extension_left, extension_right };
        }

        static std::pair<size_t, size_t> ExtendSeed(ChainLink& s, std::string const& query, std::string_view const gene, uint16_t query_left_limit=0, uint16_t query_right_limit=0) {
            size_t extension_left = 0;
            size_t extension_right = 0;
            if (query_right_limit == 0) query_right_limit = query.length();

//            std::cout << s.ToString() << std::endl;
//            std::cout << "Start: " << std::endl;
//            std::cout << "qpos:  " << s.readpos - 1 << std::endl;
//            std::cout << "rpos:  " << s.genepos - 1 << std::endl;
//            std::cout << "query_left_limit:  " << query_left_limit << std::endl;
//            std::cout << "query_right_limit:  " << query_right_limit << std::endl;
//            std::cout << "ref_left_limit:  " << 0 << std::endl;
//            std::cout << "ref_right_limit:  " << query_right_limit << std::endl;
            for (int qpos = s.readpos - 1, rpos = s.genepos - 1;
                    qpos >= query_left_limit && rpos >= 0 && query[qpos] == gene[rpos];
                    qpos--, rpos--) {
                extension_left++;
            }

            for (int qpos = s.readpos, rpos = s.genepos;
                 qpos < qpos + s.length && rpos < gene.length();
                 qpos++, rpos++) {

                auto invalid = (query[qpos] == 'A' || query[qpos] == 'C' || query[qpos] == 'G' || query[qpos] == 'T') &&
                        (gene[rpos] == 'A' || gene[rpos] == 'C' || gene[rpos] == 'G' || gene[rpos] == 'T') && query[qpos] != gene[rpos];

                if (invalid) {
                    std::cerr << "____________________________" << std::endl;
                    std::cerr << query << std::endl;
                    std::cerr << string_view(gene.data() + s.genepos, s.length) << std::endl;
                    std::cerr << "____________________________" << std::endl;
                }
            }


//            std::cout << "extension_left: " << extension_left << std::endl;
//            size_t max_extension_len = std::min(static_cast<uint32_t>(s.readpos), s.genepos);
//            for (int qpos = s.readpos - 1, rpos = s.genepos - 1;
////                 qpos >= 0 && rpos >= 0 && // Maybe remove
//                 extension_left < max_extension_len && query[qpos] == gene[rpos];
//                 qpos--, rpos--) {
//                extension_left++;
//            }
//            max_extension_len = std::min(
//                    static_cast<uint32_t>(query.length() - s.readpos - s.length),
//                    static_cast<uint32_t>(gene.length() - s.genepos - s.length));

            for (int qpos = s.readpos + s.length, rpos = s.genepos + s.length;
                    qpos < query_right_limit && rpos < gene.length() && query[qpos] == gene[rpos];
                    qpos++, rpos++) {
                extension_right++;
            }
            assert(s.ReadStart() >= extension_left);
            assert(s.GeneStart() >= extension_left);
            assert(s.ReadEnd() + extension_right <= query.length());
            assert(s.GeneEnd() + extension_right <= gene.length());

            s.ExtendLength(extension_left, extension_right);
            return { extension_left, extension_right };
        }


        static int MaxScore(double min_ani, uint32_t overlap_length, int mismatch_penalty=4) {
            return (std::ceil((1 - min_ani) * static_cast<double>(overlap_length)) * mismatch_penalty) + 1;
        }


        void ExtendAnchor(ChainAlignmentAnchor& anchor, std::string const& ref) {
            auto& gene = m_genome_loader.GetGeneOMP(anchor.taxid, anchor.geneid);
            auto const geneseq = gene.Sequence();

            for (auto i = 0; i < anchor.chain.size(); i++) {
                auto& seed = anchor.chain[i];
                auto [lefta, righta] = ExtendSeed(seed, ref, geneseq,
                                                  ((i > 0) ? anchor.chain[i-1].readpos + anchor.chain[i-1].length : 0),
                                                  ((i+1) < anchor.chain.size() ? anchor.chain[i+1].readpos : ref.length()));


                if (i > 0 && seed.OverlapsWithLeft(anchor.chain[i-1])) {
                    anchor.chain[i-1].Merge(seed.readpos, seed.length);
                    anchor.chain.erase(anchor.chain.begin() + i);
                    i--;
                }
            }
        }



        void ExtendAllAnchors(AlignmentAnchorList const& anchors, std::string const& fwd, std::string const& rev) {
            for (auto anchor : anchors) {
                bool reversed = !anchor.forward;
                ExtendAnchor(anchor, reversed ? rev : fwd);
            }
        }


        static std::tuple<size_t, size_t, size_t> GetAlignmentPositions(CAlignmentAnchor& anchor, size_t gene_length, size_t read_length) {
            int abs_pos = anchor.Front().genepos - anchor.Front().readpos;
            size_t read_start = abs_pos < 0 ? -abs_pos : 0;
            size_t gene_start = std::max(abs_pos, 0);
            size_t overlap = std::min(gene_length - gene_start, read_length - read_start);
            return { read_start, gene_start, overlap };
        }

        inline std::pair<int, std::string> Approximate(std::string &query, std::string &gene) {
//            auto [ascore, acigar] = FastAligner::FastAlign(query, gene);
//            cigar_ani = CigarANI(acigar);
//            cigar = acigar;
//            score = ascore;


            return { 0, "" };
        }

//        inline std::tuple<int, bool, std::string> Optimal(std::string_view query, std::string_view gene, int max_score) {
//            m_aligner.Alignment(query, gene, max_score);
//            return { m_aligner.GetAligner().getAlignmentScore(), m_aligner.Success(), m_aligner.GetAligner().getAlignmentCigar() };
//        }

        SamEntry sam;
        SNPList snps;
        FastxRecord record;
        // max_score_cap: a budget below the ANI floor's (LongReadAligner, --long_read_budget); an alignment that reaches it fails.
        // stretch: where fwd and rev lie in a read whose packed strands its candidates share (AlignmentScreen::ReadStretch;
        // none: the screen packs the candidate's read for it alone).
        PROTAL_CLONE_V3 bool AlignAnchor(Anchor& anchor, AlignmentResult& alignment, std::string const& fwd, std::string const& rev, bool allow_heuristic_alignment, std::string& id, int max_score_cap = INT32_MAX,
                                         AlignmentScreen::ReadStretch const& stretch = {}) {
            alignment.GetAlignmentInfo().Reset();
            alignment.Reset();
            m_aligner.Reset();
            m_attempted_alignments++;

//            auto& read = anchor.forward ? fwd : rev;
            auto& read = anchor.forward ? fwd : rev;


            // Get Resources
            auto& gene = m_genome_loader.GetGeneOMP(anchor.taxid, anchor.geneid);

//            std::cerr << "--------------- links: " << anchor.chain.size() << std::endl;
//            auto [qry, ref] = anchor.ToVisualString(read, geneseq);
//            std::cerr << qry << std::endl;
//            std::cerr << ref << std::endl;
//            auto genestart = anchor.chain.front().genepos;
//            auto readstart = anchor.chain.front().readpos;
//            std::cerr << read.substr(readstart, std::min(readstart + read.length(), read.length())) << std::endl;
//            std::cerr << geneseq.substr(genestart, std::min<size_t>(read.length(), window.ref_end > genestart ? window.ref_end - genestart : 0)) << std::endl;  // the decoded window
//            std::cerr << "---------------" << std::endl;



            // If complete anchor matches without errors dont even do alignment.
            //TODO: Update this and allow anchor.total_length + #anchors + 1 == read.length()
//            if (anchor.total_length == read.length()) {
//                auto& info = alignment.GetAlignmentInfo();
//                info.cigar = std::string(read.length(), 'M');
//                info.ResetCigarStats();
//                info.alignment_length = read.length();
//                info.compressed_cigar = std::to_string(read.length()) + 'M';
//                info.gene_alignment_start = anchor.Front().genepos;
//                info.alignment_ani = info.GetProxyANI();
//                info.matches = read.length();
//                info.UpdateScore();
//
//                alignment.Set(anchor.taxid, anchor.geneid, info.gene_alignment_start, anchor.forward, anchor.unique, anchor.unique_best_two);
//
//                bool valid = IsAlignmentValid(info, read, gene.Sequence(), true);
//
//                if (!valid) {
//                    std::cerr << "Invalid no alignment\t" << id <<  std::endl;
//                    std::cerr << "Query Name: '" << id  << "'" << std::endl;
//                    std::cerr << "Info: " << info.ToString() << std::endl;
//                    std::cerr << "---------" << std::endl;
//                    std::cerr << read << std::endl;
////                    std::cerr << gene. << std::endl;
//                    std::cerr << anchor.ToString() << std::endl;
//                    std::cerr << anchor.ToVisualString2() << std::endl;
//
//                    int abs_pos = static_cast<int>(anchor.Front().genepos) - static_cast<int>(anchor.Front().readpos);
//                    m_alignment_orientation.Update(abs_pos, read.length(), gene.GetLength(), 0);
//                    std::string reference_str = gene.Sequence().substr(m_alignment_orientation.reference_start, m_alignment_orientation.reference_len);
//                    std::cerr << reference_str << std::endl;
//
//                    return false;
//                }
//
//                return true;
//            }

            // Is indel between anchor seeds?
            int anchor_indels = SeedIndel(anchor);


            // Absolute read positioning with respect to gene
            int abs_pos = static_cast<int>(anchor.Front().genepos) - static_cast<int>(anchor.Front().readpos);
            // Long reads (chains with indels, AnchoredAligner): the read's end is placed by the last link's
            // diagonal, so that the bases past the gene's end are the free ones by the end's own diagonal
            // (AlignmentOrientation::Update). Short reads keep the one diagonal of their chain.
            int const abs_pos_right = m_anchored_aligner.IndelsAllowed() ?
                    static_cast<int>(anchor.Back().genepos) - static_cast<int>(anchor.Back().readpos) : abs_pos;

            size_t max_dove_size = 9;
            m_alignment_orientation.Update(abs_pos, abs_pos_right, read.length(), gene.GetLength(), max_dove_size);

            assert(m_alignment_orientation.query_start + m_alignment_orientation.query_len <= read.length());
            assert(m_alignment_orientation.reference_start + m_alignment_orientation.reference_len <= gene.GetLength());
            assert(m_alignment_orientation.query_start >= 0);
            assert(m_alignment_orientation.reference_start >= 0);


            bool approximate_alignment = allow_heuristic_alignment && anchor_indels == 0;
            bool dove_left_required = anchor.Front().readpos != 0 || abs_pos < 0;
            bool dove_right_required = anchor.Back().readpos + anchor.Back().length != read.length() || abs_pos_right + read.length() > gene.GetLength();

            std::string cigar = "";

            int allowed_del_left = abs_pos < 0 ? (-1 * abs_pos) + 9 : 0;
            int allowed_del_right = abs_pos_right + read.length() > gene.GetLength() ? (abs_pos_right + read.length() - gene.GetLength()) + 9 : 0;

            m_alignment_orientation.reference_start += !dove_left_required * m_alignment_orientation.reference_dove_left;
            m_alignment_orientation.reference_end -= !dove_right_required * m_alignment_orientation.reference_dove_right;
            m_alignment_orientation.reference_len = m_alignment_orientation.reference_end - m_alignment_orientation.reference_start;

//            char* reference_cstr = const_cast<char *>(gene.Sequence().c_str() +
//                                                      m_alignment_orientation.reference_start);

//            std::string_view query_view(read.c_str(), read.length());
//            std::string_view reference_view(const_cast<char *>(gene.Sequence().c_str() +
//                                                          m_alignment_orientation.reference_start),
//                                       m_alignment_orientation.reference_len);

//            std::string query_view(read.c_str(), read.length());
//            std::string reference_view(const_cast<char *>(gene.Sequence().c_str() + m_alignment_orientation.reference_start),
//                                       m_alignment_orientation.reference_len);

            AlignmentWindow window;
            window.ref_start = m_alignment_orientation.reference_start;
            window.ref_end = m_alignment_orientation.reference_end;
            window.ref_begin_free = dove_left_required ? m_alignment_orientation.reference_dove_left * 2 : 0;
            window.ref_end_free = dove_right_required ? m_alignment_orientation.reference_dove_right * 2 : 0;
            window.read_begin_free = allowed_del_left;
            window.read_end_free = allowed_del_right;
            window.max_score = std::min(MaxScore(m_max_score_ani, m_alignment_orientation.overlap), max_score_cap);

            // Too few shared k-mers for any alignment within the budget: the candidate fails here as it would
            // have in WFA2 (AlignmentScreen), at a pass over the read and the window instead of the whole budget.
            // The window's k-mers are taken from the gene's packed bytes, so a refused candidate (most of them at GTDB
            // scale) is never decoded; a gene without them (not loaded) goes through its decoded window, which is empty.
            // The read's come from its strands packed once for all its candidates (stretch), where the caller has them.
            if (m_screen_on) {
                m_bm_screen.Start();
                size_t const begin_free = static_cast<size_t>(std::max(window.read_begin_free, 0));
                size_t const end_free = static_cast<size_t>(std::max(window.read_end_free, 0));
                uint8_t const* const packed = gene.Packed();
                bool const may_align = packed == nullptr ?
                        m_screen.MayAlign(read, begin_free, end_free,
                                          gene.Window(window.ref_start, window.ref_end).substr(window.ref_start, window.ref_end - window.ref_start),
                                          window.max_score, m_aligner.Mismatch(), m_aligner.GapOpening(), m_aligner.GapExtension()) :
                        stretch.kmers != nullptr ?
                        m_screen.MayAlignPacked(*stretch.kmers, anchor.forward, anchor.forward ? stretch.fwd_offset : stretch.rev_offset,
                                                read.length(), begin_free, end_free, packed, packed::Bytes(gene.GetLength()),
                                                window.ref_start, window.ref_end, window.max_score, m_aligner.Mismatch(),
                                                m_aligner.GapOpening(), m_aligner.GapExtension()) :
                        m_screen.MayAlignPacked(read, begin_free, end_free, packed, packed::Bytes(gene.GetLength()), window.ref_start,
                                                window.ref_end, window.max_score, m_aligner.Mismatch(), m_aligner.GapOpening(),
                                                m_aligner.GapExtension());
                m_bm_screen.Stop();
                if (!may_align) {
                    m_screened_alignments++;
                    return false;
                }
            }

            // The gene, decoded where the read lies: the window, which holds every link and flank the alignment reads
            // (AnchoredAligner checks that each link is inside it) and where the alignment is checked; it lives to the
            // end of this function. The rest of the gene is not decoded.
            auto const geneseq = gene.Window(window.ref_start, window.ref_end);

            bm_alignment.Start();
            if (approximate_alignment) {
//                auto [successt, scoret, cigart] = Align(anchor, read, gene.Sequence());
//                if (!successt)  {
//                    bm_alignment.Stop();
//                    return false;
//                }
//                exit(123);
//
//                cigar = cigart;
//                score = scoret;
//                cigar_ani = CigarANI(cigart);
            } else {
                // From the anchor's exact matches (AnchoredAligner) where its chain allows, else the
                // whole read into the whole window. Both give the window's operations, free ends
                // included, which are post-processed alike.
                m_bm_alignment.Start();
                using Status = AnchoredAligner::Status;
                Status status = m_anchored ? m_anchored_aligner.Align(read, geneseq, anchor.chain, window, m_aligner, m_ops)
                                           : Status::NotApplicable;
                if (status == Status::NotApplicable) {
                    m_window.assign(geneseq, window.ref_start, window.ref_end - window.ref_start);
                    m_aligner.Reset();
                    m_aligner.Alignment(read, m_window, window.read_begin_free, window.read_end_free,
                                        window.ref_begin_free, window.ref_end_free, window.max_score);
                    status = m_aligner.Success() ? Status::Aligned : Status::Failed;
                    if (status == Status::Aligned) m_aligner.CigarInto(m_ops);
                    m_whole_window_alignments++;
                } else {
                    m_anchored_alignments++;
                }
                m_bm_alignment.Stop();
//                else {
//                    std::cout << "End2End" << std::endl;
//                    std::cout << anchor.ToString() << std::endl;
//                    std::cout << read << std::endl;
//                    std::cout << reference_str << std::endl;
//                    m_aligner.Alignment(read, reference_str,
//                                        MaxScore(m_max_score_ani, m_alignment_orientation.overlap));
//                    auto& cigar = m_aligner.GetAligner().GetAligner()->cigar;
//                    std::cout << "CIG: " << std::string(cigar->operations + cigar->begin_offset, cigar->end_offset - cigar->begin_offset) << std::endl;
//
//                    if (m_aligner.Success()) {
//                        std::cout << m_aligner.GetAlignmentCigar() << std::endl;
//                    }
//                }

                if (status != Status::Aligned) {
                    bm_alignment.Stop();
                    return false;
                }


                auto& info = alignment.GetAlignmentInfo();
                PostProcessAlignment(m_ops, info, read.length(),
                                     gene.GetLength(), m_alignment_orientation.reference_start, 0, abs_pos);

                info.UpdateScore();
                alignment.Set(anchor.taxid, anchor.geneid, info.gene_alignment_start, anchor.forward, anchor.unique, anchor.unique_best_two);

                // Safety net: reject an alignment whose CIGAR does not fit the sequences. Its counts must cover the
                // read and stay inside the gene (every alignment, from the counts alone); one alignment in 64 is also
                // walked base by base against the sequences (IsAlignmentValid), as every one was before. Lock-free
                // on the common, valid path; the rare diagnostics below are serialised.
                bool const covers = static_cast<size_t>(info.matches) + info.mismatches + info.insertions + info.softclips == read.length() &&
                                    info.gene_alignment_start >= 0 &&
                                    static_cast<size_t>(info.gene_alignment_start) + info.matches + info.mismatches + info.deletions <= geneseq.size();
                bool const valid = covers && ((++m_validity_checks & 63) != 0 || IsAlignmentValid(info, read, geneseq, 0, true));
                if (!valid) {
#pragma omp critical (invalid_align)
                    {
                        std::cerr << "Info:" << info.ToString() << std::endl;
                        std::cerr << "Record: " << id << std::endl;
                        std::cerr << "Anchor: " << anchor.ToString() << std::endl;

                        auto readstart = info.read_start_offset;
                        auto genestart = info.gene_alignment_start;
                        std::cerr << read.substr(readstart, std::min(readstart + read.length(), read.length())) << std::endl;
                        // only the window is decoded
                        std::cerr << geneseq.substr(genestart, std::min<size_t>(read.length(), window.ref_end > genestart ? window.ref_end - genestart : 0)) << std::endl;

                        std::cerr << "-----Invalid after alignment\t" << geneseq.substr(window.ref_start, window.ref_end - window.ref_start)
                                  << " " << info.ToString() << std::endl;
                    }
                    bm_alignment.Stop();
                    return false;
                }

                // (The operations cover the read: `covers` above checked it from the counts.)
            }
            bm_alignment.Stop();
            return true;
        }

        void operator() (AlignmentAnchorList& anchors, AlignmentResultList& results, std::string const& sequence, size_t align_top, std::string& header) {
            KmerUtils::ReverseComplementInto(sequence, m_reverse);
            (*this)(anchors, results, sequence, m_reverse, align_top, header);
        }

        // As above for a read whose reverse complement the caller has (the anchor finder computed it): it
        // is the same for every anchor of the read, and nothing here copies the read.
        void operator() (AlignmentAnchorList& anchors, AlignmentResultList& results, std::string const& sequence, std::string const& reverse, size_t align_top, std::string& header) {
//            std::string header = "";
            constexpr bool alignment_verbose = false;

            // The read as fwd and its reverse complement as rev
            std::string const& fwd = sequence;
            std::string const& rev = reverse;

            bool reversed = false;

            size_t read_len = sequence.length();

            int take_top = align_top;

            Anchor* last_anchor = nullptr;
            m_attempted.clear();
            m_untried.clear();
            // rev is fwd's reverse complement as KmerUtils::ReverseComplementInto writes it (here or the anchor finder's).
            m_read_kmers.Set(fwd, rev);
            AlignmentScreen::ReadStretch const stretch{ &m_read_kmers, 0, 0 };
            m_crowding = CrowdedTaxa(anchors, m_crowded_scratch);
            size_t taken = 0;  // the anchors the loop took, from the front

            for (auto& anchor : anchors) {
                if (--take_top < 0 && (last_anchor && last_anchor->total_length != anchor.total_length)) {
                    break;
                }
                taken++;
                m_attempted.emplace_back(static_cast<uint32_t>(anchor.taxid), static_cast<uint32_t>(anchor.geneid));

                auto& read = anchor.forward ? fwd : rev;

                m_alignment_result.GetAlignmentInfo().Reset();

                auto success = AlignAnchor(anchor, m_alignment_result, fwd, rev, false, header, INT32_MAX, stretch);

                if (success && m_alignment_result.GetAlignmentInfo().GetProxyANI() >= m_max_score_ani) {
                    // Moved, not copied: AlignAnchor sets every field again for the next anchor.
                    results.emplace_back(std::move(m_alignment_result));
                    total_alignments++;
                }


                last_anchor = &anchor;
            }

            if (taken < anchors.size() && m_crowding > 0) {
                size_t longest = 0;
                for (auto const& anchor : anchors) longest = std::max<size_t>(longest, anchor.total_length);
                if (m_adaptive_candidates > 0 && m_genera && !results.empty()) {
                    TryCongeners(anchors, taken, longest, results, fwd, rev, header, stretch);
                }
                // The crowd's taxa never aligned against (ZC).
                for (size_t i = taken; i < anchors.size() && m_untried.size() < kUntriedListed; i++) {
                    auto const& anchor = anchors[i];
                    if (static_cast<double>(anchor.total_length) < kCrowdedLength * static_cast<double>(longest)) continue;
                    auto const taxid = static_cast<uint32_t>(anchor.taxid);
                    if (std::any_of(m_attempted.begin(), m_attempted.end(), [taxid](FailedCandidate const& c) { return c.taxid == taxid; })) continue;
                    if (std::find(m_untried.begin(), m_untried.end(), taxid) != m_untried.end()) continue;
                    m_untried.push_back(taxid);
                }
            }

            // Sort alignment results
            std::sort(results.begin(), results.end(), [](AlignmentResult const& a, AlignmentResult const& b) {
                return a.AlignmentScore() > b.AlignmentScore();
            });
        }

    private:
        uint32_t GenusOf(uint32_t taxid) const {
            return m_genera && taxid < m_genera->size() ? (*m_genera)[taxid] : 0;
        }

        // The adaptive candidates: a read whose best alignment so far is divergent (identity below kAdaptiveIdentity)
        // also tries the anchors of its crowd beyond the `taken` ones (at least kCrowdedLength of the `longest`) whose
        // taxon is of the best alignment's genus and not tried yet, in their order, up to m_adaptive_candidates: a strain
        // whose genes lie between its own species' reference and a congener's ranks its own species below the congeners'
        // by seeds, and --align_top alone never aligns it there (docs/claude/2026-10-07-error-read-signatures, section 5).
        // Per read, so the result does not depend on other reads or on the threads.
        template<typename Anchors>
        void TryCongeners(Anchors& anchors, size_t taken, size_t longest, AlignmentResultList& results, std::string const& fwd,
                          std::string const& rev, std::string& header, AlignmentScreen::ReadStretch const& stretch) {
            auto const best = std::max_element(results.begin(), results.end(), [](AlignmentResult const& a, AlignmentResult const& b) {
                return a.AlignmentScore() < b.AlignmentScore();
            });
            if (best->GetAlignmentInfo().GetProxyANI() >= kAdaptiveIdentity) return;
            uint32_t const genus = GenusOf(static_cast<uint32_t>(best->Taxid()));
            if (genus == 0) return;
            size_t tried = 0;
            for (size_t i = taken; i < anchors.size() && tried < m_adaptive_candidates; i++) {
                auto& anchor = anchors[i];
                if (static_cast<double>(anchor.total_length) < kCrowdedLength * static_cast<double>(longest)) continue;
                auto const taxid = static_cast<uint32_t>(anchor.taxid);
                if (GenusOf(taxid) != genus) continue;
                if (std::any_of(m_attempted.begin(), m_attempted.end(), [taxid](FailedCandidate const& c) { return c.taxid == taxid; })) continue;
                m_attempted.emplace_back(taxid, static_cast<uint32_t>(anchor.geneid));
                m_alignment_result.GetAlignmentInfo().Reset();
                if (AlignAnchor(anchor, m_alignment_result, fwd, rev, false, header, INT32_MAX, stretch) &&
                    m_alignment_result.GetAlignmentInfo().GetProxyANI() >= m_max_score_ani) {
                    results.emplace_back(std::move(m_alignment_result));
                    total_alignments++;
                }
                tried++;
                m_adaptive_alignments++;
            }
        }

    public:
    };
}