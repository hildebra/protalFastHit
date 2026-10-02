#pragma once

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <memory>
#include <string>
#include <vector>
#include "ChainingStrategy.h"
#include "WFA2Wrapper2.h"
#include "TargetClones.h"

namespace protal {

    // The window a read is aligned into, as SimpleAlignmentHandler::AlignAnchor sets it up around
    // its anchor: gene positions [ref_start, ref_end), how many reference bases are free at each
    // end (the dovetails), how many read bases are free at each end (where the read runs past the
    // gene), and the penalty the alignment must stay below.
    struct AlignmentWindow {
        size_t ref_start = 0;
        size_t ref_end = 0;
        int ref_begin_free = 0;
        int ref_end_free = 0;
        int read_begin_free = 0;
        int read_end_free = 0;
        int max_score = 0;
    };

    // Aligns a read into its window from the exact matches of its anchor chain, instead of aligning
    // the whole read into the whole window. The links are taken as they are; the read left of the
    // first link is aligned to the window left of it (reversed, so that the alignment is anchored at
    // the link and its free end is the window's start), the read right of the last link likewise,
    // and between two links the bases in between. Most of a read's anchors are one link of 20-50
    // bases on a relative's gene: anchored at it, WFA gives up on such a flank as soon as it exceeds
    // the window's budget, where aligning the whole read started from every free diagonal at once.
    //
    // The result is the operations the whole-window alignment would give when its best path runs
    // through the links (M, X, I, D over the window, free ends included), so the caller
    // post-processes it the same way. The budget is shared: the pieces' penalties add up to less than
    // the window's max_score, as for a whole alignment. Chains this does not handle are left to the
    // whole-window alignment (NotApplicable): links on different diagonals (the seeds imply an
    // indel), out of order, outside the window, or not exact matches apart from Ns.
    //
    // Long reads (AllowIndels): their indels put the seeds of a gene on diagonals that drift, so the
    // links may lie on different diagonals. They are taken in read order, each cut where it overlaps
    // the one before (in the read or the gene), and the bases between two links on different
    // diagonals are aligned end to end. A long read's seeds cover a third of its gene or less (few
    // of the index's k-mers survive its errors), so the stretches before, between and after the links
    // get links of their own first (Reseed): exact matches of kSeedK bases, unique in a band of
    // diagonals around the one the neighbouring link leads to.
    class AnchoredAligner {
    public:
        enum class Status { Aligned, Failed, NotApplicable };

        AnchoredAligner() = default;
        // A copy has the settings, not the aligner of the pieces (one per thread, made when first needed).
        AnchoredAligner(AnchoredAligner const& other) : m_indels(other.m_indels) {}
        AnchoredAligner& operator=(AnchoredAligner const& other) {
            m_indels = other.m_indels;
            return *this;
        }

        void AllowIndels(bool allow) { m_indels = allow; }
        bool IndelsAllowed() const { return m_indels; }

        PROTAL_CLONE_V3 Status Align(std::string const& read, std::string_view const gene, ChainList const& input, AlignmentWindow const& w,
                     WFA2Wrapper2& aligner, std::string& ops) {
            ops.clear();
            if (input.empty() || w.ref_end > gene.size() || w.ref_start > w.ref_end) return Status::NotApplicable;
            ChainList const& chain = m_indels ? m_links : input;
            if (m_indels) {
                if (!Colinear(input, read.size(), w)) return Status::NotApplicable;
                Reseed(read, gene, w);
            } else {
                long const diagonal = static_cast<long>(chain.front().genepos) - static_cast<long>(chain.front().readpos);
                for (size_t i = 0; i < chain.size(); i++) {
                    auto const& link = chain[i];
                    if (link.length == 0 || static_cast<long>(link.genepos) - static_cast<long>(link.readpos) != diagonal) return Status::NotApplicable;
                    if (link.ReadEnd() > read.size() || link.GeneStart() < w.ref_start || link.GeneEnd() > w.ref_end) return Status::NotApplicable;
                    if (i > 0 && link.ReadStart() < chain[i - 1].ReadEnd()) return Status::NotApplicable;
                }
            }

            int used = 0;
            // Left of the first link, reversed: anchored at the link, free at the window's start.
            m_left.clear();
            auto const& first = chain.front();
            if (Flank(read, gene, 0, first.ReadStart(), w.ref_start, first.GeneStart(), true, w.ref_begin_free, w.read_begin_free,
                      w.max_score, aligner, used, m_left) != Status::Aligned) {
                return m_status;
            }
            ops.append(m_left.rbegin(), m_left.rend());

            for (size_t i = 0; i < chain.size(); i++) {
                auto const& link = chain[i];
                // The link: exact but for Ns, which count as mismatches as in any alignment.
                if (std::memcmp(read.data() + link.readpos, gene.data() + link.genepos, link.length) == 0) {
                    ops.append(link.length, 'M');  // nearly always: the seeds were exact matches
                } else {
                    for (uint32_t k = 0; k < link.length; k++) {
                        char const q = read[link.readpos + k], r = gene[link.genepos + k];
                        if (q == r) { ops += 'M'; continue; }
                        if (q != 'N' && r != 'N') return Status::NotApplicable;
                        ops += 'X';
                        if ((used += kMismatch) >= w.max_score) return Status::Failed;
                    }
                }
                if (i + 1 < chain.size()) {
                    auto const& next = chain[i + 1];
                    size_t const from = link.ReadEnd(), to = next.ReadStart();
                    if (to - from == next.GeneStart() - link.GeneEnd()) {
                        // Up to the next link on the same diagonal: as many read as gene bases.
                        if (Between(read, gene, from, to, link.GeneEnd(), w.max_score, aligner, used, ops) != Status::Aligned) return m_status;
                    } else if (Gap(read, gene, from, to, link.GeneEnd(), next.GeneStart(), w.max_score, aligner, used, ops) != Status::Aligned) {
                        return m_status;
                    }
                }
            }

            // Right of the last link: anchored at the link, free at the window's end.
            auto const& last = chain.back();
            m_right.clear();
            if (Flank(read, gene, last.ReadEnd(), read.size(), last.GeneEnd(), w.ref_end, false, w.ref_end_free, w.read_end_free,
                      w.max_score, aligner, used, m_right) != Status::Aligned) {
                return m_status;
            }
            ops += m_right;

            // The operations cover the read and the window exactly, as a whole-window alignment's do.
            size_t const read_bases = ops.size() - static_cast<size_t>(std::count(ops.begin(), ops.end(), 'D'));
            size_t const ref_bases = ops.size() - static_cast<size_t>(std::count(ops.begin(), ops.end(), 'I'));
            if (read_bases != read.size() || ref_bases != w.ref_end - w.ref_start) return Status::NotApplicable;
            return Status::Aligned;
        }

        // Tests: off, every flank goes to WFA2, which the ungapped flanks are checked against; and how many were.
        void SetUngappedFlanks(bool on) { m_ungapped_flanks = on; }
        size_t UngappedFlanks() const { return m_ungapped_count; }

    private:
        static constexpr int kMismatch = 4;  // as SimpleAlignmentHandler sets up WFA2Wrapper2 (4, 6, 2)

        static constexpr int kGapOpen = 6, kGapExtend = 2;
        // Re-seeding: exact matches of kSeedK bases, unique in a band kBand diagonals wide either side, widening by one
        // per kDrift bases from the link the search starts at (indels shift the diagonal), up to kMaxBand. A stretch
        // narrower than kMinStretch in the read or the gene is not searched.
        static constexpr int kSeedK = 12;
        static constexpr long kBand = 8, kDrift = 8, kMaxBand = 64;
        static constexpr size_t kMinStretch = 2 * kSeedK;

        std::string m_query, m_ref, m_left, m_right;
        Status m_status = Status::Aligned;
        bool m_ungapped_flanks = true;
        size_t m_ungapped_count = 0;
        bool m_indels = false;
        ChainList m_links, m_reseeded;           // long reads: the chain as aligned, and as it is re-seeded
        std::vector<int32_t> m_gene_codes, m_read_codes;
        std::unique_ptr<WFA2Wrapper2> m_pieces;  // long reads: the aligner of the pieces between links (Piece)

        // A flank: read [read_from, read_to) against gene [ref_from, ref_to), anchored at the link
        // (the flank's end for the left flank, its start for the right one) and free at the other
        // end by up to ref_free reference and read_free read bases. The left flank is aligned
        // reversed; its operations come back reversed too, for the caller to turn around.
        PROTAL_CLONE_V3 Status Flank(std::string const& read, std::string_view const gene, size_t read_from, size_t read_to, size_t ref_from,
                     size_t ref_to, bool reversed, int ref_free, int read_free, int max_score, WFA2Wrapper2& aligner,
                     int& used, std::string& out) {
            size_t const r = read_to - read_from, g = ref_to - ref_from;
            if (r == 0 && g == 0) return m_status = Status::Aligned;
            if (r == 0 || g == 0) {
                // Only reference or only read bases: free if the window allows that many at this end.
                if (r == 0 && g <= static_cast<size_t>(std::max(ref_free, 0))) { out.append(g, 'D'); return m_status = Status::Aligned; }
                if (g == 0 && r <= static_cast<size_t>(std::max(read_free, 0))) { out.append(r, 'I'); return m_status = Status::Aligned; }
                return m_status = Status::NotApplicable;
            }
            if (reversed) {
                m_query.assign(read.rbegin() + static_cast<long>(read.size() - read_to), read.rbegin() + static_cast<long>(read.size() - read_from));
                m_ref.assign(gene.rbegin() + static_cast<long>(gene.size() - ref_to), gene.rbegin() + static_cast<long>(gene.size() - ref_from));
            } else {
                m_query.assign(read, read_from, r);
                m_ref.assign(gene, ref_from, g);
            }
            int const read_end_free = std::min<int>(std::max(read_free, 0), static_cast<int>(r));
            int const ref_end_free = std::min<int>(std::max(ref_free, 0), static_cast<int>(g));
            // At most one mismatch on the link's diagonal up to the end of the read or the reference, and the bases
            // left at the other within its free end: that ungapped alignment is the only best one (0 or 4; any other
            // leaves the diagonal with a gap, 8 or more), so it is the one WFA2 returns, with the bases left over
            // written as WFA2 writes them (I for the read's, D for the reference's). WFA2 stops at a score of its
            // max_steps, so it fails where the score reaches the budget.
            size_t const m = std::min(r, g);
            if (m_ungapped_flanks && r - m <= static_cast<size_t>(read_end_free) && g - m <= static_cast<size_t>(ref_end_free)) {
                int mismatches = 0;
                for (size_t k = 0; k < m && mismatches <= 1; k++) mismatches += m_query[k] != m_ref[k];
                if (mismatches <= 1) {
                    if (kMismatch * mismatches >= max_score - used) return m_status = Status::Failed;
                    out.resize(m);
                    for (size_t k = 0; k < m; k++) out[k] = m_query[k] == m_ref[k] ? 'M' : 'X';
                    out.append(r - m, 'I');
                    out.append(g - m, 'D');
                    used += kMismatch * mismatches;
                    m_ungapped_count++;
                    return m_status = Status::Aligned;
                }
            }
            aligner.Reset();
            aligner.Alignment(m_query, m_ref, 0, read_end_free, 0, ref_end_free, max_score - used);
            if (!aligner.Success()) return m_status = Status::Failed;
            used += -aligner.GetAlignmentScore();
            aligner.CigarInto(out);
            return m_status = Status::Aligned;
        }

        // The bases between two links on one diagonal: read [from, to) against as many gene bases
        // from gene_from, both ends fixed. With up to 3 mismatches the ungapped alignment is the best
        // (12 against at least 16 for the insertion and deletion any gapped one needs); else WFA.
        Status Between(std::string const& read, std::string_view const gene, size_t from, size_t to, size_t gene_from, int max_score,
                       WFA2Wrapper2& aligner, int& used, std::string& ops) {
            size_t const n = to - from;
            if (n == 0) return m_status = Status::Aligned;
            int mismatches = 0;
            for (size_t k = 0; k < n; k++) mismatches += read[from + k] != gene[gene_from + k];
            if (mismatches <= 3) {
                size_t const old_size = ops.size();
                ops.resize(old_size + n);
                for (size_t k = 0; k < n; k++) ops[old_size + k] = read[from + k] == gene[gene_from + k] ? 'M' : 'X';
                used += kMismatch * mismatches;
                return m_status = used < max_score ? Status::Aligned : Status::Failed;
            }
            m_query.assign(read, from, n);
            m_ref.assign(gene, gene_from, n);
            return m_status = Piece(aligner, max_score, used, ops);
        }

        // The bases between two links on different diagonals: read [from, to) against gene [gene_from, gene_to),
        // both ends fixed. Only read or only gene bases are an insertion or a deletion; else WFA2 end to end.
        Status Gap(std::string const& read, std::string_view const gene, size_t from, size_t to, size_t gene_from, size_t gene_to,
                   int max_score, WFA2Wrapper2& aligner, int& used, std::string& ops) {
            size_t const r = to - from, g = gene_to - gene_from;
            if (r == 0 || g == 0) {
                size_t const n = r + g;
                ops.append(n, r == 0 ? 'D' : 'I');
                used += kGapOpen + kGapExtend * static_cast<int>(n);
                return m_status = used < max_score ? Status::Aligned : Status::Failed;
            }
            m_query.assign(read, from, r);
            m_ref.assign(gene, gene_from, g);
            return m_status = Piece(aligner, max_score, used, ops);
        }

        // m_query aligned to m_ref end to end, the operations appended to ops. Short reads on the window's aligner,
        // as before; long reads on one of their own, with WFA2's end-to-end kernels. WFA2 keeps the wavefronts of the
        // largest alignment it made and resets every one of them for each alignment, which costs a piece between
        // two links that shares the aligner with gene-long flanks more than its alignment does.
        Status Piece(WFA2Wrapper2& window_aligner, int max_score, int& used, std::string& ops) {
            if (m_indels && !m_pieces) m_pieces = std::make_unique<WFA2Wrapper2>(kMismatch, kGapOpen, kGapExtend, 0);
            WFA2Wrapper2& aligner = m_indels ? *m_pieces : window_aligner;
            aligner.Reset();
            if (m_indels) aligner.EndToEnd(m_query, m_ref, max_score - used);
            else aligner.Alignment(m_query, m_ref, 0, 0, 0, 0, max_score - used);
            if (!aligner.Success()) return Status::Failed;
            used += -aligner.GetAlignmentScore();
            aligner.CigarInto(ops, true);
            return Status::Aligned;
        }

        // The links of `input` in read order, each cut where it overlaps the one before (in the read or the gene,
        // dropped if nothing is left), into m_links. False if none is left, or one lies outside the read or window.
        bool Colinear(ChainList const& input, size_t read_length, AlignmentWindow const& w) {
            m_links.assign(input.begin(), input.end());
            std::stable_sort(m_links.begin(), m_links.end(), [](ChainLink const& a, ChainLink const& b) {
                return a.ReadStart() < b.ReadStart();
            });
            size_t kept = 0;
            for (size_t i = 0; i < m_links.size(); i++) {
                ChainLink link = m_links[i];
                if (link.length == 0) continue;
                if (link.ReadEnd() > read_length || link.GeneStart() < w.ref_start || link.GeneEnd() > w.ref_end) return false;
                if (kept > 0) {
                    auto const& prev = m_links[kept - 1];
                    long const cut = std::max<long>({ 0L, static_cast<long>(prev.ReadEnd()) - static_cast<long>(link.ReadStart()),
                                                      static_cast<long>(prev.GeneEnd()) - static_cast<long>(link.GeneStart()) });
                    if (cut >= static_cast<long>(link.length)) continue;
                    link.readpos = static_cast<uint16_t>(link.readpos + cut);
                    link.genepos = static_cast<uint32_t>(link.genepos + cut);
                    link.length = static_cast<uint16_t>(link.length - cut);
                }
                m_links[kept++] = link;
            }
            m_links.erase(m_links.begin() + static_cast<long>(kept), m_links.end());
            return kept > 0;
        }

        static int BaseCode(char c) {
            switch (c) {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
                default: return -1;
            }
        }

        // codes[i]: the 2-bit code of the k-mer of s at from + i, for every k-mer within [from, to); -1 if it holds a
        // base other than A, C, G, T.
        static void KmerCodes(char const* s, size_t from, size_t to, std::vector<int32_t>& codes) {
            codes.clear();
            if (to < from + kSeedK) return;
            uint32_t constexpr mask = (1u << (2 * kSeedK)) - 1;
            uint32_t code = 0;
            int valid = 0;
            for (size_t i = from; i < to; i++) {
                int const b = BaseCode(s[i]);
                if (b < 0) {
                    valid = 0;
                    code = 0;
                } else {
                    code = ((code << 2) | static_cast<uint32_t>(b)) & mask;
                    valid++;
                }
                if (i + 1 >= from + kSeedK) codes.push_back(valid >= kSeedK ? static_cast<int32_t>(code) : -1);
            }
        }

        // The band of diagonals searched at `distance` bases from the link the search starts at.
        static long Band(size_t distance) {
            return std::min<long>(kMaxBand, kBand + static_cast<long>(distance) / kDrift);
        }

        // The gene position of the only k-mer of m_gene_codes (from gene position g0) in [lo, hi] equal to c, or -1.
        long UniqueHit(int32_t c, long lo, long hi, size_t g0) const {
            long found = -1;
            for (long p = lo; p <= hi; p++) {
                if (m_gene_codes[static_cast<size_t>(p) - g0] != c) continue;
                if (found >= 0) return -1;
                found = p;
            }
            return found;
        }

        // Links in read [r0, r1) x gene [g0, g1), found left to right from diagonal d (gene - read position) of the
        // link before the stretch, appended to m_reseeded.
        void ScanRight(std::string const& read, std::string_view const gene, size_t r0, size_t r1, size_t g0, size_t g1, long d) {
            if (r1 < r0 + kMinStretch || g1 < g0 + kMinStretch) return;
            KmerCodes(gene.data(), g0, g1, m_gene_codes);
            KmerCodes(read.data(), r0, r1, m_read_codes);
            size_t last_end = r0, g_floor = g0;
            for (size_t q = r0; q + kSeedK <= r1; q++) {
                int32_t const c = m_read_codes[q - r0];
                if (c < 0) continue;
                long const band = Band(q - last_end);
                long const lo = std::max<long>(static_cast<long>(q) + d - band, static_cast<long>(g_floor));
                long const hi = std::min<long>(static_cast<long>(q) + d + band, static_cast<long>(g1) - kSeedK);
                long const hit = UniqueHit(c, lo, hi, g0);
                if (hit < 0) continue;
                size_t const p = static_cast<size_t>(hit);
                size_t len = kSeedK;
                while (q + len < r1 && p + len < g1 && read[q + len] == gene[p + len]) len++;
                m_reseeded.emplace_back(static_cast<uint32_t>(p), static_cast<uint16_t>(q), static_cast<uint16_t>(len));
                d = static_cast<long>(p) - static_cast<long>(q);
                last_end = q + len;
                g_floor = p + len;
                q = last_end - 1;
            }
        }

        // Links in read [r0, r1) x gene [g0, g1), found right to left from diagonal d of the link after the
        // stretch, appended to m_reseeded in read order.
        void ScanLeft(std::string const& read, std::string_view const gene, size_t r0, size_t r1, size_t g0, size_t g1, long d) {
            if (r1 < r0 + kMinStretch || g1 < g0 + kMinStretch) return;
            KmerCodes(gene.data(), g0, g1, m_gene_codes);
            KmerCodes(read.data(), r0, r1, m_read_codes);
            size_t const first = m_reseeded.size();
            size_t next_start = r1, g_ceil = g1;
            for (long q = static_cast<long>(r1) - kSeedK; q >= static_cast<long>(r0); q--) {
                int32_t const c = m_read_codes[static_cast<size_t>(q) - r0];
                if (c < 0) continue;
                long const band = Band(next_start - static_cast<size_t>(q + kSeedK));
                long const lo = std::max<long>(q + d - band, static_cast<long>(g0));
                long const hi = std::min<long>(q + d + band, static_cast<long>(g_ceil) - kSeedK);
                long const hit = UniqueHit(c, lo, hi, g0);
                if (hit < 0) continue;
                size_t qs = static_cast<size_t>(q), ps = static_cast<size_t>(hit), len = kSeedK;
                while (qs + len < next_start && ps + len < g_ceil && read[qs + len] == gene[ps + len]) len++;
                while (qs > r0 && ps > g0 && read[qs - 1] == gene[ps - 1]) { qs--; ps--; len++; }
                m_reseeded.emplace_back(static_cast<uint32_t>(ps), static_cast<uint16_t>(qs), static_cast<uint16_t>(len));
                d = static_cast<long>(ps) - static_cast<long>(qs);
                next_start = qs;
                g_ceil = ps;
                q = static_cast<long>(qs) - kSeedK + 1;
            }
            std::reverse(m_reseeded.begin() + static_cast<long>(first), m_reseeded.end());
        }

        // More links for m_links where it leaves stretches of the read and the window: before its first link,
        // between links, and after its last.
        void Reseed(std::string const& read, std::string_view const gene, AlignmentWindow const& w) {
            m_reseeded.clear();
            auto diagonal = [](ChainLink const& l) { return static_cast<long>(l.genepos) - static_cast<long>(l.readpos); };
            auto const& first = m_links.front();
            ScanLeft(read, gene, 0, first.ReadStart(), w.ref_start, first.GeneStart(), diagonal(first));
            for (size_t i = 0; i < m_links.size(); i++) {
                auto const& link = m_links[i];
                m_reseeded.push_back(link);
                size_t const r1 = i + 1 < m_links.size() ? m_links[i + 1].ReadStart() : read.size();
                size_t const g1 = i + 1 < m_links.size() ? m_links[i + 1].GeneStart() : w.ref_end;
                ScanRight(read, gene, link.ReadEnd(), r1, link.GeneEnd(), g1, diagonal(link));
            }
            std::swap(m_links, m_reseeded);
        }
    };
}
