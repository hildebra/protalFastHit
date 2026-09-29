//
// Created by fritsche on 01/09/22.
//

#pragma once


#include <algorithm>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>
#include "Constants.h"
#include "BufferedOutput.h"
#include "VarkitInterface.h"
#include "KmerLookup.h"
#include "FastxReader.h"
#include "SamHandler.h"
#include "SNP.h"
#include "AlignmentUtils.h"
#include "SNPUtils.h"
#include "GenomeLoader.h"

namespace protal {
    static bool CorrectOrientation(AlignmentResult const& a1, AlignmentResult const& a2) {
        if (!a1.IsSet() || !a2.IsSet()) {
            std::cout << "CorrectOrientation both alignments needs to be set" << std::endl;
        }
        return a1.Forward() != a2.Forward();
    }

    static int ScorePairedAlignment(AlignmentResult const& a, AlignmentResult const& b) {
        int score = 0;
        if (a.IsSet() && b.IsSet()) {
            return (a.AlignmentScore() + b.AlignmentScore()) / 2.2;
        } else if (a.IsSet()) {
            return a.AlignmentScore() - 2;
        } else if (b.IsSet()) {
            return b.AlignmentScore() - 2;
        }
        return score;
    }

    static int ScorePairedAlignment(PairedAlignment const& pair) {
        return ScorePairedAlignment(pair.first, pair.second);
    }

    static double Score(AlignmentInfo const& info1, AlignmentInfo const& info2) {
        double score_sum = std::abs(info1.alignment_score + info2.alignment_score) * 0.25;
        double length_sum = static_cast<double>(info1.alignment_length + info2.alignment_length);

        return score_sum/length_sum;
    }

    static double Score(AlignmentInfo const& info) {
        double score_sum = std::abs(info.alignment_score) * 0.25;
        double length_sum = static_cast<double>(info.alignment_length);

        return score_sum/length_sum;
    }

    static double Score(PairedAlignment const& a) {
        return Score(a.first.GetAlignmentInfo(), a.second.GetAlignmentInfo());
    }


    static bool PairedAlignmentComparator(PairedAlignment const& a, PairedAlignment const& b) {
        // Return true if a > b
        AlignmentInfo info_a_1 = a.first.GetAlignmentInfo();
        AlignmentInfo info_a_2 = a.second.GetAlignmentInfo();
        AlignmentInfo info_b_1 = b.first.GetAlignmentInfo();
        AlignmentInfo info_b_2 = b.second.GetAlignmentInfo();

        double score_a = 0.0;
        double score_b = 0.0;

        bool both = a.first.IsSet() && a.second.IsSet() &&
                b.first.IsSet() && b.second.IsSet();

        if (both) {
//            std::cout << "Compare both " << std::endl;
            score_a = Score(info_a_1, info_a_2);
            score_b = Score(info_b_1, info_b_2);
        } else if (a.first.IsSet() && b.first.IsSet()) {
//            std::cout << "Compare firsts" << std::endl;
            score_a = Score(info_a_1);
            score_b = Score(info_b_1);
        } else if (a.second.IsSet() && b.second.IsSet()) {
//            std::cout << "Compare seconds" << std::endl;
            score_a = Score(info_a_2);
            score_b = Score(info_b_2);
        } else {
//            std::cout << "Compare first/second" << std::endl;
            score_a = a.first.IsSet() ? Score(info_a_2) : Score(info_a_1);
            score_b = b.first.IsSet() ? Score(info_b_2) : Score(info_b_1);
        }



//        std::cout << "Comparator: " << score_a << " > " << score_b << " = " << (score_a > score_b) << std::endl;

        return score_a < score_b;
    }

    // SAM QNAME of a single read: its FASTQ id without a trailing "/1" or "/2" mate suffix.
    static std::string ReadQName(std::string const& id) {
        auto n = id.size();
        if (n > 2 && id[n-2] == '/' && (id[n-1] == '1' || id[n-1] == '2')) return id.substr(0, n-2);
        return id;
    }

    // SAM QNAME shared by both mates of a pair. Mate ids either match already (Casava 1.8+, SRA) or
    // differ only in a final 1/2 after a separator ("x/1" "x/2", "x.1" "x.2", "x_1" "x_2"), which is
    // dropped. Otherwise both mates take the first mate's id. The name is never truncated blindly:
    // readers tell reads apart by comparing the names of neighbouring records.
    static std::string PairQName(std::string const& id1, std::string const& id2) {
        if (id1 == id2) return id1;
        auto n = id1.size();
        bool mate_suffix = n > 2 && n == id2.size() &&
                           id1[n-1] == '1' && id2[n-1] == '2' &&
                           (id1[n-2] == '/' || id1[n-2] == '.' || id1[n-2] == '_') &&
                           id1.compare(0, n-1, id2, 0, n-1) == 0;
        return mate_suffix ? id1.substr(0, n-2) : id1;
    }

    static void ArtoSAM(SamEntry &sam, AlignmentResult const& ar, AlignmentInfo &info, FastxRecord &record, std::string const& qname) {
        sam.m_qname = qname;
        sam.m_flag = 0;
        sam.m_rname = std::to_string(ar.Taxid()) + "_" + std::to_string(ar.GeneId());
        sam.m_pos = info.gene_alignment_start + 1;
        sam.m_mapq = ar.AlignmentScore();
        sam.m_cigar = info.compressed_cigar;
        sam.m_rnext = "*";
        sam.m_pnext = 0;
        sam.m_tlen = ar.Cigar().length();
        sam.m_seq = ar.Forward() ? record.sequence : KmerUtils::ReverseComplement(record.sequence);
        sam.m_qual = record.quality;
        sam.m_uniques = ar.Uniques();
        sam.m_uniques_two = ar.UniquesTwo();
        if (!ar.Forward()) reverse(sam.m_qual.begin(), sam.m_qual.end());
    }

    using SNPList = std::vector<SNP>;

    class ProtalAlignmentDataOutputHandler {
    private:
        std::ostream& m_os;
        BufferedStringOutput m_output;
        AlignmentInfo m_info;
        double m_min_cigar_ani;

        std::string m_non_existent_alignment = ',' + std::to_string(INT32_MIN) + ',' + std::to_string(INT32_MIN) + ',' + std::to_string(INT32_MIN) + ',' + std::to_string(INT32_MIN) + ',' + std::to_string(INT32_MIN) + ",0";
    public:
        ProtalAlignmentDataOutputHandler(std::ostream& os, size_t buffer_capacity, double min_cigar_ani) :
            m_os(os),
            m_output(buffer_capacity),
            m_min_cigar_ani(min_cigar_ani) {}

        ProtalAlignmentDataOutputHandler(ProtalAlignmentDataOutputHandler const& other) :
                m_os(other.m_os),
                m_output(other.m_output.Capacity()),
                m_min_cigar_ani(other.m_min_cigar_ani) {}

        ~ProtalAlignmentDataOutputHandler() {
#pragma omp critical(adata_output)
            m_output.Write(m_os);
        }

        void operator () (PairedAlignmentResultList& alignment_results, FastxRecord& record1, FastxRecord& record2, size_t read_id=0, bool first_pair=true) {
            if (alignment_results.empty()) return;

            auto &best = alignment_results.front();

            // Mate 1, or mate 2 when only mate 2 aligned (see ProtalPairedOutputHandler).
            if (CigarANI((best.first.IsSet() ? best.first : best.second).Cigar()) < m_min_cigar_ani) {
                return;
            }



            /* ##############################################################
             * Output alignments.
             */

            bool first = true;

            std::string alignment_data_str = "";

//            std::cout << record1.id << std::endl;

            auto [taxonomic_id, gene_id] = KmerUtils::ExtractTaxIdGeneId(record1.id);

            bool global_label = false;
            size_t max_al = 5;
            size_t al_count = 0;

            auto best_taxid = alignment_results.front().first.IsSet() ? alignment_results.front().first.Taxid() : alignment_results.front().second.Taxid();
            bool best_true = best_taxid == taxonomic_id;
            bool ambiguous_best = alignment_results.size() > 1 &&
                    Score(alignment_results[0]) == Score(alignment_results[1]);


//            if (!best_true && !ambiguous_best) {
//                std::cout << record1.id << std::endl;
//                for (auto &[ar1, ar2]: alignment_results) {
//                    if (ar1.IsSet()) std::cout << "1: " << ar1.ToString() << std::endl;
//                    if (ar2.IsSet()) std::cout << "2: " << ar2.ToString() << std::endl;
//                    std::cout << "Double score: " << Score({ ar1, ar2 }) << std::endl;
//                    std::cout << "Bit score:    " << Bitscore({ ar1, ar2 }) << std::endl;
//                    std::cout << "---" << std::endl;
//
//                }
//                Utils::Input();
//            }

            for (auto &[ar1, ar2]: alignment_results) {
//                std::cout << al_count+1 << std::endl;
//                std::cerr << ar1.ToString() << "\n" << ar2.ToString() << std::endl;

                bool both = ar1.IsSet() && ar2.IsSet();

                int alignment_length = 0;
                int alignment_score = 0;
                int alignment_mismatch_quality = 0;
                int adjusted_score = 0;

                int pred = ar1.IsSet() ? ar1.Taxid() : ar2.Taxid();
                bool label = pred == taxonomic_id;
                global_label |= label;

//                if (ar1.IsSet()) {
//                    std::cout << "    " << ar1.Taxid() << " " << ar1.GeneId() << std::endl;
//                } else {
//                    std::cout << "    " << ar2.Taxid() << " " << ar2.GeneId() << std::endl;
//                }

                if (ar1.IsSet()) {
                    auto& info = ar1.GetAlignmentInfo();
                    alignment_length += info.alignment_length;
                    alignment_score += info.Score();
                    alignment_mismatch_quality += 0;

                }
                if (ar2.IsSet()) {
                    auto& info = ar2.GetAlignmentInfo();
                    alignment_length += info.alignment_length;
                    alignment_score += info.Score();
                    alignment_mismatch_quality += 0;
                }
                if (both) {

                }



                if (al_count < max_al) {
                    if (al_count > 0) alignment_data_str += ',';
                    alignment_data_str += std::to_string(pred);
                    alignment_data_str += ',' + std::to_string(alignment_score);
                    alignment_data_str += ',' + std::to_string(alignment_length);
                    alignment_data_str += ',' + std::to_string(alignment_mismatch_quality);
                    alignment_data_str += ',' + std::to_string(both);
                    alignment_data_str += ',' + std::to_string(label);
                }

                al_count++;
            }



            for (auto i = al_count; i < max_al; i++) {
                alignment_data_str += m_non_existent_alignment;
            }

            alignment_data_str += ',' + std::to_string(global_label);
            alignment_data_str += ',' + record1.id + '\n';
//            std::cout << alignment_data_str << std::endl;


            if (!m_output.Write(alignment_data_str)) {
#pragma omp critical(adata_output)
                m_output.Write(m_os);
            }

//            Utils::Input();

            // output stats

        }
     };

    /*
     * Output of single-end reads: a read's candidate alignments as unpaired records (no 0x1), the
     * best first with the read's MAPQ, the others as secondary alignments (0x100) with MAPQ 0.
     */
    template<bool DEBUG=false>
    class ProtalSingleOutputHandler {
    private:
        std::ostream& m_sam_os;
        BufferedStringOutput m_sam_output;
        SamEntry m_sam;

        GenomeLoader& m_genomes;

        size_t m_max_out = 1;
        double m_min_cigar_ani = 0;
    public:
        size_t alignments = 0;

        ProtalSingleOutputHandler(std::ostream& sam_os, size_t max_out, size_t varkit_buffer_capacity, size_t sam_buffer_capacity, GenomeLoader& genomes, double min_cigar_ani=0.0f) :
                m_sam_os(sam_os),
                m_sam_output(sam_buffer_capacity),
                m_genomes(genomes),
                m_max_out(max_out),
                m_min_cigar_ani(min_cigar_ani) {}

        ProtalSingleOutputHandler(ProtalSingleOutputHandler const& other) :
                m_sam_os(other.m_sam_os),
                m_sam_output(other.m_sam_output.Capacity()),
                m_genomes(other.m_genomes),
                m_max_out(other.m_max_out),
                m_min_cigar_ani(other.m_min_cigar_ani) {}

        ~ProtalSingleOutputHandler() {
#pragma omp critical(sam_output)
            m_sam_output.Write(m_sam_os);
        }

        // The score candidates are ranked and MAPQ is computed with, as for pairs (Bitscore).
        static int Bitscore(AlignmentResult const& ar) {
            return ar.GetAlignmentInfo().Score(2, 3, 1, 2);
        }

        // A read's candidates best first, each alignment once: anchors that lead to the same
        // alignment would otherwise count as a second, equally good hit and give the read MAPQ 0.
        static void RankCandidates(AlignmentResultList& results) {
            std::stable_sort(results.begin(), results.end(), [](AlignmentResult const& a, AlignmentResult const& b) {
                return Bitscore(a) > Bitscore(b);
            });
            auto same = [](AlignmentResult const& a, AlignmentResult const& b) {
                return a.Taxid() == b.Taxid() && a.GeneId() == b.GeneId() && a.Forward() == b.Forward() &&
                       a.GetAlignmentInfo().gene_alignment_start == b.GetAlignmentInfo().gene_alignment_start;
            };
            AlignmentResultList distinct;
            for (auto& ar : results) {
                if (std::none_of(distinct.begin(), distinct.end(), [&](AlignmentResult const& d) { return same(d, ar); })) {
                    distinct.emplace_back(std::move(ar));
                }
            }
            results = std::move(distinct);
        }

        void operator () (AlignmentResultList& alignment_results, FastxRecord& record) {
            if (alignment_results.empty()) return;
            RankCandidates(alignment_results);

            auto& best = alignment_results.front();
            if (CigarANI(best.Cigar()) < m_min_cigar_ani) {
                return;
            }

            int const best_score = Bitscore(best);
            int const second_score = alignment_results.size() > 1 ? Bitscore(alignment_results[1]) : 0;
            int const mapq = best_score > 0 ? MAPQv2(best_score, second_score) : 0;

            if constexpr (DEBUG) {
                std::vector<std::string> tokens;
                LineSplitter::Split(record.id, "-", tokens);
                auto ref = std::to_string(best.Taxid()) + "_" + std::to_string(best.GeneId());
                bool is_correct = !tokens.empty() && ref == tokens[0];
#pragma omp critical(errout)
                std::cerr << int(is_correct) << '\t' << mapq << '\t' << best_score << '\t' << second_score << '\t' << record.id << std::endl;
            }

            auto const qname = ReadQName(record.id);
            SNPList snps;
            bool first = true;
            size_t output_counter = 0;
            std::string read_records;
            for (auto& ar : alignment_results) {
                ArtoSAM(m_sam, ar, ar.GetAlignmentInfo(), record, qname);
                Flag::SetReadReverseComplement(m_sam.m_flag, !ar.Forward());
                Flag::SetNotPrimaryAlignment(m_sam.m_flag, !first);
                m_sam.m_mapq = first ? mapq : 0;
                m_sam.m_tlen = 0;

                auto const& reference = m_genomes.GetGenome(ar.Taxid()).GetGene(ar.GeneId()).Sequence();
                if (!ExtractSNPs(m_sam, reference, snps, ar.Taxid(), ar.GeneId(), 0)) {
#pragma omp critical(err_out)
                    {
                        std::cerr << record.to_string() << std::endl;
                        std::cerr << m_sam.ToString() << std::endl;
                        PrintAlignment(m_sam, reference, std::cerr);
                    }
                    // Skip only this inconsistent candidate; the read's other alignments are still written.
                    continue;
                }

                first = false;
                alignments++;
                if (!read_records.empty()) read_records += '\n';
                read_records += m_sam.ToString();
                if (++output_counter == m_max_out) break;
            }

            // All candidates of the read go into the buffer as one unit, so a flush cannot let another
            // thread's records land between them: readers take adjacent records with one name as one
            // read's candidates.
            if (!read_records.empty() && !m_sam_output.Write(std::move(read_records))) {
#pragma omp critical(sam_output)
                m_sam_output.Write(m_sam_os);
            }
        }
    };


    /*
     * Protal Output Handler
     */
    template<bool DEBUG=false>
    class ProtalPairedOutputHandler {
    private:
        std::ostream& m_sam_os;
//        ogzstream& ogz;
        BufferedStringOutput m_sam_output;
        SamEntry m_sam1;
        SamEntry m_sam2;

        AlignmentInfo m_info;

        GenomeLoader& m_genomes;

        size_t m_max_out = 1;
        double m_min_cigar_ani = 0;
    public:
        size_t alignments = 0;

        ProtalPairedOutputHandler(std::ostream& sam_os, size_t max_out, size_t varkit_buffer_capacity, size_t sam_buffer_capacity, GenomeLoader& genomes, double min_cigar_ani=0.0f) :
                m_sam_os(sam_os),
                m_sam_output(sam_buffer_capacity),
                m_min_cigar_ani(min_cigar_ani),
                m_max_out(max_out),
                m_genomes(genomes) {}

        ProtalPairedOutputHandler(ProtalPairedOutputHandler const& other) :
                m_sam_os(other.m_sam_os),
                m_sam_output(other.m_sam_output.Capacity()),
                m_min_cigar_ani(other.m_min_cigar_ani),
                m_max_out(other.m_max_out),
                m_genomes(other.m_genomes) {}

        ~ProtalPairedOutputHandler() {
#pragma omp critical(sam_output)
            m_sam_output.Write(m_sam_os);
        }


        // A mate's best alignment among a fragment's candidates: the index of a candidate holding it,
        // and the mate's MAPQ against its own alternatives. An alignment that appears in several
        // candidates (paired with different alignments of the other mate) counts once.
        struct MateBest {
            size_t index = SIZE_MAX;
            int mapq = 0;
        };

        static MateBest BestOfMate(PairedAlignmentResultList const& results, bool mate1) {
            std::vector<std::pair<int, size_t>> scored;  // (score, candidate index), one per distinct alignment
            std::vector<std::tuple<uint32_t, uint32_t, int, bool>> seen;
            for (size_t i = 0; i < results.size(); i++) {
                auto const& ar = mate1 ? results[i].first : results[i].second;
                if (!ar.IsSet()) continue;
                auto key = std::make_tuple(ar.Taxid(), ar.GeneId(), ar.GetAlignmentInfo().gene_alignment_start, ar.Forward());
                if (std::find(seen.begin(), seen.end(), key) != seen.end()) continue;
                seen.push_back(key);
                scored.emplace_back(ar.GetAlignmentInfo().Score(2, 3, 1, 2), i);
            }
            if (scored.empty()) return {};
            std::stable_sort(scored.begin(), scored.end(), [](auto const& a, auto const& b) { return a.first > b.first; });
            int const best = scored[0].first;
            int const second = scored.size() > 1 ? scored[1].first : 0;
            return { scored[0].second, best > 0 ? MAPQv2(best, second) : 0 };
        }

        // Writes both mates of a fragment whose mates align only separately: to two genes (e.g.
        // neighbouring genes of an operon; candidates pair mates only within one gene) or to one gene
        // in the wrong orientation. Each mate is written at its own best alignment with its own MAPQ,
        // flagged paired but not properly paired. Ranking the two single-mate candidates against each
        // other instead gave MAPQ ~0 and lost the fragment. Returns false, writing nothing, when only
        // one mate aligned or an alignment is unusable; the caller then writes the candidates as usual.
        // Only these two primary records are written, also with -m > 1.
        bool WriteSplitMates(PairedAlignmentResultList& results, FastxRecord& record1, FastxRecord& record2, std::string const& qname) {
            auto const best1 = BestOfMate(results, true);
            auto const best2 = BestOfMate(results, false);
            if (best1.index == SIZE_MAX || best2.index == SIZE_MAX) return false;
            auto& ar1 = results[best1.index].first;
            auto& ar2 = results[best2.index].second;
            if (CigarANI(ar1.Cigar()) < m_min_cigar_ani || CigarANI(ar2.Cigar()) < m_min_cigar_ani) return false;

            ArtoSAM(m_sam1, ar1, ar1.GetAlignmentInfo(), record1, qname);
            ArtoSAM(m_sam2, ar2, ar2.GetAlignmentInfo(), record2, qname);
            SNPList snps;
            if (!ExtractSNPs(m_sam1, m_genomes.GetGenome(ar1.Taxid()).GetGene(ar1.GeneId()).Sequence(), snps, ar1.Taxid(), ar1.GeneId(), 0) ||
                !ExtractSNPs(m_sam2, m_genomes.GetGenome(ar2.Taxid()).GetGene(ar2.GeneId()).Sequence(), snps, ar2.Taxid(), ar2.GeneId(), 0)) {
                return false;
            }
            for (auto [sam, ar, other, mapq, is_read1] : { std::make_tuple(&m_sam1, &ar1, &ar2, best1.mapq, true),
                                                            std::make_tuple(&m_sam2, &ar2, &ar1, best2.mapq, false) }) {
                Flag::SetPairedEnd(sam->m_flag, true, false, is_read1, !is_read1);
                Flag::SetReadReverseComplement(sam->m_flag, !ar->Forward());
                Flag::SetMateReverseComplement(sam->m_flag, !other->Forward());
                sam->m_mapq = mapq;
                sam->m_tlen = 0;  // the mates do not span one template on one reference
            }
            m_sam1.m_rnext = m_sam1.m_rname == m_sam2.m_rname ? "=" : m_sam2.m_rname;
            m_sam2.m_rnext = m_sam2.m_rname == m_sam1.m_rname ? "=" : m_sam1.m_rname;
            m_sam1.m_pnext = m_sam2.m_pos;
            m_sam2.m_pnext = m_sam1.m_pos;

            alignments++;
            std::string records = m_sam1.ToString() + '\n' + m_sam2.ToString();
            if (!m_sam_output.Write(std::move(records))) {
#pragma omp critical(sam_output)
                m_sam_output.Write(m_sam_os);
            }
            return true;
        }

        void operator () (PairedAlignmentResultList& alignment_results, FastxRecord& record1, FastxRecord& record2, size_t read_id=0, bool first_pair=true) {
            if (alignment_results.empty()) return;

            auto& best = alignment_results.front();

            // Identity of the mate the best candidate starts from: mate 1, or mate 2 when mate 1 did
            // not align. (Testing mate 1 alone dropped every fragment whose best alignment is mate 2's:
            // an empty CIGAR has identity 0.)
            auto const& anchor = best.first.IsSet() ? best.first : best.second;
            if (CigarANI(anchor.Cigar()) < m_min_cigar_ani) {
                return;
            }

            auto const qname = PairQName(record1.id, record2.id);
            if (!(best.first.IsSet() && best.second.IsSet()) && WriteSplitMates(alignment_results, record1, record2, qname)) {
                return;
            }

            /* ##############################################################
             * Output alignments.
             */

            bool first = true;

            std::string alignment_data_str = "";


            int mapq = MAPQv1(alignment_results);

            auto& any = best.first.IsSet() ? best.first : best.second;

            constexpr bool debug = true;

            if constexpr (DEBUG) {
                auto [mapqs, s1, s2] = MAPQv1Debug(alignment_results);
                std::vector<std::string> tokens;
                LineSplitter::Split(record1.id, "-", tokens);
                auto true_ref = tokens[0];

                auto ref = std::to_string(any.Taxid()) + "_" + std::to_string(any.GeneId());
                bool is_correct = ref == true_ref;
#pragma omp critical(errout)
                std::cerr << int(is_correct) << '\t' << mapqs << '\t' << s1 << '\t' << s2 << '\t' << record1.id << std::endl;
            }


            SNPList snps;

            size_t output_counter = 0;
            std::string read_records;
            for (auto& [ar1, ar2] :  alignment_results) {
                bool both = ar1.IsSet() && ar2.IsSet();
                // A mate without an alignment is valid; only an alignment inconsistent with its CIGAR is not.
                bool valid1 = true;
                bool valid2 = true;

                int alignment_length = 0;
                int alignment_score = 0;
                int adjusted_score = 0;

                if (ar1.IsSet()) {
                    auto& info = ar1.GetAlignmentInfo();
                    ArtoSAM(m_sam1, ar1, info, record1, qname);
                    Flag::SetPairedEnd(m_sam1.m_flag, true, both, true);
                    Flag::SetReadUnmapped(m_sam1.m_flag, false);
                    Flag::SetReadReverseComplement(m_sam1.m_flag, !ar1.Forward());
                    Flag::SetNotPrimaryAlignment(m_sam1.m_flag, !first);


                    alignment_length += info.alignment_length;
                    alignment_score += info.Score();
                    m_sam1.m_mapq = first ? mapq : 0;

                    valid1 = ExtractSNPs(m_sam1, m_genomes.GetGenome(ar1.Taxid()).GetGene(ar1.GeneId()).Sequence(), snps, ar1.Taxid(), ar1.GeneId(), 0);
                }
                if (ar2.IsSet()) {
                    auto len = std::count_if(ar2.Cigar().begin(), ar2.Cigar().end(), [](char c) {
                        return c != 'D';
                    });

                    if (len != record2.sequence.length()) {
#pragma omp critical(debug_out)
                        {
                            std::cout << " --------- Problemo --------- " << std::endl;
                            std::cout << "Len: " << len << std::endl;
                            std::cout << record2.sequence << std::endl;
                            std::cout << ar2.Cigar() << std::endl;

                        }
//                        exit(11);
                    }

                    auto& info = ar2.GetAlignmentInfo();
                    ArtoSAM(m_sam2, ar2, info, record2, qname);
                    Flag::SetPairedEnd(m_sam2.m_flag, true, both, false, true);
                    Flag::SetReadUnmapped(m_sam2.m_flag, false);
                    Flag::SetReadReverseComplement(m_sam2.m_flag, !ar2.Forward());
                    Flag::SetNotPrimaryAlignment(m_sam2.m_flag, !first);

                    alignment_length += info.alignment_length;
                    alignment_score += info.alignment_score;
                    m_sam2.m_mapq = first ? mapq : 0;

                    valid2 = ExtractSNPs(m_sam2, m_genomes.GetGenome(ar2.Taxid()).GetGene(ar2.GeneId()).Sequence(), snps, ar2.Taxid(), ar2.GeneId(), 0);
                }
                if (both) {
                    m_sam1.m_rnext = "=";
                    m_sam2.m_rnext = "=";
                    m_sam1.m_pnext = m_sam2.m_pos;
                    m_sam2.m_pnext = m_sam1.m_pos;
                    Flag::SetMateReverseComplement(m_sam1.m_flag, (FLAG_t)!ar2.Forward());
                    Flag::SetMateReverseComplement(m_sam2.m_flag, (FLAG_t)!ar1.Forward());
                    Flag::SetMateUnmapped(m_sam1.m_flag, false);
                    Flag::SetMateUnmapped(m_sam2.m_flag, false);
                    Flag::SetPairBothAlign(m_sam1.m_flag, true);
                    Flag::SetPairBothAlign(m_sam2.m_flag, true);
                } else {
                    // Only one mate aligned: flag the other as unmapped so a reader does not expect,
                    // and swallow, a mate record on the next line.
                    Flag::SetMateUnmapped(ar1.IsSet() ? m_sam1.m_flag : m_sam2.m_flag, true);
                }

                if (!valid1 || !valid2) {
                    if (!valid1) {
#pragma omp critical(err_out)
                        {
                            std::cerr << record1.to_string() << std::endl;
                            std::cerr << m_sam1.ToString() << std::endl;
                            std::string const& reference = m_genomes.GetGenome(ar1.Taxid()).GetGene(ar1.GeneId()).Sequence();
                            PrintAlignment(m_sam1, reference, std::cerr);
                        }
                    }
                    if (!valid2) {
#pragma omp critical(err_out)
                        {
                            std::cerr << record2.to_string() << std::endl;
                            std::cerr << m_sam2.ToString() << std::endl;
                            auto& reference = m_genomes.GetGenome(ar2.Taxid()).GetGene(ar2.GeneId()).Sequence();
                            PrintAlignment(m_sam2, reference, std::cerr);
                        }
                    }
                    // Skip only this inconsistent candidate; the read's other alignments are still written.
                    continue;
                }

                first = false;
                alignments++;
                if (!read_records.empty()) read_records += '\n';
                if (ar1.IsSet()) read_records += m_sam1.ToString();
                if (both) read_records += '\n';
                if (ar2.IsSet()) read_records += m_sam2.ToString();
                if (++output_counter == m_max_out) {
                    break;
                }
            }

            // All candidates of the read go into the buffer as one unit, so a flush cannot let another
            // thread's records land between them: readers pair a read1 record with the line that
            // follows, and take adjacent records with one name as one read's candidates.
            if (!read_records.empty() && !m_sam_output.Write(std::move(read_records))) {
#pragma omp critical(sam_output)
                m_sam_output.Write(m_sam_os);
            }
        }
    };
}