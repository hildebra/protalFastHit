//
// Created by fritsche on 14/08/22.
//

#pragma once

#include "FastxReader.h"
#include "omp.h"

namespace protal {

    class SeqReader {
    private:
//        const size_t m_block_size = (10 * 1024 * 1024);
        const size_t m_block_size = (1024);
        BufferedFastxReader m_reader;
        bool m_valid_fragment = false;
        bool m_valid_block = true;

        std::istream& m_is;

    public:
        SeqReader(std::istream& is) :
                m_is(is) {};

        SeqReader(SeqReader const& other) :
                m_is(other.m_is) {};

        inline void LoadBlockOMP(std::istream& is) {
#pragma omp critical(reader)
            {
                m_valid_block = m_reader.LoadBlock(is, m_block_size);
            }
        }

        bool operator() (FastxRecord &record) {
            m_valid_fragment = m_reader.NextSequence(record);
            if (m_valid_fragment) return true;

            LoadBlockOMP(m_is);
            if (!m_valid_block) return false;

            m_valid_fragment = m_reader.NextSequence(record);
            return m_valid_fragment;
        }

        bool Success() {
            return !m_reader.Error();
        }
    };

    // A read without qualities (FASTA) gets `quality` (Phred+33) for each base; 0 leaves it without.
    inline void FillMissingQuality(FastxRecord& record, char quality) {
        if (quality != 0 && record.quality.empty()) record.quality.assign(record.sequence.size(), quality);
    }

    // Single-end reads, shared by threads as SeqReaderPE: each copy takes batches of records from
    // the stream. Reads without qualities get fasta_quality for each base (FillMissingQuality).
    class SeqReaderSE {
    private:
        const size_t m_record_count = 32;
        BufferedFastxReader m_reader;
        bool m_valid_block = true;
        std::istream& m_is;
        bool m_success = true;
        char m_fasta_quality = 0;

        bool Next(FastxRecord& record) {
            if (!m_reader.NextSequence(record)) return false;
            FillMissingQuality(record, m_fasta_quality);
            return true;
        }

    public:
        explicit SeqReaderSE(std::istream& is, char fasta_quality = 0) :
                m_is(is), m_fasta_quality(fasta_quality) {};

        SeqReaderSE(SeqReaderSE const& other) :
                m_is(other.m_is), m_fasta_quality(other.m_fasta_quality) {};

        bool Success() const {
            return m_success;
        }

        void UpdateSuccess(SeqReaderSE const& other) {
            m_success &= other.m_success;
        }

        bool operator() (FastxRecord &record) {
            if (Next(record)) return true;
            if (m_reader.Error()) {
                m_success = false;
                return false;
            }
#pragma omp critical(reader)
            {
                m_valid_block = m_reader.LoadBatch(m_is, m_record_count);
            }
            if (m_valid_block && Next(record)) return true;
            m_success &= !m_reader.Error();
            return false;
        }
    };

    class SeqReaderPE {
        enum SeqReaderPEError {
            NO_ERROR,DIFF_LENGTH_SEQUENCE, DIFFERENT_LENGTH_BLOCK
        };
    private:
//        const size_t m_block_size = (10 * 1024 * 1024);
        const size_t m_block_size = (1024);
        const size_t m_record_count = 32;
        BufferedFastxReader m_reader_1;
        BufferedFastxReader m_reader_2;

        bool m_valid_fragment_1 = false;
        bool m_valid_fragment_2 = false;

        bool m_valid_block_1 = true;
        bool m_valid_block_2 = true;

        std::istream& m_is1;
        std::istream& m_is2;

        SeqReaderPEError m_error_code = NO_ERROR;
        bool m_success = true;
        char m_fasta_quality = 0;  // for reads without qualities, see FillMissingQuality

    public:
        SeqReaderPE(std::istream& is1, std::istream& is2, char fasta_quality = 0) :
                m_is1(is1),
                m_is2(is2),
                m_fasta_quality(fasta_quality) {};

        SeqReaderPE(SeqReaderPE const& other) :
                m_is1(other.m_is1),
                m_is2(other.m_is2),
                m_fasta_quality(other.m_fasta_quality) {};


        auto& GetFirstStream() {
            return m_is1;
        }

        inline void LoadBlockOMP() {
#pragma omp critical(reader)
            {
//                m_valid_block_1 = m_reader_1.LoadBlock(m_is1, m_block_size);
//                m_valid_block_2 = m_reader_2.LoadBlock(m_is2, m_block_size);
                m_valid_block_1 = m_reader_1.LoadBatch(m_is1, m_record_count);
                m_valid_block_2 = m_reader_2.LoadBatch(m_is2, m_record_count);
            }
        }

        inline void LoadBlockOMP(bool first) {
#pragma omp critical(reader)
            {
                if (first)
                    m_valid_block_1 = m_reader_1.LoadBlock(m_is1, m_block_size);
                else
                    m_valid_block_2 = m_reader_2.LoadBlock(m_is2, m_block_size);
            }
        }

        bool Success() const {
            return m_success;
        }

        void UpdateSuccess(SeqReaderPE const& other) {
            m_success &= other.m_success;
        }

        bool operator() (FastxRecord &record1, FastxRecord &record2) {
            if (!Next(record1, record2)) return false;
            FillMissingQuality(record1, m_fasta_quality);
            FillMissingQuality(record2, m_fasta_quality);
            return true;
        }

        bool Next(FastxRecord &record1, FastxRecord &record2) {
            m_valid_fragment_1 = m_reader_1.NextSequence(record1);
            m_valid_fragment_2 = m_reader_2.NextSequence(record2);

            if (m_reader_1.Error() || m_reader_2.Error()) {
                m_success = false;
                return false;
            }

            if (m_valid_fragment_1 && m_valid_fragment_2) {
                return true;
            }

            // Paired end states must always be identical.
            if (m_valid_fragment_1 != m_valid_fragment_2) {
                std::cerr << omp_get_thread_num() << "Error with next sequence. Paired end file streams are not of the same length. Abort. (R1: " << m_valid_fragment_1 << ", R2: " << m_valid_fragment_2 << ")" << std::endl;
                m_error_code = DIFF_LENGTH_SEQUENCE;
                m_success = false;
                return false;
            }

            LoadBlockOMP();
            if (!m_valid_block_1 && !m_valid_block_2) {
                return false;
            }

            // Paired end states must always be identical.
            if (m_valid_block_1 != m_valid_block_2) {
                std::cerr << omp_get_thread_num() << " Error after load block. Paired end file streams are not of the same length. Abort." << std::endl;
                m_error_code = DIFFERENT_LENGTH_BLOCK;
                m_success = false;
                return false;
            }

            m_valid_fragment_1 = m_reader_1.NextSequence(record1);
            m_valid_fragment_2 = m_reader_2.NextSequence(record2);

            if (m_reader_1.Error() || m_reader_2.Error()) {
                m_success = false;
                return false;
            }

            // Paired end states must always be identical.
            if (m_valid_fragment_1 != m_valid_fragment_2) {
                std::cerr << "2 Error with next sequence. Paired end file streams are not of the same length. Abort." << std::endl;
                m_error_code = DIFF_LENGTH_SEQUENCE;
                m_success = false;
                return false;
            }

            return m_valid_fragment_1 && m_valid_fragment_2;
        }
    };

}
