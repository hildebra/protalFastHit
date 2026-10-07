//
// Created by fritsche on 21/02/23.
//

#ifndef PROTAL_SEQUENCERANGEHANDLER_H
#define PROTAL_SEQUENCERANGEHANDLER_H

#include "SequenceRange.h"
#include <iostream>
#include <string>

using SequenceRangeList = std::vector<SequenceRange>;
class SequenceRangeHandler {
private:
    using SRIterator = const std::vector<SequenceRange>::iterator;

    // Variables
    mutable SequenceRangeList m_ranges;
    mutable CoverageVec m_cov;

    SRIterator GetOrAddSequenceRangeIterator(const size_t start, const size_t end);

public:
    // Iterator wrapper
    SRIterator begin() noexcept;
    SRIterator end() noexcept;

    // Navigation
    SRIterator FindSequenceRange(size_t start, size_t end);
    const SequenceRangeList& GetRanges() const;
    SequenceRangeList& GetRanges();

    SequenceRangeHandler Intersect(SequenceRangeHandler const& handler) const;
    SequenceRangeHandler Intersect(SequenceRangeHandler &handler, size_t min_coverage, size_t min_sequence_length);
    SequenceRangeHandler Unify(SequenceRangeHandler const& handler) const;

    const SequenceRange& GetRange(size_t pos) const;
    SequenceRange& GetRange(size_t pos);
    bool HasRange(size_t pos) const;

    size_t Size() const;
    // Drops all ranges and frees their memory.
    void Clear();
    size_t SequenceLength() const;
    // Adds a range, merged with every range it overlaps or touches, so that ranges never overlap.
    void Merge(SequenceRange range);
    void Add(size_t start, size_t end);
    void Add(SequenceRange& range);
    void Add(SequenceRange&& range);

    std::string ToString() const;
    std::string ToVerboseString() const;

    // Reads per position: all, or one strand's (SequenceRange::kForward, kReverse); with
    // max_divergence, only reads of at most that divergence from the gene.
    CoverageVec CalculateCoverageVector2(int strand = SequenceRange::kBothStrands, uint8_t max_divergence = 255) const;
    bool AreRangesValid(size_t const& reference_length) const;
    CoverageVec& GetCoverageVector() {
        return m_cov;
    }

    SequenceRangeHandler FromCoverageVector(CoverageVec const& coverage_vector, size_t min_sequence_length) const;

    static CoverageVec IntersectCoverageVectors(CoverageVec const& cov1, CoverageVec const& cov2, size_t min_coverage) {
        CoverageVec intersection(std::min(cov1.size(), cov2.size()), 0 );
        for (auto i = 0; i < std::min(cov1.size(), cov2.size()); i++) {
            intersection[i] = cov1[i] >= min_coverage && cov2[i] >= min_coverage;
        }
        return intersection;
    }

    static std::string CoverageVectorToString(CoverageVec const& cov, size_t max = 120) {
        std::string str;
        for (auto i = 0; i < cov.size(); i++) {
            str += std::to_string(i) + " " + std::string(std::min(cov[i], static_cast<uint32_t>(max)), '#') + (cov[i] >= max ? " " + std::to_string(cov[i]) : "") + '\n';
        }
        return str;
    }

    void SetCoverageVector(CoverageVec const& vector);
    size_t CoveredPortion(uint32_t min_cov=1) const;
};


#endif //PROTAL_SEQUENCERANGEHANDLER_H
