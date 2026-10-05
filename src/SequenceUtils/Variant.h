//
// Created by fritsche on 21/02/23.
//

#pragma once
#include "Constants.h"
#include <memory>
#include <string>
#include <vector>
#include <numeric>


enum VariantType {
    INS, DEL, SNP
};

using VariantID = size_t;
using Orientation = bool;
using Base = char;
using FirstRead = bool;
using ReadId = size_t;
using SSize = uint16_t;
using Qual = uint8_t;
using QualList = std::vector<Qual>;
using VariantPos = uint32_t;

// A read's divergence from its gene, 1 - identity, in steps of 0.5% (0-127; 127 = 63.5% or more).
inline uint8_t DivergenceBin(double identity) {
    double const bin = (1.0 - identity) * 200.0 + 0.5;
    return bin <= 0 ? 0 : bin >= 127 ? 127 : static_cast<uint8_t>(bin);
}

// The largest DivergenceBin of reads with at least `min_identity` (all reads for 0 or less).
inline uint8_t MaxDivergenceBin(double min_identity) {
    if (min_identity <= 0) return 127;
    double const bin = (1.0 - min_identity) * 200.0 + 1e-9;
    return bin >= 127 ? 127 : static_cast<uint8_t>(bin);
}

class Variant {

    VariantID variant_id = 0;   // size_t
    VariantType variant_type;   // enum
    VariantPos position = UINT32_MAX;        // uint32_t
    Base reference = 'X';             // char
    Base variant = 'X';               // char
    SSize structural_size = 1;  // uint16_t
    uint32_t observations_fwd = 0;
    uint32_t observations_rev = 0;
    bool is_valid = true;
    bool is_major = false;
    // Inserted/deleted bases of an INDEL; null for SNPs. Owned, and deep-copied with the Variant.
    std::unique_ptr<std::string> structural;
    QualList quals;
    // Per observation (as quals): its read's strand (bit 7) and divergence from the gene (bits 0-6,
    // DivergenceBin), so that the strain MSA can count a taxon's own reads only.
    std::vector<uint8_t> reads;

public:
    Variant(VariantPos pos, Base snp, Base ref) :
            variant_type(VariantType::SNP), position(pos), reference(ref), variant(snp) {};

    Variant(VariantType type, VariantPos pos, Base ref, std::string& structural) :
            variant_type(type), position(pos), reference(ref), structural_size(structural.length()),
            structural(std::make_unique<std::string>(structural)) {};

    Variant(Variant const& other) :
            variant_id(other.variant_id), variant_type(other.variant_type), position(other.position),
            reference(other.reference), variant(other.variant), structural_size(other.structural_size),
            observations_fwd(other.observations_fwd), observations_rev(other.observations_rev),
            is_valid(other.is_valid), is_major(other.is_major),
            structural(other.structural ? std::make_unique<std::string>(*other.structural) : nullptr),
            quals(other.quals), reads(other.reads) {};

    Variant& operator=(Variant const& other) {
        if (this != &other) *this = Variant(other);
        return *this;
    }

    Variant(Variant&&) noexcept = default;
    Variant& operator=(Variant&&) noexcept = default;

    size_t Observations() const {
        return observations_fwd + observations_rev;
    }

    auto& GetQualList() {
        return quals;
    }
    auto GetQualListCopy() const {
        return quals;
    }

    auto IsUninitialized() const {
        return variant == 'X' && reference == 'X';
    }

    size_t QualitySum() const {
        if (IsUninitialized()) {
            exit(12);
        }
        if (IsReference()) return Observations() * 40;
        return std::accumulate(quals.begin(), quals.end(), size_t{0}, [](size_t acc, const uint16_t q) {
            return acc + q;
        });
    }
    size_t MeanQuality() const {
        auto observations = Observations();
        if (observations == 0) return 0;
        return static_cast<size_t>(static_cast<double>(QualitySum())/observations);
    }


    void SetValid(bool is_valid) {
        this->is_valid = is_valid;
    }

    bool GetValid() const {
        return is_valid;
    }

    void SetMajorAllele(bool is_major) {
        this->is_major = is_major;
    }

    bool IsMajorAllele() const {
        return is_major;
    }

    std::string ToString() const {
        std::string str;
        str += '{';
        str += std::to_string(is_valid) + ' ';
        if (variant_type == VariantType::SNP) {
            str += (IsReference() ? "REF " : "SNP ");
        } else if (variant_type == VariantType::INS) {
            str += "INS ";
        } else if (variant_type == VariantType::DEL) {
            str += "DEL ";
        }

        if (variant_type == VariantType::SNP) {
            str += (IsReference() ? std::string(1, reference) : std::string(1, reference) + "->" + std::string(1, variant)) + " ";
//            str += std::to_string(quality) + '\t';
        } else if (variant_type == VariantType::INS || variant_type == VariantType::DEL) {
            str += *structural + '\t';
            str += std::to_string(structural_size);
        }

        str += std::to_string(Observations()) + " (" + std::to_string(observations_fwd) + '/' + std::to_string(observations_rev) + ") ";
        str += std::to_string(position);
        str += " " + std::to_string(MeanQuality());
        str += '}';

        return str;
    }

    std::string ToMinimalString() const {
        std::string str;

        if (variant_type == VariantType::SNP) {
            str += std::string(1, variant) + "(";
        } else if (variant_type == VariantType::INS || variant_type == VariantType::DEL) {
            str += *structural + "(";
        }

        str += std::to_string(Observations());
        str += ", ";

        if (variant_type == VariantType::SNP) {
            str += (IsReference() ? "REF)" : "SNP)");
        } else if (variant_type == VariantType::INS) {
            str += "INS)";
        } else if (variant_type == VariantType::DEL) {
            str += "DEL)";
        }

        return str;
    }

    // An observation from a read on the forward strand or not, whose alignment has this divergence
    // from the gene (DivergenceBin: 0 = identical).
    void AddObservation(Qual quality, bool from_forward, uint8_t divergence = 0) {
        observations_fwd += from_forward;
        observations_rev += !from_forward;
        quals.emplace_back(quality);
        reads.emplace_back(static_cast<uint8_t>((from_forward ? 0x80 : 0) | (divergence & 0x7f)));
    }

    // A copy with the observations of reads of at most `max_divergence` only (an inferred reference
    // allele, which has no per-read record, is copied as it is).
    Variant WithMaxDivergence(uint8_t max_divergence) const {
        Variant copy(*this);
        if (reads.size() != quals.size() || reads.size() != Observations()) return copy;
        copy.observations_fwd = copy.observations_rev = 0;
        copy.quals.clear();
        copy.reads.clear();
        for (size_t i = 0; i < reads.size(); i++) {
            if ((reads[i] & 0x7f) > max_divergence) continue;
            bool const forward = reads[i] & 0x80;
            copy.AddObservation(quals[i], forward, reads[i] & 0x7f);
        }
        return copy;
    }

    bool IsSNP() const {
        return variant_type == VariantType::SNP;
    }

    char Reference() const {
        return reference;
    }

    VariantPos Position() const {
        return position;
    }

    bool IsReference() const {
        return reference == variant;
    }

    bool HasFwdAndRev() const {
        return observations_fwd > 0 && observations_rev > 0;
    }

    uint32_t ObservationsForward() const {
        return observations_fwd;
    }

    uint32_t ObservationsReverse() const {
        return observations_rev;
    }

    // Sets inferred observation counts (the reference allele: the reads of each strand with a base
    // at the position, less those carrying another allele there).
    void SetObservations(size_t forward, size_t reverse) {
        observations_fwd = static_cast<uint32_t>(forward);
        observations_rev = static_cast<uint32_t>(reverse);
    }

    size_t GetStructuralSize() const {
        return structural_size;
    }

    std::string GetStructural() const {
        return *structural;
    }

    char GetVariant() const {
        return variant;
    }

    bool IsINDEL() const {
        return variant_type == VariantType::INS || variant_type == VariantType::DEL;
    }

    bool IsINS() const {
        return variant_type == VariantType::INS;
    }

    bool IsDEL() const {
        return variant_type == VariantType::DEL;
    }

    // Alleles of one site in an order of what they are: bases (the reference's and SNPs') before deletions before
    // insertions, then by base or by the inserted/deleted bases. Sorts break ties with it, so that a site's alleles do
    // not stay in the order the reads first showed them, which follows the order of the SAM's records.
    bool AlleleBefore(Variant const& other) const {
        auto const kind = [](Variant const& v) { return v.IsINS() ? 2 : v.IsDEL() ? 1 : 0; };
        if (kind(*this) != kind(other)) return kind(*this) < kind(other);
        if (!IsINDEL()) return variant < other.variant;
        return *structural < *other.structural;
    }

    bool Match(VariantPos pos, Base base, Base ref) const {
        return pos == position &&
            base == variant && ref == reference;
    }

    bool Match(VariantType type, VariantPos pos, std::string& structural) const {
        return type == variant_type && pos == position &&
                *this->structural == structural;
    }

    bool Match(Variant const& other) const {
        return variant_type == other.variant_type &&
               position == other.position &&
                ((variant_type == VariantType::SNP && variant == other.variant) ||
                 (IsINDEL() && *structural == *other.structural));
    };

    void SetStructural(std::string&& structural_string) {
        structural_size = structural_string.length();
        structural = std::make_unique<std::string>(std::move(structural_string));
    }
    void SetStructural(std::string& structural_string) {
        structural_size = structural_string.length();
        structural = std::make_unique<std::string>(structural_string);
    }

    Variant() {}
};