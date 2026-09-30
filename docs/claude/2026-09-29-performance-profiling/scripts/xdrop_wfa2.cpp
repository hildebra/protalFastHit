// Does WFA2's X-drop work with protal's scoring (match 0, mismatch 4, gap 6+2, ends-free)?
// Modes: wfadaptive (the library default protal runs with), xdrop alone, both.
#include <cstdio>
#include <random>
#include <string>
#include "wfa2-lib/bindings/cpp/WFAligner.hpp"
using namespace wfa;

static std::string Random(size_t n, std::mt19937& rng) {
    std::string s(n, 'A');
    for (auto& c : s) c = "ACGT"[rng() % 4];
    return s;
}
static std::string Mutate(std::string s, double rate, std::mt19937& rng) {
    std::uniform_real_distribution<double> u(0, 1);
    std::string out;
    for (size_t i = 0; i < s.size(); i++) {
        double r = u(rng);
        if (r < rate * 0.8) out += s[i] == 'A' ? 'C' : 'A';         // substitution
        else if (r < rate * 0.9) {}                                  // deletion
        else if (r < rate) { out += s[i]; out += "ACGT"[rng() % 4]; } // insertion
        else out += s[i];
    }
    return out;
}

int main() {
    std::mt19937 rng(7);
    struct Case { const char* name; size_t len; double rate; };
    for (Case c : { Case{"150 bp 3%", 150, 0.03}, Case{"150 bp 12%", 150, 0.12}, Case{"150 random", 150, 0.75},
                    Case{"2 kb 2%", 2000, 0.02}, Case{"2 kb 8%", 2000, 0.08}, Case{"8 kb 1%", 8000, 0.01} }) {
        std::string ref = Random(c.len + 18, rng);
        std::string query = Mutate(ref.substr(9, c.len), c.rate, rng);
        for (const char* mode : { "wfadaptive", "xdrop", "both" }) {
            for (int xdrop : { 1000, 200, 50, 20 }) {
                if (std::string(mode) == "wfadaptive" && xdrop != 1000) continue;
                WFAlignerGapAffine aligner(4, 6, 2, WFAligner::Alignment, WFAligner::MemoryHigh);
                if (std::string(mode) == "xdrop") { aligner.setHeuristicNone(); aligner.setHeuristicXDrop(xdrop, 1); }
                if (std::string(mode) == "both") aligner.setHeuristicXDrop(xdrop, 1);
                aligner.setMaxAlignmentSteps(INT32_MAX);
                auto status = aligner.alignEndsFree(ref, 18, 18, query, 0, 0);
                std::string cigar = aligner.getAlignment();
                size_t m = 0; for (char ch : cigar) m += ch == 'M';
                std::printf("%-11s %-10s xdrop=%-5d status=%d score=%-6d matches=%zu/%zu\n", c.name, mode,
                            std::string(mode) == "wfadaptive" ? 0 : xdrop, (int)status, aligner.getAlignmentScore(), m, query.size());
            }
        }
    }
}
