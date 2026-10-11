// repro3.cpp: the classic traceback (WFA2's high-memory mode, wf-adaptive on as protal has it) of an ends-free alignment
// with the read's end free; whether its operations are an alignment of the two sequences. Built against either tree.
#include <climits>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <string>
#include "wfa2-lib/bindings/cpp/WFAligner.hpp"

int main(int argc, char** argv) {
    std::string const ref = argv[1], read = argv[2];
    for (bool adaptive : { true, false }) {
        wfa::WFAlignerGapAffine aligner(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
        if (!adaptive) aligner.setHeuristicNone();
        int const status = aligner.alignEndsFree(ref, 0, 0, read, 0, static_cast<int>(read.size()));
        std::string const ops = aligner.getAlignment();
        size_t r = 0, q = 0, bad = SIZE_MAX;
        for (size_t i = 0; i < ops.size(); i++) {
            char const op = ops[i];
            if (op == 'M' || op == 'X') {
                if (bad == SIZE_MAX && (r >= ref.size() || q >= read.size() || (ref[r] == read[q]) != (op == 'M'))) bad = i;
                r++, q++;
            } else if (op == 'D') r++;
            else if (op == 'I') q++;
        }
        std::printf("wf-adaptive %s: status %d, score %d, %s (consumed %zu/%zu, %zu/%zu)\n", adaptive ? "on" : "off", status,
                    -aligner.getAlignmentScore(), bad == SIZE_MAX ? "valid" : ("INVALID at op " + std::to_string(bad)).c_str(),
                    r, ref.size(), q, read.size());
    }
}
