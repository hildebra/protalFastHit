"""Adds PROTAL_LR_DEBUG diagnostics to a COPY of the protal source (never the repository):
candidates, why Align() drops one, the segments' hits, and the problem of an invalid alignment."""
import sys

root = sys.argv[1]


def patch(path, pairs):
    with open(path) as fh:
        text = fh.read()
    for old, new in pairs:
        if old not in text:
            sys.exit(f"not found in {path}: {old[:80]!r}")
        text = text.replace(old, new, 1)
    with open(path, "w") as fh:
        fh.write(text)


GUARD = ("#include <cstdlib>\n#include <iostream>\n#ifndef PROTAL_LRDEBUG_DEFINED\n#define PROTAL_LRDEBUG_DEFINED\n"
         "namespace protal { inline bool LrDebug() { static bool d = std::getenv(\"PROTAL_LR_DEBUG\") != nullptr; return d; } }\n"
         "#endif\n")

lr = root + "/src/Core/LongReads.h"
patch(lr, [
    ("#include <vector>\n", "#include <vector>\n" + GUARD),
    ("            if (gene.TwiceCentre() < 2 * chunk.core.start || gene.TwiceCentre() >= 2 * chunk.core.end) return;\n",
     "            if (gene.TwiceCentre() < 2 * chunk.core.start || gene.TwiceCentre() >= 2 * chunk.core.end) {\n"
     "                if (LrDebug()) std::cerr << \"LRDBG drop-not-owned \" << anchor.taxid << '_' << anchor.geneid << \" fwd=\" << anchor.forward << \" len=\" << anchor.total_length << \" gene=[\" << gene.start << ',' << gene.end << \") chunk=\" << chunk.offset << '\\n';\n"
     "                return;\n            }\n"),
    ("            m_candidates.push_back({ anchor, chunk, gene, window });\n",
     "            m_candidates.push_back({ anchor, chunk, gene, window });\n"
     "            if (LrDebug()) std::cerr << \"LRDBG cand \" << anchor.taxid << '_' << anchor.geneid << \" fwd=\" << anchor.forward << \" len=\" << anchor.total_length << \" links=\" << anchor.chain.size() << \" gene=[\" << gene.start << ',' << gene.end << \") window=[\" << window.start << ',' << window.end << \") chain=[\" << chain.start << ',' << chain.end << \") chunk=\" << chunk.offset << '+' << chunk.length << '\\n';\n"),
    ("                if (position < 0 || position + link.length > static_cast<int64_t>(window_length)) return std::nullopt;\n",
     "                if (position < 0 || position + link.length > static_cast<int64_t>(window_length)) {\n"
     "                    if (LrDebug()) std::cerr << \"LRDBG align-fail link-outside \" << anchor.taxid << '_' << anchor.geneid << \" pos=\" << position << \" wl=\" << window_length << '\\n';\n"
     "                    return std::nullopt;\n                }\n"),
    ("            if (!m_alignment_handler.AlignAnchor(anchor, result, m_window, m_window_rev, false, m_id)) return std::nullopt;\n"
     "            auto& info = result.GetAlignmentInfo();\n"
     "            if (info.GetProxyANI() < m_min_ani) return std::nullopt;\n",
     "            if (!m_alignment_handler.AlignAnchor(anchor, result, m_window, m_window_rev, false, m_id)) {\n"
     "                if (LrDebug()) std::cerr << \"LRDBG align-fail AlignAnchor \" << anchor.taxid << '_' << anchor.geneid << \" fwd=\" << anchor.forward << \" wl=\" << window_length << \" front=\" << anchor.chain.front().readpos << '/' << anchor.chain.front().genepos << \" back=\" << anchor.chain.back().readpos << '/' << anchor.chain.back().genepos << '+' << anchor.chain.back().length << '\\n';\n"
     "                return std::nullopt;\n            }\n"
     "            auto& info = result.GetAlignmentInfo();\n"
     "            if (info.GetProxyANI() < m_min_ani) {\n"
     "                if (LrDebug()) std::cerr << \"LRDBG align-fail proxy-ani \" << anchor.taxid << '_' << anchor.geneid << \" ani=\" << info.GetProxyANI() << \" cigar=\" << info.compressed_cigar << '\\n';\n"
     "                return std::nullopt;\n            }\n"),
    ("            if (left + right >= window_length) return std::nullopt;\n",
     "            if (left + right >= window_length) {\n"
     "                if (LrDebug()) std::cerr << \"LRDBG align-fail all-clipped \" << anchor.taxid << '_' << anchor.geneid << '\\n';\n"
     "                return std::nullopt;\n            }\n"),
    ("                if (segment.hits.empty()) continue;\n                RankLongReadSegment(segment);\n",
     "                if (segment.hits.empty()) continue;\n                RankLongReadSegment(segment);\n"
     "                if (LrDebug()) { std::cerr << \"LRDBG segment mapq=\" << segment.mapq; for (auto const& h : segment.hits) std::cerr << ' ' << h.alignment.Taxid() << '_' << h.alignment.GeneId() << ':' << LongReadBitscore(h) << '@' << h.alignment.GetAlignmentInfo().gene_alignment_start; std::cerr << '\\n'; }\n"),
    ("            m_candidates.clear();\n            m_last_seeds = 0;\n",
     "            m_candidates.clear();\n            if (LrDebug()) std::cerr << \"LRDBG read \" << record.id << \" len=\" << record.sequence.size() << '\\n';\n            m_last_seeds = 0;\n"),
])

al = root + "/src/Core/AlignmentStrategy.h"
patch(al, [
    ("#include \"AnchoredAlignment.h\"\n", "#include \"AnchoredAlignment.h\"\n" + GUARD),
    ("                if (!valid) {\n#pragma omp critical (invalid_align)\n",
     "                if (!valid) {\n"
     "                    if (LrDebug()) {\n"
     "                        IsAlignmentValid(info, read, gene.Sequence(), 0, false);  // prints the problem (its own critical section)\n"
     "                        std::cerr << \"CIGAR: \" << info.compressed_cigar << \" ops=\" << m_ops.size() << \" abs_pos=\" << abs_pos << \" ref=[\" << window.ref_start << ',' << window.ref_end << \") free q=\" << window.read_begin_free << ',' << window.read_end_free << \" r=\" << window.ref_begin_free << ',' << window.ref_end_free << \" read_len=\" << read.length() << \" gene_len=\" << geneseq.length() << std::endl;\n"
     "                    }\n"
     "#pragma omp critical (invalid_align)\n"),
])
print("patched", lr, al)
