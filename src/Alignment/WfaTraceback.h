// WfaTraceback.h - how WFA2 recovers an alignment's operations (--wfa_traceback). Classic is WFA2-lib's high-memory
// mode: it keeps the M, I and D wavefronts of every score and follows them back. Singletrack (López-Villellas et al.,
// bioRxiv 10.1101/2025.10.31.685625; WFA2-lib's wavefront_memory_singletrack) keeps only the M wavefronts and
// recomputes a gap's I or D steps while tracing back, which leaves less memory to touch and was 1.2-1.4x faster in the
// paper at protal's penalties. Both give the optimal score, but among equally good alignments they may choose
// differently, so SAM records, and the features a model was trained on, can differ. A database records the traceback
// it was built with (build_settings.tsv, BuildSettings.h), and a run follows it unless --wfa_traceback says otherwise.
// --build aligns with Singletrack unless told otherwise (since 0.8.2; classic before).
#pragma once

#include <optional>
#include <string>
#include <string_view>

namespace protal {
    enum class WfaTraceback { Classic, Singletrack };

    // The name of `traceback` in --wfa_traceback and build_settings.tsv.
    inline std::string_view WfaTracebackName(WfaTraceback traceback) {
        return traceback == WfaTraceback::Singletrack ? "singletrack" : "classic";
    }

    // The traceback named `name` (classic, singletrack), or nothing if it names none.
    inline std::optional<WfaTraceback> ParseWfaTraceback(std::string_view name) {
        if (name == "classic") return WfaTraceback::Classic;
        if (name == "singletrack") return WfaTraceback::Singletrack;
        return std::nullopt;
    }

    // A run's traceback: the one given (--wfa_traceback other than auto), else the database's. A database without
    // build_settings.tsv was built before the option existed, with the classic traceback.
    inline WfaTraceback ResolveWfaTraceback(std::optional<WfaTraceback> given, std::optional<WfaTraceback> database) {
        return given.value_or(database.value_or(WfaTraceback::Classic));
    }
}
