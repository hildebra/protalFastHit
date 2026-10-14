#pragma once

// protal's run. Its stages are translation units of their own in src/Run/ (RunStages.h lists them), each including only
// the headers it uses: an edit recompiles the stages it touches, and the stages compile in parallel
// (docs/claude/2026-10-14-compile-time).

namespace protal {
    // Runs protal and returns the process exit code: 0, or 1 if any sample or output failed (see
    // RunStatus). Invalid input and fatal errors still exit directly with their own codes. In Run/Run.cpp.
    int Run(int argc, char *argv[]);
}
