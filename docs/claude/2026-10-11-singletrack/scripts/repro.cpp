// repro.cpp: which alignment of the failing PacBio pair (pattern.txt = reference piece, text.txt = read piece) makes
// WFA2's Singletrack backtrace exit. Each Singletrack alignment runs in a forked child (the library calls exit(-1));
// the classic one in the parent gives the expected outcome.
#include <climits>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include "WFA2Wrapper2.h"

using namespace protal;

static std::string Read(const char* path) {
    std::ifstream in(path);
    std::string s;
    std::getline(in, s);
    return s;
}

int main(int argc, char** argv) {
    std::string const ref = Read(argv[1]), read = Read(argv[2]);
    std::printf("ref %zu, read %zu\n", ref.size(), read.size());
    int const R = static_cast<int>(ref.size()), Q = static_cast<int>(read.size());
    int failures = 0, runs = 0;
    // Without WFA2's default wf-adaptive heuristic: the failing ends-free configuration.
    for (bool adaptive : { true, false }) {
        int const qe = Q, rb = 0, re = 0;
        WFA2Wrapper2 classic(4, 6, 2, 0, WfaTraceback::Classic);
        if (!adaptive) classic.GetAligner().setHeuristicNone();
        classic.Alignment(read, ref, 0, qe, rb, re, INT_MAX);
        fflush(stdout);
        pid_t const pid = fork();
        if (pid == 0) {
            freopen("/dev/null", "w", stderr);
            WFA2Wrapper2 single(4, 6, 2, 0, WfaTraceback::Singletrack);
            if (!adaptive) single.GetAligner().setHeuristicNone();
            single.Alignment(read, ref, 0, qe, rb, re, INT_MAX);
            _exit(single.Success() && single.GetAlignmentScore() == classic.GetAlignmentScore() ? 0 : 3);
        }
        int status = 0;
        waitpid(pid, &status, 0);
        std::printf("wf-adaptive %s: classic score %d, singletrack exit %d\n", adaptive ? "on" : "off",
                    classic.Success() ? -classic.GetAlignmentScore() : -1, WIFEXITED(status) ? WEXITSTATUS(status) : -1);
    }
    for (int mode = 0; mode < 2; mode++)
    for (size_t x_drop : { size_t{0}, size_t{1000} })
    for (int max_score : { INT_MAX, 2000, 1000, 600, 330 })
    for (int qe : { 0, Q / 2, Q })
    for (int rb : { 0, R / 2, R })
    for (int re : { 0, R / 2, R }) {
        if (mode == 0 && (qe || rb || re)) continue;  // end to end: once per budget
        WFA2Wrapper2 classic(4, 6, 2, x_drop, WfaTraceback::Classic);
        if (mode == 0) classic.EndToEnd(read, ref, max_score);
        else classic.Alignment(read, ref, 0, qe, rb, re, max_score);
        runs++;
        fflush(stdout);
        pid_t const pid = fork();
        if (pid == 0) {
            freopen("/dev/null", "w", stderr);
            WFA2Wrapper2 single(4, 6, 2, x_drop, WfaTraceback::Singletrack);
            if (mode == 0) single.EndToEnd(read, ref, max_score);
            else single.Alignment(read, ref, 0, qe, rb, re, max_score);
            bool const same = single.Success() == classic.Success() &&
                              (!single.Success() || single.GetAlignmentScore() == classic.GetAlignmentScore());
            _exit(same ? 0 : 3);
        }
        int status = 0;
        waitpid(pid, &status, 0);
        int const code = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
        if (code != 0) {
            failures++;
            std::printf("%s x_drop %zu max_score %d q_end_free %d r_begin_free %d r_end_free %d: classic %s (score %d); "
                        "singletrack %s\n", mode ? "endsfree" : "end2end", x_drop, max_score, qe, rb, re,
                        classic.Success() ? "aligned" : "failed", classic.Success() ? -classic.GetAlignmentScore() : -1,
                        code == 255 ? "EXIT(-1) in the backtrace" : code == 3 ? "other outcome or score" : "crashed");
        }
    }
    std::printf("%d of %d configurations fail\n", failures, runs);
    return 0;
}
