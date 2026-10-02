#pragma once
#include <string>

namespace protal {
    // What --version prints after the binary's name: "v0.7.3 (commit <sha>)", with ", with uncommitted changes" when
    // the source had them at the build, or "v0.7.3" for a build outside a git checkout (protal_commit.cmake).
    // build_gtdb_database.py reads it to check that protal and the simulator were built from the source its scripts
    // are at. In BuildInfo.cpp, so that a new commit recompiles only that file.
    std::string VersionText();
}
