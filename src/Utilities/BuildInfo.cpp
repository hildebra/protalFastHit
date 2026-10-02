#include "BuildInfo.h"
#include "protal_config.h"
#include "protal_commit.h"

namespace protal {
    std::string VersionText() {
        std::string text = "v" + std::to_string(protal_VERSION_MAJOR) + "." + std::to_string(protal_VERSION_MINOR) +
                           "." + std::to_string(protal_VERSION_PATCH);
        std::string const commit = PROTAL_GIT_COMMIT;
        if (!commit.empty()) {
            text += " (commit " + commit + (PROTAL_GIT_MODIFIED ? ", with uncommitted changes)" : ")");
        }
        return text;
    }
}
