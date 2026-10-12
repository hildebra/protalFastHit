#pragma once

#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>
#include <system_error>
#include <unistd.h>

namespace protal::scratch {
    // The run's folder for intermediate files (--scratch_dir DIR): DIR/protal_<host>_<pid>_<milliseconds>, made by Make and
    // removed, with everything in it, when the process ends: at the end of a run, and when protal exits early (exit() runs
    // the destructors of static objects such as this one). A killed run (SIGKILL, a cluster's time limit) leaves it; its
    // name says which machine and process made it.
    class Space {
    public:
        static Space& Get() {
            static Space space;
            return space;
        }

        Space(Space const&) = delete;
        Space& operator=(Space const&) = delete;

        ~Space() {
            Remove();
        }

        // Makes the run's folder in `dir` (once; later calls return it). Its path, or empty with `error` set.
        std::string const& Make(std::string const& dir, std::string& error) {
            if (!m_path.empty() || dir.empty()) return m_path;
            char host[256] = {};
            if (::gethostname(host, sizeof(host) - 1) != 0 || host[0] == '\0') std::snprintf(host, sizeof(host), "host");
            auto const ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
            std::filesystem::path const path = std::filesystem::path(dir) / ("protal_" + std::string(host) + "_" + std::to_string(::getpid()) + "_" +
                                                                             std::to_string(ms));
            std::error_code ec;
            std::filesystem::create_directories(dir, ec);
            if (!std::filesystem::create_directory(path, ec) || ec) {
                error = "cannot make the scratch folder " + path.string() + (ec ? ": " + ec.message() : std::string(" (it exists)"));
                return m_path;
            }
            m_path = path.string();
            return m_path;
        }

        // The run's folder, or empty without one.
        std::string const& Path() const {
            return m_path;
        }

        // Removes the folder and what is in it (once; errors ignored: a scratch folder left is no failure of the run).
        void Remove() {
            if (m_path.empty()) return;
            std::error_code ec;
            std::filesystem::remove_all(m_path, ec);
            m_path.clear();
        }

    private:
        Space() = default;

        std::string m_path;
    };
}
