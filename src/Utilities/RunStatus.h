#pragma once

#include <iostream>
#include <mutex>
#include <string>
#include <vector>

namespace protal {
    // Collects failures that do not stop a run (one sample's reads, one species' post-filter, one
    // output file) but must still make protal exit non-zero, so that workflow managers can tell a
    // partial result from a complete one. The run carries on; Finish() reports and sets the exit code.
    class RunStatus {
        std::mutex m_mutex;
        std::vector<std::string> m_failures;

    public:
        static RunStatus& Get() {
            static RunStatus status;
            return status;
        }

        void Fail(std::string message) {
            std::lock_guard<std::mutex> lock(m_mutex);
            std::cerr << "[ERROR] " << message << std::endl;
            m_failures.emplace_back(std::move(message));
        }

        bool Ok() {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_failures.empty();
        }

        // Prints a summary of all failures and returns the process exit code: 0 if none, else 1.
        int Finish(std::ostream& os = std::cerr) {
            std::lock_guard<std::mutex> lock(m_mutex);
            if (m_failures.empty()) return 0;
            os << "\nprotal finished with " << m_failures.size() << " error(s):\n";
            for (auto const& failure : m_failures) os << "  - " << failure << '\n';
            os << std::flush;
            return 1;
        }
    };
}
