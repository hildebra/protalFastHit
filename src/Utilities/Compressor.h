#pragma once

#include <cerrno>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <spawn.h>
#include <sys/wait.h>

extern char **environ;

class Compressor {
public:
    // Compresses a file in place with pigz: `inputFile` becomes `inputFile`.gz and the original is
    // removed. pigz is started directly, without a shell, so paths with spaces or shell
    // metacharacters are safe. -n leaves the file name and time stamp out of the gzip header, so the
    // same input always gives byte-identical output.
    static void compressInPlace(const std::filesystem::path& inputFile, int threads, const std::string& pigz = "pigz")
    {
        namespace fs = std::filesystem;

        // Validate inputs
        if (threads <= 0)
            throw std::invalid_argument("Thread count must be greater than zero.");

        if (inputFile.empty())
            throw std::invalid_argument("Input file path is empty.");

        if (!fs::exists(inputFile))
            throw std::invalid_argument("Input file does not exist: " + inputFile.string());

        std::vector<std::string> args = { pigz, "-f", "-n", "-p", std::to_string(threads), "--", inputFile.string() };
        std::vector<char*> argv;
        for (auto& arg : args) argv.push_back(arg.data());
        argv.push_back(nullptr);

        pid_t pid;
        int spawn_error = posix_spawnp(&pid, pigz.c_str(), nullptr, nullptr, argv.data(), environ);
        if (spawn_error != 0) {
            throw std::runtime_error("Could not start " + pigz + " (is it installed and on PATH?): " + std::strerror(spawn_error));
        }

        int status = 0;
        while (waitpid(pid, &status, 0) == -1) {
            if (errno != EINTR) throw std::runtime_error("Waiting for " + pigz + " failed: " + std::strerror(errno));
        }
        if (!WIFEXITED(status)) {
            throw std::runtime_error(pigz + " was terminated by signal " + std::to_string(WTERMSIG(status)) + " on " + inputFile.string());
        }
        if (WEXITSTATUS(status) != 0) {
            throw std::runtime_error(pigz + " exited with status " + std::to_string(WEXITSTATUS(status)) + " on " + inputFile.string());
        }

        // Verify compression result
        fs::path gzFile = inputFile;
        gzFile += ".gz";

        if (!fs::exists(gzFile)) {
            throw std::runtime_error("pigz completed but output file not found: " + gzFile.string());
        }
    }
};
