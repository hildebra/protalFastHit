// Binary writing and reading of plain values, strings and vectors of plain values, in the machine's own byte order:
// for files a run writes and reads back itself (the strain stage's spill files, --strain_spill), never moved to another
// machine.
#pragma once

#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

namespace protal::binary_io {
    template<typename T>
    inline void Write(std::ostream& os, T const& value) {
        static_assert(std::is_trivially_copyable_v<T>);
        os.write(reinterpret_cast<char const*>(&value), sizeof(T));
    }

    template<typename T>
    inline bool Read(std::istream& is, T& value) {
        static_assert(std::is_trivially_copyable_v<T>);
        return static_cast<bool>(is.read(reinterpret_cast<char*>(&value), sizeof(T)));
    }

    // A vector as its size and its elements.
    template<typename T>
    inline void WriteVector(std::ostream& os, std::vector<T> const& values) {
        static_assert(std::is_trivially_copyable_v<T> && !std::is_same_v<T, bool>);
        Write<uint64_t>(os, values.size());
        if (!values.empty()) os.write(reinterpret_cast<char const*>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(T)));
    }

    // Reads what WriteVector wrote; false if the stream ends before (or claims more than `max_size` elements).
    template<typename T>
    inline bool ReadVector(std::istream& is, std::vector<T>& values, uint64_t max_size = uint64_t{1} << 36) {
        static_assert(std::is_trivially_copyable_v<T> && !std::is_same_v<T, bool>);
        uint64_t size = 0;
        if (!Read(is, size) || size > max_size) return false;
        values.resize(size);
        if (size == 0) return true;
        return static_cast<bool>(is.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(size * sizeof(T))));
    }

    inline void WriteString(std::ostream& os, std::string const& s) {
        Write<uint64_t>(os, s.size());
        os.write(s.data(), static_cast<std::streamsize>(s.size()));
    }

    inline bool ReadString(std::istream& is, std::string& s, uint64_t max_size = uint64_t{1} << 32) {
        uint64_t size = 0;
        if (!Read(is, size) || size > max_size) return false;
        s.resize(size);
        return size == 0 || static_cast<bool>(is.read(s.data(), static_cast<std::streamsize>(size)));
    }
}
