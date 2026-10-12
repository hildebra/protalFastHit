#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <unordered_map>
#include <utility>

namespace protal {
    // A run-wide cache of values by a 64-bit key (taxid << 32 | gene), shared by the threads, in shards that each have
    // their own lock. The profiling's threads look up the gene copies of every SAM chunk they parse (ancestry sites,
    // column weights): some million lookups in an r226 paired-end sample. Behind one lock for all of them, 32 threads
    // spent most of the parse waiting (12 -> 119 thread-seconds, docs/claude/2026-10-12-feature-cost-vs-f1.md). A value
    // depends on its key only, so which thread makes it, and how often, does not change what a lookup finds.
    // A shard starts over when it holds max_entries / kShards values (the cache as a whole: about max_entries).
    template<typename Value>
    class ShardedCache {
    public:
        static constexpr size_t kShards = 64;

        explicit ShardedCache(size_t max_entries) : m_per_shard(std::max<size_t>(1, max_entries / kShards)) {}

        ShardedCache(ShardedCache const&) = delete;
        ShardedCache& operator=(ShardedCache const&) = delete;

        // The value under key: the cached one, or make()'s, made outside the lock and kept. When two threads make it
        // at once, both return the one kept first.
        template<typename Make>
        Value Lookup(uint64_t key, Make&& make) {
            Shard& shard = ShardOf(key);
            {
                std::lock_guard<std::mutex> lock(shard.mutex);
                if (auto const it = shard.values.find(key); it != shard.values.end()) return it->second;
            }
            Value value = make();
            std::lock_guard<std::mutex> lock(shard.mutex);
            if (shard.values.size() >= m_per_shard) shard.values.clear();
            return shard.values.try_emplace(key, std::move(value)).first->second;
        }

        size_t Size() const {
            size_t size = 0;
            for (auto const& shard : m_shards) {
                std::lock_guard<std::mutex> lock(shard.mutex);
                size += shard.values.size();
            }
            return size;
        }

        void Clear() {
            for (auto& shard : m_shards) {
                std::lock_guard<std::mutex> lock(shard.mutex);
                shard.values.clear();
            }
        }

    private:
        // A cache line of its own, so that two shards' locks do not share one.
        struct alignas(64) Shard {
            mutable std::mutex mutex;
            std::unordered_map<uint64_t, Value> values;
        };

        // The key's bits mixed (Fibonacci hashing), so that a species' genes and a gene's species spread over the shards.
        Shard& ShardOf(uint64_t key) {
            return m_shards[(key * 0x9E3779B97F4A7C15ull) >> 58];
        }
        static_assert(kShards == 64, "ShardOf takes the top 6 bits");

        size_t m_per_shard;
        std::array<Shard, kShards> m_shards;
    };
}
