// The database's tables and index, as protal::Run loads them (RunStages.h).

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include "RunStages.h"
#include "Seedmap.h"
#include "Stage1Probe.h"
#include "SequenceUtils/GridLevels.h"
#include "MemoryUse.h"
#include "ReferenceFingerprint.h"

namespace protal {

    ProtalDB::~ProtalDB() = default;

    void ProtalDB::SetIndexLoading(std::future<LoadedIndex> index) { m_index = std::move(index); }

    std::optional<LoadedIndex> ProtalDB::TakeIndex() {
        if (!m_index) return std::nullopt;
        auto loaded = m_index->get();
        m_index.reset();
        return loaded;
    }

