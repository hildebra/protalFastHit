// protal::Run: the run's start-up (options, the database's gene tables, the tables and models loaded side by side), then
// its stages, each in a unit of its own (RunStages.h).

#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#include "RunProtal.h"
#include "ProfileStage.h"
#include "Benchmark.h"
#include "BuildInfo.h"
#include "CpuLog.h"
#include "MemoryUse.h"
#include "RunStatus.h"
#include "ScratchSpace.h"
#include "SeqReader.h"
#include "TaxonStatisticsOutput.h"
#include "ThreadedGzStream.h"
#include "protal_config.h"

namespace protal {

