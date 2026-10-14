// Aligning the samples' reads (RunWrapper), the alignment stage of protal::Run (RunStages.h).

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include "RunStages.h"
#include "Seedmap.h"
#include "Alignment/WFA2Wrapper2.h"
#include "AlignmentStrategy.h"
#include "ChainAnchorFinder.h"
#include "Classify.h"
#include "CpuLog.h"
#include "MemoryUse.h"
#include "RunStatus.h"
#include "SamFile.h"
#include "SeedExpectation.h"
#include "Stage1Probe.h"
#include "ThreadedGzStream.h"

namespace protal {

