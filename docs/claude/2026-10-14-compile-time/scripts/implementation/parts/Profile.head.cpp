// Profiling the samples, the profiling stage of protal::Run (ProfileStage.h).

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <omp.h>
#ifdef __GLIBC__
#include <malloc.h>  // malloc_trim
#endif
#include "ProfileStage.h"
#include "Benchmark.h"
#include "RunStatus.h"

namespace protal {

