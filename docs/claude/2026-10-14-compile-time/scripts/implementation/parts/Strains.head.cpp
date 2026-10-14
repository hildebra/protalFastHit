// The strain-level MSAs over the samples' profiles, the strain stage of protal::Run (ProfileStage.h).

#include <algorithm>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#include <sys/wait.h>
#include "ProfileStage.h"
#include "Benchmark.h"
#include "ProgressBar.h"
#include "RunStatus.h"

namespace protal {

