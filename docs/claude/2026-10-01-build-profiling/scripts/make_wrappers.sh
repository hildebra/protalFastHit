#!/bin/bash
# Timing stand-ins for art_illumina and pbsim in BIN_DIR, appending to TIMES (see timed_tool.sh).
# usage: make_wrappers.sh BIN_DIR TIMES [PBSIM]
set -euo pipefail
bin=$1 times=$2 pbsim=${3:-$HOME/micromamba/envs/protal-db-build/bin/pbsim}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$bin"
for tool in art_illumina pbsim; do
  real=$([ "$tool" = pbsim ] && echo "$pbsim" || command -v art_illumina)
  printf '#!/bin/bash\nTOOL_NAME=%s TOOL_REAL=%s TOOL_TIMES=%s exec bash %s/timed_tool.sh "$@"\n' \
    "$tool" "$real" "$times" "$here" > "$bin/$tool"
  chmod +x "$bin/$tool"
done
