# Shared helpers for the I/O audit scripts (sourced).
S=$(dirname "$0")
P=${P:-$HOME/strain-build/bin/protal}
W=$HOME/audit6/io
DB=$W/db
R=$W/reads

# samtext FILE: SAM text of a plain/.gz/.zst SAM
samtext() {
  case "$1" in
    *.zst) zstd -dcq "$1" ;;
    *.gz) gzip -dc "$1" ;;
    *) cat "$1" ;;
  esac
}
# records FILE: records only, sorted
records() { samtext "$1" | grep -v '^@' | LC_ALL=C sort; }
# run protal quietly: prun LOG args...
prun() {
  local log=$1; shift
  timeout 600 "$P" "$@" > "$log" 2>&1
  local rc=$?
  echo "rc=$rc"
  return 0
}
# profiles in a dir, concatenated with names
profiles() { for f in $(find "$1" -name "*.profile*" -type f | sort); do echo "== ${f#$1/}"; cat "$f"; done; }
