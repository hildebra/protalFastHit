#!/usr/bin/env bash
# stages.sh <log>...: protal's stage timers from run logs, one column per log, in seconds.
for f in "$@"; do
  echo "== $f"
  grep -E "took|Run protal" "$f" | grep -vE "^Thread|Strain|strain|Freeing" | \
    sed -E 's/^[[:space:]]+//' | \
    awk '{
      line = $0; name = line; sub(/ took.*/, "", name)
      t = line; sub(/.* took /, "", t); sub(/ mean over.*/, "", t)
      s = 0
      if (match(t, /([0-9]+)m /, a)) s += a[1] * 60
      if (match(t, /([0-9]+)s /, a)) s += a[1]
      if (match(t, /([0-9]+)ms/, a)) s += a[1] / 1000
      if (t ~ /less than 1ms/) s = 0
      if (match(t, /^([0-9.]+)s:/, a)) s = a[1]
      printf "%-45s %8.3f\n", name, s
    }'
done
