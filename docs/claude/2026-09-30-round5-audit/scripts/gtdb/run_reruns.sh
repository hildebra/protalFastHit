#!/bin/bash
T=$(dirname "$0")
for x in "$@"; do
  echo "######## $x"
  bash $T/reruns.sh $x
  echo
done > ~/audit6/gtdb/reruns_$(echo "$@" | tr -d ' ').out 2>&1
