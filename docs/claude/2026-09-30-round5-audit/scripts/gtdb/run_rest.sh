#!/bin/bash
T=$(dirname "$0")
G=~/audit6/gtdb
for x in "$@"; do
  echo "######## $x $(date +%T)"
  case $x in
    bgfail) bash $T/bgfail.sh ;;
    term) bash $T/term.sh ;;
    s|i) bash $T/early.sh $x ;;
  esac
  echo
done > $G/rest_$(echo "$@" | tr -d ' ').out 2>&1
