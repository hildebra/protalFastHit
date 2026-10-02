#!/usr/bin/env bash
for c in pe500k ont90M; do for b in base cl; do echo "== $c $b"; cat ~/mt-work/onebtime/out.$c.$b/misc/s_runtime.tsv; done; done
