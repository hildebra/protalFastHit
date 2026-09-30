#!/usr/bin/env bash
# --add_model in bundle mode and folder mode; the stale-bundle lead.
set -u
P=~/strain-build/bin/protal
W=~/audit6/database
S=$(dirname "$0")
run(){ "$@" >/tmp/am.log 2>&1; local rc=$?; echo "  exit=$rc ($(basename $1) ${@:2})"; grep -iE 'error|fail|warn|stored|models for|already|unchanged' /tmp/am.log | head -8 | sed 's/^/    /'; return $rc; }

echo "=== A. bundle mode: add a model_se.xml (reuse model_pe.xml as the PMML)"
rm -rf $W/am_bundle; mkdir -p $W/am_bundle
cp $W/bundle_only/database.protal $W/am_bundle/database.protal
cp $W/folder/model_pe.xml $W/am_bundle/some_model.xml
before=$(stat -c%s $W/am_bundle/database.protal)
run $P --add_model $W/am_bundle/some_model.xml --read_type se --db $W/am_bundle/database.protal -t 2
after=$(stat -c%s $W/am_bundle/database.protal)
echo "  bundle size: $before -> $after"
echo "  members now:"; python3 $S/dbinfo.py $W/am_bundle/database.protal | sed 's/^/    /'
ls $W/am_bundle/*.partial 2>/dev/null && echo "    LEFTOVER .partial!" || echo "    no leftover .partial"

echo "=== B. folder mode: folder has separate files AND a stale database.protal"
rm -rf $W/am_folder; mkdir -p $W/am_folder
for f in $W/folder/*; do cp "$f" $W/am_folder/; done       # separate files (index.prx.zst etc.)
cp $W/bundle_only/database.protal $W/am_folder/database.protal   # a bundle sits next to them
b_before=$(md5sum $W/am_folder/database.protal | cut -d' ' -f1)
run $P --add_model $W/folder/model_pe.xml --read_type se --db $W/am_folder -t 2
b_after=$(md5sum $W/am_folder/database.protal | cut -d' ' -f1)
echo "  wrote separate model_se.xml? $([ -f $W/am_folder/model_se.xml ] && echo yes || echo no)"
echo "  database.protal md5 before=$b_before after=$b_after $([ "$b_before" = "$b_after" ] && echo '(bundle UNCHANGED -> now stale: its model_se differs from the folder)' || echo '(bundle updated)')"
# Does the stale bundle still hold NO model_se, while the folder now does?
python3 $S/dbinfo.py $W/am_folder/database.protal | grep -c model_se | sed 's/^/    bundle model_se members: /'

echo "=== C. add_model rejects a non-file and a bad --read_type"
run $P --add_model /no/such/file.xml --read_type se --db $W/am_bundle/database.protal -t 2
run $P --add_model $W/folder/model_pe.xml --read_type XX --db $W/am_bundle/database.protal -t 2

echo "=== D. add_model with --compress_frame_mb 0 on a bundle (frames required)"
run $P --add_model $W/folder/model_pe.xml --read_type se --db $W/am_bundle/database.protal -t 2 --compress_frame_mb 0

echo "=== E. add_model with an out-of-range --compress_level (validated?)"
run $P --add_model $W/folder/model_pe.xml --read_type se --db $W/am_bundle/database.protal -t 2 --compress_level 99
