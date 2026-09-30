#!/bin/bash
O=~/audit6/gtdb/${1:-b2}
ls -la $O/protal_db $O/training_db $O/model_logs
cat $O/training/parity/parity.txt
grep -E "PMML scored|WARNING|species held out:|out of bag:" $O/trained_model.report.txt
cat $O/protal_db/build_metadata.tsv
grep -c . $O/training/training_data.tsv
# model inside the finished database: is it the trained one?
cd /tmp && rm -rf /tmp/audit6_unpack && mkdir /tmp/audit6_unpack
~/strain-build/bin/protal --db $O/protal_db --unpack_db --unpack_dir /tmp/audit6_unpack > /tmp/audit6_unpack.log 2>&1; echo "unpack exit $?"
ls -la /tmp/audit6_unpack
cmp /tmp/audit6_unpack/model_pe.xml $O/trained_model.xml && echo "model_pe.xml in database.protal == trained_model.xml"
grep -c "protal:placeholder" /tmp/audit6_unpack/model_se.xml /tmp/audit6_unpack/model_PB.xml /tmp/audit6_unpack/model_ONT.xml
# taxa of the held-out species in the finished database vs the training database
cut -f1 $O/heldout_species.txt | head -3
