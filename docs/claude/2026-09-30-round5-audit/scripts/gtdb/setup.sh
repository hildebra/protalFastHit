#!/bin/bash
exec > ~/audit6/gtdb/setup.out 2>&1
PY=~/protal-train/bin/python
$PY -c 'import sklearn, numpy, pandas, joblib, sys; print(sys.version); print(sklearn.__version__, numpy.__version__, pandas.__version__, joblib.__version__)'
mkdir -p ~/audit6/gtdb
rsync -a --exclude /build ~/strain-build/src/ ~/audit6/gtdb/src/
ls ~/audit6/gtdb/src/scripts
# compare with the worktree copy
diff -rq ~/audit6/gtdb/src/scripts /mnt/c/Users/hildebra/Documents/locDev/protal/.claude/worktrees/strain-fixes/scripts | head
diff -q ~/audit6/gtdb/src/docs/building-a-database.md /mnt/c/Users/hildebra/Documents/locDev/protal/.claude/worktrees/strain-fixes/docs/building-a-database.md
~/strain-build/bin/protal --version 2>&1 | head -3
~/strain-build/bin/simulate_metagenomes --help 2>&1 | head -60
