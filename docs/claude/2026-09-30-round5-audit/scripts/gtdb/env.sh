#!/bin/bash
exec > ~/audit6/gtdb/env.out 2>&1
mkdir -p ~/audit6/gtdb
ls -la ~/strain-build/bin
ls ~/strain-build/src | head -40
cd ~/strain-build/src && git log --oneline -1
python3 -c 'import sklearn, numpy, pandas, joblib; print(sklearn.__version__, numpy.__version__, pandas.__version__, joblib.__version__)'
python3 -c 'import numpy; print(numpy.__version__)'
python3 -c 'import pandas; print(pandas.__version__)'
python3 -c 'import sklearn; print(sklearn.__version__)'
ls ~/.venv* ~/venv* 2>/dev/null
ls ~/miniconda3 ~/mambaforge ~/micromamba 2>/dev/null | head
which conda mamba micromamba
free -m
cat /proc/loadavg
