#!/bin/bash
# run_detached.sh SCRIPT: runs SCRIPT (a path under this directory) detached from the wsl.exe
# client, output to ~/audit6/build/<name>.out, exit code to ~/audit6/build/<name>.rc.
set -u
here=$(dirname "$0")
name=$(basename "$1" .sh)
out=~/audit6/build/$name.out
rm -f ~/audit6/build/$name.rc
nohup setsid bash -c "bash '$here/$1' > '$out' 2>&1 < /dev/null; echo \$? > ~/audit6/build/$name.rc" > /dev/null 2>&1 < /dev/null &
disown
echo "started $name (pid $!)"
