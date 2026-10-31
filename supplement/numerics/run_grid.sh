#!/bin/bash
# Usage, from the folder that contains fem.py:  bash run_grid.sh > run_grid.log 2>&1
# Interpreter: a Python with numpy and scipy; override with  PY=/path/to/python bash run_grid.sh
# Each line of cases.txt: tag L n0 rho rmin nvec
export PY=${PY:-python3}
mkdir -p out
export OMP_NUM_THREADS=5 OPENBLAS_NUM_THREADS=5 MKL_NUM_THREADS=5
grep -v '^#' cases.txt | xargs -P 2 -L 1 bash -c '
  [ -f out/$0.json ] && exit 0
  "$PY" fem.py --L $1 --n0 $2 --rho $3 --rmin $4 --nvec $5 --nev 10 --out out/$0 > out/$0.out 2>&1
  echo "done $0 $(date +%T)"' 
