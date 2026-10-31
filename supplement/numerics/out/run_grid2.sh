#!/bin/bash
# Usage, from the folder that contains fem.py (copy this file and cases2.txt there first):
#   bash run_grid2.sh > run_grid2.log 2>&1
# Interpreter: a Python with numpy and scipy; override with  PY=/path/to/python bash run_grid2.sh
# Each line of cases2.txt: tag L n0 rho rmin nvec
export PY=${PY:-python3}
mkdir -p out
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3
grep -v '^#' cases2.txt | xargs -P 1 -L 1 bash -c '
  [ -f out/$0.json ] && exit 0
  "$PY" fem.py --L $1 --n0 $2 --rho $3 --rmin $4 --nvec $5 --nev 10 --out out/$0 > out/$0.out 2>&1
  echo "done $0 $(date +%T)"' 
