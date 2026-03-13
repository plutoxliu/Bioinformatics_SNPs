#!/bin/bash -e

#SBATCH --job-name=concat
#SBATCH --time=24:00:00
#SBATCH -c 4
#SBATCH -A uoo03773
#SBATCH --output Concat.%j.out # CHANGE each run
#SBATCH --error Concat.%j.err # CHANGE each run

python /nesi/project/uoo03773/03_concat.py

echo "done"