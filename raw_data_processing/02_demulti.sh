#!/bin/bash -e

#SBATCH -A uoo03773
#SBATCH -J Demulti
#SBATCH --time 20:00:00
#SBATCH -c 4
#SBATCH --mem 32GB
#SBATCH --output demulti.%j.out
#SBATCH --error demulti.%j.err

module purge
module load Stacks/2.61-gimkl-2022a
####remember to create dictionary Demulti
echo "load stacks, start run"

process_radtags -p /nesi/project/uoo03773/NZFL2/Cut -o /nesi/project/uoo03773/NZFL2/Demulti -b /nesi/project/uoo03773/NZFL/NZFL_codes.txt -e pstI -r -c -q -P --threads 4 --inline-inline --filter-illumina

echo "done"