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
#-p: input sequence file
#-o: output dictionary
#-b: barcodes info
#-e: restriction enzyme used
#-r: rescue barcodes and RAD-Tag cut sites
#-c: clean data, remove any read with an uncalled base
#-q: discard reads with low quality scores
#-P: files contained within the directory are paired
#--threads: number of threads to run
#--inline-inline: barcode is inline with sequence, occurs on single and paired-end read
#--filter-illumina: discard reads that have been marked by Illumina’s chastity/purity filter as failing
echo "load stacks, start run"

process_radtags -p /nesi/project/uoo03773/NZFL2/Cut -o /nesi/project/uoo03773/NZFL2/Demulti -b /nesi/project/uoo03773/NZFL/NZFL_codes.txt -e pstI -r -c -q -P --threads 4 --inline-inline --filter-illumina

echo "done"
