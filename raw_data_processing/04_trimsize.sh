#!/bin/bash -e
#SBATCH -A uoo03773
#SBATCH -J TrimSize
#SBATCH --time 12:00:00
#SBATCH -c 8
#SBATCH --output TrimSize.%j.out # CHANGE each run
#SBATCH --error TrimSize.%j.err # CHANGE each run

module load cutadapt/4.4-gimkl-2022a-Python-3.11.3
###make sure to QC (fastqc and multiqc) before trimming to get a sense of the mean sequence length, and clean up samples with zero reads (i.e. delete them from the pop.txt and sample.txt file)
#-l trim length
#-o output name
#-m minimum length
#-j number of threads

cd /nesi/project/uoo03773/GBS/Trim

for sample in *.fq.gz
do
base=$(basename ${sample} .fq.gz)
echo "${base}"
cutadapt -l 100 -o "${base}".trim.fq.gz "${base}".fq.gz -m 65 -j 8
done

echo "job done"
