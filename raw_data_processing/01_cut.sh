#!/bin/bash -e

#SBATCH -A uoo03773
#SBATCH -J CutAdapt
#SBATCH --time 12:00:00
#SBATCH -c 8
#SBATCH --output 01CutAdapt.%j.out
#SBATCH --error 01CutAdapt.%j.err

module load cutadapt/4.4-gimkl-2022a-Python-3.11.3

cd /nesi/project/uoo03773/NZFL2

echo "start lane1"
echo "cut adapters"
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -j 8 -o ./cut1.fq.gz -p ./cut2.fq.gz ./9394_R1_001.fastq.gz ./9394_R2_001.fastq.gz -m 30:30 -q 25

echo "homopolymers"
cutadapt -a GGGGGGGGGGGGG -A GGGGGGGGGGGGG -m 30:30  -q 25 --overlap 5 -o Cut_R1_001.fastq -p ./Cut_R2_001.fastq ./cut1.fq.gz ./cut2.fq.gz

echo "done"
