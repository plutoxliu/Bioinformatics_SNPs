#!/bin/bash -e
#SBATCH -A uoo03773
#SBATCH -J MacroSNP
#SBATCH --time 20:00:00
#SBATCH -c 8
#SBATCH --output MacroSNP.%j.out # CHANGE each run
#SBATCH --error MacroSNP.%j.err # CHANGE each run

module load BWA SAMtools
cd /nesi/project/uoo03773/GBS/ReTrim
echo "load index"
bwa index Macro_genome.fna #load reference genome file

echo "SNP call"
#List of individual sample names, make sure one sample per line (some examples here)
files="16m1
16m10
16m2
16m3
16m4
16m5
"
#source dictionary
src=/nesi/project/uoo03773/GBS/ReTrim/

#reference genome file (if there is no reference genome check bwa manual for de novo alignment)
bwa_db=Macro_genome.fna

#align to genome
for sample in $files
do 
    bwa mem -t 8 $bwa_db $src/${sample}.trim.fq.gz  |   samtools view -b | samtools sort --threads 4 > ${sample}.bam
done

echo "align done"
