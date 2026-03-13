#!/bin/sh -e
#SBATCH --job-name=refmap
#SBATCH --time=12:00:00
#SBATCH -A uoo03773
#SBATCH --mem 4G
#SBATCH --output Refmap.%j.out # CHANGE each run
#SBATCH --error Refmap.%j.err # CHANGE each run

module purge
module load Stacks # Cluster specific

###remember to quality check all the bam files. I usually discard file sizes <50,000kb, then need to create a cleaned population file too
###create dictionary Map & Clean
###the parameters are pretty self-explanatory, and generally rule of thumb
#-p minimum number of populations
#-r minimum percentage of individuals in a population required to process a locus for that population
#--min-maf: minimum minor allele frequency

echo "start SNP calling"

ref_map.pl --samples /nesi/project/uoo03773/ALL3/Bam --popmap /nesi/project/uoo03773/ALL3/cleaned_pop.txt -T 8 -o /nesi/project/uoo03773/ALL3/Map

echo "finished calling, make population vcf"

populations -P /nesi/project/uoo03773/ALL3/Map  -M /nesi/project/uoo03773/ALL3/cleaned_pop.txt  --vcf -p 2 -r 0.2 -O /nesi/project/uoo03773/ALL3/AUClean --vcf --min-maf 0.2

echo "job done"
