#!/bin/sh -e
#SBATCH --job-name=refmap
#SBATCH --time=12:00:00
#SBATCH -A uoo03773
#SBATCH --mem 4G
#SBATCH --output Refmap.%j.out # CHANGE each run
#SBATCH --error Refmap.%j.err # CHANGE each run

echo "load Stacks, start NZ"
#remember to discard file size smaller than 50,000kb
module purge
module load Stacks # Cluster specific
ref_map.pl --samples /nesi/project/uoo03773/ALL3/Bam --popmap /nesi/project/uoo03773/ALL3/NZ_bfref.txt -T 8 -o /nesi/project/uoo03773/ALL3/NZmap

echo "finished mapping NZ, make population vcf"

populations -P /nesi/project/uoo03773/ALL3/NZmap  -M /nesi/project/uoo03773/ALL3/NZ_bfref.txt  --vcf -p 2 -r 0.2 -O /nesi/project/uoo03773/ALL3/NZClean --vcf --min-maf 0.2 

echo "done NZ, start FL"

ref_map.pl --samples /nesi/project/uoo03773/ALL3/Bam --popmap /nesi/project/uoo03773/ALL3/FL_bfref.txt -T 8 -o /nesi/project/uoo03773/ALL3/FLmap

echo "finished mapping FL, make population vcf"

populations -P /nesi/project/uoo03773/ALL3/FLmap  -M /nesi/project/uoo03773/ALL3/FL_bfref.txt  --vcf -p 2 -r 0.2 -O /nesi/project/uoo03773/ALL3/FLClean --vcf --min-maf 0.2 

echo "done FL, start AU"

ref_map.pl --samples /nesi/project/uoo03773/ALL3/Bam --popmap /nesi/project/uoo03773/ALL3/AU_bfref.txt -T 8 -o /nesi/project/uoo03773/ALL3/AUmap

echo "finished mapping AU, make population vcf"

populations -P /nesi/project/uoo03773/ALL3/AUmap  -M /nesi/project/uoo03773/ALL3/AU_bfref.txt  --vcf -p 2 -r 0.2 -O /nesi/project/uoo03773/ALL3/AUClean --vcf --min-maf 0.2

echo "done AU, job done"