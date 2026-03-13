#!/bin/sh -e
#SBATCH --job-name=filter
#SBATCH --time=2:00:00
#SBATCH -A uoo03773
#SBATCH --output filter.%j.out # CHANGE each run
#SBATCH --error filter.%j.err # CHANGE each run

cd /nesi/project/uoo03773/ALL3/AUClean
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

echo "start, 70%"

vcftools --vcf populations.snps.vcf --max-missing 0.7 --maf 0.05 --minDP 3 --recode --recode-INFO-all --out miss70maf5dp3
vcftools --vcf miss70maf5dp3.recode.vcf --missing-indv
awk '$5 > 0.80' out.imiss | cut -f1 > lowDP-80.indv
vcftools --vcf miss70maf5dp3.recode.vcf --remove lowDP-80.indv --recode --recode-INFO-all --out miss70maf5dp3INDV
vcftools --vcf miss70maf5dp3INDV.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out final70
vcftools --vcf final70.recode.vcf --missing-site

echo "done70%, do 50%"

vcftools --vcf populations.snps.vcf --max-missing 0.5 --maf 0.05 --minDP 3 --recode --recode-INFO-all --out miss50maf5dp3
vcftools --vcf miss50maf5dp3.recode.vcf --missing-indv
awk '$5 > 0.80' out.imiss | cut -f1 > lowDP-80.indv
vcftools --vcf miss50maf5dp3.recode.vcf --remove lowDP-80.indv --recode --recode-INFO-all --out miss50maf5dp3INDV
vcftools --vcf miss50maf5dp3INDV.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out final50
vcftools --vcf final50.recode.vcf --missing-site

echo "done50%, do 30%"

vcftools --vcf populations.snps.vcf --max-missing 0.3 --maf 0.05 --minDP 3 --recode --recode-INFO-all --out miss30maf5dp3
vcftools --vcf miss30maf5dp3.recode.vcf --missing-indv
awk '$5 > 0.80' out.imiss | cut -f1 > lowDP-80.indv
vcftools --vcf miss30maf5dp3.recode.vcf --remove lowDP-80.indv --recode --recode-INFO-all --out miss30maf5dp3INDV
vcftools --vcf miss30maf5dp3INDV.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out final30
vcftools --vcf final30.recode.vcf --missing-site

echo "done30%, do 15%"

vcftools --vcf populations.snps.vcf --max-missing 0.15 --maf 0.05 --minDP 3 --recode --recode-INFO-all --out miss15maf5dp3
vcftools --vcf miss15maf5dp3.recode.vcf --missing-indv
awk '$5 > 0.80' out.imiss | cut -f1 > lowDP-80.indv
vcftools --vcf miss15maf5dp3.recode.vcf --remove lowDP-80.indv --recode --recode-INFO-all --out miss15maf5dp3INDV
vcftools --vcf miss15maf5dp3INDV.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out final15
vcftools --vcf final15.recode.vcf --missing-site
echo "done"
