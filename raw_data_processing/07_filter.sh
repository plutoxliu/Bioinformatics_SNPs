#!/bin/sh -e
#SBATCH --job-name=filter
#SBATCH --time=2:00:00
#SBATCH -A uoo03773
#SBATCH --output filter.%j.out # CHANGE each run
#SBATCH --error filter.%j.err # CHANGE each run

cd /nesi/project/uoo03773/ALL3/AUClean
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

###Only example here. SNP filtering is a bit arbitrary, there's a balance betweeen filtering out uninformative SNPs, filtering out low quality SNPs, and then also ensuring that you retain enough SNPs
#--max-missing: allowed proportion of missing data
#--maf: minor allele frequency
#--minDP: minimum depth
#--missing-indv: generates a file reporting the missingness on a per-individual basis
#--minQ: can start from minimum quality of 30
#--min-meanDP: minimum mean depth of SNPs
echo "start"

vcftools --vcf populations.snps.vcf --max-missing 0.7 --maf 0.05 --minDP 3 --recode --recode-INFO-all --out miss70maf5dp3
vcftools --vcf miss70maf5dp3.recode.vcf --missing-indv
awk '$5 > 0.80' out.imiss | cut -f1 > lowDP-80.indv #make a file of the individuals that are over 80% missing
vcftools --vcf miss70maf5dp3.recode.vcf --remove lowDP-80.indv --recode --recode-INFO-all --out miss70maf5dp3INDV # filter out those individuals
vcftools --vcf miss70maf5dp3INDV.recode.vcf --min-meanDP 5 --recode --recode-INFO-all --out final70
vcftools --vcf final70.recode.vcf --missing-site #see which sites/SNPs have a lot of missingness

echo "done"
