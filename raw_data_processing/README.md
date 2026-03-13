## Process of raw data
As I used customised adapter plates, raw reads need to be cleaned and demultiplexed. <br>
All codes were writen as batch jobs submitted to New Zealand eScience Infrastructure (NeSI), therefore users need to adapt them for their own use

### Prepare: 
- a codes.txt file that matches individual adapter codes to individual samples
- a sample.txt file that containes individual sample name
- a pop.txt file that matches individuals to their respective populations
- A reference genome (if available)

### Packages used:
- Stacks2: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775. Manual: https://catchenlab.life.illinois.edu/stacks/manual/
- Cutadapt: https://journal.embnet.org/index.php/embnetjournal/article/view/200/479. Manual: https://cutadapt.readthedocs.io/en/stable/
- BWA: https://academic.oup.com/bioinformatics/article/25/14/1754/225615. Manual: https://bio-bwa.sourceforge.net/bwa.shtml
- SAMrools: https://academic.oup.com/gigascience/article/10/2/giab008/6137722. Manual: https://www.htslib.org/doc/#manual-pages
- FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- MultiQC: https://docs.seqera.io/multiqc

### Steps:
- 01_cut.sh: remove Illumina adapters and poly-G sequences
- 02_demulti.sh: Demultiplexing
- 03_concat.py: concatenate forward and reverse reads into one file for each individual sample
- 04_trimsize.sh: trim sequences to desired length, usually need to use fastqc and multiqc before to check qualities etc. before performing this step
- 05_align.sh: algin sequences to a provided reference genome, after this step needs to clean up low quality samples, the easiest way is to check file size, but can also use Stacks
- 06_refmap.sh: SNP calling
- 07_SNP_filter: there are a lot of methods to filter SNPs including using Stacks or R packages. A reference code using Stacks is included.<br>
But I mostly used "SNPfiltR" package (https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13618) as it's interactive and visualises different parameters. It has very detailed step-by-step guide on its website (https://devonderaad.github.io/SNPfiltR/). <br>
I only changes some parameters tailored to my SNP data. Therefore detailed R scripts not attached here. My filtering information can be found in the manuscipt once accepted.
