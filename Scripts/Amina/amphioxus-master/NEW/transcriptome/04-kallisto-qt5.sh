#!/bin/sh

#SBATCH --account mrobinso_default
#SBATCH --mail-user amina.echchiki@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition axiom
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 20GB
#SBATCH --job-name gtf_to_fasta
#SBATCH --export NONE
#SBATCH --time=2-00:00:00

# modules

## modules
module load Bioinformatics/Software/vital-it
module add UHTS/Analysis/kallisto/0.46.0

# idx
#kallisto index -i amphio --make-unique final_pc_genic_intergenic.fasta

# quant
# kallisto quant -i amphio -o SRR6246033 -t 10 --bias ../full_transcriptome/SRR6246033_1.fastq.gz ../full_transcriptome/SRR6246033_2.fastq.gz

kallisto quant -i amphio -o SRR6246031 -t 10 --bias ../full_transcriptome/SRR6246031_1.fastq.gz ../full_transcriptome/SRR6246031_2.fastq.gz
kallisto quant -i amphio -o SRR6246032 -t 10 --bias ../full_transcriptome/SRR6246032_1.fastq.gz ../full_transcriptome/SRR6246032_2.fastq.gz
kallisto quant -i amphio -o SRR6246030 -t 10 --bias ../full_transcriptome/SRR6246030_1.fastq.gz ../full_transcriptome/SRR6246030_2.fastq.gz
kallisto quant -i amphio -o SRR6246035 -t 10 --bias ../full_transcriptome/SRR6246035_1.fastq.gz ../full_transcriptome/SRR6246035_2.fastq.gz
kallisto quant -i amphio -o SRR6246036 -t 10 --bias ../full_transcriptome/SRR6246036_1.fastq.gz ../full_transcriptome/SRR6246036_2.fastq.gz
kallisto quant -i amphio -o SRR6246037 -t 10 --bias ../full_transcriptome/SRR6246037_1.fastq.gz ../full_transcriptome/SRR6246037_2.fastq.gz
kallisto quant -i amphio -o SRR6246038 -t 10 --bias ../full_transcriptome/SRR6246038_1.fastq.gz ../full_transcriptome/SRR6246038_2.fastq.gz
kallisto quant -i amphio -o SRR6246039 -t 10 --bias ../full_transcriptome/SRR6246039_1.fastq.gz ../full_transcriptome/SRR6246039_2.fastq.gz

