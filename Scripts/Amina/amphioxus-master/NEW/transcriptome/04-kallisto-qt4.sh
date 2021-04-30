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
# kallisto quant -i amphio -o SRR6246023 -t 10 --bias ../full_transcriptome/SRR6246023_1.fastq.gz ../full_transcriptome/SRR6246023_2.fastq.gz

kallisto quant -i amphio -o SRR6246021 -t 10 --bias ../full_transcriptome/SRR6246021_1.fastq.gz ../full_transcriptome/SRR6246021_2.fastq.gz
kallisto quant -i amphio -o SRR6246022 -t 10 --bias ../full_transcriptome/SRR6246022_1.fastq.gz ../full_transcriptome/SRR6246022_2.fastq.gz
kallisto quant -i amphio -o SRR6246023 -t 10 --bias ../full_transcriptome/SRR6246023_1.fastq.gz ../full_transcriptome/SRR6246023_2.fastq.gz
kallisto quant -i amphio -o SRR6246024 -t 10 --bias ../full_transcriptome/SRR6246024_1.fastq.gz ../full_transcriptome/SRR6246024_2.fastq.gz
kallisto quant -i amphio -o SRR6246026 -t 10 --bias ../full_transcriptome/SRR6246026_1.fastq.gz ../full_transcriptome/SRR6246026_2.fastq.gz
kallisto quant -i amphio -o SRR6246027 -t 10 --bias ../full_transcriptome/SRR6246027_1.fastq.gz ../full_transcriptome/SRR6246027_2.fastq.gz
kallisto quant -i amphio -o SRR6246028 -t 10 --bias ../full_transcriptome/SRR6246028_1.fastq.gz ../full_transcriptome/SRR6246028_2.fastq.gz
kallisto quant -i amphio -o SRR6246029 -t 10 --bias ../full_transcriptome/SRR6246029_1.fastq.gz ../full_transcriptome/SRR6246029_2.fastq.gz

