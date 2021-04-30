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
# kallisto quant -i amphio -o SRR6246003 -t 10 --bias ../full_transcriptome/SRR6246003_1.fastq.gz ../full_transcriptome/SRR6246003_2.fastq.gz

kallisto quant -i amphio -o SRR6246001 -t 10 --bias ../full_transcriptome/SRR6246001_1.fastq.gz ../full_transcriptome/SRR6246001_2.fastq.gz
kallisto quant -i amphio -o SRR6246002 -t 10 --bias ../full_transcriptome/SRR6246002_1.fastq.gz ../full_transcriptome/SRR6246002_2.fastq.gz
kallisto quant -i amphio -o SRR6246003 -t 10 --bias ../full_transcriptome/SRR6246003_1.fastq.gz ../full_transcriptome/SRR6246003_2.fastq.gz
kallisto quant -i amphio -o SRR6246004 -t 10 --bias ../full_transcriptome/SRR6246004_1.fastq.gz ../full_transcriptome/SRR6246004_2.fastq.gz
kallisto quant -i amphio -o SRR6246005 -t 10 --bias ../full_transcriptome/SRR6246005_1.fastq.gz ../full_transcriptome/SRR6246005_2.fastq.gz
kallisto quant -i amphio -o SRR6246006 -t 10 --bias ../full_transcriptome/SRR6246006_1.fastq.gz ../full_transcriptome/SRR6246006_2.fastq.gz
kallisto quant -i amphio -o SRR6246007 -t 10 --bias ../full_transcriptome/SRR6246007_1.fastq.gz ../full_transcriptome/SRR6246007_2.fastq.gz
kallisto quant -i amphio -o SRR6246008 -t 10 --bias ../full_transcriptome/SRR6246008_1.fastq.gz ../full_transcriptome/SRR6246008_2.fastq.gz
kallisto quant -i amphio -o SRR6246009 -t 10 --bias ../full_transcriptome/SRR6246009_1.fastq.gz ../full_transcriptome/SRR6246009_2.fastq.gz



