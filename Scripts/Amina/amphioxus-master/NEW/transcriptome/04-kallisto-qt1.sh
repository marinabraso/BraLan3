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
# kallisto quant -i amphio -o SRR6245993 -t 10 --bias ../full_transcriptome/SRR6245993_1.fastq.gz ../full_transcriptome/SRR6245993_2.fastq.gz
kallisto quant -i amphio -o SRR6245987 -t 10 --bias ../full_transcriptome/SRR6245987_1.fastq.gz ../full_transcriptome/SRR6245987_2.fastq.gz
kallisto quant -i amphio -o SRR6245989 -t 10 --bias ../full_transcriptome/SRR6245989_1.fastq.gz ../full_transcriptome/SRR6245989_2.fastq.gz
kallisto quant -i amphio -o SRR6245990 -t 10 --bias ../full_transcriptome/SRR6245990_1.fastq.gz ../full_transcriptome/SRR6245990_2.fastq.gz
kallisto quant -i amphio -o SRR6245991 -t 10 --bias ../full_transcriptome/SRR6245991_1.fastq.gz ../full_transcriptome/SRR6245991_2.fastq.gz
kallisto quant -i amphio -o SRR6245992 -t 10 --bias ../full_transcriptome/SRR6245992_1.fastq.gz ../full_transcriptome/SRR6245992_2.fastq.gz
kallisto quant -i amphio -o SRR6245993 -t 10 --bias ../full_transcriptome/SRR6245993_1.fastq.gz ../full_transcriptome/SRR6245993_2.fastq.gz
kallisto quant -i amphio -o SRR6245994 -t 10 --bias ../full_transcriptome/SRR6245994_1.fastq.gz ../full_transcriptome/SRR6245994_2.fastq.gz
kallisto quant -i amphio -o SRR6245995 -t 10 --bias ../full_transcriptome/SRR6245995_1.fastq.gz ../full_transcriptome/SRR6245995_2.fastq.gz
kallisto quant -i amphio -o SRR6245996 -t 10 --bias ../full_transcriptome/SRR6245996_1.fastq.gz ../full_transcriptome/SRR6245996_2.fastq.gz
kallisto quant -i amphio -o SRR6245997 -t 10 --bias ../full_transcriptome/SRR6245997_1.fastq.gz ../full_transcriptome/SRR6245997_2.fastq.gz
kallisto quant -i amphio -o SRR6245998 -t 10 --bias ../full_transcriptome/SRR6245998_1.fastq.gz ../full_transcriptome/SRR6245998_2.fastq.gz
kallisto quant -i amphio -o SRR6245999 -t 10 --bias ../full_transcriptome/SRR6245999_1.fastq.gz ../full_transcriptome/SRR6245999_2.fastq.gz
