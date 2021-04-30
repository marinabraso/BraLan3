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
# kallisto quant -i amphio -o SRR6246013 -t 10 --bias ../full_transcriptome/SRR6246013_1.fastq.gz ../full_transcriptome/SRR6246013_2.fastq.gz

kallisto quant -i amphio -o SRR6246011 -t 10 --bias ../full_transcriptome/SRR6246011_1.fastq.gz ../full_transcriptome/SRR6246011_2.fastq.gz
kallisto quant -i amphio -o SRR6246012 -t 10 --bias ../full_transcriptome/SRR6246012_1.fastq.gz ../full_transcriptome/SRR6246012_2.fastq.gz
kallisto quant -i amphio -o SRR6246013 -t 10 --bias ../full_transcriptome/SRR6246013_1.fastq.gz ../full_transcriptome/SRR6246013_2.fastq.gz
kallisto quant -i amphio -o SRR6246014 -t 10 --bias ../full_transcriptome/SRR6246014_1.fastq.gz ../full_transcriptome/SRR6246014_2.fastq.gz
kallisto quant -i amphio -o SRR6246015 -t 10 --bias ../full_transcriptome/SRR6246015_1.fastq.gz ../full_transcriptome/SRR6246015_2.fastq.gz
kallisto quant -i amphio -o SRR6246016 -t 10 --bias ../full_transcriptome/SRR6246016_1.fastq.gz ../full_transcriptome/SRR6246016_2.fastq.gz
kallisto quant -i amphio -o SRR6246017 -t 10 --bias ../full_transcriptome/SRR6246017_1.fastq.gz ../full_transcriptome/SRR6246017_2.fastq.gz
kallisto quant -i amphio -o SRR6246018 -t 10 --bias ../full_transcriptome/SRR6246018_1.fastq.gz ../full_transcriptome/SRR6246018_2.fastq.gz
kallisto quant -i amphio -o SRR6246019 -t 10 --bias ../full_transcriptome/SRR6246019_1.fastq.gz ../full_transcriptome/SRR6246019_2.fastq.gz
kallisto quant -i amphio -o SRR6246020 -t 10 --bias ../full_transcriptome/SRR6246020_1.fastq.gz ../full_transcriptome/SRR6246020_2.fastq.gz


