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
#SBATCH --time=0-02:00:00

# modules

## modules
module load Bioinformatics/Software/vital-it
module add UHTS/Analysis/kallisto/0.46.0

# idx
kallisto index -i amphio --make-unique final_pc_genic_intergenic.fasta

# quant
kallisto quant -i amphio -o SRR6245993_3 -t 10 --bias ../full_transcriptome/SRR6245993_1.fastq.gz ../full_transcriptome/SRR6245993_2.fastq.gz
