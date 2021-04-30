#!/bin/sh

#SBATCH --account mrobinso_default
#SBATCH --mail-user amina.echchiki@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition axiom
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20GB
#SBATCH --job-name gtf_to_fasta
#SBATCH --export NONE

# modules

## modules
module load Bioinformatics/Software/vital-it
module add UHTS/Assembler/cufflinks/2.2.1
module add UHTS/Analysis/kallisto/0.46.0
module add UHTS/Analysis/samtools/1.8


samtools faidx genome.fasta


# fix overhang
#perl cufftrim.pl Blan_vf3_fixname.fa.fai final_intergenic.gtf_all > final_intergenic.gtf_fixed.gtf
# convert to fasta
#gffread final_intergenic.gtf_fixed.gtf -g Blan_vf3_fixname.fa -w final_intergenic.fasta



gffread final_pc_genic_intergenic.gtf_all -g genome.fasta -w final_pc_genic_intergenic.fasta
