#!/bin/sh

#SBATCH --account mrobinso_default
#SBATCH --mail-user amina.echchiki@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition axiom
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20GB
#SBATCH --job-name addGeneLines
#SBATCH --export NONE

# modules

## modules
module load Bioinformatics/Software/vital-it
module add UHTS/Analysis/BEDTools/2.22.1

# add gene lines
gzip -d Bla_annot_preF_edited-InframeSTOP-wIRX-broken2_post.gtf.gz
cat Bla_annot_preF_edited-InframeSTOP-wIRX-broken2_post.gtf | grep -v ^# | sed 's/[;]/\t/g' | groupBy -g 1,9 -c 4,5 -o min,max > gene_lines
cat gene_lines | sed 's/gene_id/StringTie \t gene/g' | sed 's/ //g' | sed 's/gene"/gene\t"/g' | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t1000\t.\t.\tgene_id "$4}' > gene_lines.gtf
cat Bla_annot_preF_edited-InframeSTOP-wIRX-broken2_post.gtf gene_lines.gtf | sort -n -k1,1 -k4,4 > final.gtf

# prepare gtf; get intergenic regiosn
cat final.gtf | grep -v ^# | sed 's/$/ ; gene_biotype="protein_coding"/g' | sed 's/";  ;/";/g' > final_final.gtf
