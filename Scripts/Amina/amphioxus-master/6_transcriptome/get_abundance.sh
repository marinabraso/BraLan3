# tpm out 

# kallisto: index on transcripts 

# genic
# gtf to fasta 
module add UHTS/Assembler/cufflinks/2.2.1
gffread Bla_intergenic.gtf -g Bl71nemr.fa -w Bla_intergenic.fasta
# kallisto index
module add UHTS/Analysis/kallisto/0.43.0
kallisto index -i kallisto_intergenic_idx Bla_intergenic.fasta
# kallisto quantification
kallisto quant -i kallisto_intergenic_idx -o kallisto_intergenic_quant /scratch/beegfs/monthly/aechchik/amphioxus/tmp/rnaseq/rnaseq_all/all_1.fq.gz /scratch/beegfs/monthly/aechchik/amphioxus/tmp/rnaseq/rnaseq_all/all_2.fq.gz

# intergenic
module add UHTS/Assembler/cufflinks/2.2.1
gffread Bla_annot_final_allfeatures.gtf -g Bl71nemr.fa -w Bla_genic.fasta
# kallisto index
module add UHTS/Analysis/kallisto/0.43.0
kallisto index -i kallisto_genic_idx Bla_genic.fasta
# kallisto quantification 
kallisto quant -i kallisto_genic_idx -o kallisto_genic_quant /scratch/beegfs/monthly/aechchik/amphioxus/tmp/rnaseq/rnaseq_all/all_1.fq.gz /scratch/beegfs/monthly/aechchik/amphioxus/tmp/rnaseq/rnaseq_all/all_2.fq.gz
