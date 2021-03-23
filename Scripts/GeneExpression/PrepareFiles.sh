#!/bin/bash

module add R/3.5.1;
module add UHTS/Aligner/tophat/2.1.1;
module add UHTS/Analysis/kallisto/0.44.0;


# Scripts
bgee_RNA_Seq="../bgee_pipeline/pipeline/RNA_Seq"
prepare_GTF_script=${bgee_RNA_Seq}"/0Before/prepare_GTF.R"


# Files & parameters
Basename="Branchiostoma_lanceolatum.BraLan3"

ResultsFolder="Results/GeneExpression"
mkdir -p ${ResultsFolder}

gene_gtf_path=Data/Transcriptomes/${Basename}.gtf.gz
output_gtf_name=${ResultsFolder}/${Basename}
input_file_genome_file=Data/Genomes/Branchiostoma_lanceolatum.BraLan3_genome.fa

################################################
### Prepare GTF
if [[ ! -e ${ResultsFolder}/${Basename}.gene2biotype ]] || [[ ! -s ${ResultsFolder}/${Basename}.gene2transcript ]] || [[ ! -s ${ResultsFolder}/${Basename}.gtf_all ]]
then
	echo "Prepare gtf..."
	R CMD BATCH --no-save --no-restore "--args gene_gtf_path=\"$gene_gtf_path\" output_gtf_path=\"$output_gtf_name\"" $prepare_GTF_script
fi

################################################
### Generate new transcriptome
if [[ ! -s ${ResultsFolder}/${Basename}_transcriptome.fa ]] || [[ ! -s ${ResultsFolder}/${Basename}_transcriptome.fa.tlst ]] || [[ ! -s ${ResultsFolder}/${Basename}_transcriptome_final.fa ]]
then
	echo "Generate new transcriptome..."
	gtf_to_fasta $ResultsFolder/${Basename}.gtf_all $input_file_genome_file  $ResultsFolder/${Basename}_transcriptome.fa
	sed -E 's/>[0-9]* />/g' $ResultsFolder/${Basename}_transcriptome.fa > $ResultsFolder/${Basename}_transcriptome_final.fa
fi

################################################
### Index for kallisto
echo "Run kallisto index..."

if [[ ! -s $ResultsFolder/${Basename}_transcripts.idx ]]
then
	## we use by default kmer size = 31, if your read length are smaller please add -k followed by the kmer size that you intends to the command bellow
	kallisto index -i $ResultsFolder/${Basename}_transcripts.idx $ResultsFolder/${Basename}_transcriptome_final.fa
fi
















