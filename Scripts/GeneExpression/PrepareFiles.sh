#!/bin/bash

module add R/3.5.1;
module add UHTS/Aligner/tophat/2.1.1;
module add UHTS/Analysis/kallisto/0.44.0;


# Scripts

# Files & parameters
Basename="Branchiostoma_lanceolatum.BraLan3"

ResultsFolder="Results/GeneExpression"
mkdir -p ${ResultsFolder}

original_gene_gtf_path=Data/Transcriptomes/${Basename}.gtf.gz
gene_gtf_path=Data/Transcriptomes/${Basename}.strong.gtf
output_gtf_name=${ResultsFolder}/${Basename}
input_file_genome_file=Data/Genomes/Branchiostoma_lanceolatum.BraLan3_genome.fa


# Extract set of genes that are "protein coding" and have strong evidence
if [[ ! -s $gene_gtf_path ]]; then
	zcat $original_gene_gtf_path | grep 'strong' | grep 'protein_coding' | awk '{if($2 == "protein_coding"){print $0" gene_biotype \"protein_coding\";"}else{print $0}}' | sed 's/gene_type/gene_biotype/g' | sort -u > $gene_gtf_path
fi

################################################
### Generate new transcriptome
if [[ ! -s ${ResultsFolder}/${Basename}_transcriptome.fa ]] || [[ ! -s ${ResultsFolder}/${Basename}_transcriptome.fa.tlst ]] || [[ ! -s ${ResultsFolder}/${Basename}_transcriptome_final.fa ]]
then
	echo "Generate new transcriptome..."
	gtf_to_fasta $gene_gtf_path $input_file_genome_file  $ResultsFolder/${Basename}_transcriptome.fa
	sed -E 's/>[0-9]* />/g' $ResultsFolder/${Basename}_transcriptome.fa > $ResultsFolder/${Basename}_transcriptome_final.fa
fi

################################################
### Index for kallisto

if [[ ! -s $ResultsFolder/${Basename}_transcripts.idx ]]
then
	echo "Run kallisto index..."
	## we use by default kmer size = 31, if your read length are smaller please add -k followed by the kmer size that you intends to the command bellow
	kallisto index -i $ResultsFolder/${Basename}_transcripts.idx $ResultsFolder/${Basename}_transcriptome_final.fa
fi


if [[ ! -s ${ResultsFolder}/${Basename}_tx2gene.txt ]]
then
	cat $gene_gtf_path | sed 's/.*gene_id "\([A-Z0-9_]\+\)"; transcript_id "\([A-Z0-9_]\+\)".*/\1\t\2/g' | sort -u > ${ResultsFolder}/${Basename}_tx2gene.txt
fi











