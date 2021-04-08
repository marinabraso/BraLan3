#!/bin/bash

module add R/3.5.1;

# Scripts
bgee_RNA_Seq="../bgee_pipeline/pipeline/RNA_Seq"
rna_seq_analysis=${bgee_RNA_Seq}"/1Run/rna_seq_analysis.R"
rna_seq_sum=${bgee_RNA_Seq}"/1Run/rna_seq_sum_by_species.R"

# Files & parameters
Basename="Branchiostoma_lanceolatum.BraLan3"
ResultsFolder=Results/GeneExpression
sample_info="Metadata/Marletaz2018_RNAseq_bgee_sample_info.txt"
sample_excluded="Metadata/Marletaz2018_RNAseq_bgee_sample_excluded.txt"

for folder in $(ls ${ResultsFolder}/kallisto)
do
	if [[ ! -s ${ResultsFolder}/kallisto/$folder/abundance+gene_id+fpkm+intergenic.tsv ]]
	then
		echo "Analysis $folder"
		R CMD BATCH --no-save --no-restore "--args  kallisto_count_folder=\"${ResultsFolder}/kallisto/$folder\" gene2transcript_file=\"${ResultsFolder}/${Basename}.gene2transcript\" gene2biotype_file=\"${ResultsFolder}/${Basename}.gene2biotype\" library_id=\"${ResultsFolder}/kallisto/$folder\"" $rna_seq_analysis ${ResultsFolder}/kallisto/$folder/library_id.Rout
	fi
done









