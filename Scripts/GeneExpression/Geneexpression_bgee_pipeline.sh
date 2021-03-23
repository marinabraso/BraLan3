#!/bin/bash

### This script run the control quality (FASTQC) + Kallisto (QUANT) per library

module add UHTS/Quality_control/fastqc/0.11.7;
module add UHTS/Analysis/kallisto/0.44.0;

# Scripts
FASTQC_PerSample=Scripts/GeneExpression/FASTQC_PreSample.sh
KALLISTO_PerSample=Scripts/GeneExpression/Kallisto_PerSample.sh

# Files & parameters
Basename="Branchiostoma_lanceolatum.BraLan3"
ResultsFolder=Results/GeneExpression
mkdir -p ${ResultsFolder}/fastqc
mkdir -p ${ResultsFolder}/kallisto
FASTQFolder=Data/RNA-seq/Marletaz2018 # Folder fort the fastq files while computing
NASFASTQFolder=/nas/FAC/FBM/DEE/mrobinso/default/D2c/mbrasovi/Marletaz2018 # Folder for the long term storage of the fastq files (to be copied to $FASTQFolder before using them)
Metadata=Metadata/Marletaz2018_RNAseq_SRA.txt 
index_file=$ResultsFolder/${Basename}_transcripts.idx # Transcriptome index file to run Kallisto


for Sample in $(grep 'Branchiostoma lanceolatum' ${Metadata} | cut -f1 | head -1)
do
	echo $Sample
	# Get fastq file(s) form long term storage
	if [[ $(ls ${FASTQFolder}/${Sample}* 2> ~/null | wc -l) -eq 0 ]]
	then
		echo "Getting fastq file(s) from nas..."
		cp ${NASFASTQFolder}/${Sample}* ${FASTQFolder}
	fi
	# FASTQC
	for file in $(ls ${FASTQFolder}/${Sample}*)
	do
		echo ${file}
		echo Running fastqc $file .....
		sbatch -t 2:00:00 --mem=5000 -J QC$Sample -o tmp/QC$Sample.out -e tmp/QC$Sample.err ${FASTQC_PerSample} $i ${ResultsFolder}/fastqc
	done
	# KALLISTO
	if [[ $(ls ${FASTQFolder}/${Sample}* 2> ~/null | wc -l) -eq 1 ]]
	then
		file=$(ls ${FASTQFolder}/${Sample}.fastq.gz)
		echo "Running Kallisto $Sample: $file (single-end) ....."
		#sbatch -t 5:00:00 --mem=10000 -J K$Sample -o tmp/K$Sample.out -e tmp/K$Sample.err ${KALLISTO_PerSample} $Sample $index_file ${ResultsFolder}/kallisto $FASTQFolder "Single" "R"
	elif [[ $(ls ${FASTQFolder}/${Sample}* 2> ~/null | wc -l) -eq 2 ]]
	then
		file1=$(ls ${FASTQFolder}/${Sample}_1.fastq.gz)
		file2=$(ls ${FASTQFolder}/${Sample}_2.fastq.gz)
		echo "Running Kallisto $Sample: $file1 $file2 (paired-end) ....."
		#sbatch -t 5:00:00 --mem=10000 -J K$Sample -o tmp/K$Sample.out -e tmp/K$Sample.err ${KALLISTO_PerSample} $Sample $index_file ${ResultsFolder}/kallisto $FASTQFolder "Paired" "RF"
	fi
done


