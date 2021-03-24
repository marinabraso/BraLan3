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


for Sample in $(grep 'Branchiostoma lanceolatum' ${Metadata} | cut -f1 | tail -n +2)
do
	echo "#### $Sample"
	# Prevent quota problems (waii until there is les than 100G in the fastq folder)
	until [ $(du -s Data/RNA-seq/Marletaz2018 | cut -f1) -lt 100000000 ]
	do
		echo Waiting... $(du -s Data/RNA-seq/Marletaz2018 | cut -f1)
		sleep 1m
	done

	# Get fastq file(s) form long term storage
	if [[ $(ls ${FASTQFolder}/${Sample}* 2> ~/null | wc -l) -eq 0 ]]
	then
		echo "Getting fastq file(s) from nas..."
		cp ${NASFASTQFolder}/${Sample}* ${FASTQFolder}
	fi
	# FASTQC
	for file in $(ls ${FASTQFolder}/${Sample}*)
	do
		if [[ ! -s ${ResultsFolder}/fastqc/$(basename $file | sed 's/\..*/_fastqc.zip/g') ]]
		then
			echo Running fastqc $file .....
			sbatch -t 2:00:00 --mem=5000 -J QC$Sample -o tmp/QC$Sample.out -e tmp/QC$Sample.err ${FASTQC_PerSample} ${file} ${ResultsFolder}/fastqc
		else
			echo Fastqc $file done.
		fi
	done
	# KALLISTO
	if [[ ! -s ${ResultsFolder}/kallisto/${Sample}/abundance.h5 ]]
	then
		if [[ $(ls ${FASTQFolder}/${Sample}* 2> ~/null | wc -l) -eq 1 ]]
		then
			echo "Running Kallisto $Sample (single-end) ....."
			sbatch -t 5:00:00 --mem=10000 -J K$Sample -o tmp/K$Sample.out -e tmp/K$Sample.err ${KALLISTO_PerSample} $Sample $index_file ${ResultsFolder}/kallisto $FASTQFolder "Single" "R"
		elif [[ $(ls ${FASTQFolder}/${Sample}* 2> ~/null | wc -l) -eq 2 ]]
		then
			echo "Running Kallisto $Sample (paired-end) ....."
			sbatch -t 5:00:00 --mem=10000 -J K$Sample -o tmp/K$Sample.out -e tmp/K$Sample.err ${KALLISTO_PerSample} $Sample $index_file ${ResultsFolder}/kallisto $FASTQFolder "Paired" "RF"
		fi
	else
		echo Kallisto $Sample done, removing fastq files
		rm ${FASTQFolder}/${Sample}*fastq.gz 2> ~/null
	fi
done

