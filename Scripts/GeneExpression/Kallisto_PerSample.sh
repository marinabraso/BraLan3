#!/bin/bash

sample=$1
GenomeIndexFile=$2
ResultsFolder=$3
FASTQFolder=$4
SingleOrPaired=$5
Strandness=$6

if [[ ${SingleOrPaired} == "Single" ]]; then
	fastqfile=$(ls ${FASTQFolder}/${sample}*fastq.gz 2> ~/null)
	if [[ ${fastqfile} != "" ]]; then
		if [[ ${Strandness} == "R" ]]; then
			kallisto quant -i ${GenomeIndexFile} -o ${ResultsFolder}/${sample} --single --rf-stranded -l 180 -s 20 $i --bias ${fastqfile}
			err=$(echo $?)
		elif [[ ${Strandness} == "F" ]]; then
			kallisto quant -i ${GenomeIndexFile} -o ${ResultsFolder}/${sample} --single --fr-stranded -l 180 -s 20 $i --bias ${fastqfile}
			err=$(echo $?)
		else
			echo "Strandness must be either R or F for Single-end reads"
			exit 0
		fi
	fi
fi
if [[ ${SingleOrPaired} == "Paired" ]]; then
	fastqfile1=$(ls ${FASTQFolder}/${sample}*1.fastq.gz 2> ~/null)
	fastqfile2=$(ls ${FASTQFolder}/${sample}*2.fastq.gz 2> ~/null)
	if [[ ${fastqfile1} != "" && ${fastqfile2} != "" ]]; then
		if [[ ${Strandness} == "RF" ]]; then
			kallisto quant -i ${GenomeIndexFile} --rf-stranded -o ${ResultsFolder}/${sample} ${fastqfile1} ${fastqfile2}
			err=$(echo $?)
		elif [[ ${Strandness} == "FR" ]]; then
			kallisto quant -i ${GenomeIndexFile} --fr-stranded -o ${ResultsFolder}/${sample} ${fastqfile1} ${fastqfile2}
			err=$(echo $?)
		else
			echo "Strandness must be either RF or FR for Paired-end reads"
			exit 0
		fi
	fi
fi

if [[ ${err} -eq 0  ]]
then
	echo "Kallisto done, removing fastq files"
	#rm ${FASTQFolder}/${sample}*.fastq.gz
else
	echo "Error in kallisto"
fi



