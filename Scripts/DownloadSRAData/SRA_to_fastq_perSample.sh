#!/bin/bash


outdir=$1
sample=$2
		
fastq-dump --gzip --split-3 -O ${outdir} ${outdir}/${sample}/${sample}.sra

err=$(echo $?)
fastqpresent=$(ls ${outdir}/${sample}*fastq.gz 2> /dev/null | wc -l)
if [[ ${fastqpresent} -gt 0 ]] && [[ ${err} -eq 0  ]]
then
	echo "Fastq(s) present, removing sra"
	rm -f ${outdir}/${sample}/${sample}.sra
	rmdir ${outdir}/${sample}
fi
