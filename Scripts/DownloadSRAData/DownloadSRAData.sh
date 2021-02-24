#!/bin/bash

# Scripts
PerSampleDownloadScript="./Scripts/DownloadSRAData/DownloadSRAData_perSample.sh"
PerSampleSraToFastqScript="./Scripts/DownloadSRAData/SRA_to_fastq_perSample.sh"

# Files & parameters
DataFolder=$1
mkdir -p ${DataFolder}
SRAIDs=$2
errdir="Scripts/DownloadSRAData/errdir"
mkdir -p ${errdir}

while IFS= read -r sample
do
	echo "$sample"
	mkdir -p ${DataFolder}/${sample}

	jobidD=0
	jobidF=0

	# Download with prefetch
	fastqpresent=$(ls ${DataFolder}/${sample}*fastq.gz 2> ~/null | wc -l)
	if [[ ${fastqpresent} -eq 0 ]] && [[ ! -s ${DataFolder}/${sample}/${sample}.sra ]]
	then
		echo "Downloading..."
		jobidD=$(sbatch -t 10:00:00 --mem=4000 -J D${sample} -o ${errdir}/D${sample}.out -e ${errdir}/D${sample}.err ${PerSampleDownloadScript} ${DataFolder} ${sample} | rev | cut -f1 -d' ' | rev)
		echo "JobID: "${jobidD}" running"
	fi

	# SRA to fastq with fastq-dump
	if [[ ${fastqpresent} -eq 0 ]] && ([ -s ${DataFolder}/${sample}/${sample}.sra ] || [ $jobidD -gt 0 ])
	then
		if [[ $jobidD -eq 0 ]]
		then
			echo "SRA to FASTQ"
			jobidF=$(sbatch -t 10:00:00 --mem=4000 -J F${sample} -o ${errdir}/F${sample}.out -e ${errdir}/F${sample}.err ${PerSampleSraToFastqScript} ${DataFolder} ${sample} | rev | cut -f1 -d' ' | rev)
		else
			echo "SRA to FASTQ dependent on download JobID: "${jobidD}
			jobidF=$(sbatch -t 10:00:00 --mem=4000 -J F${sample} -d ${jobidD} -o ${errdir}/F${sample}.out -e ${errdir}/F${sample}.err ${PerSampleSraToFastqScript} ${DataFolder} ${sample} | rev | cut -f1 -d' ' | rev)
		fi
		echo "JobID: "${jobidF}" running"
	fi 

done < "$SRAIDs"
