#!/bin/sh

#SBATCH --account mrobinso_default
#SBATCH --mail-user amina.echchiki@unil.ch
#SBATCH --mail-type ALL
#SBATCH --partition ax-long
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 2GB
#SBATCH --job-name amphio_getRNA
#SBATCH --export NONE

## modules
module load Bioinformatics/Software/vital-it
module add UHTS/Analysis/sratoolkit/2.9.6.1

## dir
# cd /scratch/axiom/FAC/FBM/DEE/mrobinso/default/aechchik/amphioxus

## get fastq from GSM

# get GSM from GSE
# pip install --user --ignore-installed  pysradb
pysradb gse-to-gsm GSE106430 > gse-to-gsm_tmp
# fd up: select ^ GSE
cat gse-to-gsm_tmp | grep '^ GSE' > gse-to-gsm
cat gse-to-gsm | sed 's/ \+ /\t/g' | cut -f2 > list_GSM

# get SRR from GSM
for i in $(cat list_GSM); do
	pysradb gsm-to-srr $i > $i"_tmp"
done
cat *_tmp | grep ^" " |  sed 's/ \+ /\t/g' | cut -f2 > list_SRR
rm *_tmp

# get sra from SRR
for i in $(cat list_SRR); do
	prefetch --output-directory . -v $i
done

# get fastq from sra
for i in $(ls *.sra); do
	fastq-dump $i --outdir . --split-files --gzip
done

# alternative, by kamil

# for i in $(cat list_SRR); do
# 	R1_file="$i"_1.fastq.gz
# 	R2_file="$i"_2.fastq.gz
# 	URL_R1=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${i::6}/"$i"/"$i"_1.fastq.gz
# 	URL_R2=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${i::6}/"$i"/"$i"_2.fastq.gz
# 	wget $URL_R1 -O $R1_file
# 	wget $URL_R2 -O $R2_file
# done


# got for now a tmp dataset, wait for the full one


# No such directory ‘vol1/ERA645/ERA645809/fastq’.
# verified on the website, no directory :/

# SRR6245991
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR624/001/SRR6246001/

# SRR6245999
# ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR624/009/SRR6245999/

#module add UHTS/Quality_control/fastqc/0.11.7;
#module add UHTS/Aligner/bwa/0.7.17;
#module add UHTS/Analysis/samtools/1.8;
#module add UHTS/Analysis/picard-tools/2.18.11;
#module add UHTS/Analysis/GenomeAnalysisTK/4.1.3.0;
#module add UHTS/Analysis/vcftools/0.1.15;
#module add R;



