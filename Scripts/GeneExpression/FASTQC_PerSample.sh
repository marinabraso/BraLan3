#!/bin/bash

module add UHTS/Quality_control/fastqc/0.11.7;

fastq=$1
outdir=$2

fastqc $fastq -o $outdir



