#!/bin/bash

# Scripts
module add SequenceAnalysis/MultipleSequenceAlignment/T-Coffee/11.00.8cbe486
module add Blast/ncbi-blast/2.10.1+

alnFile=$1
t_coffee -infile=${alnFile}  -mode evaluate  -output=score_ascii
