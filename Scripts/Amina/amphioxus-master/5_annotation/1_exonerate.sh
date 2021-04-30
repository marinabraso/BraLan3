#!/bin/bash
# 1. argument is order of the chunk to be computed

if [[ -z $1 ]]
then
        echo "chunk was not specified (1-50)"
        exit 1
fi

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J ec_$1
#BSUB -q bgee
#BSUB -o ex_chunk$1.out
#BSUB -e ex_chunk$1.err
#BSUB -n 1
#BSUB -M 8388608
#BSUB -m cpt171

mkdir -p /scratch/local/amphio/temp$1
export TMPDIR=/scratch/local/amphio/temp$1
cd /scratch/local/amphio

time exonerate --model protein2genome --target /scratch/local/amphio/amphio_PB_scf.fa \
  --query /scratch/local/amphio/Bla_annot_final_refProteins.fa --showalignment no \
  --showcigar no --showvulgar no --bestn 1 --showquerygff no --showtargetgff yes \
  --targetchunkid $1 --targetchunktotal 50 1> amphio_exonerate_chunk$1.out

mv amphio_exonerate_chunk$1.out /scratch/beegfs/monthly/kjaron/amphioxus/data/annotation/ex_chunks
"""
