#BSUB -L /bin/bash
#BSUB -J amphio_scf
#BSUB -q long
#BSUB -o amphio_BESST.out
#BSUB -e amphio_BESST.err
#BSUB -n 1
#BSUB -M 53554432
#BSUB -R \"rusage[tmp=100000] span[ptile=1]\"

#module add UHTS/Aligner/bwa/0.7.13
#module add UHTS/Analysis/samtools/1.3

WORKING_DIR=`pwd`

# make local directory for computations
mkdir -p /scratch/local/monthly/kjaron/BESST
export TMPDIR=/scratch/local/monthly/kjaron/BESST
cd /scratch/local/monthly/kjaron/BESST

# copy required files
cp /scratch/beegfs/monthly/kjaron/amphioxus/data/asm_pb/amphio_polished.fa .
ln -s /scratch/beegfs/monthly/kjaron/amphioxus/data/mapping_mp2pb_asm/BO* $TMPDIR

# run BESST script
python ~/src/BESST/runBESST -c amphio_polished.fa -f BOSN_mp_is2k.bam BOSC_mp_is3k.bam BOSC_mp_is5k.bam BOSN_mp_is5k.bam BOSN_mp_is6k.bam BOSC_mp_is8k.bam -o $1 --no-score
