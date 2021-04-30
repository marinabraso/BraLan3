!#/bin/bash

GENOME=/scratch/beegfs/monthly/kjaron/amphioxus/data/genome/Blan_v05.fa
MAP=
READS=/scratch/beegfs/monthly/kjaron/amphioxus/data/raw_reads/reads_il/

for i in 3 5 8; do
	$MAP  BOSC_mp_is$i""k $GENOME $READS""ATQ_BOSC_IND$i""_1.fq.gz $READS""ATQ_BOSC_IND$i""_2.fq.gz
done;

for i in 2 5 6; do
        $MAP BOSN_mp_is$i""k $GENOME $READS""ATQ_BOSN_IND$i""_1.fq.gz $READS""ATQ_BOSN_IND$i""_2.fq.gz
done;

