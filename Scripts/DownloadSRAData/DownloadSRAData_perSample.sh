#!/bin/bash


outfold=$1
id=$2

echo "/repository/user/main/public/root = \"/scratch/wally/FAC/FBM/DEE/mrobinso/default/mbrasovi/tmp\"" > $HOME/.ncbi/user-settings.mkfg
rm /scratch/wally/FAC/FBM/DEE/mrobinso/default/mbrasovi/tmp/sra/${id}* 2> ~/null

cd ${outfold}
prefetch -X 50G ${id}

rm /scratch/wally/FAC/FBM/DEE/mrobinso/default/mbrasovi/tmp/sra/${id}* 2> ~/null
