
source('/Volumes/dump/scripts/generic_genomics/R_functions/getNX.R')

# read scaffold lengths of raw merged assembly

sl_mp <- read.csv('/Volumes/dump/projects/amphio/data/merging/mp_lengths.csv', header = F)$V1
sl_ref <- read.csv('/Volumes/dump/projects/amphio/data/merging/PBILref_lengths.csv', header = F)$V1
sl_mp_ref <- read.csv('/Volumes/dump/projects/amphio/data/merging/mp_ref_lengths.csv', header = F)$V1

getNX(sl_mp, 50, 1040000000)
getNX(sl_ref, 50, 520000000)
getNX(sl_mp_ref, 50, 520000000)
getNX(sl_mp_ref, 50)

sum(sl_ref)
sum(sl_mp_ref)
