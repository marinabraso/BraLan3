## Annotation of merged assembly

For synteny analysis we need to use annotation done on Bl71. To map `.gtf` file to the new assembly we will use tool `exonerate`.

Scaffolds to annotate: `data/asm_merged/amphio_PB_scf.fa`
Protein sequences file: `data/Bl71/Bla_annot_final_refProteins.fa`

```
module add SequenceAnalysis/SequenceAlignment/exonerate/2.2.0

cd data/annotation
sed -i 's/\.//g' ../Bl71/Bla_annot_final_refProteins.fa
ln -s ../Bl71/Bla_annot_final_refProteins.fa .
ln -s ../asm_merged/amphio_PB_scf.fa .
exonerate --model protein2genome --target amphio_PB_scf.fa --query Bla_annot_final_refProteins.fa --showcigar no --showvulgar no --bestn 1 --showquerygff no --showtargetgff yes > amphio_exonerate.out
```

there is a problem with big file. Divide & conquer:

For this I use script [1_exonerate](1_exonerate), that computes an alignment of 1/50 of proteins. To save a lot of I/O operations I have decided to copy files locally on disk of one host (cpt171), therefore it is submitting to one host only. So far I submited 25 jobs, let's see if they manage, computing time should be about 12 hours.

The output of exonerate is for sure blind to transctipt/gene/exon names and other information present in `gtf` of the previous assembly, therefore I need to create a script refactoring information of exonerate output.

1. sort both gtf files by query name
2. iterate trough both:
  -  if they match: print refactored gtf lines; increment both indices
  -  if not: increment index that comes first in alphabet

To compare features of the two annotation files I created a sample of the same transcript `BL00003_evm0`. Seems that PacBio assembly has very small fration of this gene. No sure why.

At some point I mingh want to do a quantitative comparison of these two annotations.
