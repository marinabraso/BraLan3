# Mediterranean amphioxus genome

This repository records progress on mediterranean amphioxus genome project. The genome is based on long read sequencing, but an assembly of the different individual was performed using Illumina, this assembly was used for comparison. The outline of project will be recorded in this document with links to individual steps.

## Data

83 PacBio SMRT cells - link to QC, stats, distribution

## Assembly

- Canu
- Falcon

## Downstream analyses

- Comparison to Illumina assembly
- HOX clusters analysis
- paraHOX cluster extraction

# Future directions

## close future

- find an edge on HOX cluster inversion (done)
- identify genes around the inversion (if it plays a role or not)
- scaffold PB assembly using mp (SOAPdenovo)
- scaffold IL assembly using PB reads (sspace)
- merge assemblies (HaploMerger - done; other tool?)

## second step

- estimate heterozygocity and structural variations in PB genome (atlas / siffels)
- estimate heterozygocity and structural variations in IL genome (atlas, ?)
- compare performance of technologies
- estimate heterozygocity and structural variations between IL and PB genomes 
- compare heterozygocity and structural variations between B. lanceloatum and B. floridae and Chinese amphioxus

## do annotation though Maker

provide all - RNAseq data, protein data (annotation from)

instead of annotation use transcript

use CEGMA (or BUSCO) gene models as well as annotation from Ferdi. gff or gtf (gene model file)

look EST for amphioxus
