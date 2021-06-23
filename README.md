# #
# #
# #
# BraLan3
# #
# #
# #

### Folder structure:
- Data > Raw data
- Metadata > Metadata files
- Scripts # Scripts to process data (detailed below)
- Results # Results of the processement of the Data by Scripts
- Analysis # Scripts generating plots/tables 
- Plots # Generated from Analysis 
- tmp # err and out files & temporary files


# #
###### ##################
## Assembly statistics
###### ##################
Extract basic assembly quality statistics from assemblies:
- Branchiostoma_lanceolatum.BraLan3
- Branchiostoma_lanceolatum.BraLan2
- Branchiostoma_belcheri.Haploidv18h27
- Branchiostoma_floridae.Bfl_VNyyK
- Danio_rerio.GRCz11
- Gallus_gallus.GRCg6a
- Homo_sapiens.GRCh38
- Mus_musculus.GRCm39
- Asterias_rubens.eAstRub1.3
- Saccoglossus_kowalevskii.Skow1.1
- Strongylocentrotus_purpuratus.Spur5.0

##### Scripts:
> Extracts basic statistics per assembly
- Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh
> Runs te above script for all assemblies and joins the results 
- Scripts/AssemblyStatistics/Run_ExtractAssemblyStatistics.sh

##### Usage:
```
sbatch -t 01:00:00 --mem=8000 -J AssemblyStatistics -o tmp/AssemblyStatistics.out -e tmp/AssemblyStatistics.err Scripts/AssemblyStatistics/Run_ExtractAssemblyStatistics.sh
```


# #
###### ##################
## Gene expression
###### ##################
Gene expression analysis of Marletaz et al. 2018 RNA seq data 

### Download SRA data

##### Scripts:
- Scripts/DownloadSRAData/DownloadSRAData.sh
- Scripts/DownloadSRAData/DownloadSRAData_perSample.sh
- Scripts/DownloadSRAData/SRA_to_fastq_perSample.sh

##### Usage:
```
./Scripts/DownloadSRAData/DownloadSRAData.sh Data/RNA-seq/Marletaz2018 <(cat Metadata/Marletaz2018_RNAseq_SRA.list)
```

### Gene expression

##### Scripts:
> Prepare GTF and kallisto index
- Scripts/GeneExpression/PrepareFiles.sh
> Run FASTQC and kallisto per sample
- Scripts/GeneExpression/Geneexpression_fastqc_kallisto.sh
- Scripts/GeneExpression/FASTQC_PerSample.sh
- Scripts/GeneExpression/Kallisto_PerSample.sh
> Join results by gene with Tximport
- Scripts/GeneExpression/JoinResultsByGene_Tximport.sh
- Scripts/GeneExpression/JoinResultsByGene_Tximport.R

##### Usage:
```
sbatch -t 01:00:00 --mem=8000 -J PrepareKallisto -o tmp/PrepareKallisto.out -e tmp/PrepareKallisto.err Scripts/GeneExpression/PrepareFiles.sh

./Scripts/GeneExpression/Geneexpression_fastqc_kallisto.sh

./Scripts/GeneExpression/JoinResultsByGene_Tximport.sh
```


# #
###### ##################
## Orthologs and paralogs analysis
###### ##################

### Filtering gene sets
Prepare gene sets for each species. 

##### Scripts:
- Scripts/FilteringGeneSets/FilteringSetGenesPerSpecies.sh
- Scripts/FilteringGeneSets/ExtractDNAseqFromGenome.pl # To extract cDNA sequences 

##### Usage:
```
sbatch -t 10:00:00 --mem=8000 -J FiltGenes -o tmp/FilteringGeneSets.out -e tmp/FilteringGeneSets.err Scripts/FilteringGeneSets/FilteringSetGenesPerSpecies.sh
```

### Find orthologs

##### Scripts:
- Scripts/FindOrthologs/FindOrthologs_brocoli.sh

##### Usage:
```
sbatch -t 10:00:00 --mem=8000 -J BAmpVertOut -o tmp/Broccoli_AmphVertebOutDeut.out -e tmp/Broccoli_AmphVertebOutDeut.err conda run -n broccoli Scripts/FindOrthologs/FindOrthologs_brocoli.sh AmphVertebOutDeut
sbatch -t 03:00:00 --mem=8000 -J BAmpVert -o tmp/Broccoli_AmphVerteb.out -e tmp/Broccoli_AmphVerteb.err conda run -n broccoli Scripts/FindOrthologs/FindOrthologs_brocoli.sh AmphVerteb
```

### dNdS between paralogs

##### Scripts:
> Extract orthologous group sequences (AA and DNA) for Amphioxus and vertebrae species
- Scripts/dNdS/ExtractOrthologousGroupSequences.sh
> MSA of the AA sequences with mafft + backtranslation to DNA (from the DNA sequences)
- Scripts/dNdS/MSA_AA_backtranslation_DNA.sh
- Scripts/dNdS/Run_MSA_AA_backtranslation_DNA.sh
> evaluate the AA MSA done with mafft t_coffe, extract score and filter mafft alignments with t_coffee score
- Scripts/dNdS/MSA_cleaning_tcoffe.sh
- Scripts/dNdS/Run_MSA_cleaning_tcoffe_PerOG.sh
> Tree from alignment with RAxML and dNdS analysis with Godon (model M8)
- Scripts/dNdS/Calculate_dDdSBetweenParalogs_Godon.sh
- Scripts/dNdS/Run_Calculate_dDdSBetweenParalogs_Godon.sh

##### Usage:
```
sbatch -t 05:00:00 --mem=8000 -J ExtrAmpVert -o tmp/ExtractSeq_AmphVerteb.out -e tmp/ExtractSeq_AmphVerteb.err Scripts/dNdS/ExtractOrthologousGroupSequences.sh AmphVerteb

./Scripts/dNdS/Run_MSA_AA_backtranslation_DNA.sh AmphVerteb

./Scripts/dNdS/Run_MSA_cleaning_tcoffe.sh AmphVerteb

./Scripts/dNdS/Run_Calculate_dDdSBetweenParalogs_Godon.sh AmphVerteb
```


```
for f in $(ls GO*); do awk -v f=$f 'BEGIN{split(f,a,"_")}{print a[1]"\t"a[2]"\t"$0}' $f;done| sed 's/.txt//g' | sed 's/^GO/GO:/g' > Human_GenesMainGOtermsMF.txt

awk '{if(NR==FNR){a[$2]=$1;next} if(a[$4]){print a[$4]"\t"$0} }' Ensembl2Uniprot_Human_GRCh38.txt Human_GenesMainGOtermsMF.txt | cut -f1,2 | sort -u > Human_GenesMainGOtermsMF_ENS.txt

tail -n +2 Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt | awk -F '\t' '{split($2, l, " ");split($5, h, " "); for(i in l){for(j in h){ print l[i]"\t"h[j]}}}' | sort -u > Results/FindOrthologs/AmphVerteb_broccoli/Blan2Hsap_genes.txt

```
















