# Gene duplication analysis in B. lanceolatum new reference genome (BraLan3) 
# #


### Folder structure:
- Data: Raw data
- Metadata: Metadata files
- Scripts: Scripts to process data (detailed below)
- Results: Results of the processement of the Data by Scripts
- ScriptsPlots: Scripts generating plots/tables 
- Plots: Generated from ScriptsPlots 
- tmp: err and out files & temporary files


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

##### Scripts:
- Scripts/AssemblyStatistics/ExtractAssemblyStatistics.sh
> Extracts basic statistics per assembly
- Scripts/AssemblyStatistics/Run_ExtractAssemblyStatistics.sh
> Runs te above script for all assemblies and joins the results 

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

### Gene expression analysis

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













