#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
#library(viridis)

######################################################################
# Files & folders

ResultsFolder <- "Plots/Synteny"
CountsFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
OhnologsSFile <- "Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt"
GTFfolder <- "Results/FilteringGeneSets/GTFs"

######################################################################
# General parameters
Species1 <- "Blan"
Species2 <- "Mmus"

ShortSpeciesNames <- c("Drer", "Ggal", "Hsap", "Mmus", "Blan")
SpeciesNames <- c("Danio_rerio", "Gallus_gallus", "Homo_sapiens", "Mus_musculus", "Branchiostoma_lanceolatum")
SpeciesGTFs <- c("Danio_rerio.GRCz11.gtf.gz", "Gallus_gallus.GRCg6a.gtf.gz", "Homo_sapiens.GRCh38.gtf.gz", "Mus_musculus.GRCm39.gtf.gz", "Branchiostoma_lanceolatum.BraLan3.gtf.gz")

SName1 <- SpeciesNames[which(ShortSpeciesNames==Species1)]
SName2 <- SpeciesNames[which(ShortSpeciesNames==Species2)]

###########################################################################
# Read data
# List of onhologs
Ohnologs <- read.table(OhnologsSFile, h=F)[,1]

# Orthologous groups counts (amphioxus + vertebrates)
CountsAV <- read.table(CountsFileAV, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFileAV, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]*_\\([a-z][a-z][a-z]\\)[A-Za-z\\.0-9_]*/\\1\\2/g'"), intern=T)
Header <- read.table(text=system_out, h=F, sep = "\t")
colnames(CountsAV) <- as.character(unlist(lapply(Header[1,], as.character)))
print(head(CountsAV))

# Read gene location information Species 1
system_out <- system(paste0("zcat ", GTFfolder, "/", SpeciesGTFs[which(ShortSpeciesNames==Species1)], " | sed 's/gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1/g' | cut -f1,4,5,9 | sort -k1,1 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2"), intern=T)
GeneCoord1 <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneCoord1) <- c("Gene1","Chr","Start","End")
print(head(GeneCoord1))

# Read gene location information Species 1
system_out <- system(paste0("zcat ", GTFfolder, "/", SpeciesGTFs[which(ShortSpeciesNames==Species2)], " | sed 's/gene_id \"\\([A-Z0-9]\\+\\)\".*/\\1/g' | cut -f1,4,5,9 | sort -k1,1 -k2,2V | awk '{if($4 == gene){end=$3}else{print gene\"\\t\"chr\"\\t\"start\"\\t\"end; gene=$4; chr =$1; start=$2; end=$3}}END{print gene\"\\t\"chr\"\\t\"start\"\\t\"end;}' | tail -n +2"), intern=T)
GeneCoord2 <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneCoord2) <- c("Gene2","Chr","Start","End")
print(head(GeneCoord2))

# OG to Sp1 gene
system_out <- system(paste("tail -n +2 ", NamesFileAV, " | awk -F '\\t' '{split($", which(colnames(CountsAV)==Species1)+1, ",a,\" \"); for(i in a){print $1\"\\t\"a[i]}}'"), intern=T)
OG2Sp1Gene <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2Sp1Gene) <- c("OG", "Gene1")
print(head(OG2Sp1Gene))

# OG to Sp2 gene
system_out <- system(paste("tail -n +2 ", NamesFileAV, " | awk -F '\\t' '{split($", which(colnames(CountsAV)==Species2)+1, ",a,\" \"); for(i in a){print $1\"\\t\"a[i]}}'"), intern=T)
OG2Sp2Gene <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2Sp2Gene) <- c("OG", "Gene2")
print(head(OG2Sp2Gene))

## To finish!
GeneData1 <- merge(OG2Sp1Gene, GeneCoord1, by="Gene1")
CountsAV$NumChrs1 <- unlist(lapply(rownames(CountsAV), function(x){length(unique(GeneData1$Chr[which(GeneData1$OG==x)]))}))
CountsAV$MaxDist1 <- unlist(lapply(rownames(CountsAV), CalcMaxDist(GeneData1)))


OGData.AV$DupType <- rep(NA, length(OGData.AV[,1]))
OGData.AV$DupType[which(OGData.AV$NumChrs==1 & OGData.AV$BlanType=="Duplicated")] <- rep("Intra",sum(OGData.AV$NumChrs==1 & OGData.AV$BlanType=="Duplicated"))
OGData.AV$DupType[which(OGData.AV$NumChrs>1 & OGData.AV$BlanType=="Duplicated")] <- rep("Inter",sum(OGData.AV$NumChrs>1 & OGData.AV$BlanType=="Duplicated"))
OGData.AV$DupType[which(OGData.AV$NumChrs==1 & OGData.AV$MaxDist<=TandemMaxDist & OGData.AV$BlanType=="Duplicated")] <- rep("Tandem",sum(OGData.AV$NumChrs==1 & OGData.AV$MaxDist<=TandemMaxDist & OGData.AV$BlanType=="Duplicated"))
table(OGData.AV$DupType)
table(OGData.AV$BlanType)
head(OGData.AV)






