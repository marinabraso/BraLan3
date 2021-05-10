#!/usr/bin/env Rscript



######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))

######################################################################
# Files & folders

ResultsFolder <- "Plots/OrthologsAnalysis"
CountsFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFileAV <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
CountsFileAVO <- "Results/FindOrthologs/AmphVertebOutDeut_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFileAVO <- "Results/FindOrthologs/AmphVertebOutDeut_broccoli/dir_step3/table_OGs_protein_names.txt"
OnhSFile <- "Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt"
RNASeqMetadataFile <- "Metadata/Marletaz2018_RNAseq_SRA.txt"
TPMFile <- "Results/GeneExpression/Gene_TPM_kallisto_tximport.tab"
dNdSFolder <- "Results/dNdSBetweenParalogs/Godon_M8"

######################################################################
# General parameters

Verteb <- c("Drer", "Ggal", "Mmus", "Hsap")
Amphi <- c("Blan", "Bflo", "Bbel")
OutDeut <- c("Spur", "Arub", "Skow")
SortedSpecies <- c(Amphi, Verteb, OutDeut)
colfunc <- colorRampPalette(c("darkred", "firebrick1"))
CNcolors <- c("grey80", "gold", colfunc(3))
TypeColors <- c("deepskyblue2", "gold", "olivedrab2", "olivedrab4")
TissueColors <- c()

###########################################################################
# Read data

# List of onhologs
OnhOG.S <- read.table(OnhSFile, h=F)[,1]

# Orthologous groups counts
CountsAV <- read.table(CountsFileAV, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFileAV, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
Header <- read.table(text=system_out, h=F, sep = "\t")
colnames(CountsAV) <- as.character(unlist(lapply(Header[1,], as.character)))
CountsAV$Sum <- rowSums(CountsAV)
CountsAV$SumVerteb <- rowSums(CountsAV[,Verteb])
CountsAV$SumAmphi <- rowSums(CountsAV[,Amphi])
CountsAV$VertebType <- rep("Missing", length(CountsAV[,1]))
CountsAV$VertebType[which(CountsAV$SumVerteb>0)] <- rep("SingleCopy", length(CountsAV$VertebType[which(CountsAV$SumVerteb>0)]))
CountsAV$VertebType[which(apply(CountsAV[,Verteb], 1, max)>1)] <- rep("Duplicated", length(CountsAV$VertebType[which(apply(CountsAV[,Verteb], 1, max)>1)]))
CountsAV$VertebType[which(rownames(CountsAV) %in% OnhOG.S)] <- rep("Ohnolog", length(CountsAV$VertebType[which(rownames(CountsAV) %in% OnhOG.S)]))
CountsAV$AmphiType <- rep("Missing", length(CountsAV[,1]))
CountsAV$AmphiType[which(CountsAV$SumAmphi>0)] <- rep("SingleCopy", length(CountsAV$AmphiType[which(CountsAV$SumAmphi>0)]))
CountsAV$AmphiType[which(apply(CountsAV[,Amphi], 1, max)>1)] <- rep("Duplicated", length(CountsAV$Amp[which(apply(CountsAV[,Amphi], 1, max)>1)]))
CountsAV$BlanType <- rep("Missing", length(CountsAV[,1]))
CountsAV$BlanType[which(CountsAV$Blan>0)] <- rep("SingleCopy", length(CountsAV$BlanType[which(CountsAV$Blan>0)]))
CountsAV$BlanType[which(CountsAV$Blan>1)] <- rep("Duplicated", length(CountsAV$BlanType[which(CountsAV$Blan>1)]))



# dNdS data
system_out <- system(paste("grep '^Final' ", dNdSFolder, "/OG_*_M8_likelihood_ratio.txt | sed 's/Results\\/dNdSBetweenParalogs\\/Godon_M8\\///g' | sed 's/_M8_likelihood_ratio.txt:Final D=/\t/g'", sep=""), intern=T)
GodonM8.D <- read.table(text=system_out, h=F, sep = "\t", row.names=1)
colnames(GodonM8.D) <- c("D")
GodonM8.D$Pval <- unlist(lapply(GodonM8.D$D, p.value.fromD))
GodonM8.D$Padj <- p.adjust(GodonM8.D$Pval, method = "bonferroni")
CountsAV$M8.D <- rep(NA, length(CountsAV[,1]))
CountsAV$M8.Pval <- rep(NA, length(CountsAV[,1]))
CountsAV$M8.Padj <- rep(NA, length(CountsAV[,1]))
CountsAV$M8.D[match(rownames(GodonM8.D), rownames(CountsAV))] <- GodonM8.D$D
CountsAV$M8.Pval[match(rownames(GodonM8.D), rownames(CountsAV))] <- GodonM8.D$Pval
CountsAV$M8.Padj[match(rownames(GodonM8.D), rownames(CountsAV))] <- GodonM8.D$Padj




# RNA-seq Metadata
system_out <- system(paste("cat ", RNASeqMetadataFile," | cut -f1,13,24", sep=""), intern=T)
MetaRNA <- read.table(text=system_out, h=F, sep = "\t")
colnames(MetaRNA) <- c("Sample","Age","Tissue")

# RNA-seq gene expression TPM data
TPM <- read.table(TPMFile, h=T)
TPM <- TPM[,MetaRNA$Sample[which(MetaRNA$Tissue!="egg" & MetaRNA$Sample%in%colnames(TPM))]]
MetaRNA <- MetaRNA[which(MetaRNA$Sample %in% colnames(TPM)),]

Tissues <- unique(MetaRNA$Tissue[which(MetaRNA$Age=="adult")])
Tissues <- c("cirri", "gills", "epidermis", "gut", "hepatic diverticulum", "muscle", "neural tube", "female gonads", "male gonads")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
TissueColors <- colfunc(length(Tissues))
EmbAges <- unique(MetaRNA$Age[which(MetaRNA$Age!="adult")])
EmbAges <- c("egg", "32cells", "Blastula", "7h", "8h", "10h", "11h", "15h", "18h", "21h", "24h", "27h", "36h", "50h", "60h", "Pre-metamorphic larvae")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
EmbAgesColors <- colfunc(length(EmbAges))

###########################################################################
# Data Processing

# TAU calculation per adult tissue
ExpData <- data.frame(row.names = rownames(TPM))
for(t in c(1:length(Tissues))){
	samples <- MetaRNA$Sample[which(MetaRNA$Tissue==Tissues[t])]
	if(length(samples) == 1){
		ExpData <- cbind(ExpData, as.numeric(TPM[, samples]))
	}else{
		ExpData <- cbind(ExpData, as.numeric(rowMeans(TPM[, samples])))	
	}
}
for(a in c(1:length(EmbAges))){
	samples <- MetaRNA$Sample[which(MetaRNA$Age==EmbAges[a])]
	if(length(samples) == 1){
		ExpData <- cbind(ExpData, as.numeric(TPM[, samples]))
	}else{
		ExpData <- cbind(ExpData, as.numeric(rowMeans(TPM[, samples])))	
	}
}
colnames(ExpData) <- c(Tissues, EmbAges)
m.ExpData <- as.matrix(ExpData[,Tissues])
Tau <- apply(log2(m.ExpData+1), 1, compute.tau)
MaxGroup <- apply(log2(m.ExpData+1), 1, max.col)
ExpData$TauTissues <- Tau
ExpData$MaxTissue <- MaxGroup
m.ExpData <- as.matrix(ExpData[,EmbAges])
Tau <- apply(log2(m.ExpData+1), 1, compute.tau)
MaxGroup <- apply(log2(m.ExpData+1), 1, max.col)
ExpData$TauEmbAge <- Tau
ExpData$MaxEmbAge <- MaxGroup

ExpData$MeanAdult <- rowMeans(TPM[rownames(ExpData),MetaRNA$Sample[which(MetaRNA$Age=="adult")]])
ExpData$MaxAdult <- apply(TPM[rownames(ExpData),MetaRNA$Sample[which(MetaRNA$Age=="adult")]],1,max)
ExpData$MeanEmbr <- rowMeans(TPM[rownames(ExpData),MetaRNA$Sample[which(MetaRNA$Age!="adult")]])
ExpData$MaxEmbr <- apply(TPM[rownames(ExpData),MetaRNA$Sample[which(MetaRNA$Age!="adult")]],1,max)


system_out <- system(paste("cat ", NamesFileAV, " | tail -n +2 | awk '{for(i=2;i<=NF;i++){print $1\"\t\"$i}}' | grep 'BLAG' | sort -u | awk '{if(a[$2]){a[$2]=\"NA\"}else{a[$2]=$1}}END{for(i in a){print a[i]\"\t\"i}}'", sep=""), intern=T)
OG2gene <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2gene) <- c("OG", "Gene")
ExpData$OG <- rep(NA, length(ExpData[,1]))
ExpData$OG[match(OG2gene$Gene, rownames(ExpData))] <- OG2gene$OG
ExpData <- ExpData[which(!is.na(ExpData$OG)),]
ExpData$BlanType <- rep(NA, length(ExpData[,1]))
ExpData$BlanType <- unlist(lapply(c(1:length(ExpData[,1])), function(x){CountsAV[which(rownames(CountsAV)==ExpData$OG[x]),"BlanType"]}))
ExpData$VertebType <- rep(NA, length(ExpData[,1]))
ExpData$VertebType <- unlist(lapply(c(1:length(ExpData[,1])), function(x){CountsAV[which(rownames(CountsAV)==ExpData$OG[x]),"VertebType"]}))
ExpData$M8.D <- rep(NA, length(ExpData[,1]))
ExpData$M8.D <- unlist(lapply(c(1:length(ExpData[,1])), function(x){CountsAV[which(rownames(CountsAV)==ExpData$OG[x]),"M8.D"]}))
ExpData$M8.Padj <- rep(NA, length(ExpData[,1]))
ExpData$M8.Padj <- unlist(lapply(c(1:length(ExpData[,1])), function(x){CountsAV[which(rownames(CountsAV)==ExpData$OG[x]),"M8.Padj"]}))
print(head(ExpData))



######################################################################
######################################################################
# Write files
cat("Write files\n")

write.table(ExpData, file = paste(ResultsFolder, "/ExpressionDataProcessed.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(CountsAV, file = paste(ResultsFolder, "/OGCountsProcessed.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)
write.table(MetaRNA, file = paste(ResultsFolder, "/RNASeqMetadataProcessed.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)





