#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.r", sep=""))
#source("ScriptsPlots/DuplicatesOntologyExpression/DuplicatesOntologyExpression_functions.r")
# source("ScriptsPlots/DuplicatesOntologyExpression/DuplicatesOntologyExpression.r")
library(ghibli)
library(MASS)
library(moments)

######################################################################
# Files & folders

ResultsFolder <- "Plots/DuplicatesOntologyExpression"
CountsFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
Ohnolog2ROG <- "Results/OhnologListing/2R_Strict_OG_pairs.txt"
Ohnolog3ROG <- "Results/OhnologListing/3R_Strict_OG_pairs.txt"

GOlistsFolder <- "Results/GeneOntology"
GOlistFile <- "Results/GeneOntology/GOlist.tbl"

BlanGeneDataFile <- paste0(ResultsFolder, "/BlanGeneData.txt")
DrerGeneDataFile <- paste0(ResultsFolder, "/DrerGeneData.txt")
DrerBgeeDataFile <- paste0(ResultsFolder, "/DrerBgeeData.txt")
RNASeqMetadataFile <- "Metadata/Marletaz2018_RNAseq_SRA.txt"
TPMFile <- "Results/GeneExpression/Gene_TPM_kallisto_tximport.tab"

######################################################################
# General parameters
PQvalThreshold <- 0.01 # for hypergeometric test
minGOsize <- 50
MaxTandemDist <- 10000
TPM.threshold <- 1

# Species
Species <- c("Blan", "Bflo", "Bbel", "Drer", "Ggal", "Mmus", "Hsap")
SpType <- c("BlanType", "AmphType", "AmphType", "VertType", "VertType", "VertType", "VertType")
Ampioxus <- c("Blan", "Bflo", "Bbel")
Vertebrates <- c("Drer", "Ggal", "Mmus", "Hsap")
SpeciesLongNames <- c("Branchiostoma_lanceolatum.BraLan3", "Branchiostoma_floridae.Bfl_VNyyK", "Branchiostoma_belcheri.Haploidv18h27", "Danio_rerio.GRCz11", "Gallus_gallus.GRCg6a", "Mus_musculus.GRCm39", "Homo_sapiens.GRCh38")
SpType <- c("AmphType", "AmphType", "AmphType", "VertType", "VertType", "VertType", "VertType")

VertTypes <- c("Single-copy", "Ohnologs", "Small-scale\nduplicates", "Missing")
VertTypes.abr <- c("SC", "O", "D", "M")
VertTypes.col <- c(ghibli_palettes$MarnieMedium2[c(3,6,7)], "royalblue2")
BlanTypes <- c("Single-copy", "Small-scale\nduplicates", "Missing")
BlanTypes.pch <- c(14, 16, 18)
BlanTypes.lty <- c(.5, 1, 2)

# Gene expression analysis
Species1 <- "Blan"
Species2 <- "Drer"

BlanMatchingTissues <- c("Blastula", "MeanEmbr", "male.gonads", "female.gonads", "muscle", "neural.tube", "gut", "gills", "hepatic.diverticulum", "epidermis")
DrerMatchingTissues <- c("blastula", "embryo", "testis", "ovary", "muscle tissue", "brain", "intestine", "pharyngeal gill", "liver", "zone of skin")
NamesMatchingTissues <- c("blastula", "embryo", "testis", "ovary", "muscle", "neural", "digestive", "gills", "hepatic", "skin")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
TissueColors <- colfunc(length(NamesMatchingTissues))
MatchingTissues <- as.data.frame(cbind(BlanMatchingTissues, DrerMatchingTissues, NamesMatchingTissues, TissueColors))
colnames(MatchingTissues) <- c("Blan", "Drer", "Name","Color")
MatchingTissues <- MatchingTissues[MatchingTissues$Name!="blastula" & MatchingTissues$Name!="ovary" & MatchingTissues$Name!="skin",]

###########################################################################
###########################################################################
# Read data

OGInfoFile <- paste(ResultsFolder, "/OGInfo.txt", sep ="")
GeneInfoFile <- paste(ResultsFolder, "/GeneInfo.txt", sep ="")
NumGenesSpeciesFile <- paste(ResultsFolder, "/NumGenesSpecies.txt", sep ="")
GenePairsFile <- paste(ResultsFolder, "/GenePairs.txt", sep ="")
GOInfoFile <- paste(ResultsFolder, "/GOInfo.txt", sep ="")

OGInfo <- read.table(OGInfoFile, h=TRUE, sep = "\t", row.names=1)
GeneInfo <- read.table(GeneInfoFile, h=TRUE, sep = "\t", row.names=1)
NumGenesSpecies <- read.table(NumGenesSpeciesFile, h=TRUE, sep = "\t", row.names=1)[,1]
GenePairs <- read.table(GenePairsFile, h=TRUE, sep = "\t", row.names=1)
GOInfo <- read.table(GOInfoFile, h=TRUE, sep = "\t", row.names=1)

DrerGeneData <- prepare_DrerGeneData(DrerGeneDataFile, DrerBgeeDataFile, MatchingTissues, ResultsFolder)
BlanGeneData <- prepare_BlanGeneData(BlanGeneDataFile, RNASeqMetadataFile, TPMFile)

OGInfo$BlanType <- gsub(" ", "\n", OGInfo$BlanType)
OGInfo$VertType <- gsub(" ", "\n", OGInfo$VertType)
OGInfo$AmphType <- gsub(" ", "\n", OGInfo$AmphType)
GeneInfo$BlanType <- gsub(" ", "\n", GeneInfo$BlanType)
GeneInfo$VertType <- gsub(" ", "\n", GeneInfo$VertType)
GeneInfo$AmphType <- gsub(" ", "\n", GeneInfo$AmphType)
GenePairs$BlanType <- gsub(" ", "\n", GenePairs$BlanType)
GenePairs$VertType <- gsub(" ", "\n", GenePairs$VertType)

# Subsets
OGInfo.BlanVert <- OGInfo[which(OGInfo$BlanType != "Missing" & OGInfo$VertType != "Missing"),]
GOInfo.MF <- GOInfo[which(GOInfo$Type=="molecular_function"),]
GOInfo.BP <- GOInfo[which(GOInfo$Type=="biological_process"),]


###########################################################################
###########################################################################
#### Printing statistics

# Amphioxus - vertebrate shared OG numbers
print(paste("Shared V-A", sum(OGInfo$VertType!="Missing" & OGInfo$AmphType!="Missing")))
print(paste("V", sum(OGInfo$VertType!="Missing")))
print(paste("A", sum(OGInfo$AmphType!="Missing")))
print(paste("Shared/V*100", sum(OGInfo$VertType!="Missing" & OGInfo$AmphType!="Missing")/sum(OGInfo$VertType!="Missing")*100))
print(paste("Shared/A*100", sum(OGInfo$VertType!="Missing" & OGInfo$AmphType!="Missing")/sum(OGInfo$AmphType!="Missing")*100))

# Amphioxus specfic OG numbers
OGInfo.A <- OGInfo[which(OGInfo$VertType=="Missing" & OGInfo$AmphType!="Missing"),]
print(paste("A specific", length(OGInfo.A[,1])))
print(paste("Shared all A", sum(unlist(apply(OGInfo.A[,Ampioxus], 1, min))>0)))
print(paste("Shared all A/A specific*100", sum(unlist(apply(OGInfo.A[,Ampioxus], 1, min))>0)/length(OGInfo.A[,1])*100))
print(paste("blan-bbel specific", sum(OGInfo.A$Blan>0 & OGInfo.A$Bflo==0 & OGInfo.A$Bbel>0)))
print(paste("blan-bbel specific/A specific*100", sum(OGInfo.A$Blan>0 & OGInfo.A$Bflo==0 & OGInfo.A$Bbel>0)/length(OGInfo.A[,1])*100))
print(paste("blan-bflo specific", sum(OGInfo.A$Blan>0 & OGInfo.A$Bflo>0 & OGInfo.A$Bbel==0)))
print(paste("blan-bflo specific/A specific*100", sum(OGInfo.A$Blan>0 & OGInfo.A$Bflo>0 & OGInfo.A$Bbel==0)/length(OGInfo.A[,1])*100))
print(paste("bflo-bbel specific", sum(OGInfo.A$Blan==0 & OGInfo.A$Bflo>0 & OGInfo.A$Bbel>0)))
print(paste("bflo-bbel specific/A specific*100", sum(OGInfo.A$Blan==0 & OGInfo.A$Bflo>0 & OGInfo.A$Bbel>0)/length(OGInfo.A[,1])*100))
BlanExist <- c("Blan0", "Blan1")[(OGInfo.A$Blan > 0)+1]
BfloExist <- c("Bflo0", "Bflo1")[(OGInfo.A$Bflo > 0)+1]
BbelExist <- c("Bbel0", "Bbel1")[(OGInfo.A$Bbel > 0)+1]
table(c(BlanExist, BfloExist, BbelExist))

SpeciesNumGenesOG <- matrix(rep(NA, 10*length(Species)),nrow=length(Species),ncol=10,byrow=T)
for(sp in c(1:length(Species))){
	GeneInfo.sp <- GeneInfo[which(GeneInfo$Species==Species[sp]),]
	SpeciesNumGenesOG[sp,1] <- Species[sp]
	SpeciesNumGenesOG[sp,2] <- sum(OGInfo[,Species[sp]]>0)
	SpeciesNumGenesOG[sp,3] <- sum(OGInfo[,Species[sp]]>1)
	SpeciesNumGenesOG[sp,4] <- NumGenesSpecies[sp]
	SpeciesNumGenesOG[sp,5] <- length(unique(GeneInfo.sp$Gene[which(GeneInfo.sp[,Species[sp]]>0)]))
	SpeciesNumGenesOG[sp,6] <- length(unique(GeneInfo.sp$Gene[which(GeneInfo.sp[,Species[sp]]>1)]))
	summary <- summary(OGInfo[which(OGInfo[,Species[sp]]>1),Species[sp]])
	SpeciesNumGenesOG[sp,7] <- summary[4] # mean num genes per duplicated OG
	SpeciesNumGenesOG[sp,8] <- summary[2] # 1st quartile
	SpeciesNumGenesOG[sp,9] <- summary[3] # median
	SpeciesNumGenesOG[sp,10] <- summary[5] # 3rd quartile
}
SpeciesNumGenesOG <- as.data.frame(SpeciesNumGenesOG)
colnames(SpeciesNumGenesOG) <- c("Species", "totalOG", "dupOG", "totalGenes", "wothologGenes", "dupGenes", "meanSizeDupOG", "1quartSizeDupOG", "medianSizeDupOG", "3quartSizeDupOG")
for(col in c(2:length(SpeciesNumGenesOG[1,]))){
	SpeciesNumGenesOG[,colnames(SpeciesNumGenesOG)[col]] <- as.numeric(SpeciesNumGenesOG[,colnames(SpeciesNumGenesOG)[col]])
}
colnames(SpeciesNumGenesOG) <- c("Species", "totalOG", "dupOG", "totalGenes", "wothologGenes", "dupGenes", "meanSizeDupOG", "1quartSizeDupOG", "medianSizeDupOG", "3quartSizeDupOG")
print(SpeciesNumGenesOG)

ContingencyTableCN(
	OGInfo$Blan[which(OGInfo$BlanType=="Small-scale\nduplicates" & OGInfo$VertType=="Small-scale\nduplicates")], 
	OGInfo$MeanVert[which(OGInfo$BlanType=="Small-scale\nduplicates" & OGInfo$VertType=="Small-scale\nduplicates")],
	paste0(ResultsFolder, "/ContingencyTableCN_SS.txt"))

ContingencyTableCN(
	OGInfo$Blan[which(OGInfo$BlanType=="Small-scale\nduplicates" & OGInfo$VertType=="Ohnologs")], 
	OGInfo$MeanVert[which(OGInfo$BlanType=="Small-scale\nduplicates" & OGInfo$VertType=="Ohnologs")],
	paste0(ResultsFolder, "/ContingencyTableCN_O.txt"))

ContingencyTableCN(
	OGInfo$Blan[which(OGInfo$BlanType=="Small-scale\nduplicates" & (OGInfo$VertType=="Small-scale\nduplicates" | OGInfo$VertType=="Ohnologs"))], 
	OGInfo$MeanVert[which(OGInfo$BlanType=="Small-scale\nduplicates" & (OGInfo$VertType=="Small-scale\nduplicates" | OGInfo$VertType=="Ohnologs"))],
	paste0(ResultsFolder, "/ContingencyTableCN_SSO.txt"))




###########################################################################
###########################################################################
### Plotting
pdf(paste(ResultsFolder, "/DuplicatesOntologyExpression.pdf", sep=""), width=15, height=10)
par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')

# Blan vs. vertebrate gene type 
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
BarPlotVertTypesInBlanGenes(OGInfo.BlanVert, BlanTypes[which(BlanTypes!="Missing")], VertTypes[which(VertTypes!="Missing")], VertTypes.col[which(VertTypes!="Missing")])
Hypergeom <- PlotHypergeomTest_VertBlanTypes(OGInfo.BlanVert, VertTypes[which(VertTypes!="Missing")], BlanTypes[which(BlanTypes!="Missing")], VertTypes.col[which(VertTypes!="Missing")], PQvalThreshold, ResultsFolder)
plot.new()
legend("bottomright", c("In vertebrates", VertTypes[which(VertTypes!="Missing")]), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VertTypes.col[which(VertTypes!="Missing")]), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)#
BarPlotVertTypesInBlanGenes(OGInfo, BlanTypes, VertTypes, VertTypes.col)
HypergeomMissing <- PlotHypergeomTest_VertBlanTypes(OGInfo, VertTypes, BlanTypes, VertTypes.col, PQvalThreshold, ResultsFolder)
plot.new()
legend("bottomright", c("In vertebrates", VertTypes), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VertTypes.col), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)


layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T), widths=c(1), heights=c(1), TRUE)
BarPlotSpeciesNumGenesOG(SpeciesNumGenesOG, Vertebrates)
BarPlotVertTypesInBlanGenes_wexpected(OGInfo.BlanVert, BlanTypes[which(BlanTypes!="Missing")], VertTypes[which(VertTypes!="Missing")], VertTypes.col[which(VertTypes!="Missing")], Hypergeom$FCmatrix)

head(GenePairs)
head(GeneInfo)


# GO term in duplicates Hsap vs. Blan
# MF
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
ScatterPlotPointSize(GOInfo.MF$HsapS/GOInfo.MF$Hsap*100, GOInfo.MF$BlanS/GOInfo.MF$Blan*100, GOInfo.MF$Hsap, "black", "% of single-copy\nH. sapiens", "B. lanceolatum\n% of single-copy", c(0,100), c(0,100))
ScatterPlotPointSize((GOInfo.MF$HsapD+GOInfo.MF$HsapO)/GOInfo.MF$Hsap*100, GOInfo.MF$BlanD/GOInfo.MF$Blan*100, GOInfo.MF$Hsap, "black", "% of small-scale duplicates + ohnologs\nH. sapiens", "B. lanceolatum\n% of single-copy", c(0,100), c(0,100))
ScatterPlotPointSize(GOInfo.MF$HsapD/GOInfo.MF$Hsap*100, GOInfo.MF$BlanD/GOInfo.MF$Blan*100, GOInfo.MF$Hsap, "black", "% of small-scale duplicates\nH. sapiens", "B. lanceolatum\n% of small-scale duplicates", c(0,100), c(0,100))
ScatterPlotPointSize(GOInfo.MF$HsapO/GOInfo.MF$Hsap*100, GOInfo.MF$BlanD/GOInfo.MF$Blan*100, GOInfo.MF$Hsap, "black", "% of ohnologs\nH. sapiens", "B. lanceolatum\n% of small-scale duplicates", c(0,100), c(0,100))
# BP
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
ScatterPlotPointSize(GOInfo.BP$HsapS/GOInfo.BP$Hsap*100, GOInfo.BP$BlanS/GOInfo.BP$Blan*100, GOInfo.BP$Hsap, "black", "% of single-copy\nH. sapiens", "B. lanceolatum\n% of single-copy", c(0,100), c(0,100))
ScatterPlotPointSize((GOInfo.BP$HsapD+GOInfo.BP$HsapO)/GOInfo.BP$Hsap*100, GOInfo.BP$BlanD/GOInfo.BP$Blan*100, GOInfo.BP$Hsap, "black", "% of small-scale duplicates + ohnologs\nH. sapiens", "B. lanceolatum\n% of single-copy", c(0,100), c(0,100))
ScatterPlotPointSize(GOInfo.BP$HsapD/GOInfo.BP$Hsap*100, GOInfo.BP$BlanD/GOInfo.BP$Blan*100, GOInfo.BP$Hsap, "black", "% of small-scale duplicates\nH. sapiens", "B. lanceolatum\n% of small-scale duplicates", c(0,100), c(0,100))
ScatterPlotPointSize(GOInfo.BP$HsapO/GOInfo.BP$Hsap*100, GOInfo.BP$BlanD/GOInfo.BP$Blan*100, GOInfo.BP$Hsap, "black", "% of ohnologs\nH. sapiens", "B. lanceolatum\n% of small-scale duplicates", c(0,100), c(0,100))


layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(15), heights=c(5), TRUE)
TandemIntraInterPerSpecies(OGInfo, Species[which(Species!="Bflo" & Species!="Bbel" & Species!="Drer")], SpType[which(Species!="Bflo" & Species!="Bbel" & Species!="Drer")])

layout(matrix(c(1,2,5,3,4,6),nrow=2,ncol=3,byrow=T), widths=c(1.2,1.2), heights=c(1.5, 1.5), TRUE)
BoxPlot_BlanTypes(BlanGeneData, GeneInfo, "MeanAdult", "Mean adult tissues expression", c(0,100), VertTypes, VertTypes.col)
BoxPlot_BlanTypes(BlanGeneData, GeneInfo, "MeanEmbr", "Mean embrionic stages expression", c(0,100), VertTypes, VertTypes.col)
BoxPlot_BlanTypes(BlanGeneData, GeneInfo, "TauTissues", "Tau among adult tissues", c(0,1), VertTypes, VertTypes.col)
BoxPlot_BlanTypes(BlanGeneData, GeneInfo, "TauEmbAge", "Tau among embrionic stages", c(0,1), VertTypes, VertTypes.col)

layout(matrix(c(1,2,5,3,4,6),nrow=2,ncol=3,byrow=T), widths=c(1.2,1.2), heights=c(1.5, 1.5), TRUE)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneInfo, "MeanAdult", "Mean adult tissues expression", c(0,100), VertTypes, VertTypes.col)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneInfo, "MeanEmbr", "Mean embrionic stages expression", c(0,100), VertTypes, VertTypes.col)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneInfo, "TauTissues", "Tau among adult tissues", c(0,1), VertTypes, VertTypes.col)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneInfo, "TauEmbAge", "Tau among embrionic stages", c(0,1), VertTypes, VertTypes.col)

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
Hist_ExpressDomanis(GenePairs, OGInfo, "Small-scale\nduplicates", "Single-copy", MatchingTissues, "Blan conditions - Drer conditions", "Relative number of\npairwise comparisons", "B. lanceolatum specific\nsmall-scale gene duplicates")
Hist_ExpressDomanis(GenePairs, OGInfo, "Single-copy", "Small-scale\nduplicates", MatchingTissues, "Blan conditions - Drer conditions", "Relative number of\npairwise comparisons", "D. rerio specific\nsmall-scale gene duplicates")
Hist_ExpressDomanis(GenePairs, OGInfo, "Single-copy", "Ohnologs", MatchingTissues, "Blan conditions - Drer conditions", "Relative number of\npairwise comparisons", "D. rerio specific\nohnolog gene duplicates")
Hist_ExpressDomanis(GenePairs, OGInfo, NA, NA, MatchingTissues, "Blan conditions - Drer conditions", "Relative frequency", "Single-copy vs. single-copy")

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
Hist_ExpressDomanis(
	GenePairs[which(GenePairs$BlanLType=="Inter" | GenePairs$BlanType=="Single-copy"),], 
	OGInfo[which(OGInfo$BlanLType=="Inter"),], "Small-scale\nduplicates", "Single-copy", MatchingTissues, "Blan conditions - Drer conditions", "Relative number of\npairwise comparisons", "B. lanceolatum specific\nmultichormosomal\ngene duplicates")
Hist_ExpressDomanis(
	GenePairs[which(GenePairs$BlanLType=="Intra" | GenePairs$BlanType=="Single-copy"),], 
	OGInfo[which(OGInfo$BlanLType=="Intra"),], "Small-scale\nduplicates", "Single-copy", MatchingTissues, "Blan conditions - Drer conditions", "Relative number of\npairwise comparisons", "B. lanceolatum specific\ndistant monochromosomal\ngene duplicates")
Hist_ExpressDomanis(
	GenePairs[which(GenePairs$BlanLType=="Tandem" | GenePairs$BlanType=="Single-copy"),], 
	OGInfo[which(OGInfo$BlanLType=="Tandem"),], "Small-scale\nduplicates", "Single-copy", MatchingTissues, "Blan conditions - Drer conditions", "Relative number of\npairwise comparisons", "B. lanceolatum specific\ntandem\ngene duplicates")

layout(matrix(c(1,2,3,4,5,6,7,8,9),nrow=3,ncol=3,byrow=T), widths=c(1.5), heights=c(1), TRUE)
Dist_ExpressDomanis(GenePairs$DiffDom[which(GenePairs$BlanType=="Single-copy" & GenePairs$VertType=="Single-copy")], "forestgreen", MatchingTissues, "Blan conditions - Drer conditions", "", "SC - SC", 2)
Dist_ExpressDomanis(GenePairs$DiffDom[which(GenePairs$BlanType=="Small-scale\nduplicates" & GenePairs$VertType=="Single-copy")], "red", MatchingTissues, "Blan conditions - Drer conditions", "", "SSD - SC", 2)
Dist_ExpressDomanis(OGInfo$DiffDom[which(OGInfo$BlanType=="Small-scale\nduplicates" & OGInfo$VertType=="Single-copy")], "darkred", MatchingTissues, "Blan conditions - Drer conditions", "", "uSSD - SC", 2)
Dist_ExpressDomanis(GenePairs$DiffDom[which(GenePairs$VertType=="Small-scale\nduplicates" & GenePairs$BlanType=="Single-copy")], "red", MatchingTissues, "Blan conditions - Drer conditions", "", "SC - SSD", 2)
Dist_ExpressDomanis(OGInfo$DiffDom[which(OGInfo$VertType=="Small-scale\nduplicates" & OGInfo$BlanType=="Single-copy")], "darkred", MatchingTissues, "Blan conditions - Drer conditions", "", "SC - uSSD", 2)
Dist_ExpressDomanis(GenePairs$DiffDom[which(GenePairs$VertType=="Ohnologs" & GenePairs$BlanType=="Single-copy")], "red", MatchingTissues, "Blan conditions - Drer conditions", "", "SC - O", 3)
Dist_ExpressDomanis(OGInfo$DiffDom[which(OGInfo$VertType=="Ohnologs" & OGInfo$BlanType=="Single-copy")], "darkred", MatchingTissues, "Blan conditions - Drer conditions", "", "SC - uO", 3)


layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
DifferenceWithSC(GenePairs, OGInfo, MatchingTissues, "Branch specific gene duplicates", "Difference with\nsingle-copy genes distribution")
DifferenceWithSC_TandemIntraInter(GenePairs, OGInfo, MatchingTissues, "Branch specific gene duplicates", "Difference with\nsingle-copy genes distribution")

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
ScatterPlot_GOExp(GOInfo.MF$BlanD.AExp, GOInfo.MF$BlanSC.AExp, GOInfo.MF$Blan, "black", "Mean adult expression\nof duplicate genes", "Mean adult expression\nof single-copy genes", c(0,120))
ScatterPlot_GOExp(GOInfo.MF$BlanD.EExp, GOInfo.MF$BlanSC.EExp, GOInfo.MF$Blan, "black", "Mean embrionic expression\nof duplicate genes", "Mean embrionic expression\nof single-copy genes", c(0,120))
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
ScatterPlot_GOExp(GOInfo.BP$BlanD.AExp, GOInfo.BP$BlanSC.AExp, GOInfo.BP$Blan, "black", "Mean adult expression\nof duplicate genes", "Mean adult expression\nof single-copy genes", c(0,120))
ScatterPlot_GOExp(GOInfo.BP$BlanD.EExp, GOInfo.BP$BlanSC.EExp, GOInfo.BP$Blan, "black", "Mean embrionic expression\nof duplicate genes", "Mean embrionic expression\nof single-copy genes", c(0,120))


dev.off()


###########################################################################
###########################################################################

