#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

#`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.r", sep=""))
# source("ScriptsPlots/OrthologsGeneOntology/OrthologsAnalysis_functions.r")
library(ghibli)
library(viridis)
library(car)
library(MASS)

######################################################################
# Files & folders

ResultsFolder <- "Plots/OrthologsGeneOntology"
CountsFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
OhnologsSFile <- "Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt"
GOMFFile <- paste("Plots/GeneOntologyMainMF_FromHsap2BlanProcessed.txt", sep ="")
GOlistsFolder <- "Results/GeneOntology"
MainGOlistFile <- "Data/GeneOntology/MainGOMF_supclass.list"
GOlistFile <- "Results/GeneOntology/GOlist_MF.list"
GTFFolder <- "Results/FilteringGeneSets/FilteredGTFs"

######################################################################
# General parameters
PQvalThreshold <- 0.01 # for hypergeometric test
minGOsizeHsap <- 50
minGOsizeBlan <- 10
MaxTandemDist <- 10000


Species <- c("Blan", "Bflo", "Bbel", "Drer", "Ggal", "Mmus", "Hsap")
SpType <- c("AmphType", "AmphType", "AmphType", "VertType", "VertType", "VertType", "VertType")
Ampioxus <- c("Blan", "Bflo", "Bbel")
Vertebrates <- c("Drer", "Ggal", "Mmus", "Hsap")
GenesBasenames <- c("BLA", "BFLO", "BBEL", "ENSDAR", "ENSGAL", "ENSMUS", "ENS")
SpeciesProteomes <- c("Branchiostoma_lanceolatum.BraLan3.fa", "Branchiostoma_floridae.Bfl_VNyyK.fa", "Branchiostoma_belcheri.Haploidv18h27.fa", "Danio_rerio.GRCz11.fa", "Gallus_gallus.GRCg6a.fa", "Mus_musculus.GRCm39.fa", "Homo_sapiens.GRCh38.fa")
SpeciesGTFs <- c("Branchiostoma_lanceolatum.BraLan3.gtf", "Branchiostoma_floridae.Bfl_VNyyK.gtf", "Branchiostoma_belcheri.Haploidv18h27.gtf", "Danio_rerio.GRCz11.gtf", "Gallus_gallus.GRCg6a.gtf", "Mus_musculus.GRCm39.gtf", "Homo_sapiens.GRCh38.gtf")
BaseColor <- ghibli_palettes$MarnieMedium2[5]
VertTypes <- c("Single-copy", "Ohnologs", "Small-scale\nduplicates", "Missing")
VertTypes.col <- c(ghibli_palettes$MarnieMedium2[c(2,4,7)], "royalblue2")
#VertTypes.col <- ghibli_palettes$MarnieMedium2[c(7,2,4,6)]
BlanTypes <- c("Single-copy", "Small-scale\nduplicates", "Missing")
BlanTypes.pch <- c(14, 16, 18)
BlanTypes.lty <- c(.5, 1, 2)

# Big molecular function terms
MFunctGO <- c("GO:0016887", "GO:0016209", "GO:0003824", "GO:0005198", "GO:0140110", "GO:0140223", "GO:0005215", "GO:0038024", "GO:0004857", "GO:0008047", "GO:0019207", "GO:0099106", "GO:0098772", "GO:0030546", "GO:0019888", "GO:0044092", "GO:0044093", "GO:0065009", "GO:0045182", "GO:0005488", "GO:0060090", "GO:0140104", "GO:0060089")
MFunctNames <- c("ATPase", "AntioxidantActivity", "CatalyticActivity", "StructuralMoleculeActivity", "TranscriptionRegulatorActivity", "GeneralTranscriptionInitiationFactorActivity", "TransporterActivity", "CargoReceptorActivity", "EnzymeInhibitorActivity", "EnzymeActivatorActivity", "KinaseRegulatorActivity", "IonChannelRegulatorActivity", "MolecularFunctionRegulator", "SignalingReceptorActivatorActivity", "ProteinPhosphataseRegulatorActivity", "NegativeRegulationOfMolecularFunction", "PositiveRegulationOfMolecularFunction", "RegulationOfMolecularFunction", "TranslationRegulatorActivity", "Binding", "MolecularAdaptorActivity", "MolecularCarrierActivity", "MolecularTransducerActivity")
GeneralMFNames <- c("Enzymes", "StructuralProteins", "TranscriptionFactors", "Transporters", "Regulators", "Binding", "Carriers", "Signaling")
GeneralMFNum <- c(1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7, 8)
GeneralMFMidp <- c(2, 4, 5.5, 7.5, 14, 20.5, 22, 23)
GeneralMFColors <- viridis(8)


###########################################################################
###########################################################################
# Read data


# List of onhologs
OnhOG.S <- read.table(OhnologsSFile, h=F)[,1]

# Orthologous groups counts
Counts <- read.table(CountsFile, h=F, sep = "\t", row.names=1)
system_out <- system(paste("head -1 ", CountsFile, " | cut -f2,3,4,5,6,7,8,9,10,11,12 | sed 's/\\([A-Z]\\)[a-z]\\+_\\([a-z][a-z][a-z]\\)[A-Za-z0-9\\._]\\+/\\1\\2/g'"), intern=T)
colnames(Counts) <- as.character(unlist(lapply(read.table(text=system_out, h=F, sep = "\t")[1,], as.character)))
Counts$OG <- rownames(Counts)
Counts$VertType <- c("Missing", "Single-copy", "Small-scale\nduplicates")[unlist(apply(Counts[,Vertebrates], 1, function(x){r <- max(x); if(r>1){r=2}; return(r)}))+1]
Counts$VertType[rownames(Counts) %in% OnhOG.S] <- "Ohnologs"
table(Counts$VertType)
Counts$AmphType <- c("Missing", "Single-copy", "Small-scale\nduplicates")[unlist(apply(Counts[,Ampioxus], 1, function(x){r <- max(x); if(r>1){r=2}; return(r)}))+1]
table(Counts$AmphType)
Counts$BlanType <- c("Missing", "Single-copy", "Small-scale\nduplicates")[unlist(lapply(Counts[,"Blan"], function(x){r <- x; if(r>1){r=2}; return(r)}))+1]
table(Counts$BlanType)
Counts$Sum <- apply(Counts[,Species], 1, sum)
Counts$MeanVert <- rowMeans(Counts[,Vertebrates])
Counts.BlanVert <- Counts[which(Counts$BlanType != "Missing" | Counts$VertType != "Missing"),]


system_out <- system(paste("cat ", NamesFile," | awk -F '\t' '{if(NR==1){split($0,h,\"\t\");next} split($0,a,\"\t\"); for(i = 2; i <= length(a); ++i){split(a[i],g,\" \"); for(j in g){print a[1]\"\t\"h[i]\"\t\"g[j]}}}'"), intern=T)
OG2Gene <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2Gene) <- c("OG", "Species", "Gene")
OG2Gene$Species <- Species[match(OG2Gene$Species, SpeciesProteomes)]

Coord <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Gene","Chr","Start", "End"))
for(s in c(1:length(Species))){
	system_out <- system(paste0("cat ",GTFFolder, "/", SpeciesGTFs[s], " | sed 's/[^\t]\\+\\?gene_id \"\\(\\S\\+\\)\";.*/\\1/g' | awk -F '\\t' '{print $9\"\t\"$1\"\t\"$4\"\t\"$5}' | sort -k1,1 -k2,2V -k3,3V | awk -F '\\t' '{if($1==gene){end=$4}else{print gene\"\t\"chr\"\t\"start\"\t\"end\"\t\"pos; if($2==chr){pos=pos+1}else{pos=1}; gene=$1; chr=$2; start=$3; end=$4}}END{print gene\"\t\"chr\"\t\"start\"\t\"end\"\t\"pos;}' | tail -n +2"), intern=T)
	tmp <- read.table(text=system_out, h=F, sep = "\t")
	print(head(tmp))
	Coord <- rbind(Coord, tmp)
}
colnames(Coord) <- c("Gene","Chr","Start", "End", "Position")
GeneInfo <- merge(Coord, OG2Gene, by = "Gene")
GeneInfo <- merge(GeneInfo, Counts, by = "OG")

for(s in c(1:length(Species))){
	supGeneInfo <- GeneInfo[which(GeneInfo$Species==Species[s]),]
	print(Species[s])
	Counts[,paste0("NumChrs",Species[s])] <- unlist(lapply(Counts$OG, function(x){length(unique(supGeneInfo$Chr[which(supGeneInfo$OG==x)]))}))
	Counts[,paste0("MaxDist",Species[s])]  <- unlist(lapply(Counts$OG, CalcMaxDist, supGeneInfo))
	Counts[,paste0("Consecutive",Species[s])]  <- unlist(lapply(Counts$OG, CalcIfConsecutive, supGeneInfo))
}




# GO terms
GO <- read.delim(GOlistFile, h=F, sep = '\t')
colnames(GO) <- c("GO", "Class", "Name")
GO <- unique(GO[,c("GO", "Name")])
GO$Hsap <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts, OG2Gene, GOlistsFolder))
GO$Blan <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Blan", Counts, OG2Gene, GOlistsFolder))
GO <- GO[which(GO$Hsap>=minGOsizeHsap & GO$Blan>=minGOsizeHsap),]
GO <- GO[order(GO$Hsap),]
GO$HsapD <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts[which(Counts$VertType=="Small-scale\nduplicates"),], OG2Gene, GOlistsFolder))
GO$HsapS <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts[which(Counts$VertType=="Single-copy"),], OG2Gene, GOlistsFolder))
GO$HsapO <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts[which(Counts$VertType=="Ohnologs"),], OG2Gene, GOlistsFolder))
GO$BlanD <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Blan", Counts[which(Counts$BlanType=="Small-scale\nduplicates"),], OG2Gene, GOlistsFolder))
GO$BlanS <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Blan", Counts[which(Counts$BlanType=="Single-copy"),], OG2Gene, GOlistsFolder))
print(head(GO))


###########################################################################
###########################################################################
# Amphioxus - vertebrate shared OG numbers
paste("Shared V-A", sum(Counts$VertType!="Missing" & Counts$AmphType!="Missing"))
paste("V", sum(Counts$VertType!="Missing"))
paste("A", sum(Counts$AmphType!="Missing"))
paste("Shared/V*100", sum(Counts$VertType!="Missing" & Counts$AmphType!="Missing")/sum(Counts$VertType!="Missing")*100)
paste("Shared/A*100", sum(Counts$VertType!="Missing" & Counts$AmphType!="Missing")/sum(Counts$AmphType!="Missing")*100)

# Amphioxus specfic OG numbers
Counts.A <- Counts[which(Counts$VertType=="Missing" & Counts$AmphType!="Missing"),]
paste("A specific", length(Counts.A[,1]))
paste("Shared all A", sum(unlist(apply(Counts.A[,Ampioxus], 1, min))>0))
paste("Shared all A/A specific*100", sum(unlist(apply(Counts.A[,Ampioxus], 1, min))>0)/length(Counts.A[,1])*100)
paste("blan-bbel specific", sum(Counts.A$Blan>0 & Counts.A$Bflo==0 & Counts.A$Bbel>0))
paste("blan-bbel specific/A specific*100", sum(Counts.A$Blan>0 & Counts.A$Bflo==0 & Counts.A$Bbel>0)/length(Counts.A[,1])*100)
paste("blan-bflo specific", sum(Counts.A$Blan>0 & Counts.A$Bflo>0 & Counts.A$Bbel==0))
paste("blan-bflo specific/A specific*100", sum(Counts.A$Blan>0 & Counts.A$Bflo>0 & Counts.A$Bbel==0)/length(Counts.A[,1])*100)

for(sp in c(1:length(Species))){
	OG2Gene.sp <- OG2Gene[which(OG2Gene$Species==Species[sp]),]
	toprint <- c(1:7)
	toprint[1] <- Species[sp]
	toprint[2] <- sum(Counts[,Species[sp]]>1)
	toprint[3] <- sum(Counts[,Species[sp]]>1)/sum(Counts[,Species[sp]]>0)*100
	toprint[4] <- length(unique(OG2Gene.sp$Gene[which(OG2Gene.sp$OG %in% rownames(Counts[which(Counts[,Species[sp]]>1),]))]))
	toprint[5] <- length(unique(OG2Gene.sp$Gene[which(OG2Gene.sp$OG %in% rownames(Counts[which(Counts[,Species[sp]]>1),]))]))/NumGenesSpecies[sp]*100
	toprint[6] <- NumGenesSpecies[sp]-length(unique(OG2Gene.sp$Gene))
	toprint[7] <- (NumGenesSpecies[sp]-length(unique(OG2Gene.sp$Gene)))/NumGenesSpecies[sp]*100
	print(paste(toprint, collapse=" "))
}


# Plotting
pdf(paste(ResultsFolder, "/OrthologsGeneOntology.pdf", sep=""), width=15, height=10)
par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)

BarPlotVertTypesInBlanGenes(Counts.BlanVert[which(Counts.BlanVert$VertType!="Missing" & Counts.BlanVert$BlanType!="Missing"),], BlanTypes[which(BlanTypes!="Missing")], VertTypes[which(VertTypes!="Missing")], VertTypes.col[which(VertTypes!="Missing")])
PlotHypergeomTest_VertBlanTypes(Counts.BlanVert[which(Counts.BlanVert$VertType!="Missing" & Counts.BlanVert$BlanType!="Missing"),], VertTypes[which(VertTypes!="Missing")], BlanTypes[which(BlanTypes!="Missing")], VertTypes.col[which(VertTypes!="Missing")], PQvalThreshold, ResultsFolder)
plot.new()
legend("bottomright", c("In vertebrates", VertTypes[which(VertTypes!="Missing")]), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VertTypes.col[which(VertTypes!="Missing")]), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)#

BarPlotVertTypesInBlanGenes(Counts.BlanVert, BlanTypes, VertTypes, VertTypes.col)
PlotHypergeomTest_VertBlanTypes(Counts.BlanVert, VertTypes, BlanTypes, VertTypes.col, PQvalThreshold, ResultsFolder)
plot.new()
legend("bottomright", c("In vertebrates", "Missing", "Single-copy", "Ohnologs", "Small-scale\nduplicates"), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VertTypes.col), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)

layout(matrix(c(1,3,5,2,4,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
ScatterPlotContour(GO$HsapD/GO$Hsap*100, GO$BlanD/GO$Blan*100, "black", "% of small-scale duplicates\nH. sapiens", "B. lanceolatum\n% of small-scale duplicates", c(0,100), c(0,100))
ScatterPlotContour(GO$HsapO/GO$Hsap*100, GO$BlanD/GO$Blan*100, "black", "% of ohnologs\nH. sapiens", "B. lanceolatum\n% of small-scale duplicates", c(0,100), c(0,100))

layout(matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=T), widths=c(1.2, 1.2), heights=c(1), TRUE)

Scatter2HeatmapPlot(
	Counts$Blan[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Small-scale\nduplicates")], 
	Counts$MeanVert[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Small-scale\nduplicates")], 
	c("white", "lightyellow", "firebrick"), 
	"Copy-number\nB. lanceolatum", 
	"Mean copy-number\nH. sapiens", 
	"Small-scale duplicates", 
	c(2,10), 
	c(2,10))

Scatter2HeatmapPlot(
	Counts$Blan[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Ohnologs")], 
	Counts$MeanVert[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Ohnologs")], 
	c("white", "lightyellow", "firebrick"), 
	"Copy-number\nB. lanceolatum", 
	"Mean copy-number\nH. sapiens", 
	"Ohnologs", 
	c(2,10), 
	c(2,10))

Scatter2HeatmapPlot(
	Counts$Blan[which(Counts$BlanType=="Small-scale\nduplicates" & (Counts$VertType=="Small-scale\nduplicates" | Counts$VertType=="Ohnologs"))], 
	Counts$MeanVert[which(Counts$BlanType=="Small-scale\nduplicates" & (Counts$VertType=="Small-scale\nduplicates" | Counts$VertType=="Ohnologs"))], 
	c("white", "lightyellow", "firebrick"), 
	"Copy-number\nB. lanceolatum", 
	"Mean copy-number\nH. sapiens", 
	"Both", 
	c(2,10), 
	c(2,10))

layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(15), heights=c(5), TRUE)
TandemIntraInterPerSpecies(Counts, Species, SpType, MaxTandemDist, "MaxDistance")
layout(matrix(c(1,2),nrow=2,ncol=1,byrow=T), widths=c(15), heights=c(5), TRUE)
TandemIntraInterPerSpecies(Counts, Species, SpType, MaxTandemDist, "Consecutive")


dev.off()

ContingencyTableCN(
	Counts$Blan[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Small-scale\nduplicates")], 
	Counts$MeanVert[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Small-scale\nduplicates")],
	paste0(ResultsFolder, "/ContingencyTableCN_SS.txt"))

ContingencyTableCN(
	Counts$Blan[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Ohnologs")], 
	Counts$MeanVert[which(Counts$BlanType=="Small-scale\nduplicates" & Counts$VertType=="Ohnologs")],
	paste0(ResultsFolder, "/ContingencyTableCN_O.txt"))

ContingencyTableCN(
	Counts$Blan[which(Counts$BlanType=="Small-scale\nduplicates" & (Counts$VertType=="Small-scale\nduplicates" | Counts$VertType=="Ohnologs"))], 
	Counts$MeanVert[which(Counts$BlanType=="Small-scale\nduplicates" & (Counts$VertType=="Small-scale\nduplicates" | Counts$VertType=="Ohnologs"))],
	paste0(ResultsFolder, "/ContingencyTableCN_SSO.txt"))

GO$PercentHsapD = GO$HsapD/GO$Hsap*100
GO$PercentHsapO = GO$HsapO/GO$Hsap*100
GO$PercentBlanD = GO$BlanD/GO$Blan*100
write.table(GO[order(GO$PercentBlanD),], file = paste(ResultsFolder, "/GOterms_PropData.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = FALSE)



