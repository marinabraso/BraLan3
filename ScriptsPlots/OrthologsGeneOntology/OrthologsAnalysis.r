#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

#`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.r", sep=""))
library(ghibli)
library(viridis)
library(car)

######################################################################
# Files & folders

ResultsFolder <- "Plots/OrthologsGeneOntology"
CountsFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt"
NamesFile <- "Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt"
OhnologsSFile <- "Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt"
GOMFFile <- paste("Plots/GeneOntologyMainMF_FromHsap2BlanProcessed.txt", sep ="")
ProteomesFolder <- "Results/FilteringGeneSets/Proteomes"
GOlistsFolder <- "Results/GeneOntology"
MainGOlistFile <- "Data/GeneOntology/MainGOMF_supclass.list"
GOlistFile <- "Results/GeneOntology/GOlist_MF.list"

######################################################################
# General parameters
PQvalThreshold <- 0.01 # for hypergeometric test
minGOsize <- 200

Species <- c("Blan", "Bflo", "Bbel", "Drer", "Ggal", "Mmus", "Hsap")
Ampioxus <- c("Blan", "Bflo", "Bbel")
Vertebrates <- c("Drer", "Ggal", "Mmus", "Hsap")
GenesBasenames <- c("BLA", "BFLO", "BBEL", "ENSDAR", "ENSGAL", "ENSMUS", "ENS")
SpeciesProteomes <- c("Branchiostoma_lanceolatum.BraLan3.fa", "Branchiostoma_floridae.Bfl_VNyyK.fa", "Branchiostoma_belcheri.Haploidv18h27.fa", "Danio_rerio.GRCz11.fa", "Gallus_gallus.GRCg6a.fa", "Mus_musculus.GRCm39.fa", "Homo_sapiens.GRCh38.fa")

BaseColor <- ghibli_palettes$MarnieMedium2[5]
VertTypes <- c("Missing", "SingleCopy", "Ohnolog", "Duplicated")
VertTypes.col <- c("royalblue2", ghibli_palettes$MarnieMedium2[c(2,4,7)])
#VertTypes.col <- ghibli_palettes$MarnieMedium2[c(7,2,4,6)]
BlanTypes <- c("Missing", "SingleCopy", "Duplicated")
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
Counts$VertType <- c("Missing", "SingleCopy", "Duplicated")[unlist(apply(Counts[,Vertebrates], 1, function(x){r <- max(x); if(r>1){r=2}; return(r)}))+1]
Counts$VertType[rownames(Counts) %in% OnhOG.S] <- "Ohnolog"
table(Counts$VertType)
Counts$AmphType <- c("Missing", "SingleCopy", "Duplicated")[unlist(apply(Counts[,Ampioxus], 1, function(x){r <- max(x); if(r>1){r=2}; return(r)}))+1]
table(Counts$AmphType)
Counts$BlanType <- c("Missing", "SingleCopy", "Duplicated")[unlist(lapply(Counts[,"Blan"], function(x){r <- x; if(r>1){r=2}; return(r)}))+1]
table(Counts$BlanType)
Counts$Sum <- apply(Counts[,Species], 1, sum)

system_out <- system(paste("cat ", NamesFile," | awk -F '\t' '{if(NR==1){split($0,h,\"\t\");next} split($0,a,\"\t\"); for(i = 2; i <= length(a); ++i){split(a[i],g,\" \"); for(j in g){print a[1]\"\t\"h[i]\"\t\"g[j]}}}'"), intern=T)
OG2Gene <- read.table(text=system_out, h=F, sep = "\t")
colnames(OG2Gene) <- c("OG", "Species", "Gene")
OG2Gene$Species <- Species[match(OG2Gene$Species, SpeciesProteomes)]

NumGenesSpecies <- rep(0, length(Species))
for(sp in c(1:length(Species))){
	system_out <- system(paste0("cat ", ProteomesFolder, "/", SpeciesProteomes[sp], " | grep '>' | sort | uniq | wc -l"), intern=T)
	NumGenesSpecies[sp] <- read.table(text=system_out, h=F, sep = "\t")
}
NumGenesSpecies <- unlist(NumGenesSpecies)


# GO terms
mGO <- read.table(MainGOlistFile, h=F, sep = '\t')
colnames(mGO) <- c("GO", "Class", "Name")
colfunc <- colorRampPalette(c("forestgreen", "gold2", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
GOClass <- as.data.frame(cbind(unique(mGO$Class), sample(colfunc(length(unique(mGO$Class))))))
colnames(GOClass) <- c("Class", "Color")
mGO$ClassColors <- GOClass$Color[match(mGO$Class, GOClass$Class)]
GO <- read.delim(GOlistFile, h=F, sep = '\t')
colnames(GO) <- c("GO", "Class", "Name")
GO$ClassColors <- GOClass$Color[match(GO$Class, GOClass$Class)]
print(head(GO))

GO$Hsap <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts, OG2Gene, GOlistsFolder, minGOsize))
GO <- GO[which(GO$Hsap>0),]
GO <- GO[order(GO$Hsap),]
GO$HsapD <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts[which(Counts$VertType=="Duplicated"),], OG2Gene, GOlistsFolder, minGOsize))
GO$HsapS <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts[which(Counts$VertType=="SingleCopy"),], OG2Gene, GOlistsFolder, minGOsize))
GO$HsapO <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Hsap", Counts[which(Counts$VertType=="Ohnolog"),], OG2Gene, GOlistsFolder, minGOsize))
GO$Blan <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Blan", Counts, OG2Gene, GOlistsFolder, minGOsize))
GO$BlanD <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Blan", Counts[which(Counts$BlanType=="Duplicated"),], OG2Gene, GOlistsFolder, minGOsize))
GO$BlanS <- unlist(lapply(GO$GO, GetNumberOfGOtermGenes, "Blan", Counts[which(Counts$BlanType=="SingleCopy"),], OG2Gene, GOlistsFolder, minGOsize))

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

BarPlotVertTypesInBlanGenes(Counts[which(Counts$VertType!="Missing" & Counts$BlanType!="Missing"),], BlanTypes[which(BlanTypes!="Missing")], VertTypes[which(VertTypes!="Missing")], VertTypes.col[which(VertTypes!="Missing")])
PlotHypergeomTest_VertBlanTypes(Counts[which(Counts$VertType!="Missing" & Counts$BlanType!="Missing"),], VertTypes[which(VertTypes!="Missing")], BlanTypes[which(BlanTypes!="Missing")], VertTypes.col[which(VertTypes!="Missing")], PQvalThreshold)
plot.new()
legend("bottomright", c("In vertebrates", VertTypes[which(VertTypes!="Missing")]), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VertTypes.col[which(VertTypes!="Missing")]), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)

BarPlotVertTypesInBlanGenes(Counts, BlanTypes, VertTypes, VertTypes.col)
PlotHypergeomTest_VertBlanTypes(Counts, VertTypes, BlanTypes, VertTypes.col, PQvalThreshold)
plot.new()
legend("bottomright", c("In vertebrates", "Missing", "Single copy", "Ohnolog", "Duplicated"), pch=c(NA,15,15,15,15), text.col="black", col=c(NA, VertTypes.col), bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)

ScatterPercentagePlot(GO$HsapD/GO$Hsap*100, GO$BlanD/GO$Blan*100, GO$ClassColors, "% of small scale duplicated genes\nH. sapiens", "B. lanceolatum\n% of duplicated genes", c(20,100), c(0,100))
ScatterPercentagePlot(GO$HsapO/GO$Hsap*100, GO$BlanD/GO$Blan*100, GO$ClassColors, "% of ohnolog duplicated genes\nH. sapiens", "B. lanceolatum\n% of duplicated genes", c(0,40), c(0,100))
plot.new()
legendvec <- unique(GO[,c("Class", "ClassColors")])
legend("center", legendvec[,1], pch=15, col=legendvec[,2], bty = "n", pt.cex=2, cex=1.5, xjust = 0, yjust = 0)


dev.off()


write.table(GO[order(GO$BlanD/GO$Blan*100),], file = paste(ResultsFolder, "/GOterms_PropData.txt", sep =""), quote = F, sep="\t", col.names = TRUE, row.names = TRUE)





