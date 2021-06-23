#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
source(paste(dirname(script), "/Amphioxus_Vertebrate_categories.R", sep=""))
source(paste(dirname(script), "/GeneExpression.R", sep=""))
source(paste(dirname(script), "/PositiveSelection.R", sep=""))
source(paste(dirname(script), "/Synteny.R", sep=""))
source(paste(dirname(script), "/GeneFunction.R", sep=""))
library(qvalue)
library(ghibli)
library(alluvial)
library(viridis)

######################################################################
# Files & folders

RDataFolder <- "Plots"
ResultsFolder <- "Plots/OrthologsAnalysis"
GeneDataFile <- paste(RDataFolder, "/GeneDataProcessed.txt", sep ="")
OGCountsFile <- paste(RDataFolder, "/OGCountsProcessed.txt", sep ="")
RNASeqMetadataFile <- paste(RDataFolder, "/RNASeqMetadataProcessed.txt", sep ="")
GOMFFile <- paste(RDataFolder, "/GeneOntologyMainMF_FromHsap2BlanProcessed.txt", sep ="")
ChrLengthFile <- "Results/AssemblyStatistics/Branchiostoma_lanceolatum.BraLan3/Branchiostoma_lanceolatum.BraLan3_lengths.txt"

######################################################################
# General parameters
QvalT <- 0.1 # q-value threshold for positive selection
PQvalThreshold <- 0.01 # for hypergeometric test
TandemMaxDist <- 1000 

Verteb <- c("Drer", "Ggal", "Mmus", "Hsap")
Amphi <- c("Blan", "Bflo", "Bbel")
OutDeut <- c("Spur", "Arub", "Skow")
SortedSpecies <- c(Amphi, Verteb, OutDeut)

Tissues <- c("cirri", "gills", "epidermis", "gut", "hepatic.diverticulum", "muscle", "neural.tube", "female.gonads", "male.gonads")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
TissueColors <- colfunc(length(Tissues))

EmbAges <- c("egg", "32cells", "Blastula", "7h", "8h", "10h", "11h", "15h", "18h", "21h", "24h", "27h", "36h", "50h", "60h", "Pre.metamorphic.larvae")
#EmbAges <- c("32cells.Blastula", "7h.8h", "10h.11h", "15h.18h", "21h.24h", "27h.36h", "50h.60h")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
EmbAgesColors <- colfunc(length(EmbAges))

VTypes <- c("Missing", "SingleCopy", "Ohnolog", "Duplicated")
BlanTypes <- c("SingleCopy", "Duplicated")
VTypeColors <- c("deepskyblue2", "gold", "olivedrab2", "olivedrab4")
VTypeColors <- c("#481567FF", "#FDE725FF", "#1F968BFF", "#39568CFF")
VTypeColors <- ghibli_palettes$MarnieMedium2[c(7,2,4,6)]
BaseColor <- ghibli_palettes$MarnieMedium2[5]
BlanType.pch <- c(16,18)
BlanType.lty <- c(1,2)

# Big molecular function terms
MFunctGO <- c("GO:0016887", "GO:0016209", "GO:0003824", "GO:0005198", "GO:0140110", "GO:0140223", "GO:0005215", "GO:0038024", "GO:0004857", "GO:0008047", "GO:0019207", "GO:0099106", "GO:0098772", "GO:0030546", "GO:0019888", "GO:0044092", "GO:0044093", "GO:0065009", "GO:0045182", "GO:0005488", "GO:0060090", "GO:0140104", "GO:0060089")
MFunctNames <- c("ATPase", "AntioxidantActivity", "CatalyticActivity", "StructuralMoleculeActivity", "TranscriptionRegulatorActivity", "GeneralTranscriptionInitiationFactorActivity", "TransporterActivity", "CargoReceptorActivity", "EnzymeInhibitorActivity", "EnzymeActivatorActivity", "KinaseRegulatorActivity", "IonChannelRegulatorActivity", "MolecularFunctionRegulator", "SignalingReceptorActivatorActivity", "ProteinPhosphataseRegulatorActivity", "NegativeRegulationOfMolecularFunction", "PositiveRegulationOfMolecularFunction", "RegulationOfMolecularFunction", "TranslationRegulatorActivity", "Binding", "MolecularAdaptorActivity", "MolecularCarrierActivity", "MolecularTransducerActivity")
GeneralMFNames <- c("Enzymes", "StructuralProteins", "TranscriptionFactors", "Transporters", "Regulators", "Binding", "Carriers", "Signaling")
GeneralMFNum <- c(1, 1, 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7, 8)
GeneralMFMidp <- c(2, 4, 5.5, 7.5, 14, 20.5, 22, 23)
GeneralMFColors <- sample(c(ghibli_palettes$MononokeMedium, "darkred"))
GeneralMFColors <- viridis(8)


###########################################################################
###########################################################################
# Read data

GeneData <- read.table(GeneDataFile, h=T, sep = "\t", row.names=1)
colnames(GeneData) <- sub("^X", "", colnames(GeneData))

OGData.AV <- read.table(OGCountsFile, h=T, sep = "\t", row.names=1)
OGData.AV$DupType <- rep(NA, length(OGData.AV[,1]))
OGData.AV$DupType[which(OGData.AV$NumChrs==1 & OGData.AV$BlanType=="Duplicated")] <- rep("Intra",sum(OGData.AV$NumChrs==1 & OGData.AV$BlanType=="Duplicated"))
OGData.AV$DupType[which(OGData.AV$NumChrs>1 & OGData.AV$BlanType=="Duplicated")] <- rep("Inter",sum(OGData.AV$NumChrs>1 & OGData.AV$BlanType=="Duplicated"))
OGData.AV$DupType[which(OGData.AV$NumChrs==1 & OGData.AV$MaxDist<=TandemMaxDist & OGData.AV$BlanType=="Duplicated")] <- rep("Tandem",sum(OGData.AV$NumChrs==1 & OGData.AV$MaxDist<=TandemMaxDist & OGData.AV$BlanType=="Duplicated"))
table(OGData.AV$DupType)
table(OGData.AV$BlanType)
head(OGData.AV)

MetaRNA <- read.table(RNASeqMetadataFile, h=T, sep = "\t", row.names=1)

ChrLengths <- read.table(ChrLengthFile, h=F, sep = "\t", row.names=1)
colnames(ChrLengths) <- c("Length", "Chr", "AccLength")
rownames(ChrLengths) <- sub("^>", "", rownames(ChrLengths))

GOMFData <- read.table(GOMFFile, h=T, sep = "\t")
print(head(GOMFData))

###########################################################################
###########################################################################
#Amphioxus_Vertebrate_categories(OGData.AV)

#GeneExpression(GeneData)

#PositiveSelection(GeneData, OGData.AV)

#Synteny(GeneData, OGData.AV, ChrLengths)

GeneFunction(GeneData, OGData.AV, GOMFData)

print("Nodal (Z linked?)")
print(GeneData[which(GeneData$Gene=="BLAG15000766"),])


