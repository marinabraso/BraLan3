#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions
script <- sub(".*=", "", commandArgs()[4])
source(paste0(dirname(script), "/", substr(basename(script),1, nchar(basename(script))-2), "_functions.r", sep=""))
# source("./ScriptsPlots/ComparativeGeneExpression/ComparativeGeneExpression_functions.r")
library(ghibli)

######################################################################
# Files & folders
RDataFolder <- "Plots"
ResultsFolder <- "Plots/ComparativeGeneExpression"
#OrganOrthBgeeFile <- paste0(pwd, "/Data/AnatomicalOntology/Julien_AnatOrht_Bgee15.tsv")
#CountsFile <- paste0(pwd, "/Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_counts.txt")
NamesFile <- paste0("Results/FindOrthologs/AmphVerteb_broccoli/dir_step3/table_OGs_protein_names.txt")
BlanGeneDataFile <- paste0(ResultsFolder, "/BlanGeneData.txt")
OhnologsSFile <- paste0("Results/FindOrthologs/OGOnhologs/OG_w_Onhologs.Strict.DR_GG_MM_HS.txt")
DrerGeneDataFile <- paste0(ResultsFolder, "/DrerGeneData.txt")
DrerBgeeDataFile <- paste0(ResultsFolder, "/DrerBgeeData.txt")
RNASeqMetadataFile <- "Metadata/Marletaz2018_RNAseq_SRA.txt"
TPMFile <- "Results/GeneExpression/Gene_TPM_kallisto_tximport.tab"

######################################################################
# General parameters
TPM.threshold <- 1
Species1 <- "Blan"
Species2 <- "Drer"

ShortSpeciesNames <- c("Drer", "Ggal", "Hsap", "Mmus", "Blan")
SpeciesNames <- c("Danio_rerio", "Gallus_gallus", "Homo_sapiens", "Mus_musculus", "Branchiostoma_lanceolatum")
SpeciesGRefs <- c("Danio_rerio.GRCz11", "Gallus_gallus.GRCg6a", "Homo_sapiens.GRCh38", "Mus_musculus.GRCm39", "Branchiostoma_lanceolatum.BraLan3", "Branchiostoma_floridae.Bfl_VNyyK")
SName1 <- SpeciesNames[which(ShortSpeciesNames==Species1)]
SName2 <- SpeciesNames[which(ShortSpeciesNames==Species2)]

BlanTissues <- c("cirri", "gills", "epidermis", "gut", "hepatic.diverticulum", "muscle", "neural.tube", "female.gonads", "male.gonads", "egg", "32cells", "Blastula", "MeanEmbr")
DrerTissues <- c("blastula", "embryo", "head", "tail", "granulocyte", "brain", "heart", "muscle tissue", "mesonephros", "intestine", "pharyngeal gill", "testis", "ovary", "bone element", "liver", "mature ovarian follicle", "zone of skin", "swim bladder", "head kidney", "spleen")
DrerTissuesIDs <- c("UBERON:0000307", "UBERON:0000922", "UBERON:0000033", "UBERON:0002415", "CL:0000094", "UBERON:0000955", "UBERON:0000948", "UBERON:0002385", "UBERON:0000080", "UBERON:0000160", "UBERON:0000206", "UBERON:0000473", "UBERON:0000992", "UBERON:0001474", "UBERON:0002107", "UBERON:0003982", "UBERON:0000014", "UBERON:0006860", "UBERON:0007132", "UBERON:0002106")

BlanMatchingTissues <- c("Blastula", "MeanEmbr", "male.gonads", "female.gonads", "muscle", "neural.tube", "gut", "gills", "hepatic.diverticulum", "epidermis")
DrerMatchingTissues <- c("blastula", "embryo", "testis", "ovary", "muscle tissue", "brain", "intestine", "pharyngeal gill", "liver", "zone of skin")
NamesMatchingTissues <- c("blastula", "embryo", "testis", "ovary", "muscle", "neural", "digestive", "gills", "hepatic", "skin")
colfunc <- colorRampPalette(c("forestgreen", "gold", "darkorange", "firebrick", "darkslateblue", "deepskyblue3"))
TissueColors <- colfunc(length(NamesMatchingTissues))

MatchingTissues <- as.data.frame(cbind(BlanMatchingTissues, DrerMatchingTissues, NamesMatchingTissues, TissueColors))
colnames(MatchingTissues) <- c("Blan", "Drer", "Name","Color")
MatchingTissues <- MatchingTissues[MatchingTissues$Name!="blastula" & MatchingTissues$Name!="ovary" & MatchingTissues$Name!="skin",]

VertTypes <- c("Missing", "SingleCopy", "Ohnolog", "Duplicated")
VertTypes.col <- c("royalblue2", ghibli_palettes$MarnieMedium2[c(2,4,7)])

###########################################################################
###########################################################################
# Read data

# List of onhologs
Ohnologs <- read.table(OhnologsSFile, h=F)[,1]

# Read Drer gene information previously processed or process it and write itinto a file
DrerGeneData <- prepare_DrerGeneData(DrerGeneDataFile, DrerBgeeDataFile, SName2, MatchingTissues, ResultsFolder, pwd)
print(head(DrerGeneData))

BlanGeneData <- prepare_BlanGeneData(BlanGeneDataFile, RNASeqMetadataFile, TPMFile)
print(head(BlanGeneData))

# Read gene pairs and the gene copy number in each species
head <- unlist(strsplit(system(paste0("head -1 ", NamesFile, " | sed 's/#//g'"), intern=T), "\t"))
columnS1 <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames==Species1)],".fa"))
columnS2 <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames==Species2)],".fa"))
columnDrer <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames=="Drer")],".fa"))
columnGgal <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames=="Ggal")],".fa"))
columnMmus <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames=="Mmus")],".fa"))
columnHsap <- which(head==paste0(SpeciesGRefs[which(ShortSpeciesNames=="Hsap")],".fa"))

system_out <- system(paste0("tail -n +2 ", NamesFile, " | awk -F '\t' '{print $1\"\t\"$",columnS1,"\"\t\"$",columnDrer,"\"\t\"$",columnGgal,"\"\t\"$",columnMmus,"\"\t\"$",columnHsap,"}' | awk -F '\t' '{split($2,a,\" \"); vcn=0; for(v = 3; v <= 6; ++v){split($v,b,\" \"); if(length(b)>vcn){vcn=length(b)}}; for(i in a){print $1\"\t\"a[i]\"\t\"length(a)\"\t\"vcn}}'"), intern=T)
GeneCN <- read.table(text=system_out, h=F, sep = "\t")
colnames(GeneCN) <- c("OG","Gene","CNBlan","CNVert")
GeneCN$TypeBlan <- c("SC", "D")[as.numeric(GeneCN$CNBlan>1)+1]
GeneCN$TypeVert <- c("SC", "D")[as.numeric(GeneCN$CNVert>1)+1]
GeneCN$TypeVert[which(GeneCN$CNVert==0)] <- "M"
GeneCN$TypeVert[which(GeneCN$OG %in% Ohnologs)] <- "O"

system_out <- system(paste0("tail -n +2 ", NamesFile, " | awk -F '\t' '{print $1\"\t\"$", columnS1, "\"\t\"$", columnS2,"}' | awk -F '\t' '{split($2,a,\" \");split($3,b,\" \"); for(i in a){for(j in b){print $1\"\t\"a[i]\"\t\"b[j]\"\t\"length(a)\"\t\"length(b)}}}'"), intern=T)
GenePairs <- read.table(text=system_out, h=F, sep = "\t")
colnames(GenePairs) <- c("OG","Gene1","Gene2","CopyNumber1","CopyNumber2")
GenePairs$Type1 <- c("SC", "D")[as.numeric(GenePairs$CopyNumber1>1)+1]
GenePairs$Type2 <- c("SC", "D")[as.numeric(GenePairs$CopyNumber2>1)+1]
GenePairs$Type2[which(GenePairs$OG %in% Ohnologs)] <- "O"
print(head(GenePairs))

OGProfile <- unique(GenePairs[,c("OG","Type1","Type2")])
for(t in c(1:length(MatchingTissues[,1]))){
	GenePairs[,paste0(MatchingTissues$Name[t],1)] <- BlanGeneData[match(GenePairs$Gene1, BlanGeneData$Gene),MatchingTissues$Blan[t]]>TPM.threshold
	GenePairs[,paste0(MatchingTissues$Name[t],2)] <- DrerGeneData[match(GenePairs$Gene2, DrerGeneData$Gene),paste0(MatchingTissues$Name[t],"Presence")]=="present"
	OGProfile[,paste0(MatchingTissues$Name[t],1)] <- unlist(lapply(OGProfile[,1], function(x){sum(GenePairs[which(GenePairs$OG==x),paste0(MatchingTissues$Name[t],1)])>0}))
	OGProfile[,paste0(MatchingTissues$Name[t],2)] <- unlist(lapply(OGProfile[,1], function(x){sum(GenePairs[which(GenePairs$OG==x),paste0(MatchingTissues$Name[t],2)])>0}))
}

GenePairs$BlanSum <- rowSums(GenePairs[,paste0(MatchingTissues$Name,1)])
GenePairs$DrerSum <- rowSums(GenePairs[,paste0(MatchingTissues$Name,2)])
GenePairs$DiffDom <- GenePairs$BlanSum-GenePairs$DrerSum
OGProfile$BlanSum <- rowSums(OGProfile[,paste0(MatchingTissues$Name,1)])
OGProfile$DrerSum <- rowSums(OGProfile[,paste0(MatchingTissues$Name,2)])
OGProfile$DiffDom <- OGProfile$BlanSum-OGProfile$DrerSum


###########################################################################
###########################################################################
# Plotting
pdf(paste0(ResultsFolder, "/ComparativeGeneExpression.pdf"), width=15, height=10)
par(mar=c(10,10,5,5),oma=c(1,1,1,1), yaxs='i', xaxs='i')
layout(matrix(c(1,2,3,4),nrow=1,ncol=4,byrow=T), widths=c(1,1,1,1), heights=c(2), TRUE)

BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneCN, "MeanAdult", "Mean adult TPM", c(0,100), VertTypes.col)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneCN, "MeanEmbr", "Mean embrionic TPM", c(0,100), VertTypes.col)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneCN, "TauTissues", "Tau among tissues", c(0,1), VertTypes.col)
BoxPlot_BlanTypes_VertTypes(BlanGeneData, GeneCN, "TauEmbAge", "Tau among developmental stages", c(0,1), VertTypes.col)

layout(matrix(c(1,2,3,4,5,6),nrow=2,ncol=3,byrow=T), widths=c(1), heights=c(1), TRUE)
ScatterPlot_ExpressDomanis(GenePairs[which(GenePairs$Type1=="SC" & GenePairs$Type2=="SC"),], OGProfile[which(OGProfile$Type1=="D" & OGProfile$Type2=="SC"),], MatchingTissues, "Number of B. lanceolatum domains", "Numbre of D. rerio domains", c(0,length(MatchingTissues[,1])))
ScatterPlot_ExpressDomanis(GenePairs[which(GenePairs$Type1=="D" & GenePairs$Type2=="SC"),], OGProfile[which(OGProfile$Type1=="D" & OGProfile$Type2=="SC"),], MatchingTissues, "Number of B. lanceolatum domains", "Numbre of D. rerio domains", c(0,length(MatchingTissues[,1])))
ScatterPlot_ExpressDomanis(GenePairs[which(GenePairs$Type1=="SC" & GenePairs$Type2=="D"),], OGProfile[which(OGProfile$Type1=="D" & OGProfile$Type2=="SC"),], MatchingTissues, "Number of B. lanceolatum domains", "Numbre of D. rerio domains", c(0,length(MatchingTissues[,1])))

Hist_ExpressDomanis(GenePairs, OGProfile, "SC", "D", MatchingTissues, "B. lanceolatum donains - D. rerio domains", "Relative number of\npairwise comparisons", "one-to-many small scale")
Hist_ExpressDomanis(GenePairs, OGProfile, "SC", "O", MatchingTissues, "B. lanceolatum donains - D. rerio domains", "Relative number of\npairwise comparisons", "one-to-many ohnologs")
Hist_ExpressDomanis(GenePairs, OGProfile, "D", "SC", MatchingTissues, "B. lanceolatum donains - D. rerio domains", "Relative number of\npairwise comparisons", "many-to-one small scale")

ExprProfileLines(GenePairs, OGProfile, MatchingTissues, "lab1", "lab2", VertTypes, VertTypes.col)

dev.off()


