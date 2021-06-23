#!/usr/bin/env Rscript

Sys.setenv(LANG = "en")

######################################################################
# Libraries & functions

`%out%` <- Negate(`%in%`)
script <- sub(".*=", "", commandArgs()[4])
#source(paste(substr(script,1, nchar(script)-2), "_functions.R", sep=""))
#library(BgeeCall)
#library(devtools)
library(BgeeDB)

######################################################################
# Files & folders
pwd <- system("pwd", intern=TRUE)
OrganOrthBgeeFile <- paste0(pwd, "/Data/AnatomicalOntology/Julien_AnatOrht_Bgee15.tsv")

######################################################################
# General parameters
args <- commandArgs(trailingOnly=TRUE)
ResultsFolder <- args[1]
setwd(paste0("./", ResultsFolder))
Species1 <- args[2]
Species2 <- args[3]

ShortSpeciesNames <- c("Drer", "Ggal", "Hsap", "Mmus", "Blan")
SpeciesNames <- c("Danio_rerio", "Gallus_gallus", "Homo_sapiens", "Mus_musculus", "Branchiostoma_lanceolatum")

SName1 <- SpeciesNames[which(ShortSpeciesNames==Species1)]
SName2 <- SpeciesNames[which(ShortSpeciesNames==Species2)]

###########################################################################
###########################################################################
# Read data

BgeeS1 <- Bgee$new(species=SName1, dataType="rna_seq")
BgeeS2 <- Bgee$new(species=SName2, dataType="rna_seq")
aBgeeS1 <- getAnnotation(BgeeS1)
aBgeeS2 <- getAnnotation(BgeeS2)
dBgeeS1 <- getData(BgeeS1) # SRP123447
dBgeeS2 <- getData(BgeeS2, experimentId= "GSE30352") # The evolution of gene expression levels in mammalian organs


system_out <- system(paste0("cut -f1,2,3,4 ", OrganOrthBgeeFile), intern=T)
OrganOrthBgee <- read.table(text=system_out, h=T, sep = "\t")
colnames(OrganOrthBgee) <- c("Anat.IDs", "Anat.Names", "Source.Anat.IDs", "Source.Anat.Names")

OrganOrthBgee[which(OrganOrthBgee$Anat.IDs %in% unique(dBgeeS1$Anatomical.entity.ID) & OrganOrthBgee$Source.Anat.IDs %in% unique(dBgeeS2$Anatomical.entity.ID)),]


















