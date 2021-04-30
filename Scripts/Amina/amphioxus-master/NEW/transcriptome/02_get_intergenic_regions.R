# followed steps from prepare_GTF.R, Bgee dev RNAseq pipeline, a bit modified to our needs
# TODO verify with Marc if this makes sense

# call packages

setwd('/Users/aechchik/Desktop/amphio_0520/')

#BiocManager::install("IRanges")
#BiocManager::install(c("GenomicFeatures", "chipseq"))
library(IRanges)
library(GenomicFeatures)

library(chipseq)
# custom functions
# get_annot_value
# Function for obtaining the part of the annotation field from gtf file
# Input: whole field, already splitted, for example: gene_id "FBgn0264003"; gene_name "mir-5613"; gene_source "FlyBase"; gene_biotype "pre_miRNA"
# Output: could be for example gene_id: FBgn0264003
get_annot_value <- function(split_annotation, field_name){
	## find the right field
	field_all <- split_annotation[grep(field_name, split_annotation, fixed=T)];
	## split the field
	field_value <- strsplit(field_all, ' ', fixed=T)[[1]][2];
	## remove the last ';' if necessairy
    field_value <- sub(';', '', field_value,fixed=T)
    return(field_value)
}

# user-defined variables
# take the subset of the transcriptome assembly: cufflinks on rna-seq mapping to PacBio genome by Kamil. Illumina genome annotation is weird. they just had the exon and CDs features, aka not consistent with the cufflinks standard output. e-mailed Ferdi, not helping.
# path & fullname for the gtf file
gene_gtf_path='/Users/aechchik/Desktop/amphio_0520/final_final.gtf'
# path & basename for the output gtf file
output_gtf_path='/Users/aechchik/Desktop/amphio_0520/final_intergenic'
output_genic_gtf_path='/Users/aechchik/Desktop/amphio_0520/final_genic'


# read gtf
gene_gtf <- as.matrix(read.table(gene_gtf_path, sep='\t', strip.white=TRUE, as.is=TRUE, colClasses="character", comment.char='#'))

## selecting exon lines
gene_gtf_exon <- gene_gtf[gene_gtf[,3]=="exon",]
## selecting gene lines
gene_gtf <- gene_gtf[gene_gtf[,3]=="gene",]

##  selecting only genes from assembled chromosome sequence (potentially useful code to select only fully asembled genome sequence)
##  chromosomes_all <- c(1:100, "X", "Y", "MT")
##  gene_gtf_exon <- gene_gtf_exon[gene_gtf_exon[,1] %in% chromosomes_all,]

## getting gene start, stop, chromosome, strand, biotype for each gene (from the exon GTF)
## GTF format reports 1-based coordinates (start and end)
cat("Extracting gene informations...\n")
## splitting the annotation field using single space "; " as a pattern
split_annotation_list <- strsplit(gene_gtf_exon[,9], "; ",  fixed =T)

## getting the vector of the gene IDs (1 for every exon)
gene_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') })

## getting the vector of the transcript IDs (1 for every exon)
transcript_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'transcript_id') })

## getting the table with mappings between transcript and gene IDs (for export)
gene_transcript_ids <- unique(cbind(gene_ids, transcript_ids), MARGIN=1)

## getting the vector of gene_biotypes: splitting gene gtf
split_annotation_list <- strsplit(gene_gtf[,9], "; ",  fixed =T)
gene_biotypes <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_biotype') })
names(gene_biotypes) <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') })

## getting the chromosome, start and the end of the gene
## For start, take the minimum exon start position
## For end, take the maximum exon end position
gene_start <- sapply(split(as.numeric(gene_gtf_exon[,4]), gene_ids), function(x){ sort(as.numeric(x))[1] })
gene_stop <- sapply(split(as.numeric(gene_gtf_exon[,5]), gene_ids), function(x){ rev(sort(as.numeric(x)))[1] })
gene_chr <- sapply(split(gene_gtf_exon[,1], gene_ids), function(x){ x[1] })

## chromosome/contig names from given gtf files
chromosomes <- unique(gene_gtf_exon[,1])
## removing patch contigs before selecting intergenic regions
chromosomes <- chromosomes[grep('PATCH', chromosomes, invert=TRUE, ignore.case=TRUE)]

## This object will include the coordinates of the selected intergenic regions
intergenic_regions <- matrix(ncol=3, nrow=0)
colnames(intergenic_regions) <- c("chr", "start", "end")

for(chr in chromosomes){
    cat(chr, " ")
    ## keeping genes from selected chromosome
    ## skip chromosome/contig if 1 or less gene
    if(( sum(gene_chr==as.character(chr)) <= 1 )){ next }

    ## constructing coverage map of the chromosomes using start and stop coordinates of the genes and selecting regions with 0 coverage (intergenic)
    gene_IR <- IRanges(as.numeric(gene_start[gene_chr==as.character(chr)]), as.numeric(gene_stop[gene_chr==as.character(chr)]))
    inter_IR <- slice(coverage(gene_IR), lower=0, upper=0, rangesOnly=TRUE)
    inter_gene_data <- as.data.frame(cbind(inter_IR@start, end(inter_IR), inter_IR@width))
    colnames(inter_gene_data) <- c("start", "end", "width")

    ## if selected chromosome/contig has only no intergenic regions with length min 2000nt then skip this contig
    if( sum(inter_gene_data[,3] > 2000) == 0 ){ next }

    ## restrict to regions larger than 2000nt
    inter_gene_data <- inter_gene_data[inter_gene_data[,3] > 2000,]

    ## finding the center of the intergenic region
    inter_gene_data$center <- round((inter_gene_data[,1]+inter_gene_data[,2])/2)

    ## getting the size of "usable" intergenic regions around the center point
    ## width - 500nt around each flank
    inter_gene_data$size <- inter_gene_data[,3] - 1000
    inter_gene_data$size[inter_gene_data$size > 20000] <- 20000 ## limit the size of usable regions to max 20000

    ## storing information (chr, start, stop, size) for selected intergenic regions on this chromosome
    intergenic_regions <- rbind(intergenic_regions, cbind(chr, inter_gene_data$center - ceiling(inter_gene_data$size/2) + 1 , inter_gene_data$center + floor(inter_gene_data$size/2)))
}

## preparing intergenic gtf data
cat("\nPreparing intergenic GTF data...\n")
intergenic_regions_gtf <- matrix(ncol=9, nrow=nrow(intergenic_regions))
intergenic_regions_gtf[,1] <- intergenic_regions[,1]
intergenic_regions_gtf[,2] <- "intergenic"
intergenic_regions_gtf[,3] <- "exon"
intergenic_regions_gtf[,c(4,5)] <- intergenic_regions[,c(2,3)]
intergenic_regions_gtf[,c(6,8)] <- "."
## Strand is chosen randomly
set.seed(12) ## setting seed for random number generator
intergenic_regions_gtf[,7] <- sample(c("+", "-"), nrow(intergenic_regions), replace =TRUE)

## intergenic_id - chr "_" start "_" stop
intergenic_id <- apply(intergenic_regions[,1:3], 1, function(x){ paste(x, collapse="_") })
gene_id <- paste("gene_id ", paste(intergenic_id, ";", sep=""), sep="")
transcript_id <- paste("transcript_id ", paste(intergenic_id, ";", sep=""), sep="")
intergenic_regions_gtf[,9] <- apply(cbind(gene_id, transcript_id), 1, function(x){ paste(x, collapse=" ") })

# write to output: will only include intergenic regions

# a priori, can now run the kallisto index on this gtf once transformed to fasta


write.table(intergenic_regions_gtf, file=paste(output_intergenic_gtf_path, ".gtf", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(gene_gtf_exon, file=paste(output_genic_gtf_path, ".gtf", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)


## Table with mappings between transcript and gene IDs
write.table(gene_transcript_ids, file=paste(output_gtf_path, ".gene2transcript", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Table with biotype for all gene IDs
write.table(gene_biotypes, file=paste(output_gtf_path, ".gene2biotype", sep=""), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)



## GTF file with both genic exons and intergenic regions
write.table(rbind(gene_gtf_exon, intergenic_regions_gtf), file=paste(output_gtf_path, ".gtf_all", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Table with mappings between transcript and gene IDs
write.table(gene_transcript_ids, file=paste(output_gtf_path, ".gene2transcript", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Table with biotype for all gene IDs
write.table(gene_biotypes, file=paste(output_gtf_path, ".gene2biotype", sep=""), sep="\t", row.names=TRUE, col.names=FALSE, quote=FALSE)




## run kallisto on the fasta from gtf: intergenic+genic+small scaffolds(pretending it's the transcripts)
# e.g. in: /scratch/axiom/FAC/FBM/DEE/mrobinso/default/aechchik/amphioxus/test_bl6/kallisto_interg-scaffolds


kallisto_count_folder = "/Users/aechchik/Desktop/amphioxus_wrk"
gene2transcript_file = "final_intergenic.gene2transcript"
    gene2biotype_file = "final_intergenic.gene2biotype"

    ## Session info
    print(sessionInfo())


## checking if all necessary arguments were passed in command line
command_arg <- c("kallisto_count_folder", "gene2transcript_file", "gene2biotype_file", "library_id")
for( c_arg in command_arg ){
    if( !exists(c_arg) ){
        stop( paste(c_arg,"command line argument not provided\n") )
    }
}

## reading kallisto's output. If file not exists, script stops
kallisto_count_file <- paste0(kallisto_count_folder, "/abundance.tsv")
if( file.exists(kallisto_count_file) ){
    kallisto_count <- read.table(kallisto_count_file, h=T, sep="\t")
} else {
    stop( paste("Kallisto results file not found [", kallisto_count_file, "]\n"))
}


## reading gene to transcript file. If file not exists, script stops
if( file.exists(gene2transcript_file) ){
    gene2transcript <- read.table(gene2transcript_file, h=F, sep="\t")
    names(gene2transcript) <- c("gene_id", "transcript_id")
} else {
    stop( paste("Gene to transcript file not found [", gene2transcript_file, "]\n"))
}

## reading gene biotype file. If file not exists, script stops
if( file.exists(gene2biotype_file) ){
    gene2biotype <- read.table(gene2biotype_file, h=F, sep="\t")
    names(gene2biotype) <- c("gene_id", "biotype")
    gene2biotype <- gene2biotype[order(gene2biotype$gene_id), ] ## order by gene Id
} else {
    stop( paste("Gene biotype file not found [", gene2biotype_file, "]\n"))
}

###############################################################################
## Add gene ids to the kallisto output. Effectively, this removes all intergenic regions,
genic_count <- merge(kallisto_count, gene2transcript, by.x=1, by.y=2)[, c(1,6,2,3,4,5)]
names(genic_count)[2] <- "gene_id"
## Add biotype for all genes
genic_count <- merge(genic_count, gene2biotype, by.x=2, by.y=1)[, c(2,1,7,3,4,5,6)]
## Resort the table by transcript_id
genic_count <- genic_count[order(genic_count$target_id), ]

## Add gene id column to kallisto's output (genic regions only)
kallisto_count$gene_id <- rep(NA, times=length(kallisto_count$target_id))
kallisto_count <- kallisto_count[order(kallisto_count$target_id),] ## sort by transcript_id
## Check transcripts are matching
## summary(kallisto_count$target_id[kallisto_count$target_id %in% genic_count$target_id] == genic_count$target_id)
kallisto_count$gene_id[kallisto_count$target_id %in% genic_count$target_id] <- as.character(genic_count$gene_id)
## Add region type
kallisto_count$type <- rep("intergenic", times=length(kallisto_count$target_id))
kallisto_count$type[kallisto_count$target_id %in% genic_count$target_id] <- "genic"
## Add biotype for genic regions
kallisto_count$biotype <- rep(NA, times=length(kallisto_count$target_id))
kallisto_count$biotype[kallisto_count$target_id %in% genic_count$target_id] <- as.character(genic_count$biotype)

## Calculate FPKMs for all transcripts + intergenic regions, from their TPM values
## TPM = RPKM * 10^6 / sum(RPKM)
fpkmToTpm <- function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
## RPKM = TPM * (sum(rg / flg)/ R) * 10^3
##      = TPM * (sum(estimated counts/effective length) / sum(estimated counts)) * 10^3
tpmToFpkm <- function(tpm, counts, effLen){
    exp(log(tpm) + log(sum(counts/effLen)) - log(sum(counts)) + log(1e3))
}
kallisto_count$fpkm <- tpmToFpkm(kallisto_count$tpm, kallisto_count$est_counts, kallisto_count$eff_length)
## reorder columns and rows before exporting
kallisto_count <- kallisto_count[order(kallisto_count$gene_id), c(1, 6, 2:5, 9, 7, 8)]
## Saving to Kallisto output folder
write.table(kallisto_count, file = paste0(kallisto_count_folder, "/abundance+gene_id+fpkm+intergenic.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

###############################################################################
## All reads mapped to both genic and intergenic regions are used for the total number of reads, and for the length of the transcriptome. This is the only way to have comparable FPKM/TPM values between genic and intergenic regions, for example to be able to compute a cutoff for presence/absence of expression. But for users, these TPMs or FPKMs are not correct. It is easy to recalculate them using only genic regions!

## Recalculate TPMs and calculate FPKMs, using functions from https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
countToTpm <- function(counts, effLen){
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}
countToFpkm <- function(counts, effLen){
    N <- sum(counts)
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
genic_count$tpm <- countToTpm(genic_count$est_counts, genic_count$eff_length)
genic_count$fpkm <- countToFpkm(genic_count$est_counts, genic_count$eff_length)
## reorder columns and rows before exporting
genic_count <- genic_count[order(genic_count$gene_id), c(1:2, 4:8, 3)]
## Saving to Kallisto output folder
write.table(genic_count, file = paste0(kallisto_count_folder, "/abundance+gene_id+new_genic_tpm+new_genic_fpkm.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

###############################################################################
## Gene-level expression
## Sum TPMs, FPKMs and counts of all transcripts of each gene
gene_count <- aggregate(genic_count[,5:7], list(genic_count$gene_id), sum)
names(gene_count)[1] <- "gene_id"
## These are the values that will be inserted into the database
## Add biotype for all genes. First, check:
## summary(gene_count$gene_id == gene2biotype$gene_id)
gene_count$biotype <- gene2biotype$biotype

## Gene-level expression for kallisto's output:
## Nothing is summed for intergenic regions
kallisto_gene_count <- kallisto_count[kallisto_count$type == "intergenic", c(1, 5, 6, 7, 8, 9)]
names(kallisto_gene_count)[1] <- "gene_id"
## For genic regions sum read counts, TPMs and FPKMs
temp <- kallisto_count[kallisto_count$type == "genic", c(2, 5, 6, 7)]
temp <- aggregate(temp[,2:4], list(temp$gene_id), sum)
names(temp)[1] <- "gene_id"
temp$type <-  rep("genic", times=length(temp$gene_id))
## summary(temp$gene_id == gene2biotype$gene_id)
temp$biotype <- gene2biotype$biotype
## Make final table with both genic and intergenic regions
kallisto_gene_count <- rbind(temp, kallisto_gene_count)
## This will be used for presence / absence cutoff calculation

## Export gene level files
write.table(gene_count, file = paste0(kallisto_count_folder, "/abundance_gene_level+new_tpm+new_fpkm.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(kallisto_gene_count, file = paste0(kallisto_count_folder, "/abundance_gene_level+fpkm+intergenic.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


###############################################################################
## Plotting of the distribution of TPMs for genic and intergenic regions
pdf(file = paste0(kallisto_count_folder, "/distribution_TPM_genic_intergenic.pdf"), width = 6, height = 5)
par(mar=c(5,6,1,1)) ## bottom, left, top and right margins

dens <- density(log2(na.omit(kallisto_gene_count$tpm) + 10^-6))
## Subgroups densities. Visualization trick: we add an invisible set of points at x=-30, to make densities comparable
## genic regions
dens_genic <- density(c(rep(-30, times=sum(kallisto_gene_count$type != "genic")), log2(kallisto_gene_count$tpm[kallisto_gene_count$type == "genic"] + 10^-6)))
## protein-coding genes only (had to take care of NAs strange behavior)
dens_coding <- density(c(rep(-30, times=sum(!kallisto_gene_count$biotype %in% "protein_coding")), log2(kallisto_gene_count$tpm[kallisto_gene_count$biotype %in% "protein_coding"] + 10^-6)))
## intergenic
dens_intergenic <- density(c(rep(-30, times=sum(kallisto_gene_count$type != "intergenic")), log2(kallisto_gene_count$tpm[kallisto_gene_count$type == "intergenic"] + 10^-6)))

## Plot whole distribution
plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-23, 21), lwd=2, main="", bty="n", axes=F, xlab="")
axis(2, las=1)
## Add 2 x-axes: TPMs and RPKMs
## See http://stackoverflow.com/questions/8443820/r-multiple-x-axis-with-annotations
axis(1, at=seq(-30 , 30, by=10), line=0, mgp = c(3, 0.5, 0), cex.axis=0.8)
mtext(expression(log[2]('TPM'+10^-6)), 1,  adj = 1, padj = 0, line=0.2, at=par("usr")[1], col="black", cex=0.8)
## To make FPKM scale, we need to know what log2(TPM + 10^-6) value corresponds to any log2(FPKM + 10^-6) value. We know that for any given gene, TPMg/FPKMg = coef
##    log2(FPKM + 10^-6) = x
## <=>              FPKM = exp(x*log(2)) - 10^-6
## <=>               TPM = coef * (exp(x*log(2)) - 10^-6)
## <=> log2(TPM + 10^-6) = log2( coef * (exp(x*log(2)) - 10^-6) + 10^-6)
coef <- na.omit(kallisto_gene_count$tpm / kallisto_gene_count$fpkm)[1]
## We generate scale from -20 to 100 FPKMs
axis(1, at=log2( coef * (exp(seq(-20 , 100, by=10)*log(2)) - 10^-6) + 10^-6), labels=seq(-20 , 100, by=10), line=2, mgp = c(3, 0.5, 0), cex.axis=0.8)
mtext(expression(log[2]('FPKM'+10^-6)), 1,  adj = 1, padj = 0, line=2.2, at=par("usr")[1], col="black", cex=0.8)

## Add subgroups distributions (genic, intergenic, etc):
## genic
lines(dens_genic, col="firebrick3", lwd=2)
## protein-coding genes
lines(dens_coding, col="firebrick3", lwd=2, lty=2)
## intergenic
lines(dens_intergenic, col="dodgerblue3", lwd=2)

## legend
legend("topright", c("all", "genic", "protein-coding genes", "intergenic"), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
dev.off()


file <- "/abundance_gene_level+fpkm+intergenic.tsv"

kallisto_gene_counts <- read.table("abundance_gene_level+fpkm+intergenic.tsv", h=T, sep="\t")

summed <- kallisto_gene_counts

dens <- density(log2(na.omit(summed$tpm) + 10^-6))
## Subgroups densities. Visualization trick: we add an invisible set of points at x=-30, to make densities comparable
## genic regions
dens_genic <- density(c(rep(-30, times=sum(summed$type != "genic")), log2(summed$tpm[summed$type == "genic"] + 10^-6)))
## protein-coding genes only (had to take care of NAs strange behavior)
dens_coding <- density(c(rep(-30, times=sum(!summed$biotype %in% "protein_coding")), log2(summed$tpm[summed$biotype %in% "protein_coding"] + 10^-6)))
## intergenic
dens_intergenic <- density(c(rep(-30, times=sum(summed$type != "intergenic")), log2(summed$tpm[summed$type == "intergenic"] + 10^-6)))
## Plot whole distribution
plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-23, 21), lwd=2, main=paste0(as.character(unique(sampleInfo$organism[sampleInfo$speciesId == species])), " (", numLibs, " libraries)"), bty="n", axes=T, xlab="log2(TPM + 10^-6)")
## Add subgroups distributions (genic, intergenic, etc):
## genic
lines(dens_genic, col="firebrick3", lwd=2)
## protein-coding genes
lines(dens_coding, col="firebrick3", lwd=2, lty=2)
## intergenic
lines(dens_intergenic, col="dodgerblue3", lwd=2)
## legend
legend("topleft", c(paste0("all (", length(summed[,1]),")"), paste0("genic (", sum(summed$type == "genic"), ")"), paste0("coding (", sum(summed$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
dev.off()

subgroup.densitie


