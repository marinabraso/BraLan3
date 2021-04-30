# followed steps from prepare_GTF.R, Bgee dev RNAseq pipeline, a bit modified to our needs
# TODO verify with Marc if this makes sense 

# call packages
# source("https://bioconductor.org/biocLite.R")
# biocLite('chipseq')
# biocLite('GenomicFeatures')
library(chipseq)
library(GenomicFeatures)

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
gene_gtf_path='/Users/aechchik/amphio_test_data/Bla_annot_final_allfeatures.gtf'
# path & basename for the output gtf file
output_gtf_path='/Users/aechchik/amphio_test_data/Bla_intergenic'

# read gtf
gene_gtf <- as.matrix(read.table(gene_gtf_path, sep='\t', strip.white=TRUE, as.is=TRUE, colClasses="character", comment.char='#'))

# subset feature by type, aka exon and transcript 
# note: according to the official pipeline, would have subset the features 'gene', conceptually not different.
gene_gtf_exon <- gene_gtf[gene_gtf[,3]=="exon",]
gene_gtf <- gene_gtf[gene_gtf[,3]=="gene",]

# splitting the annotation field using single space "; " as a pattern
split_annotation_list <- strsplit(gene_gtf_exon[,9], "; ",  fixed =T)
# getting the vector of the gene IDs (1 for every exon)
# note: changed from original to fetch transcript_id field instead of gene_id
gene_ids <- sapply(split_annotation_list, function(x){ get_annot_value(x, 'gene_id') })

# getting the chromosome, start and the end of the gene
# For start, take the minimum exon start position
# For end, take the maximum exon end position
gene_start <- sapply(split(as.numeric(gene_gtf_exon[,4]), gene_ids), function(x){ sort(as.numeric(x))[1] })
gene_stop <- sapply(split(as.numeric(gene_gtf_exon[,5]), gene_ids), function(x){ rev(sort(as.numeric(x)))[1] })
gene_chr <- sapply(split(gene_gtf_exon[,1], gene_ids), function(x){ x[1] })
# chromosome/contig names from given gtf files
chromosomes <- unique(gene_gtf_exon[,1])
# removing patch contigs before selecting intergenic regions
chromosomes <- chromosomes[grep('PATCH', chromosomes, invert=TRUE, ignore.case=TRUE)]

# This object will include the coordinates of the selected intergenic regions
intergenic_regions <- matrix(ncol=3, nrow=0)
colnames(intergenic_regions) <- c("chr", "start", "end")

for(chr in chromosomes){
	cat(chr, " ")
    if(( sum(gene_chr==as.character(chr)) <= 1 )){ next }
	gene_IR <- IRanges(as.numeric(gene_start[gene_chr==as.character(chr)]), as.numeric(gene_stop[gene_chr==as.character(chr)]))
	inter_IR <- slice(coverage(gene_IR), lower=0, upper=0, rangesOnly=TRUE)
	inter_gene_data <- as.data.frame(cbind(inter_IR@start, end(inter_IR), inter_IR@width))
	colnames(inter_gene_data) <- c("start", "end", "width")
	if( sum(inter_gene_data[,3] > 2000) == 0 ){ next }
	inter_gene_data <- inter_gene_data[inter_gene_data[,3] > 2000,]
	inter_gene_data$center <- round((inter_gene_data[,1]+inter_gene_data[,2])/2)
	inter_gene_data$size <- inter_gene_data[,3] - 1000
	inter_gene_data$size[inter_gene_data$size > 20000] <- 20000 ## limit the size of usable regions to max 20000
	intergenic_regions <- rbind(intergenic_regions, cbind(chr, inter_gene_data$center - ceiling(inter_gene_data$size/2) + 1 , inter_gene_data$center + floor(inter_gene_data$size/2)))    
}

# preparing intergenic gtf data
intergenic_regions_gtf <- matrix(ncol=9, nrow=nrow(intergenic_regions))
intergenic_regions_gtf[,1] <- intergenic_regions[,1]
intergenic_regions_gtf[,2] <- "intergenic"
intergenic_regions_gtf[,3] <- "exon"
intergenic_regions_gtf[,c(4,5)] <- intergenic_regions[,c(2,3)]
intergenic_regions_gtf[,c(6,8)] <- "."
# Strand is chosen randomly
set.seed(12) ## setting seed for random number generator
intergenic_regions_gtf[,7] <- sample(c("+", "-"), nrow(intergenic_regions), replace =TRUE)

# intergenic_id - chr "_" start "_" stop
intergenic_id <- apply(intergenic_regions[,1:3], 1, function(x){ paste(x, collapse="_") })
gene_id <- paste("gene_id ", paste(intergenic_id, ";", sep=""), sep="")
transcript_id <- paste("transcript_id ", paste(intergenic_id, ";", sep=""), sep="")
intergenic_regions_gtf[,9] <- apply(cbind(gene_id, transcript_id), 1, function(x){ paste(x, collapse=" ") })

# write to output: will only include intergenic regions 
write.table(intergenic_regions_gtf, file=paste(output_gtf_path, ".gtf", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# a priori, can now run the kallisto index on this gtf once transformed to fasta
