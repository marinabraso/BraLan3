# rnaseq analysis: to run for each kallisto abundance file 
setwd("/Users/aechchik/Desktop/amphio_0520/")
# import biotype and gene2transcripts
gene2biotype <- read.table("final_pc_genic_intergenic.gene2biotype", h=F)
gene2biotype[is.na(gene2biotype)] <- "protein_coding"
colnames(gene2biotype)[2]  <- "biotype"
gene2transcript <- read.table("final_pc_genic_intergenic.gene2transcript", h=F)

setwd("./abundances/")
args <- commandArgs(trailing = TRUE)
fn <- args[1]
kallisto_count <- read.table(fn, h=T) # replace with: args[1]
outpath <- paste0(fn, "_dir")

# add genes to kallisto output
# rm intergenic regions too
genic_count <- merge(kallisto_count, gene2transcript, by.x=1, by.y=2)[, c(1,6,2,3,4,5)]
names(genic_count)[2] <- "gene_id"
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

write.table(kallisto_count, file = paste0(outpath, "/abundance+gene_id+fpkm+intergenic.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
# replace with  # replace with: args[2]


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
write.table(genic_count, file=paste0(outpath, "/abundance+gene_id+new_genic_tpm+new_genic_fpkm.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
# replace with args[3]


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
write.table(gene_count, file=paste0(outpath, "/abundance_gene_level+new_tpm+new_fpkm.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(kallisto_gene_count, file=paste0(outpath, "/abundance_gene_level+fpkm+intergenic.tsv"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


#plot distribution

pdf(file = paste0(outpath, "/distribution_TPM_genic_intergenic.pdf"), width = 6, height = 5)
dens <- density(log2(na.omit(kallisto_gene_count$tpm) + 10^-6))
## Subgroups densities. Visualization trick: we add an invisible set of points at x=-30, to make densities comparable
## genic regions
dens_genic <- density(c(rep(-30, times=sum(kallisto_gene_count$type != "genic")), log2(kallisto_gene_count$tpm[kallisto_gene_count$type == "genic"] + 10^-6)))
## protein-coding genes only (had to take care of NAs strange behavior)
dens_coding <- density(c(rep(-30, times=sum(!kallisto_gene_count$biotype %in% "protein_coding")), log2(kallisto_gene_count$tpm[kallisto_gene_count$biotype %in% "protein_coding"] + 10^-6)))
## intergenic
dens_intergenic <- density(c(rep(-30, times=sum(kallisto_gene_count$type != "intergenic")), log2(kallisto_gene_count$tpm[kallisto_gene_count$type == "intergenic"] + 10^-6)))

## Plot whole distribution
plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-23, 21), lwd=2, main="", bty="n", axes=F, xlab="log2(TPM+10^-6)")
axis(2, las=1)
## Add 2 x-axes: TPMs and RPKMs
## See http://stackoverflow.com/questions/8443820/r-multiple-x-axis-with-annotations
axis(1, at=seq(-30 , 30, by=10), mgp = c(3, 0.5, 0), cex.axis=0.8)
## To make FPKM scale, we need to know what log2(TPM + 10^-6) value corresponds to any log2(FPKM + 10^-6) value. We know that for any given gene, TPMg/FPKMg = coef
##    log2(FPKM + 10^-6) = x
## <=>              FPKM = exp(x*log(2)) - 10^-6
## <=>               TPM = coef * (exp(x*log(2)) - 10^-6)
## <=> log2(TPM + 10^-6) = log2( coef * (exp(x*log(2)) - 10^-6) + 10^-6)
coef <- na.omit(kallisto_gene_count$tpm / kallisto_gene_count$fpkm)[1]
## We generate scale from -20 to 100 FPKMs
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

