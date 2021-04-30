setwd("/Users/aechchik/Desktop/amphio_0520/abundances/")

args <- commandArgs(trailing = TRUE)

# args call & read

title <- args[1]

outpath <- paste0(title, ".tsv_dir")

sum_by_species <- read.table("sum_abundance_gene_level+intergenic+classification.tsv", h=T)

# kallisto_gene <- args[2] # read.table("SRR6246037.tsv_dir/abundance_gene_level+fpkm+intergenic.tsv", h=T) 
kallisto_gene_counts <- read.table(paste0(outpath, "/abundance_gene_level+fpkm+intergenic.tsv"), h=T)

# gene <-  args[3] #read.table("SRR6246037.tsv_dir/abundance_gene_level+new_tpm+new_fpkm.tsv", h=T)
gene_counts <- read.table(paste0(outpath, "/abundance_gene_level+new_tpm+new_fpkm.tsv"), h=T)


intergenic_choice <-  read.table("../intergenic1.ids", h=F)

pdf(file = paste0(outpath, "/", title, "_distribution_TPM_genic_intergenic+cutoff.pdf"), width = 6, height = 5)


# funtions from bgeecall

plot_distributions <- function(counts, selected_coding, selected_intergenic, cutoff, title){
    ## Plotting of the distribution of TPMs for coding and intergenic regions + cutoff
    ## Note: this code is largely similar to plotting section in rna_seq_analysis.R
    
    par(mar=c(5,6,1,1)) ## bottom, left, top and right margins
    dens <- density(log2(na.omit(counts$tpm) + 10^-6))
    
    ## protein-coding genes only (had to take care of NAs strange behavior)
    dens_coding <- density(log2(counts$tpm[selected_coding] + 10^-6))
    ## Normalize density for number of observations
    dens_coding$y <- dens_coding$y * sum(selected_coding) / length(counts$tpm)
    
    ## intergenic
    dens_intergenic <- density(log2(counts$tpm[selected_intergenic] + 10^-6))
    dens_intergenic$y <- dens_intergenic$y * sum(selected_intergenic) / length(counts$tpm)
    
    ## Plot whole distribution
    plot(dens, ylim=c(0, max(dens$y)*1.1), xlim=c(-23, 21), lwd=2, main= title, bty="n", axes=F, xlab="")
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
    ## We generate scale from -20 to 100 FPKMs

    ## Plot the TPM cutoff
    ## abline(v=cutoff, col="gray", lty=1, lwd=2)
    arrows(log2(cutoff + 10e-6), par("usr")[3], log2(cutoff + 10e-6), par("usr")[4]/2, col="gray", lty=1, lwd=2, angle=160, length=0.1)
    
    ## Add subgroups distributions (coding, intergenic, etc):
    ## protein-coding genes
    lines(dens_coding, col="firebrick3", lwd=2, lty=2)
    ## intergenic
    lines(dens_intergenic, col="dodgerblue3", lwd=2, lty=2)
    
    ## legend
    legend("topright", c("all", "selected protein-coding genes", "selected intergenic regions"), lwd=2, col=c("black", "firebrick3", "dodgerblue3"), lty=c(1, 2, 2), bty="n")
    return();
}


calculate_and_plot_r <- function(counts, selected_coding, selected_intergenic, desired_r_cutoff, title){
    ## r = (number of intergenic regions with TPM values higher than x * number of coding regions) /
    ##     (number of coding regions with TPM values higher than x * number of intergenic regions)
    ##   = 0.05
    ## What is value of x (cutoff)? calculate the distribution of r for a range of TPMs, then select the closest value to 0.05
    
    ## Counting how many intergenic regions have equal or higher value of TPM for every value of TPM
    ## For each gene's TPM (sorted), calculate r
    
    ##   r <- sapply(sort(unique(counts$tpm[selected_coding])), function(x){
    ##     return(
    ##            ( sum(counts$tpm[selected_intergenic] >= x) / sum(selected_intergenic) ) /
    ##            ( sum(counts$tpm[selected_coding] >= x) / sum(selected_coding) )
    ##            )
    ##   })
    ##   ## This is too long! Takes >30s for each library
    
    ## For each unique value of coding TPM, let's calculate the number of intergenic TPMs that are larger:
    ## ptm <- proc.time()
    summed_intergenic <- sapply(unique(sort(counts$tpm[selected_coding])), function(x){
        return( sum(counts$tpm[selected_intergenic] >= x) )
    })
    ## It is not necessary to do the same for coding regions: for a sorted vector the number of greater elements is equal to lenght(vector) - position + 1. Here, it is a bit trickier since we do not consider all coding TPM values, but the unique ones, so we use the rle function to know the lengths of the runs of similar values and sum them
    summed_coding <- c(0, cumsum(rle(sort(counts$tpm[selected_coding]))$lengths))
    summed_coding <- summed_coding[-(length(summed_coding))]
    summed_coding <- sum(selected_coding) - summed_coding
    
    ## Now we can calculate r
    r <- ( summed_intergenic / sum(selected_intergenic) ) /
        ( summed_coding / sum(selected_coding) )
    ## This is twice faster as code above!
    
    percent <- (1-desired_r_cutoff)*100
    
    ## Select the minimal value of TPM for which the ratio of genes and intergenic regions is equal to 0.05 or lower (first test if at least 1 TPM value has this property):
    if (sum(r < desired_r_cutoff) == 0){
        TPM_cutoff <- sort(unique(counts$tpm[selected_coding]))[which(r == min(r))[1]]
        r_cutoff <- min(r)
        cat(paste0("    There is no TPM cutoff for which " , percent,"%", " of the expressed genes would be coding. TPM cutoff is fixed at the first value with maximum coding/intergenic ratio. r=", r_cutoff, " at TPM=", TPM_cutoff,"\n"))
    } else {
        TPM_cutoff <- sort(unique(counts$tpm[selected_coding]))[which(r < desired_r_cutoff)[1]]
        r_cutoff <- desired_r_cutoff
    }
    
    plot(log2(sort(unique(counts$tpm[selected_coding]))+10e-6), r, pch=16, xlab="log2(TPM + 10^-6)", type="l")
    abline(h=desired_r_cutoff, lty=2, col="gray")
    if (r_cutoff > desired_r_cutoff){
        abline(h=desired_r_cutoff, lty=3, col="gray")
    }
    arrows(log2(TPM_cutoff + 10e-6), par("usr")[3], log2(TPM_cutoff + 10e-6), par("usr")[4]/2, col="gray", lty=1, lwd=2, angle=160, length=0.1)
    return(c(TPM_cutoff, r_cutoff))
}


cutoff_info <- function(title, counts, column, max_intergenic, TPM_cutoff, r_cutoff){
    ## Calculate summary statistics to export in cutoff info file
    genic_present <- sum(counts[[column]][counts$type == "genic"] == "present")/sum(counts$type == "genic") * 100
    number_genic_present <- sum(counts[[column]][counts$type == "genic"] == "present")
    
    coding_present <- sum(counts[[column]][counts$biotype %in% "protein_coding"] == "present")/sum(counts$biotype %in% "protein_coding") * 100
    number_coding_present <- sum(counts[[column]][counts$biotype %in% "protein_coding"] == "present")
    
    intergenic_present <- sum(counts[[column]][counts$type == "intergenic"] == "present")/sum(counts$type == "intergenic") * 100
    number_intergenic_present <- sum(counts[[column]][counts$type == "intergenic"] == "present")
    
    ## Export cutoff_info_file
    to_export <- c(title,
                   max_intergenic,
                   TPM_cutoff,
                   genic_present, number_genic_present, sum(counts$type == "genic"),
                   coding_present,  number_coding_present, sum(counts$biotype %in% "protein_coding"),
                   intergenic_present, number_intergenic_present, sum(counts$type == "intergenic"),
                   r_cutoff
    )
    names(to_export) <- c("libraryId",
                          "max_intergenic",
                          "cutoffTPM",
                          "proportionGenicPresent", "numberGenicPresent", "numberGenic",
                          "proportionCodingPresent", "numberPresentCoding", "numberCoding",
                          "proportionIntergenicPresent", "numberIntergenicPresent", "numberIntergenic",
                          "ratioIntergenicCodingPresent")
    return(to_export)
}


selected_coding <- kallisto_gene_counts$biotype %in% "protein_coding"
max_intergenic <- max(sum_by_species$tpm[sum_by_species$classification %in% paste0("intergenic_", 1) ])
selected_intergenic <- kallisto_gene_counts$gene_id %in% intergenic_choice$V1

results <- calculate_and_plot_r(kallisto_gene_counts, selected_coding, selected_intergenic, 0.05, title)
TPM_cutoff <- results[1]
r_cutoff <- results[2]

plot_distributions(kallisto_gene_counts, selected_coding, selected_intergenic, TPM_cutoff, title)


cutoff_info_file <- cutoff_info(title, kallisto_gene_counts, "call", max_intergenic, TPM_cutoff, r_cutoff)

dev.off()


write.table(kallisto_gene_counts,
            file = paste0(outpath, "/abundance_gene_level+fpkm+intergenic+calls.tsv"),
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(t(t(cutoff_info_file)),
            file = paste0(outpath, "/cutoff_info_file.tsv"),
            quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE)
