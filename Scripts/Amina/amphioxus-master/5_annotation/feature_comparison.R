ex_gff <- read.table('PB_BL00003_evm0.gff', header = F, sep = '\t')
cf_gtf <- read.table('IL_BL00003_evm0.gtf', header = F, sep = '\t')

ex_gff$V3 #features
ex_gff$V4 #from
ex_gff$V5 #to
# for cf_gtf - the same

# total length of transctipt annotation on reference
max(ex_gff$V5) - min(ex_gff$V4)
max(cf_gtf$V5) - min(cf_gtf$V4)
# wow, it is very different, I guess I would like to record this statistics.

# draw a vertical line; x coordinates are give by from/to gtf entries (4,5)
# vertical position is given by factor "features"
draw_feature <- function(gtfline){
  lines(gtfline[4:5], rep(gtfline[3],2), lwd = 5)
}

plot_features <- function(gtf){
  labels <- levels(gtf$V3)
  # conversion to matrix is not very successful for factors
  gtf$V3 <- as.numeric(gtf$V3)
  # for every letter in the longest label the plot will be wider by 0.1
  par(mar=c(5, 4 + (max(nchar(labels)) * 0.1), 4, 2) + 0.1)
  plot(as.numeric(0),
       xlim = c(min(gtf$V4), max(gtf$V5)), xlab = 'position on reference',
       ylim = c(0.9, length(labels) + 0.1 ), ylab = '', yaxt="n", main = 'features')
  axis(2, at=1:length(labels), labels=labels, las=2)
  apply(gtf, 1, draw_feature)
}

plot_features(cf_gtf)
plot_features(ex_gff)

sum_of_all_exons <- function(gtf){
  gtf <- gtf[gtf$V3 == 'exon',]
  sum(gtf$V5 - gtf$V4)
}

sum_of_all_exons(cf_gtf)
sum_of_all_exons(ex_gff)
