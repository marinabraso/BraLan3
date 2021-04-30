setwd('/Users/aechchik/amphio_test_data')

remove <- c('0')

genic=read.table('./kallisto_genic_quant/abundance.tsv', header=TRUE)
tpm_genic=genic$tpm
tpm_genic_nonzero <- tpm_genic [! tpm_genic %in% remove]
log_tpm_genic_nonzero <- log2(tpm_genic_nonzero)

intergenic=read.table('./kallisto_intergenic_quant/abundance.tsv', header=TRUE)
tpm_intergenic=intergenic$tpm
tpm_intergenic_nonzero <- tpm_intergenic [! tpm_intergenic %in% remove]
log_tpm_intergenic_nonzero <- log2(tpm_intergenic_nonzero)

plot(density(log_tpm_intergenic_nonzero), xlim=c(-10, 10), ylim=c(0, 0.15), col='red', main='amphio - distribution of log2(tpm)', sub='blue=genic, red=intergenic', xlab='')
lines(density(log_tpm_genic_nonzero), xlim=c(-10, 10), ylim=c(0, 0.15), col='blue')
