a20x <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_cov_tests/amphio20_lengths.csv', header = F)$V1
a40x <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_cov_tests/amphio40_lengths.csv', header = F)$V1
a60x <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_cov_tests/amphio60_lengths.csv', header = F)$V1
a80x <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_cov_tests/amphio80_lengths.csv', header = F)$V1

source('/Volumes/dump/scripts/R/genomics/getNX.R')
source('/Volumes/dump/scripts/R/genomics/getN50.R')

full <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_e025/amphio_contig_lengths.csv', header = F)$V1
full_ref <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_e025/amphio_contigs_lengths_ref_haplotype.csv', header = F)$V1
full_alt <- read.csv('/Volumes/dump/data/sequences/draft_genomes/Amphioxus/canu_e025/amphio_contigs_lengths_alt_haplotype.csv', header = F)$V1

library(RColorBrewer)
SUBpal <- brewer.pal(4, 'BrBG')
FULpal <- brewer.pal(4, 'YlOrRd')[2:4]


sum(a20x)
length(a20x)
max(a20x)
# computed N50
getN50(a20x)

png('/Volumes/dump/pictures/pool/ngs/assembly/amphio_cov_NX_graph.png')
X <- seq(0,100, by = 0.02)
NX <- getNX(a20x, X, 520000000)
plot(NX ~ X, type = 'l', main = 'amphio', ylim = c(1, max(full_alt)), col = SUBpal[1])
NX <- getNX(a40x, X, 520000000)
lines(NX ~ X, col = SUBpal[2])
NX <- getNX(a60x, X, 520000000)
lines(NX ~ X, col = SUBpal[3])
NX <- getNX(a80x, X, 520000000)
lines(NX ~ X, col = SUBpal[4])

NX <- getNX(full, X, 520000000)
lines(NX ~ X, col = FULpal[1])
NX <- getNX(full_ref, X, 520000000)
lines(NX ~ X, col = FULpal[2])
NX <- getNX(full_alt, X, 520000000)
lines(NX ~ X, col = FULpal[3])

legend('topright', pch = 20, col = c(SUBpal, FULpal), legend = c('20x','40x','60x', '80x', '153x', '153x ref', '153x alt'))

dev.off()
