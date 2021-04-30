rm(list = ls())

library(data.table)

pb2pb <- read.csv('/Volumes/dump/data/mapping/amphio/hox/pb2pb/PBreads2PBalthox.csv', sep = ' ', header = F)
colnames(pb2pb) <- colnames(read.csv('/Volumes/dump/data/mapping/amphio/hox/variables', sep = ' '))
head(pb2pb)

# extracting a flowcell and zmw for every subread
pb2pb$zmw <- factor(unlist(lapply(strsplit(as.character(pb2pb$qName), split = '/'), "[[", 2)))
pb2pb$cells <- factor(unlist(lapply(strsplit(as.character(pb2pb$qName), split = '/'), "[[", 1)))

# I would like to capture data of sequencing depth
depth <- rep(0, pb2pb$tLength[1])

# and links of reads connecting distinct regions
scaffold_name <- pb2pb$tName[1]
links <- data.frame(ch1 = rep(scaffold_name,500000), s1 = 0, e1 = 0,
                    ch2 = scaffold_name, s2 = 0, e2 = 0)
row_to_write = 1

# I cut useless stuff form the table
pb2pb <- subset(pb2pb, (((pb2pb$qEnd - pb2pb$qStart) / pb2pb$qLength) > 0.15) & pb2pb$percentSimilarity > 70)
pb2pb <- pb2pb[,c(7:8,9:12,14,15)]

# for every SMART cell
bad = 0
oneplace = 0
several = 0

list_of_cells <- levels(pb2pb$cells)



for(cell in 1:length(list_of_cells)){
  SMRTcell <- subset(pb2pb, pb2pb$cells == list_of_cells[cell])
  list_of_zmw <- levels(factor(SMRTcell$zmw))
  # and all zmw wells  - i.e. every iteration is one set of subreads of one modelucle
  for(well in 1:length(list_of_zmw)){
    subreads <- subset(SMRTcell, SMRTcell$zmw == list_of_zmw[well])
    if(((max(subreads$qEnd) - min(subreads$qStart)) / max(subreads$qLength)) < 0.99){
      bad = bad + 1;
      next;
    }
    if(any(((subreads$qEnd - subreads$qStart) / subreads$qLength) > 0.98)){
#     chose the best mapping and increase coverage where it mappes;
      the_best <- which.max((subreads$qEnd - subreads$qStart) / subreads$qLength)
      depth[subreads$tStart[the_best]:subreads$tEnd[the_best]] <- depth[subreads$tStart[the_best]:subreads$tEnd[the_best]] + 1
      
      oneplace = oneplace + 1
      next;
    }
    several = several + 1

    #ordering might not be necessery
    subreads <- subreads[order(subreads$qStart),]
    
    for(i in 1:(nrow(subreads)-1)) {
      set(links,as.integer(row_to_write),2L, subreads$tStart[i])
      set(links,as.integer(row_to_write),3L, subreads$tEnd[i])
      set(links,as.integer(row_to_write),5L, subreads$tStart[i+1])
      set(links,as.integer(row_to_write),6L, subreads$tEnd[i+1])
      row_to_write = row_to_write + 1;
    }
  }
}
bad
oneplace
several


# to get approx coverage
oneplace * mean(pb2pb$qLength) / 1500000

length(depth)

dm <- matrix(depth, ncol = 100)
circus_depth <- data.frame(chr = scaffold_name, 
                           from = seq(from = 1, by = 100, length.out = nrow(dm)), 
                           to = seq(from = 100, by = 100, length.out = nrow(dm)), 
                           depth = rowSums(dm) /100)

circus_depth <- format(circus_depth, scientific = FALSE)
links <- format(links, scientific = FALSE)

write.table(x = circus_depth, file = '/Volumes/dump/pictures/source/circos/pb2pb_depth', quote = F, row.names = F, col.names = F, sep = ' ')
write.table(x = links, file = '/Volumes/dump/pictures/source/circos/pb2pb_links', quote = F, row.names = F, col.names = F, sep = ' ')

links <- read.csv('/Volumes/dump/pictures/source/circos/pb2pb_links', header = F, sep = ' ')
head(links)

max(circus_depth$depth)
head(links)
links <- links[links$e1 != links$s1[500000], ]



plot(circus_depth$from, as.numeric(circus_depth$depth) * 1000, type = 'l')

pdf('/Volumes/dump/pictures/pool/ngs/ampiho/pb2pb_coverage.pdf')
  plot(depth, type = 'l', xlab = 'Position', ylab = 'Coverage')
dev.off()

file_name <- '/Volumes/dump/data/mapping/amphio/hox/pb2il/PBreads2ILhox.csv'
header_file <- '/Volumes/dump/data/mapping/amphio/hox/variables'
out_name_pattern <- 'pb2il'

MappingCsv2Circos(file_name, header_file, out_name_pattern)

# EXPLORATORY
# oneplace
# head(pb2pb)
# 
# sum(mapping_lengths > 0.98)
# 
# hist(mapped_per_subread)
# #hist(mapping_lengths, breaks = 320, xlim = c(0.95,1))
# hist(mapping_lengths, breaks = 20000, xlim = c(0,50))
# sum(mapping_lengths < 50)
# hist(pb2pb$qLength)
# hist((pb2pb$qEnd - pb2pb$qStart) / pb2pb$qLength, breaks = 320)
