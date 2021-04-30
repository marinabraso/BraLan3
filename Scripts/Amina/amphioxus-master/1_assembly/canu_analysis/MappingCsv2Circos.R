library(data.table)

# So far it works with mapping to one scaffold only, but I guess it wont be difficult to extend it to more

MappingCsv2Circos <- function(file_name, header_file, out_name_pattern){
  print('loading data')
  map_csv <- read.csv(file_name, sep = ' ', header = F)
  colnames(map_csv) <- colnames(read.csv(header_file, sep = ' '))

  # extracting a flowcell and zmw for every subread
  map_csv$zmw <- factor(unlist(lapply(strsplit(as.character(map_csv$qName), split = '/'), "[[", 2)))
  map_csv$cells <- factor(unlist(lapply(strsplit(as.character(map_csv$qName), split = '/'), "[[", 1)))
  
  # I would like to capture data of sequencing depth
  depth <- rep(0, map_csv$tLength[1])
  
  # and links of reads connecting distinct regions
  scaffold_name <- map_csv$tName[1]
  links <- data.frame(ch1 = rep(scaffold_name,500000), s1 = 0, e1 = 0,
                      ch2 = scaffold_name, s2 = 0, e2 = 0)
  row_to_write = 1
  
  # I cut useless stuff form the table
  map_csv <- subset(map_csv, (((map_csv$qEnd - map_csv$qStart) / map_csv$qLength) > 0.15) & map_csv$percentSimilarity > 70)
  map_csv <- map_csv[,c(7:8,9:12,14,15)]
  
  # for every SMART cell
  bad = 0
  oneplace = 0
  several = 0
  
  list_of_cells <- levels(map_csv$cells)
  print('computing depth and links')
  
  for(cell in 1:length(list_of_cells)){
    SMRTcell <- subset(map_csv, map_csv$cells == list_of_cells[cell])
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
  print('Computation finished, numbers of bad, oneone and multiple mapping')
  print(bad)
  print(oneplace)
  print(several)

  dm <- matrix(depth, ncol = 100)
  circus_depth <- data.frame(chr = scaffold_name, 
                             from = seq(from = 1, by = 100, length.out = nrow(dm)), 
                             to = seq(from = 100, by = 100, length.out = nrow(dm)), 
                             depth = rowSums(dm) /100)
  
  circus_depth <- format(circus_depth, scientific = FALSE)
  
  links <- links[links$e1 != 0,]
  links <- format(links, scientific = FALSE)
  
  print(c('saving data to...',paste0('/Volumes/dump/pictures/source/circos/')))
  write.table(x = circus_depth, file = paste0('/Volumes/dump/pictures/source/circos/', out_name_pattern ,'_depth'), quote = F, row.names = F, col.names = F, sep = ' ')
  write.table(x = links, file = paste0('/Volumes/dump/pictures/source/circos/', out_name_pattern ,'_links'), quote = F, row.names = F, col.names = F, sep = ' ')  

  print(c('saving coverage plot to...','/Volumes/dump/pictures/pool/ngs/ampiho/hox/'))
  pdf(paste0('/Volumes/dump/pictures/pool/ngs/ampiho/hox',out_name_pattern,'_coverage.pdf'))
    plot(depth, type = 'l', xlab = 'Position', ylab = 'Coverage')
  dev.off()
}
