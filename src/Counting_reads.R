library(Rsubread)

list_bam <- list.files(path='/maindisk/camilo/Projects', pattern = "*.bam")
#bams <- readline(prompt = 'inpt list')
#list_bam <- readLines(bams)
typeof(list_bam)
group <- list_bam
########### Counting for single end ################

counts_table <- featureCounts(group, 
                  annot.ext = "New_Annotation.gtf",
                  isGTFAnnotationFile = T,
                  GTF.featureType = 'CDS',
                  useMetaFeatures = F,
                  isPairedEnd = F,
                  nthreads = 25,
                  countMultiMappingReads = F,
                  allowMultiOverlap = T)
  
porcent_assigned <- round((counts_table$stat[1,-1]*100)/colSums(counts_table$stat[list_bam]), digits = 1)
porcent_assigned <- cbind(Status ='porcent_assigned', porcent_assigned)
counts_table$stat <- rbind(counts_table$stat, porcent_assigned)
counts_table$stat


############## Data Frame with the counts ##################

df_single <- data.frame(counts_table$annotation[,c("GeneID","Strand","Length")], counts_table$counts, stringsAsFactors=FALSE)

colnames(df_single) <- gsub(".bam", "", colnames(df_single))

write.table(df_single, paste(bams,"_single.txt", sep=""), row.names =F, quote =F, sep ='\t')
write.table(counts_table$stat, 'counts_single_stat.txt', row.names =F, quote =F, sep ='\t')

annotation <- read.table(file = 'New_Annotation.gtf')
