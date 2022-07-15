######## Libraries required
library(foreach)
library(doParallel)## Two parallelization packages
library(edgeR)
args <- commandArgs(trailingOnly = TRUE)

########### Import Data##############

counts <- read.table(file='..', header = T, sep = '\t')
counts_stat <- read.table(file='..', head=TRUE)
genes <- counts[,1]
Samples <- colnames(counts_stat)[-1]

######### Filtering Samples by count percents assigned ###################

percents <- counts_stat[15,]
caja <- boxplot(x=as.numeric(percents[1,-1])) # set a box plot to get the outliersQ1 <- caja$stat[1] ## First quantile 
std <- sd(x=percents[-1]) ## standard deviation
mu <- mean(as.numeric(percents[-1])) # Mean

# It's gonna remove the outliers
removed_samples <- percents[, percents >= Q1] # It chose Q1 but but could be use mu-2*std

#### Rename removed_samples columns: Removing the ".bam" pattern
colnames(removed_samples) <- gsub(pattern='*.bam', replacement = "", colnames(removed_samples))


### Getting the filtered data
new_counts <- subset(counts, select = c(colnames(removed_samples)[-1])) ### Counts to be normalized
new_counts <- counts[[colnames(removed_samples)[-1]]]
########################## MA-plots before normalization #####

maPlot(new_counts[,1], new_counts[,200], deCol = 'green')
abline(h=0, col = 'red')


########################## Normaization ################################

library_size <- counts_stat[1,][-1] ## Total Reads Assigned by sample vector
colnames(library_size) <- gsub(pattern='*.bam', replacement = "", colnames(library_size)) ## Remove ".bam"
library_size <- library_size[, colnames(removed_samples)[-1]] ## Library size for samples after QC on the counts

######## TPM #########
L <- sum(counts[["Length"]])/1000 ## a number
length <- counts[["Length"]]/1000 ## a vector

TPM <- function(sample){
  
  Sample <- new_counts[[sample]]
  Counts <- sum(Sample)
  A <- Counts/L
  #tpm <- mapply((Sample*10^6)/length)
  tpm <- mapply(`/`, Sample*10/6, length)
  tpm <- tpm/A
  
  return(tpm)
}

######## CTF ##########

### Calculate the scale factors for each sample with calcNormFactor 

CTF <- function(Counts){
  ### Calculate the normalization factor for every sample
  factors_norm <- calcNormFactors(Counts, library_size = library_size, method = 'TMM')
  ### Normalize the counts
  counts_norm <- sweep(Counts, 2, factors_norm, FUN = '/')
  return(counts_norm)
}

######## CTF ##########

### Calculate the scale factors for each sample with calcNormFactor 

TMM <- function(Counts){
  ### Calculate the normalization factor for every sample
  factors_norm <- calcNormFactors(Counts, library_size = library_size, method = 'TMM')
  ### Effective size
  effective_size <- factors_norm*as.numeric(library_size)
  ### Normalize the counts
  counts_norm <- sweep(Counts, 2, effective_size, FUN = '/')
  return(counts_norm)
}

## Execute TPM
system.time( tpm <- mapply(TPM, colnames(new_counts)))

## Paralleization option, much faster than mapply
numCores <- detectCores()
registerDoParallel(numCores)

l <- length(ncol(new_counts))
system.time(norm_tpm_counts <- foreach (i = 1:l, .combine = cbind) %dopar% {
  TPM(colnames(new_counts)[i])
})
## Execute CTF

ctf <- CTF(new_counts)

## Execute TMM
tmm <- TMM(new_counts)


maPlot(tmm[,1], tmm[,200],deCol = 'green') ## Plot mat plot with ctf normalization
abline(h=0, col = 'red')
