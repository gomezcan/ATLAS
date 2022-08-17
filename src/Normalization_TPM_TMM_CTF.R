######## Libraries required
library(foreach)
library(doParallel)## Two parallelization packages
library(edgeR)

############################
args <- commandArgs( TRUE) ## Input Counts and stats
methods <- c('TMM','CTF','TPM')  ## Select the method
method <- methods[2] ## SELECT ONE OF THE THREE METHODS
########### Import Data##############

counts <- read.table(file = args[1], header = T, sep = '\t')
counts_stat <- read.table(file = args[2], head=TRUE)
genes <- counts[,1]
Samples <- colnames(counts_stat)[-1]

######### Filtering Samples by count percents assigned ###################

percents <- counts_stat[15,]
caja <- boxplot(x=as.numeric(percents[1,-1])) # set a box plot to get the outliersQ1 <- caja$stat[1] ## First quantile 
std <- sd(x=percents[-1]) ## standard deviation
mu <- mean(as.numeric(percents[-1])) # Mean

# It's gonna remove the outliers
removed_samples <- percents[, percents >= Q1] # It chose Q1 but could be use mu-2*std

#### Rename removed_samples columns: Removing the ".bam" pattern
colnames(removed_samples) <- gsub(pattern='*.bam', replacement = "", colnames(removed_samples))


### Getting the filtered data
new_counts <- subset(counts, select = c(colnames(removed_samples)[-1])) ### Counts to be normalized
new_counts <- counts[[colnames(removed_samples)[-1]]]

########################## MA-plots before normalization #####

maPlot(new_counts[,1], new_counts[,200], deCol = 'green')
abline(h=0, col = 'red')


########################## Normalization ################################

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
  factors_norm <- calcNormFactors(Counts, lib_size = library_size, method = 'TMM')
  ### Normalize the counts
  counts_norm <- sweep(Counts, 2, factors_norm, FUN = '/')
  return(counts_norm)
}

######## CTF ##########

### Calculate the scale factors for each sample with calcNormFactor 

TMM <- function(Counts){
  ### Calculate the normalization factor for every sample
  factors_norm <- calcNormFactors(Counts, lib_size = library_size, method = 'TMM')
  ### Effective size
  effective_size <- factors_norm*as.numeric(library_size)/10**6
  ### Normalize the counts
  counts_norm <- sweep(Counts, 2, effective_size, FUN = '/')
  return(counts_norm)
}

##################### DEFINE THE NORMALIZATIONFUNCTION ##############

Normalization <- function(method, X){
  if(method=='TMM'){
    X_normed <- TMM(X)
  }
  
  else if (method=='CTF'){
    X_normed <- CTF(X)
  }
  
  else{
    ## Parallelization option, much faster than mapply
    numCores <- detectCores()
    registerDoParallel(numCores)
    
    l <- length(ncol(new_counts))
    
      X_normed <- foreach (i = 1:l, .combine = cbind) %dopar% {
      TPM(colnames(new_counts)[i])
    }
  }
  return(X_normed)
}
###################### EXECUTING ##########################

Counts_normalized <- Normalization(new_counts)


maPlot(Counts_normalized[,1], Counts_normalized[,200],deCol = 'green') ## Plot mat plot with tmm normalization
abline(h=0, col = 'red')
