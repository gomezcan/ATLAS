library(foreach)
library(doParallel)

data <- read.table(file = '../Projects/counts_single.txt', header = T, sep='\t')
data <- subset(data, select = -c(Strand, Length))
genes <- data[,1][!duplicated(data[,1])]
l_bim <- 10/1000

TPM <- function(g){ ## g is the gen name
  gen <- data[data$GeneID==g,] ## select all the gen set from data counts
  counts <- colSums(gen[,-1]) ## count sums vector for each sample
  Length <- as.numeric(nrow(gen))*l_bim ## Length of gen
  A <- counts/Length 
  
  tpm <- mapply(`/`, gen[,-1], A) ## Divide each sample by its corresponding A value
  tpm <- tpm*(10^6/l_bim)
  return(tpm)
}

#system.time(norm_tpm_counts <- mapply(TPM, genes[1:3]))

## Parallelization over all the genes
numCores <- detectCores()
registerDoParallel(numCores)

l <- length(genes)
system.time(norm_tpm_counts <- foreach (i = 1:l, .combine = rbind) %dopar% {
  TPM(genes[i])
})
