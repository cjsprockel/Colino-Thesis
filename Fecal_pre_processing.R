
## This script is used only to prepare the combined and transformed datasets to be used during analyses

#######################
#### Fecal Data L5 ####
#######################
library(tidyverse)

workingDir = "C:/Users/cooli/Documents/School/Scriptie";
setwd(workingDir);

data <- read.table("Feces/roos_fecal_feature-table-dada2-maxee4_trunc280-L5.tsv", sep = "\t", skip=-2)


sample <- read.csv("Feces/mapping_file_feces.csv", sep = ";")

# remove samples of which we have no microbiome data.
sample <- sample[1:174,]




# check for later
mean(rowMeans(data[,-1]))

## Select only the bacteria with a mean count of over 10
goodbacteria <- rowMeans(data[,-1])>10
table(goodbacteria)
filtered <- data[goodbacteria,]
bad <- data[!goodbacteria,]

# Check to see if selection went through
mean(rowMeans(filtered[,-1]))
mean(rowMeans(bad[,-1]))

# create vector of bacteria names
names <- filtered[,1]

# split bacteria names by ; which separates names by taxonomic rank
nam <- strsplit(names, ";")
for (i in 1:length(names)){
  if (nchar(nam[[i]][5]) < 4) {
    x <- paste(nam[[i]][4], " ; " ,nam[[i]][5], sep="")
  nam[[i]] <- x
  }
  else {
    nam[[i]] <- nam[[i]][5]
  }
}


# set rownames to bacteria L4;L5 names
rownames(filtered) <- nam
# remove column 1 (contains rownames)
fecal_dat <- filtered[,-1]
fecal_dat <- t(fecal_dat)

######################################################
# We now have clean bacteria count data in the dataframe "fecal_dat"
# And we have clean patient information in the dataframe "sample"
######################################################



#####################################################
# Preparation for network analysis.
# Create correlation matrices ready to use in network construction
#######################################################


 # Compositional + log transform (Pearson correlation)


# function to make data compositional
comp <- function(x) {
  x <- apply(x,1,function(x) (x/sum(x)))
  return(t(x))
}

# clr1: add + 1 to count data for possible missing reads and for maths
# clr1: make data compositional
# clr1: take logarithm
# clr1: minus rowMeans to show relative abundance.
clr1 <- function(x){
  x <- x + 1
  x <- comp(x)
  x <- log(x)
  x <- x - rowMeans(x)
  return(x)
}
fecal_clr1 <- clr1(fecal_dat)
fecal_clr1_pearson <- cor(fecal_clr1, method="pearson")



  # Build correlation matrix using SparCC


library(SpiecEasi)

fecal_sparcc <- sparcc(fecal_dat, iter = 20, inner_iter = 10, th = 0.1)
# fecal_sparcc contains 2 matrices, the first is the covariance matrix, the second is the correlation matrix.
fecal_sparcc_cor <- fecal_sparcc[[2]]

###############################################
# Correlation matrices have been constructed
# We will now save the data and construct the networks in the Fecal_networks R file
###############################################


workingDir = "C:/Users/cooli/Documents/School/Scriptie";
setwd(workingDir);


save(fecal_sparcc_cor, fecal_clr1_pearson, sample, fecal_dat, fecal_clr1, file="Fecal_pre_processing_completed.RData")




