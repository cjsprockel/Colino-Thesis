library(philr)
library(ape)
library(phyloseq)
library(ggplot2)

source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
biocLite("philr")

setwd("T:/microbiologie/Research/Microbiota/2018.02 Project Roos_CVID feces samples/Results biostastics_HW/DADA2-HW/philr")

## first, fix the data
data <- read.delim("feces_L8_counts.txt", header = T, row.names = 1)
data.1<- data[1:174]

# they recommend lumping some OTUs together that are too infrequent
count <- 100 #this is the chosen cutoff
data.2 <- data.frame(data.1[which(apply(data.1, 1, sum) > count),],check.names=F)
data.3 <- data.frame(data.1[which(apply(data.1, 1, sum) < count),],check.names=F)
data.4 <- apply(data.3, 2, sum)
data.5 <- rbind(data.2, data.4)
# then they say to add a pseudocount of 1 to everything
# or do some alternative but let's just stick to this for now
data.6 <- apply(data.5,1, function(x) x+1)

## next, fix the phylogenic information
tree<- read.tree("tree.nwk")
# make sure tree is rooted and binary
tree.1 <- multi2di(tree)

# label the nodes?
tree.2 <- makeNodeLabel(phy_tree(tree.1), method="number", prefix='n')

## let's do the philr
data.philr <- philr(data.6, tree.1)
## Error in sbp[colnames(df), ] : subscript out of bounds
## apparently this means that the column names of the data file don't match the phylo name? 
## hang on is this because I lumped some OTUs together?

## let's try:
data.6 <- apply(data.1,1, function(x) x+1)
data.philr <- philr(data.6, tree.1)
# YES 
# IT WORKS 

## can we do pca with this?
data.philr.df <- as.data.frame(data.philr)





