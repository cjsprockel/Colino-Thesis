
## This script is used only to prepare the combined and transformed datasets to be used during analyses

#######################
#### Fecal Data L5 ####
#######################






workingDir = "C:/Users/cooli/Documents/School/Scriptie";
setwd(workingDir);







data <- read.table("Feces/roos_fecal_feature-table-dada2-maxee4_trunc280-L5.tsv", sep = "\t", skip=-2)

dataL7 <- read.table("Feces/roos_fecal_feature-table-dada2-maxee4_trunc280-L7.tsv", sep = "\t", skip=-2)

sample <- read.csv("Feces/mapping_file_feces.csv", sep = ";")







