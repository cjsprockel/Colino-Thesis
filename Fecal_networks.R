####################################################
# Network construction from Fecal L5 count data
# 
####################################################

rm(list = ls())

workingDir = "C:/Users/cooli/Documents/School/Scriptie";
setwd(workingDir);

library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


# Loading correlation matrices
lnames = load(file = "Fecal_pre_processing_completed.RData")
lnames



nSamples = 174
nGenes = 24
nEdges = nGenes * (nGenes-1) / 2

gamma = log(nEdges, base=(sqrt(nSamples)))
print(gamma)

# A gamma of 2.18 indicates that any soft power scaling of over 2.18 is okay.
# Power should not be a lot higher than gamma.
# 


##  Looking at Scale independence and mean connectivity of clr1 correlation matrix
## Results are quite poor, A power of 5 seems most appropriate bt still a very poor fit.
## picking a power of 3 based on the gamma calculated above might be better.

sim_matrix <- abs(fecal_clr1_pearson)


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(sim_matrix, powerVector = powers, verbose = 5)


# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")






##########################################################
##########################################################
##########################################################
## TOM construction
## 
##########################################################
##########################################################
##########################################################

similarity_clr1 <- abs(fecal_clr1_pearson)

## Picking soft power of 3 based on calculation above
softPower = 3
adjacency_clr1 <- adjacency.fromSimilarity(similarity_clr1, power = softPower)
TOM_clr1 = TOMsimilarity(adjacency_clr1, TOMType = "unsigned")
dissTOM_clr1 = 1 - TOM_clr1

# Call the hierarchical clustering function
geneTree_clr1 = hclust(as.dist(dissTOM_clr1), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree_clr1, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


# Setting module size is important and has great impact on the modules formed.
# We pick a module size of 3, keeping in mind the relatively low amount of bacteria (24)
minModuleSize = 3;
# Module identification using dynamic tree cut:
dynamicMods_clr1 = cutreeDynamic(dendro = geneTree_clr1, distM = dissTOM_clr1,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods_clr1)





# Convert numeric lables into colors
dynamicColors_clr1 = labels2colors(dynamicMods_clr1)
table(dynamicColors_clr1)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_clr1, dynamicColors_clr1, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")








datExpr = fecal_clr1

# Calculate eigengenes


MEList_clr1 = moduleEigengenes(datExpr, colors = dynamicColors_clr1)
MEs_clr1 = MEList_clr1$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss_clr1 = 1-cor(MEs_clr1);
# Cluster module eigengenes
METree_clr1 = hclust(as.dist(MEDiss_clr1), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree_clr1, main = "Clustering of module eigengenes, based of clr1 tranformed data",
     xlab = "", sub = "")

####
# Plot shows that no merging is necessarry
###
moduleColors_clr1 <- dynamicColors_clr1



##########################################################
##########################################################
##########################################################
## Basic Network construction
## 
##########################################################
##########################################################
##########################################################

library(igraph)

Probes=colnames(fecal_dat)

dimnames(TOM_clr1) = list(Probes, Probes)


net_clr1 <- graph_from_adjacency_matrix(TOM_clr1, mode = "undirected", weighted = TRUE, diag=F)

# Calculate network density:
sum(E(net_clr1)$weight/(gsize(net_clr1)))


net_clr1 <- set_vertex_attr(net_clr1, "color", value = moduleColors_clr1)
plot(net_clr1)


# calculate node strength: the sum of all edge weights of a node.
# stronger nodes are better connected and persumably of more biological importance.
V(net_clr1)$strength <- strength(net_clr1, mode="all")

plot(net_clr1, 
     vertex.size= V(net_clr1)$strength * 90)

# Delete all edges with a value under 0.01 in the TOM, only the strongest edges remain
trimmed_clr1 <- delete.edges(net_clr1, which(E(net_clr1)$weight <0.01))


## Contruct the finished network 

plot(trimmed_clr1, 
     vertex.size= V(net_clr1)$strength * 90,
     vertex.label.color="black",
     vertex.color = adjustcolor(V(net_clr1)$color, alpha.f=0.6), # make nodes slightly opaque to improve readability
     frame.color=NA, #trying to remove the black frame surrounding the vertices (unsuccesfully)
     main="Network construction based on clr + 1 transformed count data",
     sub="edges>0.01 in the TOM are shown. Node size indicates strength of node")




















































## Looking at scale independence and Mean connectivity of Sparcc correlation matrix
## Results aren't great but better than those of the CLR1 correlation matrix.
## A power of 4 seems most appropriate

sim_matrix <- abs(fecal_sparcc_cor)


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(sim_matrix, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



