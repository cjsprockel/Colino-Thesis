## This script is used to analyse the differencesbetween reguar log ratio transform followed by pearson 
# #compared to sparCC when cnstructing networks using WGCNA

#######################
#### Oral Data L7 ####
#######################

library(igraph)

library(WGCNA)
options(stringsAsFactors = FALSE);


myDensity <- function(adjmatrix) {
  
  net <- graph_from_adjacency_matrix(adjmatrix, mode = "undirected", weighted = TRUE, diag=F)
  
  nGenes = ncol(adjmatrix)
  
  totalWeigths = sum(E(net)$weight)
  totalEdges = (nGenes * (nGenes - 1))
  
  density = (totalWeigths/totalEdges)
  
  return(density)
}


myCentralization <- function(adjmatrix){
  
  
  net <- graph_from_adjacency_matrix(adjmatrix, mode = "undirected", weighted = TRUE, diag=F)
  
  strengths <- strength(net, mode="all")
  density = myDensity(adjmatrix)
  
  centralization = (nGenes/nGenes - 2) * ((max(strengths)/nGenes-1) - density)
  
  return(centralization)
}

workingDir = "C:/Users/cooli/Documents/School/Scriptie";
setwd(workingDir);


dat <- readRDS('L6_8000.RDS')
sample <- readRDS('sample_8000.RDS')
###?dat_merge <- readRDS('merged_L6_8000.RDS')
# dat <- readRDS('L6_8000.RDS')
# sample <- readRDS('sample.RDS')
# all.equal(rownames(dat), rownames(dat_HW))
all.equal(rownames(dat), rownames(sample))
# all.equal(dat, dat_HW[,1:170])
###?rm(dat_merge)

# replace troublesome symbols in colnames (these symbols will mess up in RF)
colnames(dat) <- gsub(colnames(dat), pattern = "\\[|\\]", replacement = "")
colnames(dat) <- gsub(colnames(dat), pattern = " ", replacement = "_")









## Removing Bacteria with a mean count under 10. very low abundant bacteria are an unreliable source of information.


mean(colMeans(dat))


goodbacteria <- colMeans(dat)>10
table(goodbacteria)
bacteria_names = colnames(dat[goodbacteria])
dat <- dat[goodbacteria]
colnames(dat) = bacteria_names
mean(colMeans(dat))


comp <- function(x) {
  x <- apply(x,1,function(x) (x/sum(x)))
  return(t(x))
}

# clr1
clr1 <- function(x){
  x <- x + 1
  x <- comp(x)
  x <- log(x)
  x <- x - rowMeans(x)
  return(x)
}
dat_clr1 <- clr1(dat)



dat_clr1_spear <- cor(dat_clr1, method="spearman")

## Pearson will be used since it is advised by the WGCNA tutorial
dat_clr1_pearson <- cor(dat_clr1, method="pearson")


# To determine the power used to turn our similarity matrix into an adjacency matrix we will use the formula used 
# by Bartzis et al (2017)

nSamples = 152
nGenes = 164
nEdges = nGenes * (nGenes-1) / 2

gamma = log(nEdges, base=(sqrt(nSamples)))
print(gamma)

# power must be at least as large as gamma but not much larger


# Plotting scale free topology shows very high variance, This shows that this similarity matrix is unreliable.



similarity <- abs(dat_clr1_pearson)


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(similarity, powerVector = powers, verbose = 5)


# Plot the results:
sizeGrWindow(5,3)
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



# A power of 4 will be chosen


softPower = 4;


# Construct adjacency matrix
adjacency <- adjacency.fromSimilarity(similarity, power = 4)

# Construct Topological Overlap Matrix
TOM = TOMsimilarity(adjacency, TOMType = "unsigned")

TOM_clr1 <- TOM

dissTOM = 1 - TOM


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);




# We like large modules, so we set the minimum module size relatively high:
# Set in such a way that it results in 4 modules.
minModuleSize = 3;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
table(dynamicColors)



# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")








datExpr = dat_clr1

# Calculate eigengenes

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# No merging necessary

colnames(TOM) = bacteria_names
net <- graph_from_adjacency_matrix(TOM, mode = "undirected", weighted = TRUE, diag=F)

clr1_colors <- dynamicColors

myDensity(TOM)
myCentralization(TOM)

net <- set_vertex_attr(net, "color", value = dynamicColors)
plot(net, vertex.label = "")
trimmed <- delete.edges(net, which(E(net)$weight <0.02))
plot(trimmed, vertex.label = "")


weightdeg <- strength(net, mode="all")
V(net)$size <- weightdeg * 20

trimmed <- delete.edges(net, which(E(net)$weight <0.03))
plot(trimmed, vertex.label = "")

plot(trimmed)

############################################################################################
############################################################################################
############################################################################################
############################################################################################
# Now network construction using SparCC correlation matrix
# mostly same code as above, now using SparCC
############################################################################################
############################################################################################
############################################################################################
############################################################################################

library(SpiecEasi)


sparccmatrix <- sparcc(dat, iter = 20, inner_iter = 10, th = 0.3)

cor_matrix <- data.frame(sparccmatrix[2])
cor_matrix <- data.matrix(cor_matrix, rownames.force = NA)

colnames(cor_matrix) = bacteria_names
rownames(cor_matrix) = bacteria_names




similarity_sparcc <- (abs(cor_matrix))

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(similarity_sparcc, powerVector = powers, verbose = 5)


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



# A power of 4 will be chosen


softPower = 4;


# Construct adjacency matrix
adjacency <- adjacency.fromSimilarity(similarity_sparcc, power = 4)

# Construct Topological Overlap Matrix
TOM = TOMsimilarity(adjacency, TOMType = "unsigned")

dissTOM = 1 - TOM


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);




# We like large modules, so we set the minimum module size relatively high:
# Set in such a way that it results in 4 modules.
minModuleSize = 3;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
table(dynamicColors)



# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")








datExpr = dat_clr1

# Calculate eigengenes

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# No merging necessary

colnames(TOM) = bacteria_names
net <- graph_from_adjacency_matrix(TOM, mode = "undirected", weighted = TRUE, diag=F)

sparcc_colors <- dynamicColors

myDensity(TOM)
myCentralization(TOM)

net <- set_vertex_attr(net, "color", value = dynamicColors)
plot(net, vertex.label = "")
trimmed <- delete.edges(net, which(E(net)$weight <0.02))
plot(trimmed, vertex.label = "")


weightdeg <- strength(net, mode="all")
V(net)$size <- weightdeg * 20

trimmed <- delete.edges(net, which(E(net)$weight <0.03))
plot(trimmed, vertex.label = "")

plot(trimmed)




plot(trimmed, vertex.frame.color= clr1_colors, vertex.label.color = clr1_colors)



cor(c(similarity), c(similarity_sparcc))
cor(c(TOM), c(TOM_clr1))










#############################################################################
## Consensus analysis
#############################################################################