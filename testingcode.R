
## test for github


library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)



BiocManager::install("RCy3")

workingDir = "C:/Users/cooli/Documents/School/Scriptie";
setwd(workingDir);

dat <- readRDS('L6_8000.RDS')

sparccmatrix <- sparcc(dat, iter = 20, inner_iter = 10, th = 0.3)

matrix <- sparccmatrix[[2]]
sparccmax = max(sparccmatrix[[1]])

n1 = sparccmatrix[[1]]
diag(n1) <- 1
n1 <- n1/max(n1)
diag(n1) <- 1






df = dat[,!sapply(dat, function(x) mean(x))<100]





goodgenes <- colMeans(dat)>10
table(goodgenes)
filtered <- dat[goodgenes]

####################
# 43 bacteria with mean counts > 10 
# The rest has been discarded
####################






plot(net)



plot(net,
  
     vertex.size = V(trimmed)$strength * 5,
     edge.width= (edge_attr(trimmed)$weight * 12))






l <- layout_with_fr(net)
plot(net, layout=l)


plot(net, main="Blue Module.   node size = sum(edge weights)",
     vertex.size = V(net)$strength * 5,
     edge.width= (edge_attr(net)$weight * 12))


plot(net, main="Blue module.    node size = sum(edge weights)",
     vertex.size = V(net)$strength * 5,
     edge.width= (edge_attr(trimmed)$weight * 12))



V(net)$eigenvector <- evcent(net)$vector





plot(net,
     main="Brown module network.    Node size = eigenvvector centrality",
     vertex.size = V(net)$eigenvector/max(V(net)$eigenvector) * 15)


plot(net,
     main="Brown module network.    Node size = eigenvvector centrality",
     sub="smaller edges removed, edge size related to weight",
     edge.width= (edge_attr(trimmed)$weight * 20)^1.3,
     vertex.size = V(net)$eigenvector/max(V(net)$eigenvector) * 15,
     layout=l)




plot(trimmed, vertex.color="orange",
     main="Brown Module network TOM \n node size = sum(edge weights) \n edges > 0.05 shown",
     vertex.size = V(trimmed)$strength * 15,
     edge.width=2)





plot(trimmed,
     main="Brown module network TOM.  \n  Node size = eigenvector centrality \n edges >0.05 shown",
     vertex.size = V(net)$eigenvector/max(V(net)$eigenvector) * 30,
     edge.width=2)





plot(trimmed, main="Blue module TOM \n  Node size = sum(edge weights) \n edges >0.1 shown",
     vertex.size = V(trimmed)$strength * 5,
     edge.width=2)





plot(trimmed,
     main="Blue module network.  \n  Node size = eigenvvector centrality \n edges > 0.1 shown",
     vertex.size = V(net)$eigenvector/max(V(net)$eigenvector) * 20,
     edge.width=2)


plot(net, vertex.label = "")


V(net)$strength <- strength(net, mode="all")
plot(net,
     vertex.label = "", 
     vertex.size = V(net)$strength * 1)


plot(trimmed, main="Network of TOM \n  Node size = sum(edge weights) \n edges >0.05 shown",
     vertex.label="",
     vertex.size = V(trimmed)$strength * 2,
     edge.width=1)

plot(trimmed,
     main="Full TOM network  \n  Node size = eigenvector centrality \n edges >0.05 shown",
     vertex.label="",
     vertex.size = V(net)$eigenvector/max(V(net)$eigenvector) * 18,
     edge.width=1)




sim_matrix <- abs(n1)


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
