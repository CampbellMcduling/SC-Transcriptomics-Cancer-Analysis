#=========================== TRANSCRIPTOMICS CLUSTERING ============
rm(list=ls())
require(tidyverse)
require(Seurat)
library(matrixStats)
library(caret)
library(M3C)
library(factoextra)
#require(cowplot)
#require(scatr)
#require(scran)
#require(igraph)

dat = sc_pro_red

### investigate true labels
table(sc_info$index) #first level
table(sc_info[,3:4]) #second level
table(sc_info[,4:5]) #third level
table(sc_info[,c(3,5)])

#====================== DIMENSION REDUCTION ============
### PCA
pca.1 = M3C::pca(dat)
pca.1
pca.1$data
plot(pca.1$data[,1:3])
plot(pca.1$data[,1:3], col=as.factor(sc_info$index)) #coloured by tumor/non-tumor
plot(pca.1$data[,1:3], col=as.factor(sc_info$index2)) #coloured by type2
plot(pca.1$data[,1:3], col=as.factor(sc_info$index3)) #coloured by type3

#check- should be same as prcomp
tst = prcomp(x=t(dat), center = F, scale. = F)
plot(tst$sdev, type="l", xlab = "Principle Component", ylab="Standard deviation") + abline(v=30, lty=2, col="red")  # "elbow" plot
plot(tst$x[,1:3]) 
propvar = (cumsum(tst$sdev^2))/sum(tst$sdev^2) #prop var explained
plot(propvar, type="l", xlab = "Principle Component", ylab="Proportion Variance Explained") + abline(v=30, lty=2, col="red") 
propvar[30]


#====================== CLUSTERING ============
#=======determine k=====
#ELBOW METHOD
# within SS for k=2:20
k.max <- 8 # Max no. clusters
set.seed(200)
WSS <- sapply(1:k.max,
              function(k){kmeans(pca.1$data[,1:30], k, nstart=10)$tot.withinss})
plot(1:k.max, WSS, type="b", pch = 19, frame = FALSE, xlab="Number of clusters K", ylab="Total within-clusters sum of squares")



#GAP
library(cluster)
set.seed(200)
gap_stat <- clusGap(pca.1$data[,1:30], kmeans, nstart=10, K.max=8, B=5)
plot(gap_stat, frame = FALSE, xlab = "Number of clusters k", main = "")


###--------------- K-MEANS on PCA

### K=2
#run k-means on first X pca scores
set.seed(200)
km.pca1 = kmeans(pca.1$data[,1:30], centers=2, nstart = 10)

#compare k-means clustering to actual labels
table(km.pca1$cluster)
table(sc_info$index)
pred= as.factor(km.pca1$cluster)
table(pred)
true= sc_info$index
true = ifelse(true=="Tumor", 1, 2)
true=as.factor(true)
table(km.pca1$cluster)
table(true)

cfm1 = confusionMatrix(pred, true)
cfm1

plot(pca.1$data[,1:3], col=as.factor(-1*km.pca1$cluster))
plot(pca.1$data[,1:3], col=as.factor(sc_info$index)) 


#find mean expressions per cluster
#gene expressions
index.tum = which(km.pca1$cluster==1)
index.nontum = which(km.pca1$cluster!=1)
mean.tum = rowMeans(dat[,index.tum])
mean.nontum = rowMeans(dat[,index.nontum])
plot(mean.nontum, type='l') + lines(1:35000, mean.tum, col="red")

#PCA scores
clust_means <- km.pca1$centers
plot(clust_means[1,], type="l", xlab="Principle Component", ylab="Cluster mean", ylim=c(-23, 11), col="red") + lines(1:30, clust_means[2,], col="black")

###-------------- HEIRARCHICHAL 
###----determine K
temp = dat[,sample(1:ncol(dat), 0.5*ncol(dat), replace = F)] #random subset for computational efficiency and user-sanity maintenance
fviz_nbclust(t(temp), FUNcluster = hcut, method = "wss")
fviz_nbclust(t(temp), FUN = hcut, method = "silhouette")


### run clustering algorithm
# create dissimilarity matrix - euclidean
d.eucl = dist(t(dat), method="euclidean")

#perform clustering with ward linkage
set.seed(200)
hc.eucl.ward = hclust(d.eucl, method="ward.D2")

plot(hc.eucl.ward) 
#Cutting the cluster tree to make
# 2 groups
hc.cut2 <- cutree(hc.eucl.ward,k =2)
# 5 groups
hc.cut5 = cutree(hc.eucl.ward, k=5)

#compare HC clustering to actual labels
hc.cut2 = as.factor(hc.cut2)
table(hc.cut2)
table(true)
str(true); str(hc.cut2)

cfm3 = confusionMatrix(hc.cut2, true)
cfm3 #worse than kmeans
plot(pca.1$data[,1:3], col=hc.cut2)
plot(pca.1$data[,1:3], col=true)

#compare expression profiles of clusters
# Load necessary packages
library(ggplot2)
library(dplyr)
library(pheatmap)
library(gplots)

# Create a dataframe with labels
df <- data.frame(t(dat), hc.cut2) #with expression matrix
df.pc <- data.frame(pca.1$data, hc.cut2) #with PCs

# Rename column with labels to "Type"
colnames(df)[ncol(df)] <- "Type"
colnames(df.pc)[ncol(df.pc)] <- "Type"

#generate PC means for each cluster
hc.tum = df.pc[which(df.pc$Type==1),]
hc.nontum = df.pc[which(df.pc$Type==2),]
pcmean.tum = colMeans(hc.tum[,-ncol(hc.tum)])
pcmean.nontum = colMeans(hc.nontum[,-ncol(hc.nontum)])

# 1. Mean profiles
plot(pcmean.nontum[1:30], type="l", xlab="Principle Component", ylab="Cluster mean", ylim=c(-17, 12)) + lines(1:30, pcmean.tum[1:30], col="red")


# 4. Density plots
ggplot(df.pc, aes(x=df.pc[,1], fill=Type)) +
  geom_density(alpha=.5) +
  labs(title="Density plot", x="Expression", y="Density", fill="Type")
ggplot(df.pc, aes(x=df.pc[,2], fill=Type)) +
  geom_density(alpha=.5) +
  labs(title="Density plot", x="Expression", y="Density", fill="Type")
ggplot(df.pc, aes(x=df.pc[,3], fill=Type)) +
  geom_density(alpha=.5) +
  labs(title="Density plot", x="Expression", y="Density", fill="Type")
ggplot(df.pc, aes(x=df.pc[,4], fill=Type)) +
  geom_density(alpha=.5) +
  labs(title="Density plot", x="Expression", y="Density", fill="Type")
ggplot(df.pc, aes(x=df.pc[,5], fill=Type)) +
  geom_density(alpha=.5) +
  labs(title="Density plot", x="Expression", y="Density", fill="Type")
ggplot(df.pc, aes(x=df.pc[,6], fill=Type)) +
  geom_density(alpha=.5) +
  labs(title="Density plot", x="Expression", y="Density", fill="Type")
ggplot(df.pc, aes(x=df.pc[,7], fill=Type)) +
  geom_density(alpha=.5) +
  labs(title="Density plot", x="Expression", y="Density", fill="Type")
ggplot(df.pc, aes(x=df.pc[,8], fill=Type)) +
  geom_density(alpha=.5) +
  labs(title="Density plot", x="Expression", y="Density", fill="Type")


