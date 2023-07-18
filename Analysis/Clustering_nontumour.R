#=========================== NON-TUMOUR TRANSCRIPTOMICS CLUSTERING ============
rm(list=ls())
require(tidyverse)
require(Seurat)
library(matrixStats)
library(caret)
library(M3C)
library(elasticnet)
#require(cowplot)
#require(scatr)
#require(scran)
#require(igraph)
### ================ DATA MANAGEMENT ======
setwd("~/OneDrive - University of Cape Town/2023/MSc Biostatistics/Coursework Semester 1/Multivariate/Assignment/Datat/Breast Cancer ")
### DATA PREPROCESSING ###

# read in data on expression
bc_dat <- read.table(gzfile("GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt.gz"), header = TRUE)
# read in data on sample info
bc_info <- read.table(gzfile("GSE75688_final_sample_information.txt.gz"),header = TRUE)
str(bc_info)

# remove pooled samples and only keep single cell data from sample info 
sc_info <- bc_info[which(bc_info$type == "SC"),]
# remove pooled samples and only keep single cell data from expression data
sc_dat <- bc_dat[-(57913:57915),-(4:17)]

# REMOVE TUMOUR CELLS
sc_nontum_info <- sc_info[which(sc_info$index != "Tumor"),]

# order sample info to match order of expression data (alphabetical)
sc_nontum_info <- sc_nontum_info[order(sc_nontum_info$sample),]
rownames(sc_nontum_info) <- 1:nrow(sc_nontum_info)

# the sample info only contains high quality samples where the expression matrix
# contains all samples so match the expression data samples with sample info
# by removing poor quality samples.
sc_nontum_match <- sc_dat[,colnames(sc_dat) %in% sc_nontum_info$sample]
sc_nontum_match <- cbind(sc_dat[,1:3], sc_nontum_match)
rownames(sc_nontum_match) <- sc_nontum_match[,1]
sc_nontum_match <- sc_nontum_match[,-(1:3)]

## Pre-process data using Seurat
library(Seurat)
# create seurat object
seurat_dat <- CreateSeuratObject(counts = sc_nontum_match)
# reduce dataset to only 10000 most differentially expressed genes
var_gene <- FindVariableFeatures(seurat_dat, selection.method = 'vst', nfeatures = 20000)
# scale data
scale_sc_red <- ScaleData(object = var_gene)
# extract scaled and reduced expression matric from Seurat object
sc_pro_red <- GetAssayData(object = scale_sc_red, slot = "scale.data")
sc_pro_red <- as.data.frame(sc_pro_red)

dat = sc_pro_red



### investigate true labels
table(sc_nontum_info$index) #first level
table(sc_nontum_info[,3:4]) #second level
table(sc_nontum_info[,4:5]) #third level
table(sc_nontum_info[,c(3,5)])

#====================== DIMENSION REDUCTION ============
### PCA
pca.1 = M3C::pca(dat)
pca.1
pca.1$data
plot(pca.1$data[,1:3])
plot(pca.1$data[,1:3], col=as.factor(sc_nontum_info$index2)) #coloured by tumor/non-tumor
plot(pca.1$data[,1:3], col=as.factor(sc_nontum_info$index3)) #coloured by type2
plot(pca.1$data[,1:3], col=as.factor(sc_nontum_info$sample)) #coloured by sample

#check- should be same as prcomp
tst = prcomp(x=t(dat), center = F, scale. = F)
plot(tst$sdev, type="l") # "elbow" plot
plot(tst$x[,1:3]) 
propvar = (cumsum(tst$sdev^2))/sum(tst$sdev^2) #prop var explained
plot(propvar, type="l", xlab = "Principle Component", ylab="Proportion Variance Explained") + abline(v=200, lty=2, col="red") 
propvar[50]

### tSNE
library(Rtsne)
tsne.1 = Rtsne(t(dat), perplexity = 10)
tsne.2 = Rtsne(t(dat), perplexity = 30)
tsne.3 = Rtsne(t(dat), perplexity = 40)

#visualize tsne scores
#coloured by tumour/non-tumour
plot(tsne.1$Y[,1], tsne.1$Y[,2], col=as.factor(sc_nontum_info$index2), main="tSNE (perp=10)")
plot(tsne.2$Y[,1], tsne.2$Y[,2], col=as.factor(sc_nontum_info$index2), main="tSNE (perp=30)")
plot(tsne.3$Y[,1], tsne.3$Y[,2], col=as.factor(sc_nontum_info$index2), main="tSNE (perp=40)")
#coloured by 5 types
plot(tsne.1$Y[,1], tsne.1$Y[,2], col=as.factor(sc_nontum_info$index3), main="tSNE (perp=10)")
plot(tsne.2$Y[,1], tsne.2$Y[,2], col=as.factor(sc_nontum_info$index3), main="tSNE (perp=30)")
plot(tsne.3$Y[,1], tsne.3$Y[,2], col=as.factor(sc_nontum_info$index3), main="tSNE (perp=40)")


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
set.seed(2023)
gap_stat <- clusGap(pca.1$data[,1:30], kmeans, nstart=10, K.max=8, B=5)
plot(gap_stat, frame = FALSE, xlab = "Number of clusters k", main = "")


###--------------- K-MEANS on PCA

### K=2
#run k-means on first X pca scores
set.seed(200)
km.pca1 = kmeans(pca.1$data[,1:30], centers=2, nstart = 10)

#compare k-means clustering to actual labels
table(km.pca1$cluster)
table(sc_nontum_info$index2)
pred= as.factor(km.pca1$cluster)
table(pred)
true= sc_nontum_info$index2
true = ifelse(true=="Immune", 2, 1)
true=as.factor(true)
table(km.pca1$cluster)
table(true)

cfm1 = confusionMatrix(pred, true)
cfm1

plot(pca.1$data[,1:3], col=as.factor(-1*km.pca1$cluster))
plot(pca.1$data[,1:3], col=as.factor(sc_nontum_info$index2)) 


### K=4
#run k-means on first X pca scores
set.seed(200)
km.pca1.1 = kmeans(pca.1$data[,1:50], centers=4, nstart = 10)

table(km.pca1.1$cluster)
table(sc_nontum_info$index3)
plot(pca.1$data[,1:3], col=as.factor(-1*km.pca1.1$cluster))
plot(pca.1$data[,1:3], col=as.factor(sc_nontum_info$index3)) 


###--------------- K-MEANS on tSNE
#run k-means on first X pca scores
set.seed(200)
km.tsne3 = kmeans(tsne.3$Y[,1:2], centers=2)

#compare k-means clustering to actual labels
table(km.tsne3$cluster)
table(sc_nontum_info$index)
pred= as.factor(km.tsne3$cluster)
table(pred)
true= sc_nontum_info$index2
true = ifelse(true=="Tumor", 1, 2)
true=as.factor(true)
table(km.tsne3$cluster)
table(sc_nontum_info$index)

cfm2 = confusionMatrix(pred, true)
cfm2

plot(tsne.3$Y, col=as.factor(-1*km.tsne3$cluster))
plot(tsne.3$Y, col=as.factor(sc_nontum_info$index)) 

###-------------- HEIRARCHICHAL 
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

# 5 groups - allocation of cell types must be done visually: results look poor
plot(pca.1$data[,1:3], col=hc.cut5)
plot(pca.1$data[,1:3], col=as.factor(sc_nontum_info$index3))


