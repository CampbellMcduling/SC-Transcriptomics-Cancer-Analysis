
### PCA ###

# perform PCA on data
PCA <- prcomp(t(sc_pro_red))

# plot first 3 PCs and colour by tumour and non-tumour
pairs(PCA$x[,1:3], col = as.factor(sc_info$index))
# plot first 3 PCs and colour by 5 cell types
pairs(PCA$x[,1:3], col = as.factor(sc_info$index3))
# plot first 3 PCs and colour by sample
pairs(PCA$x[,1:3], col = samp)

# plot proportion of variance explained by PC
plot(PCA$sdev/sum(PCA$sdev), type = "n",
     ylab = "Proportion variance explained",xlab = "Principal component")
lines(PCA$sdev/sum(PCA$sdev))
abline(v=30, col = "red")

sum(PCA$sdev[1:30]/sum(PCA$sdev)) # total prop var explained by first 30 PCs


### t-SNE ###

library(Rtsne)
# perform tSNE with perplexity of 30
tsne <- Rtsne(t(sc_pro_red), dims = 2, perplexity=30,
              verbose=TRUE, max_iter = 500)

plot(tsne$Y[,1], tsne$Y[,2], col = as.factor(sc_info$index),
     xlab = "tSNE 1", ylab = "tSNE 2")
legend("bottomleft",legend = c("non-tumor","tumor"), col = pal, pch = 16)
plot(tsne$Y[,1], tsne$Y[,2], col = as.factor(sc_info$index3),
     xlab = "tSNE 1", ylab = "tSNE 2")
legend("bottomleft",legend = sort(unique(sc_info$index3)), col = pal, pch = 16)
plot(tsne$Y[,1], tsne$Y[,2], col = samp, xlab = "tSNE 1", ylab = "tSNE 2")
legend("bottomleft",legend = unique(samp), col = pal, pch = 16)

### Kernel PCA ###

library(kernlab)

# Get datset with 10000 genes instead of 35000 because kernlab can't handle 35000
var_gene2 <- FindVariableFeatures(seurat_dat, selection.method = 'vst',
                                  nfeatures = 10000)
# scale data
scale_sc_red2 <- ScaleData(object = var_gene2)
# extract scaled and reduced expression matric from Seurat object
sc_pro_red2 <- GetAssayData(object = scale_sc_red2, slot = "scale.data")
sc_pro_red2 <- as.data.frame(sc_pro_red2)
t_sc_pro <- as.data.frame(t(sc_pro_red2))

# do kPCA using polydot
kPCA1 <- kpca(~., data = t_sc_pro, kernel = "polydot", 
              kpar = list(degree = 1),features = 10)

plot(rotated(kPCA1), col = as.factor(sc_info$index))
plot(rotated(kPCA1), col = as.factor(sc_info$index3))
plot(rotated(kPCA1), col = samp)

# do kPCA using laplacedot
kPCA2 <- kpca(~., data = t_sc_pro, kernel = "laplacedot", kpar = list(sigma = 0.001),
              features = 10)
plot(rotated(kPCA2), col = as.factor(sc_info$index))
plot(rotated(kPCA2), col = as.factor(sc_info$index3))
plot(rotated(kPCA2), col = samp)

# do kPCA using anovadot
kPCA3 <- kpca(~., data = t_sc_pro, kernel = "anovadot", kpar = list(sigma = 0.001, degree = 2),
              features = 10)
plot(rotated(kPCA3), col = as.factor(sc_info$index))
plot(rotated(kPCA3), col = as.factor(sc_info$index3))
plot(rotated(kPCA3), col = samp)

# do kPCA using rbfdot
kPCA4 <- kpca(~., data = t_sc_pro, kernel = "rbfdot", 
              kpar = list(sigma = 0.0001), features = 10)
plot(rotated(kPCA4), col = as.factor(sc_info$index),
     xlab = "PC1", ylab = "PC2")
plot(rotated(kPCA4), col = as.factor(sc_info$index3))
plot(rotated(kPCA4), col = samp)
