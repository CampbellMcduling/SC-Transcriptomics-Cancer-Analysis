
### k-means (no dim reduction) ###

kmax <- 10
wss <- sapply(1:kmax,function(k){kmeans(t(sc_pro_red),k,nstart = 10)$tot.withinss})
plot(1:kmax,wss,type = 'b',xlab = 'Number of clusters',ylab = 'Total within-cluster sum of squares')
abline(v=2, col = "red")

set.seed(2002)
km <- kmeans(t(sc_pro_red),2,nstart = 10)
misclass <- length(which((km$cluster-as.numeric(factor(sc_info$index, levels = c("Tumor","nonTumor") , labels = c(1,2))))!=0))/515

clust_means <- km$centers
par(mar=c(7, 4, 4, 2))
plot(c(1,ncol(t(sc_pro_red))),range(clust_means), type = 'n',
     xaxt='n',ylab = 'Average scaled expression',xlab = '')
axis (side=1, 1:ncol(t(sc_pro_red)), red_gene_names, las=2)
colvec <- c("black","red")
lines (1:ncol(t(sc_pro_red)),clust_means[2,],col=colvec[2])
lines (1:ncol(t(sc_pro_red)),clust_means[1,],col=colvec[1])
par(mar=c(5, 4, 4, 2))

### Graph-Based (using Seurat) ###

S_obj <- RunPCA(object = scale_sc_red)
S_obj <- FindNeighbors(S_obj, reduction = "pca", dims = 1:30)
S_obj <- FindClusters(S_obj, resolution = 0.8)
S_obj <- RunTSNE(object = S_obj)

DimPlot(object = S_obj, reduction = "tsne")
DimPlot(object = S_obj, reduction = "pca")


### Biclustering ###

library(biclust)

#set.seed(2000)
#bi_res <- biclust(as.matrix(sc_match),method = BCPlaid(),cluster = "b", fit.model = y~m+a+b)
set.seed(2000)
bi_res <- biclust(as.matrix(sc_match),method = BCPlaid(),cluster = "b", fit.model = y~m+a+b,
                  row.release = 0.5, col.release = 0.5, iter.startup = 15, 
                  iter.layer = 20)

summary(bi_res)

parallelCoordinates(as.matrix(sc_match),bi_res,number = 1)
parallelCoordinates(as.matrix(sc_match),bi_res,number = 2)
parallelCoordinates(as.matrix(sc_match),bi_res,number = 3)
parallelCoordinates(as.matrix(sc_match),bi_res,number = 4)
parallelCoordinates(as.matrix(sc_match),bi_res,number = 5)
parallelCoordinates(as.matrix(sc_match),bi_res,number = 6)

drawHeatmap(as.matrix(sc_match),bi_res,number = 1)
drawHeatmap(as.matrix(sc_match),bi_res,number = 2)
drawHeatmap(as.matrix(sc_match),bi_res,number = 3)
drawHeatmap(as.matrix(sc_match),bi_res,number = 4)
drawHeatmap(as.matrix(sc_match),bi_res,number = 5)
drawHeatmap(as.matrix(sc_match),bi_res,number = 6)


bi_res2 <- biclust(as.matrix(sc_match),method = BCCC(), delta = 0.1, alpha = 1.5, number = 3)
parallelCoordinates(as.matrix(sc_pro_red),bi_res2,number = 1)
parallelCoordinates(as.matrix(sc_pro_red),bi_res2,number = 2)
parallelCoordinates(as.matrix(sc_pro_red),bi_res2,number = 3)

drawHeatmap(as.matrix(sc_pro_red),bi_res2,number = 1)
drawHeatmap(as.matrix(sc_pro_red),bi_res2,number = 2)
drawHeatmap(as.matrix(sc_pro_red),bi_res2,number = 3)


