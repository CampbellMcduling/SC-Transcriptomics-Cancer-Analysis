#=========================== TRANSCRIPTOMICS - ENSEMBLE METHODS ============
rm(list=ls())
require(tidyverse)
require(Seurat)
library(matrixStats)
library(caret)
library(M3C)
library(factoextra)
library(dplyr)
library(caret)
library(xgboost)
library(randomForest)
### ================ DATA MANAGEMENT ======
setwd("~/OneDrive - University of Cape Town/2023/MSc Biostatistics/Coursework Semester 1/Multivariate/Assignment/Datat/Breast Cancer ")
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


###============= RANDOM FORESTS ========
pca.1 = M3C::pca(dat)
dat.pc = pca.1$data[,1:120]
i2 = as.factor(sc_nontum_info$index2)
i3 = as.factor(sc_nontum_info$index3)
dat.pc = cbind(dat.pc, i3)


#split into train and validation sets
set.seed(2022)
train_index <- createDataPartition(y = dat.pc$i3, p = 0.6, list = FALSE)
dat.train <- dat.pc[train_index,]
dat.valid <- dat.pc[-train_index,]

#============================= RANDOM FORESTS ==============================
#tuning params = B, mtry

# Bagging - mtry = 18
set.seed(2022)
 bag_dat1 <- randomForest(i3 ~ ., data = dat.train,
                             mtry = ncol(dat.train) -2,
                             ntree = 1000,
                             importance = F,
                             na.action = na.exclude,
                             do.trace = 100)
# Bagging - mtry = 14
set.seed(2022)
 bag_dat2 <- randomForest(i3 ~ ., data = dat.train,
                          mtry = ncol(dat.train) - 6,
                          ntree = 1000,
                          importance = F,
                          na.action = na.exclude,
                          do.trace = 100)

# Bagging - mtry = 10
set.seed(2022)
 bag_dat3 <- randomForest(i3 ~ ., data = dat.train,
                          mtry = ncol(dat.train) - 9,
                          ntree = 1000,
                          importance = F,
                          na.action = na.exclude,
                          do.trace = 100)


# CONFUSION
bag_dat1$confusion
bag_dat2$confusion
bag_dat3$confusion

which.min(bag_dat1$err.rate[,1]); min(bag_dat1$err.rate[,1]) #minimum err rate and corresponding ntree
which.min(bag_dat2$err.rate[,1]); min(bag_dat2$err.rate[,1]) #minimum err rate and corresponding ntree
which.min(bag_dat3$err.rate[,1]); min(bag_dat3$err.rate[,1]) #minimum err rate and corresponding ntree


## Compare OOB Errors
par(mfrow = c(1,1))
plot(bag_dat1$err.rate[,1], type = 'l', xlab = 'Number of trees', ylab = 'OOB Error Rate', 
     col = 'blue', lwd = 2, ylim = c(0, max(bag_dat1$err.rate[,1], bag_dat2$err.rate[,1])))
lines(bag_dat2$err.rate[,1], col = 'darkgreen', lwd = 2)
lines(bag_dat3$err.rate[,1], col = 'pink', lwd = 2)
abline(v=which.min(bag_dat3$err.rate[,1]), lty=2, col="pink")
legend('topright', legend = c('mtry=100', 'mtry=96', 'mtry=92'), 
       col = c('blue', 'darkgreen', 'yellow'), lwd = 2, lty = c('solid', 'solid', 'solid'))

#Run optimal forest
nt = which.min(bag_dat2$err.rate[,1])
mt = bag_dat2$mtry
set.seed(2022)
bag_dat <- randomForest(i3 ~ ., data = dat.train, 
                         mtry = mt,
                         ntree = nt, 
                         importance = TRUE, 
                         na.action = na.exclude, 
                         do.trace = 100)
plot(bag_dat)

#model performance - INSAMPLE
bag_dat$confusion
min(bag_dat$err.rate[,1])
plot(bag_dat$err.rate[,1], type="l")

#model performance - out of sample
bag_dat.pred = predict(bag_dat, dat.valid)
table(bag_dat.pred)
confusionMatrix(bag_dat.pred, dat.valid$i3)

#============================= GRADIENT BOOSTING ==============================

#Parallelisation initialisation 
library(doParallel)
cores <- detectCores() - 2 #Don't use all your cores
cl <- makePSOCKcluster(cores)
registerDoParallel(cl)

#grid search to determine hyperparams
xgb_grid1 <- expand.grid(nrounds = c( 30, 50, 70, 100),  #B 
                         max_depth = c(2:8),      #d 
                         eta = c(0.01, 0.005, 0.001),       #lambda 
                         gamma = 0.001,            #mindev
                         colsample_bytree = 1,     #proportion random features per tree
                         min_child_weight = 1,     #also controls tree depth
                         subsample = 1             #bootstrap proportion
)
ctrl <-  trainControl(method = 'cv', number = 5, verboseIter = T)

set.seed(2022)
 xgb_dat1 <- train(i3 ~ ., data = dat.train,
                    method = 'xgbTree',
                    trControl = ctrl,
                    verbose = F,
                    tuneGrid = xgb_grid1, allowParallel=TRUE)
 stopCluster(cl)
# save(xgb_dat1, file = 'xgb_dat1.1.Rdata')
# load('xgb_dat1.1.Rdata')

#compare xgb configurations 
plot(xgb_dat1)
xgb_dat1$bestTune

#model performance of best configuration
xgb_dat1$results["66",]

#compare RF and XGB
# INSAMPLE
bag_dat
xgb_dat1$results["66",]
#accuracy of XGB is slightly better

# OUTOFSAMPLE
xgb.pred = predict(xgb_dat1, dat.valid)
table(xgb.pred)
confusionMatrix(xgb.pred, dat.valid$i3)
confusionMatrix(bag_dat.pred, dat.valid$i3)
table(bag_dat.pred)
#OOS accuracy is better for RF

#====================== VARIABLE IMPORTANCE =========
impToPlot.rf <- importance(bag_dat, scale=FALSE)

varImpPlot(bag_dat, scale = TRUE, type = 1, main="")
#gradient boosted tree

plot(varImp(xgb_dat1, scale = TRUE), main="", xlab="Variable Importance (%)") 




