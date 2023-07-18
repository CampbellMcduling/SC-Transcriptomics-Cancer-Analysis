
dat = read.csv("dat.csv", row.names = 1)
pca.1 = M3C::pca(dat)
dat.pc = pca.1$data[,1:100]
i2 = as.factor(sc_info$index2)
i3 = as.factor(sc_info$index3)
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
                         mtry = ncol(dat.train) - 1,
                         ntree = 1000,
                         importance = F,
                         na.action = na.exclude,
                         do.trace = 100)
# Bagging - mtry = 14
set.seed(2022)
bag_dat2 <- randomForest(i3 ~ ., data = dat.train,
                         mtry = ncol(dat.train) - 5,
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
legend('topright', legend = c('mtry=101', 'mtry=97', 'mtry=93'), 
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

#model performance - INSAMPLE
bag_dat$confusion
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
save(xgb_dat1, file = 'xgb_dat1.1.Rdata')
load('xgb_dat1.1.Rdata')

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
#OOS accuracy is better for RF
