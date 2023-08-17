# Martin King. Created 10 Aug 2023.
# Linear Discriminant Analysis for classifying stratified methane emission.
# 1. LDA.
# 2. PLSR-DA.

# 1. LDA.------------------------------------------------------------------------
rm(list=ls())

set.seed(2023)
library(MASS)

data.df = read.csv("/Users/martinpeterking/ucc_courseworks.dir/semester3.dir/data_work.dir/PredMethEMk_DMKComp_cowsinsingleyears_mpkedited.csv", header=T)

# Calculate quantile thresholds
quantiles <- quantile(data.df$avg_ch4_smilk, probs = c(0, 0.5, 1))
#quantiles <- quantile(data.df$avg_ch4_smilk, probs = c(0, 0.33, 0.66, 1))
#quantiles <- quantile(data.df$avg_ch4_smilk, probs = c(0, 0.20, 0.4, 0.6, 0.8, 1))
# Add small tolerance so that all values are within range.
quantiles <- quantiles + c(-0.1, 0, 0.1)
#quantiles <- quantiles + c(-0.1, 0, 0, 0.1)
#quantiles <- quantiles + c(-0.1, 0, 0, 0, 0, 0.1)
# Create a new column with categorical levels
data.df$methcat <- cut(data.df$avg_ch4_smilk, breaks = quantiles, labels = c("lo", "hi"))
#data.df$methcat <- cut(data.df$avg_ch4_smilk, breaks = quantiles, labels = c("lo", "me", "hi"))
#data.df$methcat <- cut(data.df$avg_ch4_smilk, breaks = quantiles, labels = c("lo", "ml", "me", "mh", "hi"))

data_y = as.factor(data.df$methcat)

grouping = data.df$TB_NUM
# Number of folds.
K = 4
# Bootstrap.
B = 100

accuracy.train = numeric(K*B)
accuracy.test = numeric(K*B)
misclassification.train = matrix(nrow=K*B, ncol=2)
misclassification.test = matrix(nrow=K*B, ncol=2)
#misclassification.train = matrix(nrow=K*B, ncol=3)
#misclassification.test = matrix(nrow=K*B, ncol=3)
#misclassification.train = matrix(nrow=K*B, ncol=5)
#misclassification.test = matrix(nrow=K*B, ncol=5)


for (b in 1:B)
{
  folds = groupKFold(grouping, k = K)
  for (k in 1:K)
  {
   itrain = as.numeric(unlist(folds[k]))
   
   data.spec.df = data.frame(data.df[itrain, 19:549])
   mu = apply(data.spec.df, 2, mean)
   Xminusmu = data.spec.df - matrix(mu, nrow = nrow(data.spec.df), ncol = ncol(data.spec.df), byrow = TRUE)
   pc.mir = prcomp(Xminusmu, center=FALSE, scale=FALSE)
   data_x_train = data.frame(pc.mir$x[,1:20], data.df$dim[itrain], data.df$yield[itrain])
   
   data.spec.df = data.frame(data.df[-itrain, 19:549])
   Xminusmu = data.spec.df - matrix(mu, nrow = nrow(data.spec.df), ncol = ncol(data.spec.df), byrow = TRUE)
   pc.test = as.matrix(Xminusmu) %*%  pc.mir$rotation
   data_x_test = data.frame(pc.test[,1:20], data.df$dim[-itrain], data.df$yield[-itrain])
  
   data_x_test <- data_x_test %>%
     rename(data.df.dim.itrain. = data.df.dim..itrain.)
   data_x_test <- data_x_test %>%
     rename(data.df.yield.itrain. = data.df.yield..itrain.)
   
   ldaout = lda(x=data_x_train, grouping=data_y[itrain])
  
   ldaout.train = predict(ldaout, newdata=data_x_train)
   confmat.train = table(ldaout.train$class, data_y[itrain])
   ldaout.test = predict(ldaout, newdata=data_x_test)
   confmat.test = table(ldaout.test$class, data_y[-itrain])
   
   indx = (b-1)*K+k
  
   # Training accuracy. 
   accuracy.train[indx] = sum(diag(confmat.train))/sum(confmat.train)
   # Confusion matrix in sum of all rows in each col to 1.
   confmat.train = round(sweep(confmat.train, 2, apply(confmat.train,2,sum), FUN="/" ), 3)
   # Training misclassification
   misclassification.train[indx,] = 1-diag(confmat.train)
   
   # Test accuracy. 
   accuracy.test[indx] = sum(diag(confmat.test))/sum(confmat.test)
   # Confusion matrix in sum of all rows in each col to 1.
   confmat.test = round(sweep(confmat.test, 2, apply(confmat.test,2,sum), FUN="/" ), 3)
   # Training misclassification
   misclassification.test[indx,] = 1-diag(confmat.test)
  }
}  

# Just some plotting.

plot(ldaout$svd^2, main="F-statistic")
par(mfrow=c(2,2))
par(cex.main = 1.5)
par(cex.axis = 1.4) 
par(cex.lab = 1.4)
par(mar = c(5, 5, 4, 2)) #Bottom, then clockwise.

Z = as.matrix(data_x_train) %*% ldaout$scaling
boxplot(Z[,1]~data_y[itrain], xlab="", ylab="", main="Discrim Var 1")
boxplot(Z[,2]~data_y[itrain], xlab="", ylab="", main="Discrim Var 2")
boxplot(Z[,3]~data_y[itrain], xlab="", ylab="", main="Discrim Var 3")
boxplot(Z[,4]~data_y[itrain], xlab="", ylab="", main="Discrim Var 4")

# 2. PLSR-DA.------------------------------------------------------------------------

rm(list=ls())

set.seed(2023)
library(caret)

data.df = read.csv("/Users/martinpeterking/ucc_courseworks.dir/semester3.dir/data_work.dir/PredMethEMk_DMKComp_cowsinsingleyears_mpkedited.csv", header=T)

# Calculate quantile thresholds
#quantiles <- quantile(data.df$avg_ch4_smilk, probs = c(0, 0.33, 0.66, 1))
quantiles <- quantile(data.df$avg_ch4_smilk, probs = c(0, 0.20, 0.4, 0.6, 0.8, 1))
# Add small tolerance so that all values are within range.
#quantiles <- quantiles + c(-0.1, 0, 0, 0.1)
quantiles <- quantiles + c(-0.1, 0, 0, 0, 0, 0.1)
# Create a new column with categorical levels
#data.df$methcat <- cut(data.df$avg_ch4_smilk, breaks = quantiles, labels = c("lo", "ne", "hi"))
data.df$methcat <- cut(data.df$avg_ch4_smilk, breaks = quantiles, labels = c("lo", "ml", "ne", "mh", "hi"))

data_y = as.factor(data.df$methcat)

grouping = data.df$TB_NUM
# Number of folds.
K = 4
# Bootstrap.
B = 100

accuracy.train = numeric(K*B)
accuracy.test = numeric(K*B)
#misclassification.train = matrix(nrow=K*B, ncol=3)
#misclassification.test = matrix(nrow=K*B, ncol=3)
misclassification.train = matrix(nrow=K*B, ncol=5)
misclassification.test = matrix(nrow=K*B, ncol=5)


for (b in 1:B)
{
  folds = groupKFold(grouping, k = K)
  for (k in 1:K)
  {
    itrain = as.numeric(unlist(folds[k]))
    
    data_x_train = data.frame(data.df[itrain, 19:549], data.df$yield[itrain], data.df$dim[itrain])
    data_x_test = data.frame(data.df[-itrain, 19:549], data.df$yield[-itrain], data.df$dim[-itrain])
  
    ata_x_test <- data_x_test %>%
      rename(data.df.dim.itrain. = data.df.dim..itrain.)
    data_x_test <- data_x_test %>%
      rename(data.df.yield.itrain. = data.df.yield..itrain.)
    
    plsda.out = plsda(x = data_x_train, y = data_y[itrain], ncomp=40)

    ldaout.train = predict(plsda.out, newdata=data_x_train)
    confmat.train = table(ldaout.train, data_y[itrain])
    ldaout.test = predict(plsda.out, newdata=data_x_test)
    confmat.test = table(ldaout.test, data_y[-itrain])
    
    indx = (b-1)*K+k
    
    # Training accuracy. 
    accuracy.train[indx] = sum(diag(confmat.train))/sum(confmat.train)
    # Confusion matrix in sum of all rows in each col to 1.
    confmat.train = round(sweep(confmat.train, 2, apply(confmat.train,2,sum), FUN="/" ), 3)
    # Training misclassification
    misclassification.train[indx,] = 1-diag(confmat.train)
    
    # Test accuracy. 
    accuracy.test[indx] = sum(diag(confmat.test))/sum(confmat.test)
    # Confusion matrix in sum of all rows in each col to 1.
    confmat.test = round(sweep(confmat.test, 2, apply(confmat.test,2,sum), FUN="/" ), 3)
    # Training misclassification
    misclassification.test[indx,] = 1-diag(confmat.test)
  }
}  
