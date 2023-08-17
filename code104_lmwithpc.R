# Martin King. Created 4 Aug 2023.
# Linear Regression Model using PCs.

rm(list=ls())

library(caret)
library(dplyr)

set.seed(2023)

data.df = read.csv("/Users/martinpeterking/ucc_courseworks.dir/semester3.dir/data_work.dir/PredMethEMk_DMKComp_cowsinsingleyears_mpkedited.csv", header=T)

# Number of folds.
K = 4
# Bootstrap.
B = 100

rmse.train.k4 = numeric(K*B)
rmse.test.k4 = numeric(K*B)
bias.train.k4 = numeric(K*B)
bias.test.k4 = numeric(K*B)
cor.train.k4 = numeric(K*B)
cor.test.k4 = numeric(K*B)
reg.train.k4 = numeric(K*B)
reg.test.k4 = numeric(K*B)
rpiq.train.k4 = numeric(K*B)
rpiq.test.k4 = numeric(K*B)

grouping = data.df$TB_NUM

for (b in 1:B)
{
  folds = groupKFold(grouping, k = K)
  
  for (k in 1:K)
  {
    itrain = as.numeric(unlist(folds[k]))
    
    data.spec.df = data.frame(data.df[itrain, 19:ncol(data.df)])
    mu = apply(data.spec.df, 2, mean)
    Xminusmu = data.spec.df - matrix(mu, nrow = nrow(data.spec.df), ncol = ncol(data.spec.df), byrow = TRUE)
    
    pc.mir = prcomp(Xminusmu, center=FALSE, scale=FALSE)
    data.train.df = data.frame(data.df$avg_ch4_smilk[itrain], pc.mir$x[,1:20])
    
    # Test data are centred to mu and projected on pcs of training data.
    data.spec.df = data.frame(data.df[-itrain, 19:ncol(data.df)])
    Xminusmu = data.spec.df - matrix(mu, nrow = nrow(data.spec.df), ncol = ncol(data.spec.df), byrow = TRUE)
    
    pc.test = as.matrix(Xminusmu) %*%  pc.mir$rotation
    data.test.df = data.frame(data.df$avg_ch4_smilk[-itrain], pc.test[,1:20])
    
    # Just renaming so predict can work.
#   data.test.df <- data.test.df %>%
#      rename(data.df.dim.itrain. = data.df.dim..itrain.)
#    data.test.df <- data.test.df %>%
#      rename(data.df.yield.itrain. = data.df.yield..itrain.)
#    data.test.df <- data.test.df %>%
#      rename(data.df.fat.itrain. = data.df.fat..itrain.)
#    data.test.df <- data.test.df %>%
#      rename(data.df.protein.itrain. = data.df.protein..itrain.)
#    data.test.df <- data.test.df %>%
#      rename(data.df.lactose.itrain. = data.df.lactose..itrain.)
    
    lm.model = lm(data.df.avg_ch4_smilk.itrain. ~ ., data=data.train.df)
    lm.train = predict(lm.model, newdata=data.train.df)
    lm.pred = predict(lm.model, newdata=data.test.df)
    indx = (b-1)*K+k
    
    # RMSE.
    rmse.train.k4[indx] =  sqrt(mean((data.train.df[,1] - lm.train)^2.0))
    rmse.test.k4[indx] =  sqrt(mean((data.test.df[,1] - lm.pred)^2.0))
    # Bias.
    bias.train.k4[indx] = mean(lm.train) - mean(data.train.df[,1])
    bias.test.k4[indx] = mean(lm.pred) - mean(data.test.df[,1])
    # Correlation.
    cor.train.k4[indx] = cor(lm.train, data.train.df[, 1])
    cor.test.k4[indx] = cor(lm.pred, data.test.df[, 1])
    # Regression coefficient.
    # reg.train.k4[k] = lm(datanew.df[itrain,1] ~ pls.train[,,20])
    # reg.test.k4[k] = lm(datanew.df[-itrain,1] ~ pls.test[,,20])
    # RPIQ
    rpiq.train.k4[indx] = IQR(data.train.df[, 1])/rmse.train.k4[indx]
    rpiq.test.k4[indx] = IQR(data.test.df[, 1])/rmse.test.k4[indx]
  }
}

#par(cex.axis = 1.5) 
#par(cex.main = 1.5)
#par(cex.lab = 1.5)
#plot(cumsum(pc.mir$sdev**2)/sum(pc.mir$sdev**2
#), log="x", xlab="PC No.", ylab="", 
#main="Cumulative fraction of variance explained")

