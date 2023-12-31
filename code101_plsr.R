# Martin King. Created 4 Aug 2023.
# PLSR. Repeating McParland et al. 2023.
# Redone with 218 cows that appeared in single years only.
# PLSR.
# 1. CV. Group measurements by same cows together.
# 1.5 Model diagnostics.
# 2. Leave one treatment out validation.

# 1. CV. Group measurements by same cows together.---------------

library(caret)
library(pls)

rm(list=ls())

# Read edited file. 
data.df = read.csv("/Users/martinpeterking/ucc_courseworks.dir/semester3.dir/data_work.dir/PredMethEMk_DMKComp_cowsinsingleyears_mpkedited.csv", header=T)

set.seed(2023)

# Defining data to use.
# Select AM and/or PM.
#selrow = which(data.df$milking_time==1 | data.df$milking_time==2)
#selrow = which(data.df$TRT!="Gra")
#date.df = data.frame(as.Date(data.df$milk_date))
#selrow = which(format(date.df, "%Y")==2022 & data.df$TRT!="Gra")
selrow = which(data.df$recs_smilk > 30)
datanew.df = data.frame(data.df$avg_ch4_smilk[selrow], data.df$dim[selrow], 
                        data.df$yield[selrow], data.df[selrow,19:ncol(data.df)])
#datanew.df = data.frame(data.df$avg_ch4_smilk[selrow], data.df$yield[selrow], data.df$dim[selrow])
#datanew.df = data.frame(data.df$avg_ch4_smilk[selrow], data.df[selrow,19:ncol(data.df)])

# Group measurements by the same cow together.
grouping = data.df$TB_NUM[selrow]

# Number of folds.
K = 4
# Bootstrap.
B = 100

# Defining variables to store metrics
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


for (b in 1:B)
{
  folds = groupKFold(grouping, k = K)
  for (k in 1:K)
  {
    itrain = as.numeric(unlist(folds[k]))
    # Calibration/training.
    pls.model <- plsr(data.df.avg_ch4_smilk.selrow. ~ ., data = datanew.df[itrain, ], 
                     ncomp = 20, scale = FALSE, validation="CV")
    #pls.model <- lm(data.df.avg_ch4_smilk.selrow. ~ ., data = datanew.df[itrain, ])
    # Prediction.
    pls.train <- predict(pls.model, newdata=datanew.df[itrain, ])
    pls.test <- predict(pls.model, newdata=datanew.df[-itrain, ])
    
    indx = (b-1)*K+k
    
    # RMSE.
    rmse.train.k4[indx] =  sqrt(mean((datanew.df[itrain, 1] - pls.train[,,20])^2.0))
    rmse.test.k4[indx] =  sqrt(mean((datanew.df[-itrain, 1] - pls.test[,,20])^2.0))
    # Bias.
    bias.train.k4[indx] = mean(pls.train[,,20]) - mean(datanew.df[itrain, 1])
    bias.test.k4[indx] = mean(pls.test[,,20]) - mean(datanew.df[-itrain, 1])
    # Correlation.
    cor.train.k4[indx] = cor(pls.train[,,20], datanew.df[itrain, 1])
    cor.test.k4[indx] = cor(pls.test[,,20], datanew.df[-itrain, 1])
    # Regression coefficient.
    # reg.train.k4[k] = lm(datanew.df[itrain,1] ~ pls.train[,,20])
    # reg.test.k4[k] = lm(datanew.df[-itrain,1] ~ pls.test[,,20])
    # RPIQ
    rpiq.train.k4[indx] = IQR(datanew.df[itrain, 1])/rmse.train.k4[indx]
    rpiq.test.k4[indx] = IQR(datanew.df[-itrain, 1])/rmse.test.k4[indx]
  }
}
mean(rmse.test.k4)
sd(rmse.test.k4)
mean(bias.test.k4)
sd(bias.test.k4)
mean(cor.test.k4)
sd(cor.test.k4)
mean(rpiq.test.k4)
sd(rpiq.test.k4)

# 1.5 Model diagnostics.--------------------------------------------------------

library(pls)

rm(list=ls())

# Read edited file. 
data.df = read.csv("/Users/martinpeterking/ucc_courseworks.dir/semester3.dir/data_work.dir/PredMethEMk_DMKComp_cowsinsingleyears_mpkedited.csv", header=T)

datanew.df = data.frame(data.df$avg_ch4_smilk, data.df$yield, data.df$dim, data.df[,19:ncol(data.df)])
pls.model <- plsr(data.df.avg_ch4_smilk ~ ., data = datanew.df, ncomp = 20, scale = FALSE)

datanew.df = data.frame(data.df$avg_ch4_smilk, data.df$yield, data.df$dim)
lm.model <- lm(data.df.avg_ch4_smilk ~ ., data = datanew.df)

par(mfrow=c(2,2))
par(cex.main = 1.5)
par(cex.axis = 1.4) 
par(cex.lab = 1.4)
par(mar = c(5, 5, 4, 2)) #Bottom, then clockwise.

library(scales)

# Residuals vs. Fitted values
plot(lm.model, which=1, col = alpha("black", 0.3), main="(a)")

# QQ plot
plot(lm.model, which=2, col = alpha("black", 0.3), main="(b)")

#Cook's Distance
plot(lm.model, which=4, col = alpha("black", 0.3), main="(c)")

#Cook's Distance vs. Leverage (h/1-h)
plot(lm.model, which=6, col = alpha("black", 0.3), main="(d)")

# Residuals vs. fitted values.
#plot(pls.model$fitted.values[,,20], pls.model$residuals[,,20], main="Residuals vs. fitted values", xlab="Fitted values", ylab="Residuals", col = alpha("black", 0.3))
#plot(lm.model$fitted.values, lm.model$residuals, main="Residuals vs. fitted values", xlab="Fitted values", ylab="Residuals", col = alpha("black", 0.3))
#abline(h = 0, col = "coral1", lwd = 2)

#qqnorm(pls.model$residuals[,,20], pch=16, col = alpha("black", 0.3))
#qqline(pls.model$residuals[,,20])
#qqnorm(lm.model$residuals, pch=16, col = alpha("black", 0.3))
#qqline(lm.model$residuals)


#plot(data.df$yield, pls.model$residuals[,,20], main="Residuals vs. yield", xlab="yield", ylab="Residuals", col = alpha("black", 0.3))
#plot(data.df$yield, lm.model$residuals, main="Residuals vs. yield", xlab="yield", ylab="Residuals", col = alpha("black", 0.3))
#abline(h = 0, col = "coral1", lwd = 2)

#plot(data.df$dim, pls.model$residuals[,,20], main="Residuals vs. dim", xlab="dim", ylab="Residuals", col = alpha("black", 0.3))
#plot(data.df$dim, lm.model$residuals, main="Residuals vs. dim", xlab="dim", ylab="Residuals", col = alpha("black", 0.3))
#abline(h = 0, col = "coral1", lwd = 2)

# Partial residuals plot
#library(car)
#datanew.df = data.frame(data.df$avg_ch4_smilk, data.df$yield, data.df$dim)
#pls.model <- lm(data.df.avg_ch4_smilk ~ ., data = datanew.df)
#crPlots(pls.model, terms = ~ data.df.yield)
#crPlots(pls.model, terms = ~ data.df.dim)


# 2. Leave one treatment out.---------------------------------------------------
rm(list=ls())

data.df = read.csv("/Users/martinpeterking/ucc_courseworks.dir/semester3.dir/data_work.dir/PredMethEMk_DMKComp_cowsinsingleyears_mpkedited.csv", header=T)

library(caret)
library(pls)

#set.seed(2023)

# Defining data to use.
#selrow = which(data.df$TRT=="Gra")
#data.df = data.df[-selrow, ]
datanew.df = data.frame(data.df$avg_ch4_smilk, 
                        data.df$dim, 
                        data.df$yield,
                        data.df[,19:ncol(data.df)])
#datanew.df = data.frame(data.df$avg_ch4_smilk,
#                        data.df[,19:ncol(data.df)])

treatment=c("C15", "G25", "Gra", "HHG", "HLG", "LHG", "LLG")
#treatment= c("C15", "G25", "HHG", "HLG", "LHG", "LLG")

# 7 treatments.
K=7
#K=6

rmse.train.k4 = numeric(K)
rmse.test.k4 = numeric(K)
bias.train.k4 = numeric(K)
bias.test.k4 = numeric(K)
cor.train.k4 = numeric(K)
cor.test.k4 = numeric(K)
reg.train.k4 = numeric(K)
reg.test.k4 = numeric(K)
rpiq.train.k4 = numeric(K)
rpiq.test.k4 = numeric(K)

for (k in 1:K)
{
  itest = which(data.df$TRT==treatment[k])
  
  # Calibration/training.
  pls.model <- plsr(data.df.avg_ch4_smilk ~ ., data = datanew.df[-itest, ], ncomp = 20, scale = TRUE)
  # Prediction.
  pls.train <- predict(pls.model, newdata=datanew.df[-itest, ])
  pls.test <- predict(pls.model, newdata=datanew.df[itest, ])
  
  # RMSE.
  rmse.train.k4[k] =  sqrt(mean((datanew.df[-itest, 1] - pls.train[,,20])^2.0))
  rmse.test.k4[k] =  sqrt(mean((datanew.df[itest, 1] - pls.test[,,20])^2.0))
  # Bias.
  bias.train.k4[k] = mean(pls.train[,,20]) - mean(datanew.df[-itest, 1])
  bias.test.k4[k] = mean(pls.test[,,20]) - mean(datanew.df[itest, 1])
  # Correlation.
  cor.train.k4[k] = cor(pls.train[,,20], datanew.df[-itest, 1])
  cor.test.k4[k] = cor(pls.test[,,20], datanew.df[itest, 1])
  # Regression coefficient.
  # reg.train.k4[k] = lm(datanew.df[-itest,1] ~ pls.train[,,20])
  #  reg.test.k4[k] = lm(datanew.df[itest,1] ~ pls.test[,,20])
  # RPIQ
  rpiq.train.k4[k] = IQR(datanew.df[-itest, 1])/rmse.train.k4[k]
  rpiq.test.k4[k] = IQR(datanew.df[itest, 1])/rmse.test.k4[k]
}

