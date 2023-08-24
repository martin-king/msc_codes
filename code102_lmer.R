# Martin King. Created 4 Aug 2023.
# Linear Mixed Effects Model.
# 0. Just fit LMM using whole data.
# 1. lmer fixed AND random effects predictions.

# 0. Just fit LMM using whole data.---------------------------------------------
rm(list=ls())

data.df = read.csv("/Users/martinpeterking/ucc_courseworks.dir/semester3.dir/data_work.dir/PredMethEMk_DMKComp_cowsinsingleyears_mpkedited.csv", header=T)

# Using LMER
# Random effects for measurements among individual cows.
data.plus.df = data.frame(data.df$avg_ch4_smilk, 
                          data.df$yield, data.df$dim)
data.plus.df$id = as.factor(data.df$TB_NUM)
data.plus.df$trt = as.factor(data.df$TRT)
library(lme4)
lmer.1 <- lmer(data.df.avg_ch4_smilk ~ . - id - trt +
                 (1 | trt), data=data.plus.df)
sqrt(mean(residuals(lmer.1)^2))

# 1. lmer predictions.--------------------------------------------
# NB. Mean of random effect is zero in our case. Therefore, 
# commented out below in prediction.

# Read original file. 
rm(list=ls())
library(lme4)
data.df = read.csv("/Users/martinpeterking/ucc_courseworks.dir/semester3.dir/data_work.dir/PredMethEMk_DMKComp_cowsinsingleyears_mpkedited.csv", header=T)

set.seed(2023)

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

data.new.df = data.frame(data.df$avg_ch4_smilk, data.df$dim, data.df$yield)
#, data.df$fat, data.df$protein, data.df$lactose)
data.new.df$trt = as.factor(data.df$TRT)
data.new.df$tb_num = as.factor(data.df$TB_NUM)
names(data.new.df) = c("avg_ch4_smilk", "dim", "yield", "trt", "tb_num") 
#"fat", "protein", "lactose", "trt", "tb_num")

grouping = data.new.df$tb_num

for (b in 1:B)
{
  folds = groupKFold(grouping, k = K)
  for (k in 1:K)
  {
    itrain = as.numeric(unlist(folds[k]))
    
    # PC for train data.
    data.spec.df = data.frame(data.df[itrain, 19:ncol(data.df)])
    mu = apply(data.spec.df, 2, mean)
    Xminusmu = data.spec.df - matrix(mu, nrow = nrow(data.spec.df), ncol = ncol(data.spec.df), byrow = TRUE)
    
    pc.mir = prcomp(Xminusmu, center=FALSE, scale=FALSE)
    #data.train.df = data.frame(data.new.df[itrain,], pc.mir$x[,1:40])
    data.train.df = data.frame(data.new.df[itrain,])
    
    #data.spec.df = data.frame(data.df[-itrain, 19:ncol(data.df)])
    #Xminusmu = data.spec.df - matrix(mu, nrow = nrow(data.spec.df), ncol = ncol(data.spec.df), byrow = TRUE)
    
    # Projecting test data on train pcs.
    #pc.test = as.matrix(Xminusmu) %*%  pc.mir$rotation
    #data.test.df = data.frame(data.new.df[-itrain,], pc.test[,1:40])
    #data.test.df = data.frame(data.new.df[-itrain,])
    
    # Fitting model.
    lme.model <- lmer(avg_ch4_smilk ~ . - trt - tb_num + (1 | tb_num), data=data.train.df)
    
    #random_effects <- ranef(lme.model, condVar = TRUE)[[1]]
    #predicted_random_effects <- mean(random_effects[,1])
    #fixed_effects <- fixef(lme.model)
    #trained_responses <- rep(fixed_effects["(Intercept)"], length(data.new.df$tb_num[itrain]))
    #predicted_responses <- rep(fixed_effects["(Intercept)"], length(data.new.df$tb_num[-itrain]))
    #for (i in 2:length(fixed_effects)) 
    #{
    #  predictor <- names(fixed_effects)[i]
    #  trained_responses <- trained_responses + fixed_effects[predictor] * data.new.df[itrain,predictor]
    #  predicted_responses <- predicted_responses + fixed_effects[predictor] * data.new.df[-itrain,predictor]
    #}
    
    # Predicting.
    
    # Forcing to use only fixed effects.
    #lm.train <- trained_responses
    lm.train <- predict(lme.model, newdata = data.train.df)
    
    # Forcing to use only fixed effects.
    #lm.pred <- predicted_responses #+ predicted_random_effects
    
    # Predicting would use only the fixed effects for new individuals.
    #lm.pred <- predict(lme.model, newdata=data.test.df, allow.new.levels = TRUE )
    
    indx = (b-1)*K+k
    
    # RMSE.
    rmse.train.k4[indx] =  sqrt(mean((data.new.df[itrain, 1] - lm.train)^2.0))
    #rmse.test.k4[indx] =  sqrt(mean((data.new.df[-itrain, 1] - lm.pred)^2.0))
    # Bias.
    bias.train.k4[indx] = mean(lm.train) - mean(data.new.df[itrain, 1])
    #bias.test.k4[indx] = mean(lm.pred) - mean(data.new.df[-itrain, 1])
    # Correlation.
    cor.train.k4[indx] = cor(lm.train, data.new.df[itrain, 1])
    #cor.test.k4[indx] = cor(lm.pred, data.new.df[-itrain, 1])
    # Regression coefficient.
    # reg.train.k4[k] = lm(datanew.df[itrain,1] ~ pls.train[,,20])
    # reg.test.k4[k] = lm(datanew.df[-itrain,1] ~ pls.test[,,20])
    # RPIQ
    rpiq.train.k4[indx] = IQR(data.new.df[itrain, 1])/rmse.train.k4[indx]
    #rpiq.test.k4[indx] = IQR(data.new.df[-itrain, 1])/rmse.test.k4[indx]
  }
}  
  
}