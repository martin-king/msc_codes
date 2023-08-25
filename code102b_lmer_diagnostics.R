# Martin King. Created 23 Aug 2023.
# Linear Mixed Effects Model.
# 1. Diagnostics

# Read original file. 
rm(list=ls())
library(lme4)
data.df = read.csv("/Users/martinpeterking/ucc_courseworks.dir/semester3.dir/data_work.dir/PredMethEMk_DMKComp_cowsinsingleyears_mpkedited.csv", header=T)

set.seed(2023)

data.new.df = data.frame(data.df$avg_ch4_smilk, data.df$dim, data.df$yield, data.df$sd_ch4_smilk)
#, data.df$fat, data.df$protein, data.df$lactose)
data.new.df$trt = as.factor(data.df$TRT)
data.new.df$tb_num = as.factor(data.df$TB_NUM)
names(data.new.df) = c("avg_ch4_smilk", "dim", "yield", "sd", "trt", "tb_num") 
#"fat", "protein", "lactose", "trt", "tb_num")

# PC for train data. (Not used for now)
#data.spec.df = data.frame(data.df[, 19:ncol(data.df)])
#mu = apply(data.spec.df, 2, mean)
#Xminusmu = data.spec.df - matrix(mu, nrow = nrow(data.spec.df), ncol = ncol(data.spec.df), byrow = TRUE)
#pc.mir = prcomp(Xminusmu, center=FALSE, scale=FALSE)
#data.train.df = data.frame(data.new.df, pc.mir$x[,1:40])

data.train.df = data.new.df

# The 'standard' linear regression.
lm.model <- lm(avg_ch4_smilk ~ yield + dim, data=data.train.df)
lme.model.reduced <- lmer(avg_ch4_smilk ~ 1 + (1 | tb_num), data=data.train.df, )
lme.model <- lmer(avg_ch4_smilk ~ yield + dim + (1 | tb_num), data=data.train.df)

sqrt(mean(residuals(lme.model)^2))

summary(lm.model)
summary(lme.model)

anova(lme.model.reduced, lme.model)

par(mfrow=c(2,2))
par(cex.main = 1.5)
par(cex.axis = 1.4) 
par(cex.lab = 1.4)
par(mar = c(5, 5, 4, 2)) #Bottom, then clockwise.

library(scales)
# Residuals vs. fitted values.
plot(fitted(lme.model), residuals(lme.model), main="Residuals vs. fitted values", xlab="Fitted values", ylab="Residuals", col = alpha("black", 0.3))
abline(h = 0, col = "coral1", lwd = 2)

qqnorm(residuals(lme.model), pch=16, col = alpha("black", 0.3))
qqline(residuals(lme.model))

plot(data.new.df$yield, residuals(lme.model), main="Residuals vs. yield", xlab="yield", ylab="Residuals", col = alpha("black", 0.3))
abline(h = 0, col = "coral1", lwd = 2)

plot(data.new.df$dim, residuals(lme.model), main="Residuals vs. dim", xlab="dim", ylab="Residuals", col = alpha("black", 0.3))
abline(h = 0, col = "coral1", lwd = 2)


