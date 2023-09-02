# Martin King. Created 8 Aug 2023.
# Redoing data exploration with new extracted data.
# 1. Distribution of contribution to data (observations).
# 2. Bar charts for categorical variables.
# 3. Histograms.
# 4. Correlation matrix of the 531 spectral values.
# 5. Plotting some basic statistics.
# 6. Cows and their TRT.

rm(list=ls())

data.df = read.csv("/Users/martinpeterking/ucc_courseworks.dir/semester3.dir/data_work.dir/PredMethEMk_DMKComp_cowsinsingleyears_mpkedited.csv", header=T)

# 1. Distribution of contribution to data (observations).--------------------------

cowids = unique(data.df$TB_NUM) 
print(length(cowids))

cowcontrinum = c(rep(0,length(cowids)))

for (i in 1:length(cowids))
{
  cowcontrinum[i] = length(data.df$TB_NUM[data.df$TB_NUM == cowids[i]])
}
par(cex.main = 1.8)
par(cex.axis = 1.8) 
par(cex.lab = 1.8)
par(mar = c(5, 5, 4, 2)) #Bottom, then clockwise.

hist(cowcontrinum, freq = TRUE, breaks = seq(0,30,1), ylab="No. of cows", 
     xlab="No. of observations", main="Number of cows vs number of observations")

# 2. Bar charts for categorical variables.-----------------------------------------

par(cex.main = 1.8)
par(cex.axis = 1.4) 
par(cex.lab = 1.8)
par(mar = c(5, 5, 4, 2)) #Bottom, then clockwise.

par(mfrow=c(1,2))
categories <- c(data.df$TRT)
category_counts <- table(categories)
# Add count labels
text(x = barplot(category_counts, xlab="Treatment", ylab="Count", main="", ylim=c(0,1400)),
     y = category_counts,
     label = category_counts,
     pos = 3, cex=1.3)

#categories <- c(data.df$milking_time)
#category_counts <- table(categories)
# Add count labels
#text(x = barplot(category_counts, xlab="milking_time", ylab="Count", main="milking_time counts", ylim=c(0,2000)),
#     y = category_counts,
#     label = category_counts,
#     pos = 3, cex=1.6)

categories <- c(data.df$parity)
category_counts <- table(categories)
# Add count labels
text(x = barplot(category_counts, xlab="Parity", ylab="Count", main="", ylim=c(0,1000)),
     y = category_counts,
     label = category_counts,
     pos = 3, cex=1.3)

# 3. Histograms.-------------------------------------------------------------------

par(mfrow=c(2,2))
par(cex.main = 1.8)
par(cex.axis = 1.8) 
par(cex.lab = 1.8)
par(mar = c(5, 5, 4, 2)) #Bottom, then clockwise.

hist(data.df$avg_ch4_smilk, breaks=30, xlab="Average methane, g/day", main="", ylab="Count")
hist(data.df$sd_ch4_smilk, breaks=20, xlab="Methane SD, g/day", main="", ylab="Count")
hist(data.df$recs_smilk, breaks=20, xlab="No. of measurements", main="", ylab="Count")
hist(data.df$time_in_machine_s, breaks=20, xlab="Time in machine, sec", main="", ylab="Count")

par(mfrow=c(3,2))
par(cex.main = 1.8)
par(cex.axis = 1.8) 
par(cex.lab = 1.8)
par(mar = c(5, 5, 4, 2)) #Bottom, then clockwise.

hist(data.df$yield, breaks=20, xlab="Yield", main="", ylab="Count")
hist(data.df$fat, breaks=20, xlab="Fat, %", main="", ylab="Count")
hist(data.df$protein, breaks=23, xlab="Protein, %", main="", ylab="Count")
hist(data.df$lactose[-496], breaks=20, xlab="Lactose, %", main="", ylab="Count")
hist(data.df$dim, breaks=20, xlab="Days in milk", main="", ylab="Count")

# 4. Correlation matrix of the 531 spectral values.--------------------------------

corspec = cor(data.df[,19:ncol(data.df)])
count = 0.
for (col in 1:530)
{
  rowstart = col + 1
  rowsel = which(abs(corspec[rowstart:531,col])>=0.95)
  count = count + length(corspec[rowsel,col])
}
print(count)
#count = 22282 out of (531x531-531)/2 = 140715
par(mfrow=c(1,1))
par(cex.main = 1.8)
par(cex.axis = 1.8) 
par(cex.lab = 1.8)
library(corrplot)
corrplot(corspec, method = "color", tl.pos="n", cl.axis=3.0)

# 5. Plotting some basic statistics.-----------------------------------------------

datanew.df = data.frame(data.df[,19:ncol(data.df)])

# Remove the word index from column names.
colnames(datanew.df) <- gsub("index", "", colnames(datanew.df))
#colnames(dataspec.df) <- gsub("index", "", colnames(dataspec.df))
#yzerodataspec = c(rep(0,ncol(dataspec.df)))

par(mfrow=c(1,1))
par(cex.main = 1.8)
par(cex.axis = 1.8) 
par(cex.lab = 1.8)
plot(colnames(datanew.df[,1:531]), apply(datanew.df[,1:531], 2, mean), col="azure4", ylab="", xlab="Index")
points(colnames(datanew.df[,1:531]), apply(datanew.df[,1:531], 2, max), col="coral1")
points(colnames(datanew.df[,1:531]), apply(datanew.df[,1:531], 2, min), col="deepskyblue")
#points(colnames(dataspec.df[,1:ncol(dataspec.df)]), yzerodataspec, pch=15, cex=1.2, col="gray")
title(main="Max, mean, min MIR spectral values")

# 6. Cows and their TRT.--------------------------------------------------------

cowids = unique(data.df$TB_NUM)  
cowuniquetrtnum = c(rep(0,length(cowids)))
count = 1
for (i in 1:length(cowids))
{
  # All obs for a particular cow.
  obsforacow = which(data.df$TB_NUM == cowids[i])
  # How many unique treatments for this cow.
  cowuniquetrtnum[i] = length(unique(data.df$TRT[obsforacow]))
  if (cowuniquetrtnum[i]==3)
  {
    print(count)
    print(cowids[i]) 
    print(unique(data.df$TRT[obsforacow]))
    count = count + 1
  }
}

# 218 on just one TRT.
# 0 on two TRT.
# 0 on three TRT.

# 7. Cows and their parities.---------------------------------------------------

cowids = unique(data.df$TB_NUM)  
cowuniqueparitynum = c(rep(0,length(cowids)))
count = 1
for (i in 1:length(cowids))
{
  # All obs for a particular cow.
  obsforacow = which(data.df$TB_NUM == cowids[i])
  # How many unique parities for this cow.
  cowuniqueparitynum[i] = length(unique(data.df$parity[obsforacow]))
  if (cowuniqueparitynum[i]==1)
  {
    print(count)
    print(cowids[i]) 
    print(unique(data.df$parity[obsforacow]))
    count = count + 1
  }
}





