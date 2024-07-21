# question 1 

library(snpStats)
library(broom)

data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
snp3 = as.numeric(snpdata[,3])
snp3[snp3==0] = NA 
glm1 = glm(status ~ snp3,family="binomial") # to make sure it will calcualte as a logistic regression 
tidy(glm1)

glm2 = glm(status ~ snp3) # to make sure it will calcualte as a linear regression 
tidy(glm2)
# linear = -0.039 
# log = -1.158
#Both models are fit on the additive scale. So in the linear model case, the coefficient is the decrease in probability associated with each additional copy of the minor allele. In the logistic regression case, it is the decrease in the log odds ratio associated with each additional copy of the minor allele.

# question 2 

par(mfrow = c(1,2))

plot(status ~ snp3, pch = 19 )
abline(glm2 , col = "blue", lwd= 5 )
plot(glm1$residuals)

# look at the linear graph = first one. 
# more values -> means negative values 

# question 3 
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
# 10th snp 
snp10 = as.numeric(snpdata[,10])
snp10[snp10==0] = NA 
# log regression 
glm1 = glm(status ~ snp10,family="binomial") # to make sure it will calcualte as a logistic regression 
tidy(glm1)

# dominant one 
snp10_dom = (snp10 == 2)
glm_dom = glm(status ~ snp10_dom,family="binomial") # to make sure it will calcualte as a logistic regression 
tidy(glm_dom)

plot(glm1$residuals)
plot(glm_dom$residuals)

# no differences

# question 4 
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
# create variable to keep the values of effects  
results = rep(NA, dim(snpdata)[2])

# create the for-loop 

for(i in 1:ncol(snpdata)){
  snpdata_loc = as.numeric(snpdata[,i])
  snpdata_loc[snpdata_loc ==0] = NA ## replace nulls
  glm_loc = glm(status ~ snpdata_loc, family = 'binomial')
  results[i] = tidy(glm_loc)$statistic[2]
}
mean(results)
# 0.007155377
max(results)
# 3.900891
min(results)
# -4.251469 

# quesiton 5
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

# create variable to keep the values of effects  
results = rep(NA, dim(snpdata)[2])

# create the for-loop 
for(i in 1:ncol(snpdata)){
  snpdata_loc = as.numeric(snpdata[,i])
  snpdata_loc[snpdata_loc ==0] = NA ## replace nulls
  glm_loc = glm(status ~ snpdata_loc, family = 'binomial')
  results[i] = tidy(glm_loc)$statistic[2]
}

results_squar = results^2

# correlation

glm_all = snp.rhs.tests(status ~1, snp.data = sub.10) # default value is binomial -> 
cor(results_squar, chi.squared(glm_all))
# answer: 0.99. They are both testing for the same association using the same additive regression model on the logistic scale but using slightly different tests. 

# question 6 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
edata = log2(as.matrix(edata) + 1)

# Ftest using gene filter 
fstats_obj = rowFtests(edata,as.factor(pdata$population))
names(fstats_obj)
hist(fstats_obj$statistic,col=2)
tidy(fstats_obj)
# ttest 
tstats_obj = rowttests(edata,pdata$population)
names(tstats_obj)
hist(tstats_obj$statistic,col=2)
tidy(tstats_obj)
# same p values, different statistics 

# question 7 
library(DESeq2)
library(limma)
library(edge)
library(genefilter)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
edata = edata[rowMeans(edata) > 100,]
fdata = fData(mp)

# using the DesSeq 
de = DESeqDataSetFromMatrix(edata, pdata, ~study)
glm_de = DESeq(de)
result_de = results(glm_de)

#limma 
edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~ as.factor(pdata$study))
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma) 
top = topTable(ebayes_limma,number=dim(edata)[1], sort.by="none")

cor(result_de$stat, top$t)
# 0.9278568

# make an MA-plot
y = cbind(result_de$stat, top$t)
limma::plotMA(y)

# question 9 
# desseq
fp_bh = p.adjust(result_de$pvalue, method="BH")
sum(fp_bh < 0.05)
#1995
#limma 
fp_bh = p.adjust(top$P.Value, method="BH")
sum(fp_bh < 0.05)
#2807




