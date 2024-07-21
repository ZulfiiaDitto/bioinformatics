library(devtools)
library(Biobase)
library(snpStats)
library(MASS)
library(broom)
library(DESeq2)

tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=18) 


data(for.exercise)

use = seq(1, ncol(snps.10),10)
sub.10 = snps.10[, use]
xxmat <-- xxt(sub.10, correct.for.missing = FALSE)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5]


snpdata = sub.10@.Data
status = subject.support$cc
snp1 = as.numeric(snpdata[,1])
snp1[snp1==0] = NA # missing values was as 0,
#replace with NA so R will take care of it during calculation 

glm1 = glm(status ~ snp1,family="binomial") # to make sure it will calcualte as a logistic regression 
tidy(glm1)

# we can fit the dominant model as well 
snp1_dom = (snp1 == 1) # dummy variable 
glm1_dom = glm(status ~ snp1_dom,family="binomial")
tidy(glm1_dom)
tidy(glm1) 


# adjusting for principal component 
glm2 = glm(status ~ snp1 + pcs[,1:5],family="binomial")
tidy(glm2)
# adjusting for many variables 
glm_all = snp.rhs.tests(status ~ 1,snp.data=sub.10)
slotNames(glm_all)
qq.chisq(chi.squared(glm_all),df=1)


# adjusted for variables and pcs
glm_all_adj = snp.rhs.tests(status ~ pcs,snp.data=sub.10)
qq.chisq(chi.squared(glm_all_adj),df=1)


# for poison distribution 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)
ls()


## ------------------------------------------------------------------------
edata = edata[rowMeans(edata) > 10, ]

## ------------------------------------------------------------------------
glm3 = glm(edata[1, ] ~ pdata$strain,family="poisson")
tidy(glm3) # tells us on log scales 

# negative binomial 
glm.nb1 = glm.nb(edata[1, ] ~ pdata$strain)
tidy(glm.nb1)

# many negative binomial regressions 
de = DESeqDataSetFromMatrix(edata, pdata, ~strain)
glm_all_nb = DESeq(de) # estimates smooth relationship between the variance and mean 
result_nb = results(glm_all_nb)
hist(result_nb$stat)

# comparing statistics 

library(limma)
library(edge)
library(genefilter)
 

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)

# normalize data 
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

# tstats for 2 group comparisons 
tstats_obj = rowttests(edata,pdata$strain)
names(tstats_obj)
hist(tstats_obj$statistic,col=2)
# for  multigroup comparison 
fstats_obj = rowFtests(edata,as.factor(pdata$lane.number))
names(fstats_obj)
hist(fstats_obj$statistic,col=2)

# no adjustment version of moderated statistics useng limma 
mod = model.matrix(~ pdata$strain)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
head(ebayes_limma$t)

plot(ebayes_limma$t[,2],-tstats_obj$statistic,col=4, xlab="Moderated T-stat",ylab="T-stat")
abline(c(0,1),col="darkgrey",lwd=3)

# adjasted model 
mod_adj = model.matrix(~ pdata$strain + as.factor(pdata$lane.number))
fit_limma_adj = lmFit(edata,mod_adj)
ebayes_limma_adj = eBayes(fit_limma_adj)
head(ebayes_limma_adj$t)


plot(ebayes_limma_adj$t[,2],-tstats_obj$statistic,col=3,
     xlab="Moderated T-stat",ylab="T-stat")
abline(c(0,1),lwd=3,col="darkgrey")

# another option to model multifactor statistical models 
mod_lane = model.matrix(~ as.factor(pdata$lane.number))
fit_limma_lane = lmFit(edata,mod_lane)
ebayes_limma_lane = eBayes(fit_limma_lane) # emperical way to estimate statistics 
head(ebayes_limma_lane$t)


top_lane = topTable(ebayes_limma_lane, coef=2:7,
                    number=dim(edata)[1],sort.by="none")
head(top_lane)

plot(top_lane$F,fstats_obj$statistic,
     xlab="Moderated F-statistic",ylab="F-statistic",col=3)

# same thing in edge 
edge_study = build_study(edata, grp = as.factor(pdata$lane.number))
de_obj = lrt(edge_study) # likelihood ration test
qval = qvalueObj(de_obj)
plot(qval$stat,fstats_obj$statistic,col=4,
     xlab="F-stat from edge",ylab="F-stat from genefilter") # genefilter and edge


# adjusting with edge 
edge_study2 = build_study(edata, grp = as.factor(pdata$lane.number),
                          adj.var=pdata$strain)
de_obj2 = lrt(edge_study2)
qval2 = qvalueObj(de_obj2)
plot(qval2$stat,fstats_obj$statistic,col=4,
     xlab="F-stat from edge",ylab="F-stat from genefilter")


## permutation 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)

edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

tstats_obj = rowttests(edata,pdata$strain)
hist(tstats_obj$statistic,col=2,xlim=c(-5,2))


set.seed(135)
strain = pdata$strain
strain0 = sample(strain)
tstats_obj0 = rowttests(edata,strain0)
hist(tstats_obj0$statistic,col=2,xlim=c(-5,2))

quantile(tstats_obj0$statistic)
quantile(tstats_obj$statistic)

# calculation of P values ad correction for different error rates 

library(qvalue)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)

edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

fstats_obj = rowFtests(edata,as.factor(pdata$strain))
hist(fstats_obj$p.value,col=2)

# same with the edge package 
edge_study = build_study(edata, grp = pdata$strain, 
                         adj.var = as.factor(pdata$lane.number))
de_obj = lrt(edge_study)
qval = qvalueObj(de_obj)
hist(qval$pvalues,col=3) # model is not quite right based on the grafic 

# moderated  p values 
mod = model.matrix(~ pdata$strain + pdata$lane.number)
fit_limma = lmFit(edata,mod)
ebayes_limma = eBayes(fit_limma)
limma_pvals = topTable(ebayes_limma,number=dim(edata)[1])$P.Value
hist(limma_pvals,col=4)
# permutation -> recalculation of T test with the permutation 
set.seed(3333)
B = 1000
tstats_obj = rowttests(edata,pdata$strain)
tstat0 = matrix(NA,nrow=dim(edata)[1],ncol=B)
tstat = tstats_obj$statistic
strain = pdata$strain
for(i in 1:B){
  strain0 = sample(strain)
  tstat0[,i] = rowttests(edata,strain0)$statistic
}

emp_pvals = empPvals(tstat,tstat0)
hist(emp_pvals,col=2)

# eroor correction for Bonferony = family wise discovery rate 

fp_bonf = p.adjust(fstats_obj$p.value,method="bonferroni")
hist(fp_bonf,col=3)
quantile(fp_bonf)
sum(fp_bonf < 0.05) # = zero -> no statistically significant 

# error correction for the BH - false discovery rate  
fp_bh = p.adjust(fstats_obj$p.value,method="BH")
hist(fp_bh,col=3)
quantile(fp_bh)
sum(fp_bh <0.05) # nothing is significant with the BH adj 

# same with limma package 
limma_pvals_adj = topTable(ebayes_limma,number=dim(edata)[1])$adj.P.Val
hist(limma_pvals_adj,col=2)
quantile(limma_pvals_adj)

## ------------------------------------------------------------------------
qval_limma = qvalue(limma_pvals)
summary(qval_limma)
qval$pi0

## ------------------------------------------------------------------------
qval = qvalueObj(de_obj)
summary(qval)





