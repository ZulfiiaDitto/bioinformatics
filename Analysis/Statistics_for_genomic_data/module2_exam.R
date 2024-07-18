library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)
library(broom)
library(limma)

# question 1
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# no transformation 
svd1 = svd(edata)
names(svd1) 
pc = svd1$d^2/sum(svd1$d^2)
pc[1] # 0.8873421

#log transformation 
edata_log = log2(edata+1)
svd2 = svd(edata_log)
pc_log = svd2$d^2/sum(svd2$d^2)
pc_log[1] # 0.9737781

#log2(data + 1) transform and subtract row means
edata_cent = edata_log - rowMeans(edata_log)
svd3 = svd(edata_cent)
pc_cent = svd3$d^2/sum(svd3$d^2)
pc_cent[1] #0.3463729

# question 2
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
set.seed(333) # set seed 
edata_log = log2(edata +1) # log data 
edata_cent = edata_log- rowMeans(edata_log) # centr data 
# clustering 
kmeanC = kmeans(t(edata_cent), centers = 2)
names(kmeanC)
svd1 = svd(edata_cent)

cor.test(svd1$v[,1], kmeanC$cluster) #0.8678247 

#Question 3 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
# we have very unbalanced dataset,
#meaning majority of the replicates are 2, but we have 5, 6 and count of those is small

lm1 = lm(edata[1,] ~ as.factor(pdata_bm$num.tech.reps))
broom::tidy(lm1)

# plot the data
plot(pdata_bm$num.tech.reps,edata[1,])
abline(lm1$coefficients[1], lm1$coefficients[3] )

# question 4
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

lm2 = lm(edata[1,] ~ pdata_bm$age + pdata_bm$gender)
broom::tidy(lm2)

# question 5 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
# log transform 
edata_log = log(edata+1)
# fitting the regression model for each sample based on the population outcome 
mod1 = model.matrix(~ pdata$population)
fit1 = lm.fit(mod1, t(edata_log))

dim(fit1$residuals) # 129 52580

dim(fit1$coefficients) #2 52580

dim(fit1$effects) # 129 52580

# question6 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
dim(pdata)
dim(edata)
edata_log = log2(edata+1)

mod1 = model.matrix(~ pdata$population)
fit1 = lm.fit(mod1, t(edata_log))

dim(fit1$residuals) # 129 52580

dim(fit1$coefficients) #2 52580

dim(fit1$effects) # 129 52580 - samples genes

fit1$effects[,1]

#question 7 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

# cleaning the Na in age
pdata_bm = na.omit(pdata_bm)
edata = edata[, rownames(pdata_bm), drop = FALSE]

# many regressions

model_adj = model.matrix(~ pdata_bm$age)
fit_l = lmFit(edata, model_adj)
fit_l$coefficients[1000] 

summary(fit_l)

intercept = fit_l$coefficients[1000,][1]
slope = fit_l$coefficients[1000,][2]
x = edata[1000,]*slope+intercept

plot(x,pdata_bm$age) # not a great fit 

# Question 8 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

# cleaning the Na in age
pdata_bm = na.omit(pdata_bm)
edata = edata[, rownames(pdata_bm), drop = FALSE]

mod_adj = model.matrix(~pdata_bm$age + pdata_bm$tissue.type)

fit_l = lmFit(edata, mod_adj)

fit_l$cov.coefficients

pdata_bm$tissue.type # we have 17 different tiddue types 
dim(edata) # 16 samples -> we have too many co-vaents in tissue type for 16 samples 

# question9 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

pdata$population

pdata$study
# see the lecture for answer 

# question 10 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

set.seed(33353)
# removing age NA
pheno = na.omit(pdata_bm)
edata = edata[,rownames(pheno), drop=FALSE]
# log transform
edata_log = log2(edata+1)
#rowMEans >1 
edata_cent = edata_log[rowMeans(edata_log) >1,]
#sva
mod = model.matrix(~age, data=pheno)
mod0 = model.matrix(~1, data=pheno)
sva1 = sva(edata_cent, mod,mod0, n.sv=2)
names(sva1)

cor(sva1$sv, pheno$age)

# correlation with gender 
cor(sva1$sv, as.numeric(pheno$gender))

#cor with race
cor(sva1$sv, as.numeric(pheno$race))


