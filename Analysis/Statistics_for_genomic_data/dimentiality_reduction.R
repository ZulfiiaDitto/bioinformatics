#millions measurements(features) per sample,
# to better communicate patterns and identify relationship -> best way to do it 
# to reduce the dimentiality 
library(devtools)
library(Biobase)
library(RSkittleBrewer)

# color schema from RSkittleBrewer
tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=18) # set up filled circles

# loading the data 
fileUrl = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file = fileUrl)
close(fileUrl)

mp = montpick.eset
pdata = pData(mp)
fdata = fData(mp)
edata = as.data.frame(exprs(mp))

# filtering data 
edata = edata[rowMeans(edata) > 100, ]
# normalizing data 
edata = log2(edata +1)
# centering data for single value decomposition 
# if this is not done than first single value vector always will be the most variation one and will represent the mean 
edata_centered = edata - rowMeans(edata)

# svd = single value decomposition 

svd1 = svd(edata_centered)

names(svd1) # returns 3 names of matrixs: 'd' = diagonal, 'u' - variation across samples,
#'v'- variation across genes

# plotting

plot(svd1$d, ylab='Singular values', col = 4)

plot(svd1$d^2/sum(svd1$d^2), ylab = 'Percent varience explained', col = 2)
# the first value in a plot -> explained 50% of variant. it is highly explanatory varieble  

# ploting the samples 
par(mfrow = c(1,2))
plot(svd1$v[,1], col = 2, ylab = '1nd PC')
plot(svd1$v[,2], col = 2, ylab = '2st PC')

par(mfrow = c(1,1))
plot(svd1$v[,1], svd1$v[,2],
     ylab = '1nd PC', xlab = '2nd PC', col = as.numeric(pdata$study) )
# the difference very depends on the what study they come from 

boxplot(svd1$v[,1] ~ pdata$study,  boarder = c(1,2))
points(svd1$v[,1] ~jitter(as.numeric(pdata$study)), col = as.numeric(pdata$study))
# you can also utilize the build in function 

pc1 = prcomp(edata)
plot(pc1$rotation[,1], svd1$v[,1])

# scaling data 
edata_centered2 = t(t(edata)- colMeans(edata))

svd2= svd(edata_centered2)
plot(pc1$rotation[,1], svd2$v[,1], col = 3) # you can see that pc1 and svd from the columns centred is same. cuz pc is calcualted based on the column values 

# outlier also can infere pc 

edata_outlier = edata_centered
edata_outlier[6,] = edata_centered[6,] * 10000

svd3 = svd(edata_outlier)
plot(svd1$v[,1], svd3$v[,1], xlab = "without outlier", ylab = 'with outlier')

# normalization 

library(preprocessCore)

# reload the file 
fileUrl = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file = fileUrl)
close(fileUrl)

mp = montpick.eset
pdata = pData(mp)
fdata = fData(mp)
edata = as.data.frame(exprs(mp))
# preprocess
edata = log2(edata +1)
edata = edata[rowMeans(edata) > 3, ]
dim(edata)

# lets see the distribution of the first sample
colramp = colorRampPalette(c(3,'white',2))(20)
plot(density(edata[,1]), col = colramp[1], lwy = 3, ylim = c(0, .3))
for(i in 2:20){lines(density(edata[,i]), lwd = 3,col = colramp[i])}

# the graph shows us that there are some differences between
# samples and it is possible it is due to technology -> so we are doing the 
# quantile normalization 

norm_edata = normalize.quantiles(as.matrix(edata))
plot(density(norm_edata[,1]),
     col = colramp[1],lwd = 3, ylim = c(0, 0.2) )
for(i in 2:20){lines(density(norm_edata[,i]), lwd = 3,col = colramp[i])}

# there is small variability in the lower left of the graph it is 
# due to the small values it is very hard to match the procentile 
# quantile normalization removes the bulk of the differences but not the 
# gene to gene variations

# by ploting the norm_data by the study you can see th egene vs gene variation 
plot(norm_edata[1,], col = as.numeric(pdata$study))

# let do the svd
svd1 = svd(norm_edata- rowMeans(norm_edata))
plot(svd1$v[,1], svd1$v[,2], xlab = 'PC1', ylab= 'PC2', 
     col = as.numeric(pdata$study))
# analyzing the graph -> still have gene to gene variation between studies 
# meaning even if we normalize data we still have batch effect, or other type of artifacts


########## linear regression 
library(broom)
library(tidyr)

# loading the data 
fileUrl = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file = fileUrl)
close(fileUrl)

# expretion data 
exp_data = exprs(bodymap.eset)
dim(exp_data) # 52580 different features = genes = rows, 19 samples/patients = columns
head(exp_data)
# phenotype data 
pheno = pData(bodymap.eset)
dim(pheno) # rows = 19 should match the columns in expression table
head(pheno)

# feature table 
feature = fData(bodymap.eset)
dim(feature) # the row should match the row count of the exression data 
head(feature)

exp_data = as.matrix(exp_data)

lm1 = lm(exp_data[1,] ~pheno$age)
lm1

broom::tidy(lm1) 
# ploting the gene expression of first sample 
plot(pheno$age , exp_data[1,], col = 1)
abline(lm1, col = 2, lmd = 3)

# let review the gender covarient 

table(pheno$gender)

boxplot(exp_data[1,] ~ pheno$gender)
points(exp_data[1,] ~ jitter(as.numeric(pheno$gender)), 
       col = as.numeric(pheno$gender))
# looking into boxplot - there are variability btween points in gender groups
# but how we can quantify this? by creating the dummy data , 
#R in linear regretiiion does it for you 

dummy_m = pheno$gender == "M"
dummy_m*1

dummy_f = pheno$gender == "F"
dummy_f*1

lm2 = lm(exp_data[1,] ~pheno$gender)
broom::tidy(lm2)
# to see what happend under the hood. 
# the model converts M gender to the 1 and F to 0 
mod2 = model.matrix(~ pheno$gender)
mod2

# lets see the tissue variables 
table(pheno$tissue.type)

pheno$tissue.type == "adipose"
pheno$tissue.type == 'adrenal'

# if you tissue type into model 
#- intercept will actually show you the 
broom::tidy(lm(exp_data[1,] ~ pheno$tissue.type))

##lets expend the model formula 
lm3 = lm(exp_data[1,] ~pheno$age + pheno$gender)

broom::tidy(lm3)

# to apply interaction model you ise the * within the covariants

lm4 =  lm(exp_data[1,] ~pheno$age*pheno$gender)
broom::tidy(lm4)

# lets plot our data points

lm4 = lm(exp_data[6,] ~ pheno$age)
plot(pheno$age, exp_data[6,], col =2)
abline(lm4, col = 1, lwd=3)
# in this case th eputlier does not pull the line/ interfere with the line 

# in order to see where putliers are interfere with the regresion line 
# lets do the following
index = 1:19 
lm5 = lm(exp_data[6,] ~index)
plot(index, exp_data[6,], col = 2)
abline(lm5, col=1, lwd = 3)
# line is slightly pulled up -> lets remove the point = outlier and remodel it 
lm6 = lm(exp_data[6,-19] ~ index[-19])
abline(lm6, col = 3, lwd = 3)

legend(5,1000, c('with outlier', 'without outlier'), col = c(1,3), lwd = 3)

par(mfrow = c(1, 2))
hist(lm6$residuals)
hist(lm5$residuals)
# the hist of not very pretty aligned - we need to do the normalization 
gene1 = log2(exp_data[1,] + 1)
lm7 = lm(gene1 ~ index)
hist(lm7$residuals, col = 4)

# lets fill the models with as many co-efficients as we want 

lm8=lm(gene1~pheno$tissue.type+pheno$age)
broom::tidy(lm8)
# the stats for this model doe snot calcualted cuz you have 16 points and 18 covarients 
# too many co-varients for small amount of points
dim(model.matrix(~pheno$tissue.type+pheno$age))

lm9 = lm(exp_data[2,]~ pheno$age)
plot(lm9$residuals, col = colramp[as.numeric(pheno$tissue.type)])


##### many linear models 
library(limma)
library(edge)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata=pData(bot)
edata=as.matrix(exprs(bot))
fdata = fData(bot)
ls()
# removing the low genes 
edata = log2(as.matrix(edata) + 1)
edata = edata[rowMeans(edata) > 10, ]

mod = model.matrix( ~ pdata$strain)
fit = lm.fit(mod, t(edata))
names(fit)

fit$coefficients[,1]
head(fit$coefficients)
broom::tidy(lm(as.numeric(edata[1, ]) ~ pdata$strain))

par(mfrow=c(1,2))
hist(fit$coefficients[1,],breaks=100,col=2,xlab="Intercept")
hist(fit$coefficients[2,],breaks=100,col=2,xlab="Strain")
abline(v=0,lwd=3,col=1)

mod_adj = model.matrix(~pdata$strain + as.factor(pdata$lane.number))
fit_adj = lm.fit(mod_adj, t(edata))
fit_adj$coefficients[,1]

### doing it with limma package
fit_limma = lmFit(edata, mod_adj)
names(fit_limma)
fit_limma$coefficients[1,]
fit_adj$coefficients[,1]

edge_stady = build_study(data = edata, grp = pdata$strain, 
                         adj.var = as.factor(pdata$lane.number))

fit_edge = fit_models(edge_stady)

summary(fit_edge)

fit_edge@beta.coef[1,]

####### batch effects 
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)


data(bladderdata)

pheno = pData(bladderEset)
edata = exprs(bladderEset)

mod = model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
fit = lm.fit(mod,t(edata))
hist(fit$coefficients[2,],col=2,breaks=100)

table(pheno$cancer,pheno$batch)

# how to adjast for the batch effect 
# if we know the batch effect reasons do the Combat 
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
modcancer = model.matrix(~cancer, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
combat_fit = lm.fit(modcancer,t(combat_edata))
hist(combat_fit$coefficients[2,],col=2,breaks=100)

plot(fit$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="Linear Model",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)


# of we dont know the batch effect 
mod = model.matrix(~cancer,data=pheno)
mod0 = model.matrix(~1, data=pheno)
sva1 = sva(edata,mod,mod0,n.sv=2)
names(sva1)
dim(sva1$sv) # it is found potential batche effect 
summary(lm(sva1$sv ~ pheno$batch))
boxplot(sva1$sv[,2] ~ pheno$batch)
points(sva1$sv[,2] ~ jitter(as.numeric(pheno$batch)),col=as.numeric(pheno$batch))

# with sva package we actually found new co-variants in our data, we did not cleant them . 
# so lets add them into our model 


modsv = cbind(mod,sva1$sv)
fitsv = lm.fit(modsv,t(edata))

par(mfrow=c(1,2))
plot(fitsv$coefficients[2,],combat_fit$coefficients[2,],col=2,
     xlab="SVA",ylab="Combat",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)
plot(fitsv$coefficients[2,], fit$coefficients[2,],col=2,
     xlab="SVA",ylab="linear model",xlim=c(-5,5),ylim=c(-5,5))
abline(c(0,1),col=1,lwd=3)


# biological varibility 

data(for.exercise)
print(data(for.exercise))
controls <- rownames(subject.support)[subject.support$cc==0]
controls
use <- seq(1, ncol(snps.10), 10)
ctl.10 <- snps.10[controls,use]
# calculation of principal component 
xxmat <- xxt(ctl.10, correct.for.missing=FALSE)
evv <- eigen(xxmat, symmetric=TRUE)
pcs <- evv$vectors[,1:5]


pop <- subject.support[controls,"stratum"]
plot(pcs[,1],pcs[,2],col=as.numeric(pop),
     xlab="PC1",ylab="PC2")
legend(0,0.15,legend=levels(pop),pch=19,col=1:2)
