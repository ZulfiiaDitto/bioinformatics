#setwd("~/Documents/bioinformatics/Analysis/Statistics_for_genomic_data")

# coding along 
library(devtools)
library(Biobase)
library(gplots)
library(dplyr)
library(RSkittleBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dendextend)

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

# color schema from RSkittleBrewer
tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=18) # set up filled circles

# tables 
table(pheno$gender)
table(pheno$gender, pheno$race) # unbalanced 

summary(exp_data) # gives you descriptive stats for every columns ->
#almost all of them skewed 

table(pheno$age) 
table(pheno$age,useNA = 'ifany') # 3 missing values 

sum(pheno$age == ' ', na.rm = TRUE) # no values with ' ' values 

sum(is.na(exp_data)) # check if there any missing values in expression dataset 
# no missing values in expressio dataset 

# lets check is there missing values by rows= genes in expretion data
gene_na = rowSums(is.na(exp_data))
table(gene_na)
# lets check is there missing values by columns= sample in expretion data
sample_na = colSums(is.na(exp_data))
table(sample_na)

# plotting

boxplot(exp_data[,1])
boxplot(log2(exp_data[,1]+1)) # first columns 

boxplot(log2(exp_data+1), col= 2, range= 0) # data is skewed 

# to look at the each individual sample you can do below 
par(mfrow =c(1,2)) # setting up the screen so we will have 1 row and 2 columns
hist(log2(exp_data[,1]+1), col = 2) # hist for 1st sample
hist(log2(exp_data[,2]+1), col = 2) # hist for 2nd sample

# to look at the all of the samples at once 
#-> density plot to overlay the samples and see if distribution is all similar
par(mfrow = c(1,1))
plot(density(log2(exp_data[,1]+1)),col=3)
lines(density(log2(exp_data[,2]+1)),col=4) # to overay the second sample 
# comparison samples with the qqplot = Quantile-Quantile Plots
qqplot(log2(exp_data[,1]+1),log2(exp_data[,2]+1), col = 3 )
abline(c(0,1))
# another plot is MA or a Bland-Altman plot
ma = log2(exp_data[,1]+1) - log2(exp_data[,2]+1) # y axis
aa = log2(exp_data[,1]+1) + log2(exp_data[,2]+1) # x axis 

plot(aa, ma, col = 2) # so if there is not difference between the samples is will align in zero. 
# this plot is useful to check the replicas difference 

# And so, one thing that you can often do to make the plots a little bit better, 
#especially for count-based data. 
#You often need to remove the low expression or the low count features 
#to be able to really be able to see the distribution of the data.
edata = as.data.frame(exp_data)
filt_edata = filter(edata, rowMeans(edata) >1)
dim(filt_edata) # less in dimentionality 
# lets plot new data 
boxplot(as.matrix(log2(filt_edata +1)), col=2)

# third part of the analysis is to check for consistency 
# comparing metadata with the annotation from something else such as public database of genoms

# lets get the idies of features 
aeid = as.character(feature[,1])
chr = AnnotationDbi::select(org.Hs.eg.db, keys=aeid,
                            keytype = "ENSEMBL",
                            columns = "CHR")
head(chr)
dim(chr) # problem the dim of the chr should be == dim (feature). 
#It is not. dim(chr)> dim(feature) -> issue we have duplicates in chr table

chr = chr[!duplicated(chr[,1]),] # after we removed duplicates -> dim is equals 

all(chr[,1]==rownames(edata))

# select the chromosome Y samples 
edata = as.data.frame(edata)
dim(edata)
filter(edata,chr$CHR=="Y") -> edataY
dim(edataY)

boxplot(colSums(edataY)~pheno$gender) 
points(colSums(edataY) ~ jitter(as.numeric(pheno$gender)),
       col = as.numeric(pheno$gender), 
       pch =18)

# heatplot 
# first let pull only features with the high row means -> high expression 
ematrix = as.matrix(edata)[rowMeans(edata) > 100000,]
heatmap(ematrix)
# lets define new color palette 
colramp = colorRampPalette(c(3, 'white',2))(9)
heatmap(ematrix,col=colramp)
# remove the clustering part of the heatmap
heatmap(ematrix, col = colramp, Rowv = NA, Colv = NA)

heatmap.2(ematrix, col = colramp, Rowv = NA, Colv = NA, 
          dendrogram = 'none', scale = 'row', 
          trace = 'none')

# data tansforms of the small gene expression values 

# this is how norm distribution is look slike 
hist(rnorm(1000), col=1)
hist(edata[,1], col = 2, breaks = 100) # almost all values as zeros :( 

#in our case with highly skewed data, use of the log data is the one way to viz it 
hist(log(edata[,1]), col = 2, breaks = 100) # little bit better But
# min(log(edata)) = -inf
min(log(edata))
quantile(log(edata[,1])) # can not calculate -> a lot of -inf values
# that is why we adding +1 -> it wil not drastically affect any small values 
hist(log(edata[,1]+1), col=1)
# we can also do the log2 transformation 
hist(log2(edata[,1]+1), col = 3)
# setting up x limits and y limits in order to ignore the zeros 
hist(log2(edata[,1]+1), col = 3, breaks = 100, xlim = c(1,15),
     ylim = c(0,400))

# lets see how may is egenes expresisons = zeros 
hist(rowSums(edata ==0), col = 1)

#lets remove the low genes expression values by mean
low_genes = rowMeans((edata)) < 6

table(low_genes)
head(low_genes)
# removing the low genes expresttion by median 
fil_edata = filter(as.data.frame(edata),!low_genes)
dim(fil_edata)

low_genes2 = rowMedians(as.matrix(edata)) <6
table(low_genes, low_genes2)

fil_edata2 = filter(edata, !low_genes2)
dim(fil_edata2)

hist(log2(fil_edata[,1]+1), col =3) # hist filtered by the mean 
hist(log2(fil_edata2[,1]+1), col = 4) # hist filtered by median 


# clustering 
# main idea is to identify data points that are close to
# each other and other cluster those togethergreoup somehow:
#two ways to do it for genomic data:
#1. hierachical clustering
#2. k-mean ( randomly assign position of the clusters center and re-iterate it over)
# filtering out by the row mean (feature mean)
edataFinal = edata[rowMeans(edata) > 5000, ]
dim(edataFinal)
edataFinal = log2(edataFinal +1)

# lets calculate the euclidian distance between rows 
dist = dist(t(edataFinal))
head(dist)

heatmap(as.matrix(dist), col = colramp, Colv = NA, Rowv = NA)

hclustE = hclust(dist)
plot(hclustE, hang = -1) # we see what **samples** are close to each other 

# lets define the clusters using the colors as well 
dend = as.dendrogram(hclustE)
dend = color_labels(hclustE, 4, 1:4) # here we define how many clusters
plot(dend) # you see 4 different color == 4 different clusters 

# k-mean clustering 

kmeanC = kmeans(edataFinal, centers = 3)
names(kmeanC)
matplot(t(kmeanC$centers), col = 1:3, type='l', lwd = 3) # review the centers in data 
table(kmeanC$cluster)
kmeanC$cluster[1:10]
# re-order the data based on the clusters 
newdata = as.matrix(edataFinal)[order(kmeanC$cluster),]
dim(newdata)

heatmap(newdata, col = colramp, Colv = NA, 
        Rowv = NA)






