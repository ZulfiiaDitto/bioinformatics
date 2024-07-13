library(Biobase)
library(GenomicRanges)
library(SummarizedExperiment)
# quesiton 2
data(sample.ExpressionSet, package = "Biobase")
head(data)

se = makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)
genome = assay(se)
phenotype = colData(se)
feature = rowData(se)
var = rowRanges(se)

dim(genome)
dim(phenotype)
dim(feature)
var$`AFFX-BioC-5_st`
# question 5
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)

# question6 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)

library(plotrix)
pie3D(pdata_bm$num.tech.reps,labels=pdata_bm$tissue.type)

# question7 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
# this one 
row_sums = rowSums(edata)
edata = edata[order(-row_sums),]
head(edata)
index = 1:500
heatmap(edata[index,],Rowv=NA,Colv=NA)

#question 8 
library(DESeq2)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata = pData(bm)
edata = exprs(bm)
# ma plot 
ma = log2(edata[,1]+1) - log2(edata[,2]+1) # y axis
aa = log2(edata[,1]+1) + log2(edata[,2]+1) # x axis
plot(aa, ma, col = 2) 

# ma using the Deseq2 package 

rlogData = rlog(edata)
y = rlogData[,1] - rlogData[,2]
x = rlogData[,1] + rlogData[,2]
plot(x, y, col = 3)

#question 9 
library('rafalib')
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# original 
# lets calculate the euclidian distance between rows 
distO = dist(t(edata))
head(distO)
heatmap(as.matrix(distO), col = colramp, Colv = NA, Rowv = NA)
hclustO = hclust(distO)
plot(hclustO, hang = -1, labels = FALSE, title = 'Original') # total mess

# After filtering all genes with rowMeans
edataF = edata[rowMeans(edata) < 100, ]
dim(edataF)
distF = dist(t(edataF))
heatmap(as.matrix(distF), col = colramp, Colv = NA, Rowv = NA)
hclustF = hclust(distF)
plot(hclustF, hang = -1, labels = FALSE)

# log transformation 
log_edata = log2(edata + 1)
l_dist1 = dist(t(log_edata))
l_hclust1 = hclust(l_dist1)

par(mar=c(0, 4, 4, 2))
plot(l_hclust1, hang=-1, main="perform log2 transform", labels=FALSE)


# question 10 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

edata = log2(edata + 1)

# perfrom k-means clustering
set.seed(1235)
k2 = kmeans(edata,centers=2)
matplot(t(k2$centers),col=1:2,type="l",lwd=3)

dist1 = dist(t(edata))
hclust1 = hclust(dist1)
tree = cutree(hclust1, 2)

par(mar=c(0, 4, 4, 2))
plot(hclust1, tree, main="cutree")