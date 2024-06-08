# INSTALATION OF THE LIBRARIES 
library(preprocessCore)

# review the file, first 66 line is metdata, the table starts from line 67, 
# it is tab separated file 
microarray <- read.csv('GSE7765-GPL96-series-matrix 2.txt', skip = 66, 
                       header = TRUE, sep = "\t", row.names = 1)

# shape of the table 
dim(microarray)

#draw boxplot 
boxplot(microarray, las = 2, cex.axis = 0.7 )

#analysis of box plot display pretty wide range of the 
#values from close to zero to all the way to the 100000 
# see the max an min in whole table, 
#appending the [n], where n = 1, 6 -> choces the columns

max(microarray)
min(microarray)

# lets do log transfornation 
microarrayLog <- apply(microarray, 2, log2) 
# box plot looks so much better 
boxplot(microarrayLog, las = 2, cex.axis = 0.7 )

# normalization 
microarrayLogNorm <- normalize.quantiles(microarrayLog)

# renaming the columns 
colnames(microarrayLogNorm) <- colnames(microarrayLog)
rownames(microarrayLogNorm) <-rownames(microarrayLog)

#box plot on the log and norm data -> revealed beatiful data with the 
# [3,17] 
boxplot(microarrayLogNorm, las = 2, cex.axis = 0.7 )

#clusteirng to find out the groups

clusters <- hclust(dist(t(microarrayLogNorm)))
plot(clusters)

# grouping of the samples 
scaledVar <- apply(microarrayLogNorm,1,function(x){
  var(x)/mean(x)
})

# highly variable genes 
highVar <- names(which(scaledVar >1))

#viz with the heatmap 
heatmap(microarrayLogNorm[highVar,], cexRow = 0.7, cexCol = 0.7)





