library(ALL)
library(hgu95av2.db)
data(ALL)
ALL
experimentData(ALL)

# exploration analysis 

exprs(ALL)[1:4, 1:4]

ALL$age

featureData(ALL)

ids = featureNames(ALL)
ids[1:5]

as.list(hgu95av2ENTREZID[ids])

library(airway)
data(airway)
airway
colData(airway)
airway$SampleName

head(rownames(airway))

airway
assay(airway, 'counts')[1:4, 1:4]

rowRanges(airway)

sum(elementLengths(rowRanges(airway)))

library(GEOquery)
elist = getGEO('GSE11675')
edata = elist[[1]]
pData(edata)
eList2 = getGEOSuppFiles('GSE11675')
eList2

library(biomaRt)
head(listMarts())
mart = useMart('ensembl')
ensembl = useDataset('hsapiens_gene_ensembl', mart)

ensembl







