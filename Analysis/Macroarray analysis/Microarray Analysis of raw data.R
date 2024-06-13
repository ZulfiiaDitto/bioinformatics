# setwd("~/Documents/bioinformatics/Analysis/Macroarray analysis")

# installing library 
library(affy)
library(GEOquery)
library(tidyverse)

# we are going to pull data utilizing the Geoquery 
# 1 pull supplemental file

getGEOSuppFiles("GSE69078") # upload the file in your working directory 

# 2 unzip the file 
untar("GSE69078/GSE69078_RAW.tar", exdir = 'data/') # creates folder with your files in it

# 3 read cel files 
raw_data <- ReadAffy(celfile.path = "data/") # pulls all CEL reads into affyBatch file 

# 4 lets normalize files 
# rma function AffyBatch and utilize robust multi-array average (rma)
# expression normalize it, by default it also correct for background 

normalized <- rma(raw_data)
# lets return the dataframe of the matrices 
normalized_exp <- as.data.frame(exprs(normalized))

# 5 extracting the gene names from supporting file

support_table <- getGEO("GSE69078", GSEMatrix = TRUE)

gene_vector <- support_table$GSE69078_series_matrix.txt.gz@featureData@data

# 6 attaching the column

final <- merge(normalized_exp, gene_vector, by = 0, all.x = TRUE)

head(final)

# 7 cleaning the columns names 
support_table <- getGEO("GSE69078", GSEMatrix = TRUE)
column_vector <- support_table$GSE69078_series_matrix.txt.gz@phenoData@data

test <- final[, c(2:13, 23:25, 27)]

test2 <- as.data.frame(t(test))
columns <- as.character(test2[nrow(test2) -2, ])
colnames(test2) <- columns

test2 <- head(test2, nrow(test2) -4)

test2 <- cbind(column_vector$`cell line:ch1`, test2)

test2$`column_vector$\`cell line:ch1\``[1:6] <- "Cfz"

test2$`column_vector$\`cell line:ch1\``[7:12] <- "Kms"

rownames(test2) <- NULL

# 8 stat analysis 

# TODO: 
#1 loop thought the df 
#2 perfomr T-test that genes really expressed differently between group -> t -test 
# you will have to do some data wrangling 

# for further analysis we need to divide samples into two groups: KMS and KMS/Cfz,
# thus rename the columns to group

for( i in colnames(2:ncol(test2))) {
  print(i)
  local <- test2[c(1, i)]
  print(local)
  break
}

head(test2)
