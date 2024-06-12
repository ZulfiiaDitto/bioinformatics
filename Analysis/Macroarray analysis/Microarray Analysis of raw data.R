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

# TODO: 
  #1 columns -> better namign convention
  #2 perfomr T-test that genes really expressed differently between group -> t -test 
   # you will have to do some data wrangling 

  
