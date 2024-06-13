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

# 8. statistical analysis

# for further analysis we need to divide samples into two groups: KMS and KMS/Cfz,
# thus rename the columns to group

# we do not have much of the samples to perform statistics (min 10 samples in a group, we have 6 in each).
#However, as a practice lets do T-test 

p_values <- c()
gene_names <- c()

# Loop through columns from 2 to ncol(test2)
for (i in 2:ncol(test2)) {
  print(paste("Processing column:", i, "(", colnames(test2)[i], ")"))
  
  # Extract the relevant column data
  local <- test2[, c(i)]
  
  # Extract sample and control values
  sample <- as.numeric(local[1:6])
  control <- as.numeric(local[7:12])
  
  # Perform t-tests
  t_test<- t.test(control, sample)
  # Extract p-values
  p_value<- t_test$p.value
    
   # Check if the p-value is significant (e.g., < 0.05)
  if (p_value < 0.05) {
      p_values <- c(p_values, p_value)
      gene_names <- c(gene_names, colnames(test2)[i]) 
    }
}

# Create a data frame with significant p-values and gene names
significant_results <- data.frame(
  Gene = gene_names,
  p_value = p_values
)

# we identify 2117 genes has different gene expression 
# Note: we did not check if those genes express less or greater. 
# Just a difference in a mean

# 9. visualize 
# lets visualize some of the gene expression with statistical difference in expression 
# lets vizualize gene expression of INSIG1 (11075 index)
#col_index <- grep("INSIG1", colnames(test2))
ins_gene <- test2[c(1, 11075)]
print(ins_gene)
boxplot(as.numeric(ins_gene$INSIG1[1:6]), 
        as.numeric(ins_gene$INSIG1[7:12]), names = c("KMS-Cfz", "KMS"), main = "Gene Expression Levels", ylab = "Expression")

# 10. further steps in analysis will be to review biological process of those genes.