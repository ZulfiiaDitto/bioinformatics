# setwd("~/Documents/bioinformatics/Analysis/Macroarray analysis")
# analysis used publically available data from EMBL-EBI
# R analysis happend on Jun 11 2024

# getting data from the the database
# data pubmed ID 16860862  - pre-eclampsia data 
# Includes 47,000 genes screening, two chanel microarray 
# expression of the genes from 
# placenta of 13 normal and 14 pathological pregnancy, 


library(preprocessCore)

# 1. read the table with ref data 
# downloaded files have isue and need to be pre-processed: 
# run the following bash command, whitch replaced "#" to the tab delimeter and output cleanRefFile.txt 
#  sed 's/#/\t/g' "E-GEOD-4707.sdrf.txt" > 'CleanRefFile.txt' and then 
#  sed 's/ftp:/https:/g' CleanRefFile.txt > Reference.txt


df <- read.table('Reference.txt', sep="\t", header = TRUE)

url <- df$Comment..Derived.ArrayExpress.FTP.file.[1]

# 2. get data from the ftp.file

temp <- tempfile()
download.file(url,temp)
con <- unz(temp, "E-GEOD-4707-processed-data-1618157355.txt")

data <- read.table(con, header = TRUE,  sep="\t")
unlink(temp)

head(data, 7) # some cleaning of the columns names are needed

# the first row in data and column names need to be combined
first_row = as.character(data[1, ])
new_columns = paste(names(data),first_row, sep ="_")
print(new_columns)

final <- data[-1,]
colnames(final) <- new_columns

# 3. selecting only columns with mean green and red fluorescent light 

final <- final[, grep("Mean|Scan", names(final), value = TRUE)]


# Step 3. 1: Identify columns matching the pattern
pattern <- "GEO:AGILENT_gMeanSignal"
colnames(final) <- gsub(pattern, "greenMeanSignal", colnames(final))

pattern2 <- "GEO:AGILENT_rMeanSignal"
colnames(final) <- gsub(pattern2, "redMeanSignal", colnames(final))

# after reviwing the table -> certain rows need to be drop since they are not the probes

final <- final[1:41675,]

# 4 normalization 

oligo <- final$`Scan.REF_Reporter REF`

#4a converting everything to numeric 
for (i in 2:ncol(final)) {
  final[, i] <- as.numeric(final[, i])
}
#4b converting to matcix
matrix <- as.matrix(final[,-1])
#4c finally normalization 
normalized_data <- normalize.quantiles(matrix)

# Combine normalized data with gene identifiers

normalized_final <- cbind(Gene = oligo, as.data.frame(normalized_data))
colnames(normalized_final)[-1] <- colnames(final)[-1]

# Check the structure of the normalized data
str(normalized_final)

boxplot(as.data.frame(normalized_data), main = "Boxplot of Normalized Data", las = 2)

# Log2 transformation
log2_normalized_data <- log2(normalized_data + 1)

# box polot looks little bit better 
boxplot(as.data.frame(log2_normalized_data), main = "Boxplot of Normalized Data", las = 2)

# T test 
control_samples = normalized_final[, c(6,7, 10,11, 14, 15,18,19 )]
sample = normalized_final[, -c(6,7, 10,11, 14, 15,18,19 )]

p_values <- apply(normalized_data, 1, function(row) {
  control_samples <- row[c(6,7, 10,11, 14, 15,18,19)]
  sample <- row[-c(6,7, 10,11, 14, 15,18,19 )]
  t.test(control_samples, sample)$p.value
})

results <- data.frame(Gene = oligo, p_value = p_values)

# Filter significant genes
significant_genes <- results[results$p_value < 0.05, ]

# next step is to pull what each oligo is represent, gene location 
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 

annotated_ids <- c(significant_genes$Gene) 
# TODO: figurate out the biomart API 
genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
               filters = "", values = annotated_ids,
               mart = ensembl)

# Display the retrieved gene annotations
print(listFilters(ensembl))
