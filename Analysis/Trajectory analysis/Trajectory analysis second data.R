library(tidyverse) 
library(monocle3) 
library(ggplot2)


expression_matrix <- read.delim2('ABC_umi_matrix_7551_cells.csv', header = TRUE, sep= ',')
cell_metadata <- read.delim('ABC_Meta.txt', header = TRUE)
gene_annotation <- read.delim('ABC_Marker.txt', header = TRUE)

expression_matrix <- t(expression_matrix)
dim(expression_matrix)
#NOTE according the documentation in monocle3 package the gene_annotation file 
# shoudl have column names 'gene_short_name'
gene_annotation <- rename(gene_annotation, Gene = 'gene_short_name')

# step 1 Make the Cell Data Sell object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
