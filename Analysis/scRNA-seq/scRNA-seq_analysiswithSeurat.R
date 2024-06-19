# analysis on scRNa of publicly available data 
# (the data we took is a sample from the following link 
# https://www.10xgenomics.com/datasets/20k-mixture-of-nsclc-dtcs-from-7-donors-3-v3-1-with-intronic-reads-3-1-standard),
# we uploaded cell matrix HDF5 (raw) only for this project
# utilizing Seurat package and others 
# you can install all the packages using following command install.packages('Seurat')

# setwd("~/Documents/bioinformatics/Analysis/scRNA-seq")

# Note if during the scaling process your allocated memory for R process will exost itsefl utilize the following command to 
# library(usethis) 
# usethis::edit_r_environ()
# those will prompt the new tab to open and you can type there R_MAX_VSIZE=100Gb (or how much you want to allocate)
# this will also promt to restarting the R and all variables will be lost. 
# Thus, I am recomending to see how much is allocated and add small amount prior the project 

library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(ggplot2)
library(clustree)

# step 1 reading data 
# creating matrix
matrix_data <- Seurat::Read10X_h5(filename = '20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_count_raw_feature_bc_matrix.h5')

# pulling out only the gene expression parameter
gexp <- matrix_data$`Gene Expression`

# first lets create SeuratObj 

non_normalized_data <- CreateSeuratObject(counts = gexp, 
                                          project = "lung_cancer", 
                                          min.cell = 5 , # we want at lest expression in 5 cells 
                                          min.features = 200) # we want at least 200 genes to be present 
str(non_normalized_data)

# step 2 QC 
# lets look into mitochondrial genes -> if we have those in cell, 
# means the cell is dying, due to cell membrane disruption

# the step bellow will add the pecent-mt_gene column based on which we going to filter data 
non_normalized_data[['percent_mt_genes']] <- PercentageFeatureSet(non_normalized_data, pattern = '^MT-')

# adding the ribosomal genes 
non_normalized_data[['percent_rib_genes']] <- PercentageFeatureSet(non_normalized_data, pattern = '^RP')

# adding the hemoglobin genes 
non_normalized_data[['percent_hgb_genes']] <- PercentageFeatureSet(non_normalized_data, pattern = '^HB')

dim(non_normalized_data@meta.data) # we have 71880 rows 

str(non_normalized_data)

# plot 
VlnPlot(non_normalized_data, features = c("percent_mt_genes", 
                                          "percent_rib_genes", 
                                          "percent_hgb_genes"))

# checking correlation b/w mt_genes and hgb_genes. = -0.01
FeatureScatter(non_normalized_data, feature1 = 'percent_hgb_genes', feature2 = 'percent_mt_genes')

FeatureScatter(non_normalized_data, feature1 = 'nCount_RNA', 
               feature2 = 'nFeature_RNA') + geom_smooth(method = 'lm')
# now let filter out the rows which are < 5% mitochondrial genes, 
# and gene # should be above 200 and less than 2500

non_normalized_data <- subset(non_normalized_data, subset = nFeature_RNA > 200 & 
                                nFeature_RNA < 2500 & percent_mt_genes < 5 & percent_rib_genes < 10)

dim(non_normalized_data@meta.data) # now we have 14681 rows 


# step 3 Normalization 
# P.S I would like to create the new alias for this purpose 

normalized_data <- NormalizeData(non_normalized_data)

# step 4 feature selection 
# we are going to use the default features, method and nfeatures
normalized_data <- FindVariableFeatures(normalized_data)
dim(normalized_data)

most_variable <- head(VariableFeatures(normalized_data), 15)
str(most_variable)

most_var_gene_plot <- VariableFeaturePlot(normalized_data)
LabelPoints(plot = most_var_gene_plot, points = most_variable, repel = TRUE)

# step 5 scaling 
genes <- rownames(normalized_data)
normalized_data <- ScaleData(normalized_data, features = genes)  

# dimensionality reduction 
# first we will do the linear PCA dimentionality reduction and after the non-linear UMAP

normalized_data <- RunPCA(normalized_data) # takes time 
str(normalized_data)
print(normalized_data@reductions$pca, dims = 1:15, nfeatures = 10)

# elbow plot to define what PCs we need 
ElbowPlot(normalized_data) # looks like 15 Pcs is best value from the albow plot 

DimHeatmap(normalized_data, dims = 1:15, cells = 500, balanced = TRUE)

DimPlot(normalized_data, reduction = "pca")

# running the UMAP ( non-linear dimentionality reduction algorithm)
normalized_data <- RunUMAP(normalized_data, dims = 1:15)

# plotting the above
DimPlot(normalized_data, reduction = "umap")

#step 6 Clustering 
normalized_data <- FindNeighbors(normalized_data, dims = 1:15)

normalized_data <- FindClusters(normalized_data, resolution = seq(0.1, 0.8, by=0.1))

head(normalized_data@meta.data)

clustree(normalized_data,  return= c('plot', 'graph','layout'))
# looks like more reasonable resolution is 0.3
DimPlot(normalized_data, group.by = "RNA_snn_res.0.3")

# setting up the identity
# first checking the original identity 
Idents(normalized_data)

#setting up based on choised resolution 0.3 -> where we define 8 communities(clusters)
Idents(normalized_data) <- "RNA_snn_res.0.3"

# Summary: We identify at least 2 big clusters in our dataset 


# step 10 cell type annotation 

# saving counts into the reference 
counts <- GetAssayData(normalized_data, assay="RNA", layer="counts")

# reading the reference data 
library(celldex)
reference <- HumanPrimaryCellAtlasData()
View(as.data.frame(colData(reference)))

library(SingleR)
attempt <- SingleR(test = counts, 
        ref = reference, 
        labels = reference$label.main)

# saving labels into Seurat object 

normalized_data[['singleR_labels']] <- attempt$labels[match(rownames(normalized_data@meta.data), rownames(attempt))]
DimPlot(normalized_data, reduction = 'umap', group.by = 'singleR_labels')
# after viz of the clusters with annotated values we can see that we have cluster of T cells, B cell and biggest cluster with different gene expression
# I am not really happy with the way it assigned the third biggest cluster. 
# lets see how well it's annotated 

attempt$scores 
library(pheatmap)
plotScoreHeatmap(attempt)

# another check with the unsupervised clustering
library(tidyverse)
data <- table(assigned = attempt$labels, clusters = normalized_data$seurat_clusters)
View(data)
# review of heatmap gives us review of the clusters and assigned cells 
pheatmap(log10(data+10), color = colorRampPalette(c('white', 'green'))(10))






