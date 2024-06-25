# setwd("~/Documents/bioinformatics/Analysis/Trajectory analysis") 

# trajectory analysis -> utilizing the Monocle3 library 

library(tidyverse) 
library(monocle3) 
library(ggplot2)


expression_matrix <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("https://depts.washington.edu:/trapnell-lab/software/monocle3/celegans/data/cao_l2_rowData.rds"))

# step 1 Make the Cell Data Sell object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# step 2 preprocess the data (usually same steps as scRNA analysis)
# monocle3 utilize PCA for dimentinality reduction -> linear approach 
cds <- preprocess_cds(cds, num_dim = 80)
plot_pc_variance_explained(cds)

# reduce dimentinality utilizing default UMAP 
cds <- reduce_dimension(cds)
# plotting by cell type 
plot_cells(cds,color_cells_by="cao_cell_type")
# plotting by gene expression 
plot_cells(cds, genes=c("cpna-2", "egl-21", "ram-2", "inos-1"))
# reduction by tSNE
cds <- reduce_dimension(cds, reduction_method="tSNE")
# plotting by cells 
plot_cells(cds, reduction_method="tSNE", color_cells_by="cao_cell_type")

# step 3 removing batch effect 
# first lets check the batch effect in this dataset 
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)
# if there was only one plate -> well we suffer the batch effect 
# in this case the bacth effect is minimal 

# we can also remove the batch effect by following function, 
# knn algorithm utilized in the case of removing the batch effect 
cds <- align_cds(cds, num_dim = 100, alignment_group = "plate")
cds <- reduce_dimension(cds)
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)

# step 4 clusting -> clusting based on the PAGA algorithm
cds <- cluster_cells(cds, resolution=1e-5)
plot_cells(cds)

plot_cells(cds, color_cells_by="cao_cell_type")

# step 5 finding market gene expression 
marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)



top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(4, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))


plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=2)

# step 6 cell annotation 
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$assigned_cell_type,
                                                 "1"="Body wall muscle",
                                                 "2"="Germline",
                                                 "3"="Motor neurons",
                                                 "4"="Seam cells",
                                                 "5"="Sex myoblasts",
                                                 "6"="Socket cells",
                                                 "7"="Marginal_cell",
                                                 "8"="Coelomocyte",
                                                 "9"="Am/PH sheath cells",
                                                 "10"="Ciliated neurons",
                                                 "11"="Intestinal/rectal muscle",
                                                 "12"="Excretory gland",
                                                 "13"="Chemosensory neurons",
                                                 "14"="Interneurons",
                                                 "15"="Unclassified eurons",
                                                 "16"="Ciliated neurons",
                                                 "17"="Pharyngeal gland cells",
                                                 "18"="Unclassified neurons",
                                                 "19"="Chemosensory neurons",
                                                 "20"="Ciliated neurons",
                                                 "21"="Ciliated neurons",
                                                 "22"="Inner labial neuron",
                                                 "23"="Ciliated neurons",
                                                 "24"="Ciliated neurons",
                                                 "25"="Ciliated neurons",
                                                 "26"="Hypodermal cells",
                                                 "27"="Mesodermal cells",
                                                 "28"="Motor neurons",
                                                 "29"="Pharyngeal gland cells",
                                                 "30"="Ciliated neurons",
                                                 "31"="Excretory cells",
                                                 "32"="Amphid neuron",
                                                 "33"="Pharyngeal muscle")


plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type")

# choosing small cluster of cells for deep dive 
cds_subset <- choose_cells(cds) # nice pop-up window interface for cluster selection 

# NOte: annotation cells can alco happend autolatically based on 
# the Garnet alg -> not provided here, additinal dependency will be needed 

# step7 travectory analysis 

cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)


cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

