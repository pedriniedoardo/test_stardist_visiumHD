# libraries ---------------------------------------------------------------
library(Seurat)
library(spacexr)
library(rhdf5)
library(tidyverse)

# -------------------------------------------------------------------------
# read in the sample reference
reference_obj <- readRDS("../../../../data/test/tabula_muris.rds")


# full dataset
# make sure there are more than 25 cells in the reference
id_cell_keep <- reference_obj$cell_type %>% table() %>% .[.>=25] %>% names()
reference_obj_filter <- subset(reference_obj, subset = cell_type %in% id_cell_keep)

# make sure the name of the features is mathching the two objects
new_features <- reference_obj_filter@assays$RNA@meta.features$feature_name %>% str_remove_all("_ENSM.*") %>% make.unique()

ref_raw_mat <- reference_obj_filter@assays$RNA@counts
rownames(ref_raw_mat) <- new_features

cell_types <- as.factor(reference_obj_filter$cell_type) # in che variabile se trovano i celltypes? è standard?
names(cell_types) <- colnames(reference_obj_filter)
reference <- Reference(counts = ref_raw_mat, cell_types = cell_types)

saveRDS(reference, file = "/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/test/reference_tabula_muris_full.rds")

# reduce test file
# reduce the size of the reference just of testing
set.seed(2144)
cell_keep <- reference_obj@meta.data %>%
  sample_n(10000,replace = F) %>%
  rownames()

# subset the object
reference_obj2 <- subset(reference_obj, cells = cell_keep)

# make sure there are more than 25 cells in the reference
id_cell_keep2 <- reference_obj2$cell_type %>% table() %>% .[.>=25] %>% names()
reference_obj_filter2 <- subset(reference_obj2, subset = cell_type %in% id_cell_keep2)

# make sure the name of the features is mathching the two objects
new_features2 <- reference_obj_filter2@assays$RNA@meta.features$feature_name %>% str_remove_all("_ENSM.*") %>% make.unique()

ref_raw_mat2 <- reference_obj_filter2@assays$RNA@counts
rownames(ref_raw_mat2) <- new_features2

cell_types2 <- as.factor(reference_obj_filter2$cell_type) # in che variabile se trovano i celltypes? è standard?
names(cell_types2) <- colnames(reference_obj_filter2)
reference2 <- Reference(counts = ref_raw_mat2, cell_types = cell_types2)

saveRDS(reference2, file = "/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/test/reference_tabula_muris_tiny.rds")

