# libraries ---------------------------------------------------------------
library(Seurat)
library(spacexr)
library(rhdf5)
library(tidyverse)

# functions ---------------------------------------------------------------

# this can be swapped with schard implementation
h5ad_to_seurat = function(file_in) {
  
  stopifnot(require(rhdf5))
  stopifnot(require(Matrix))
  stopifnot(require(Seurat))
  
  cell_counts = h5read(file = file_in, 'X')
  cell_obs = h5read(file = file_in, 'obs')
  cell_var = h5read(file = file_in, 'var')
  
  n_rows <- length(cell_counts$indptr) - 1  # Number of rows in the matrix
  n_cols <- max(cell_counts$indices) + 1    # Number of columns (assuming 0-based indices)
  
  i = rep(1:n_rows, times = diff(cell_counts$indptr))
  j = cell_counts$indices + 1
  x = as.numeric(cell_counts$data)
  
  dcg_mat = sparseMatrix(
    i = i, j = j, x = as.numeric(x), dims = c(n_rows, n_cols),
    dimnames = list(
      as.character(cell_obs$`_index`),
      as.character(cell_var$`_index`)
    ))
  
  CreateSeuratObject(CreateAssayObject(counts = t(dcg_mat)), meta.data = cell_obs)
  
}

classify_type = function(RCTD_results_nuclei, RCTD_results_cells) {
  cg = intersect(rownames(RCTD_results_nuclei),
                 rownames(RCTD_results_cells))
  df_nuclei = RCTD_results_nuclei[cg,]
  df_cells = RCTD_results_cells[cg,]
  feature_type = ifelse(df_nuclei$first_type == df_cells$first_type, "cells", "nuclei")
  feature_type[df_nuclei$spot_class != "singlet" & df_cells$spot_class != "singlet"] = "rejected"
  feature_type[df_nuclei$spot_class == "singlet" & df_cells$spot_class != "singlet"] = "nuclei"
  feature_type[df_nuclei$spot_class != "singlet" & df_cells$spot_class == "singlet"] = "cells"
  feature_class = ifelse(feature_type=="cells", df_cells$first_type, df_nuclei$first_type)
  feature_class[feature_type=="rejected"] = "rejected"
  data.frame(type=feature_type, class=feature_class, row.names = cg)
}

# processing --------------------------------------------------------------
# INPUT

# REFERENCE 
# for the generation of the reference see 00_generate_reference.R
reference <- readRDS("../../../../data/test/reference_tabula_muris_tiny.rds")

# fai stessa rule con i due input diversi?
# ADATA: NUCLEI
# seu_obj <- h5ad_to_seurat("../../../../data/Mouse_Embryo/Mouse_Embryo_nuclei_grouped.h5ad")
seu_obj <- readRDS("../../../../data/Mouse_Embryo/Mouse_Embryo_nuclei_grouped_tiny.rds")
meta_obs <- seu_obj@meta.data[,c("array_row", "array_col")]

# rouding is needed as the stardist produce non integer assignament
# puck <- spacexr::SpatialRNA(meta_obs, round(seu_obj@assays$RNA$counts, 0))
puck <- spacexr::SpatialRNA(meta_obs, round(seu_obj@assays$Spatial$counts, 0))
# alternatively if not filtering, it is possible to change the default value of the CELL_MIN_INSTANCE parameter
myRCTD <- spacexr::create.RCTD(puck,reference,max_cores = 15,test_mode = F, UMI_min = 20,CELL_MIN_INSTANCE = 25)
myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = "doublet")
RCTD_nuclei <- myRCTD@results

write.csv(RCTD_nuclei$results_df,"../../../../data/Mouse_Embryo/Mouse_Embryo_RCTD_nuclei_tiny.csv")


# ADATA: CELLS
# seu_obj2 <- h5ad_to_seurat("../../../../data/Mouse_Embryo/Mouse_Embryo_expanded_nuclei.h5ad")
seu_obj2 <- readRDS("../../../../data/Mouse_Embryo/Mouse_Embryo_expanded_nuclei_tiny.rds")
meta_obs2 <- seu_obj2@meta.data[,c("array_row", "array_col")]

# rouding is needed as the stardist produce non integer assignament
# puck2 <- spacexr::SpatialRNA(meta_obs2, round(seu_obj2@assays$RNA$counts, 0))
puck2 <- spacexr::SpatialRNA(meta_obs2, round(seu_obj2@assays$Spatial$counts, 0))
# alternatively if not filtering, it is possible to change the default value of the CELL_MIN_INSTANCE parameter
myRCTD2 <- spacexr::create.RCTD(puck2,reference,max_cores = 15,test_mode = F, UMI_min = 20,CELL_MIN_INSTANCE = 25)
myRCTD2 <- spacexr::run.RCTD(myRCTD2, doublet_mode = "doublet")
RCTD_cells <- myRCTD2@results

write.csv(RCTD_cells$results_df,"../../../../data/Mouse_Embryo/Mouse_Embryo_RCTD_cells_tiny.csv")

# PARTE 2: FILTRI ON RCTD RESULTS

# INPUT: 
# RCTD_cells = read.csv(paste0(path, "RCTD_cells.csv"), row.names = "X")
# RCTD_nuclei = read.csv(paste0(path, "RCTD_nuclei.csv"), row.names = "X")

df_class_RCTD = classify_type(RCTD_results_nuclei = RCTD_nuclei$results_df, RCTD_results_cells = RCTD_cells$results_df) 
table(df_class_RCTD$type)    # vedo quante ok e quante rejected da filtri RCTD
table(df_class_RCTD$class)   # number cell types detected

write.csv(df_class_RCTD, "../../../../data/Mouse_Embryo/Mouse_Embryo_RCTD_filters.csv")
