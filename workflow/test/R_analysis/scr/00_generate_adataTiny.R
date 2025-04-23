# libraries ---------------------------------------------------------------
library(Seurat)
library(spacexr)
library(rhdf5)
library(tidyverse)
library(schard)

# -------------------------------------------------------------------------
# read in the sample reference
seu_obj_nuclei <- schard::h5ad2seurat_spatial("../../../../data/Mouse_Embryo/Mouse_Embryo_nuclei_grouped.h5ad",use.raw = F)

# reduce test file
# reduce the size of the reference just of testing
set.seed(2144)
keep_nuclei <- seu_obj_nuclei@meta.data %>%
  sample_n(10000,replace = F) %>%
  rownames()

# subset the object
seu_obj_nuclei2 <- subset(seu_obj_nuclei, cells = keep_nuclei)
SpatialDimPlot(seu_obj_nuclei2)

saveRDS(seu_obj_nuclei2,"../../../../data/Mouse_Embryo/Mouse_Embryo_nuclei_grouped_tiny.rds")

# read in the sample reference
seu_obj_cell <- schard::h5ad2seurat_spatial("../../../../data/Mouse_Embryo/Mouse_Embryo_expanded_nuclei.h5ad",use.raw = F)

# reduce test file
# reduce the size of the reference just of testing
# use the same cells
keep_cell <- keep_nuclei

# subset the object
seu_obj_cell2 <- subset(seu_obj_cell, cells = keep_cell)
SpatialDimPlot(seu_obj_cell2)

saveRDS(seu_obj_cell2,"../../../../data/Mouse_Embryo/Mouse_Embryo_expanded_nuclei_tiny.rds")
