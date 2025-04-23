library(schard)
library(Seurat)

tabula_muris <- schard::h5ad2seurat('/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/test/4df3d323-25b2-4749-b475-81a1be98e12b.h5ad',use.raw = T)

tabula_muris@meta.data

saveRDS(tabula_muris, file = '/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/test/tabula_muris.rds')

# https://datasets.cellxgene.cziscience.com/ae0c62a1-a30f-4033-97d2-0edb2e146c53.rds
# https://datasets.cellxgene.cziscience.com/98770743-77c9-4747-948e-940f77e3ed0b.rds

test_01 <- readRDS("/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/test/ae0c62a1-a30f-4033-97d2-0edb2e146c53.rds")
test_01@meta.data

test_02 <- readRDS("/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/test/98770743-77c9-4747-948e-940f77e3ed0b.rds")
test_02@meta.data
