RCTD_nuclei <- read.csv("../../../../data/Mouse_Embryo/Mouse_Embryo_RCTD_nuclei_tiny.csv",row.names = "X")
RCTD_cells <- read.csv("../../../../data/Mouse_Embryo/Mouse_Embryo_RCTD_cells_tiny.csv",row.names = "X")

df_class_RCTD <- classify_type(RCTD_results_nuclei = RCTD_nuclei, RCTD_results_cells = RCTD_cells) 
table(df_class_RCTD$type)    # vedo quante ok e quante rejected da filtri RCTD
table(df_class_RCTD$class)   # number cell types detected

write.csv(df_class_RCTD, "../../../../data/Mouse_Embryo/Mouse_Embryo_RCTD_filters.csv")

str(RCTD_nuclei$results_df)

RCTD_results_nuclei = RCTD_nuclei$results_df
RCTD_results_cells = RCTD_nuclei$results_df

cg = intersect(rownames(RCTD_results_nuclei),
               rownames(RCTD_results_cells))

  cg = intersect(rownames(RCTD_results_nuclei),
                 rownames(RCTD_results_cells))
  df_nuclei = RCTD_results_nuclei[cg,]
  df_cells = RCTD_results_cells[cg,]
  
  df_cells$first_type <- as.character(df_cells$first_type)
  
  feature_type = ifelse(df_nuclei$first_type == df_cells$first_type, "cells", "nuclei")
  feature_type[df_nuclei$spot_class != "singlet" & df_cells$spot_class != "singlet"] = "rejected"
  feature_type[df_nuclei$spot_class == "singlet" & df_cells$spot_class != "singlet"] = "nuclei"
  feature_type[df_nuclei$spot_class != "singlet" & df_cells$spot_class == "singlet"] = "cells"
  feature_class = ifelse(feature_type=="cells", df_cells$first_type, df_nuclei$first_type)
  feature_class[feature_type=="rejected"] = "rejected"
  data.frame(type=feature_type, class=feature_class, row.names = cg)