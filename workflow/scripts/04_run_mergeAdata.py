# ======================================================================
# == load libraries ==
# ======================================================================

import pandas as pd
import scanpy as sc
import numpy as np
from scipy.sparse import vstack
import anndata as ad


# ======================================================================
# == Snakemake integation ==
# ======================================================================

# Redirect stderr/stdout to Snakemake log file
log_file = snakemake.log[0]
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

print(f"--- Running adata merge for sample {snakemake.wildcards.sample_id} ---")

# 1. Get input from Snakemake object
df_class_RCTD_file = snakemake.input.csv_filters
print(f"Input RCTD filter file: {df_class_RCTD_file}")

nuclei_grouped_file = snakemake.input.nuclei_grouped_adata
print(f"Input adata nuclei: {nuclei_grouped_file}")

nuclei_expanded_file = snakemake.input.nuclei_expanded_adata
print(f"Input adata cells: {nuclei_expanded_file}")

# 2. Get output from Snakemake object
merged_adata_file = snakemake.output.merged_adata
print(f"Output merged adata path: {merged_adata_file}")

# ======================================================================
# == load input ==
# ======================================================================

# INPUT
df_class_RCTD = pd.read_csv(df_class_RCTD_file,index_col=0)
print("Loaded RCTD filter file")

# NUCLEI GROUPED
nuclei_grouped = sc.read_h5ad(nuclei_grouped_file)
print("Loaded adata nuclei")

# EXPANDED NUCLEI
nuclei_expanded = sc.read_h5ad(nuclei_expanded_file)
print("Loaded adata cells")


# ======================================================================
# == processing ==
# ======================================================================

# --- filter the nuclei ---

nuclei_grouped.obs_names = nuclei_grouped.obs_names.astype(str)
nuclei_selected = (np.array(df_class_RCTD.index[df_class_RCTD['type'] == 'nuclei'])).astype(str)
nuclei_grouped_filtered = nuclei_grouped[nuclei_selected, :]
nuclei_grouped_filtered = nuclei_grouped_filtered.copy()
nuclei_grouped_filtered.obs["cell_nucleus"] = "nucleus"

# --- filter the cells ---

nuclei_expanded.obs_names = nuclei_expanded.obs_names.astype(str)
cells_selected = (np.array(df_class_RCTD.index[df_class_RCTD['type'] == 'cells'])).astype(str)
nuclei_expanded_filtered = nuclei_expanded[cells_selected, :]
nuclei_expanded_filtered = nuclei_expanded_filtered.copy()
nuclei_expanded_filtered.obs["cell_nucleus"] = "cell"

# --- merge the adata ---

if np.all(nuclei_grouped_filtered.var_names == nuclei_expanded_filtered.var_names):
    combined_matrix = vstack([nuclei_grouped_filtered.X, nuclei_expanded_filtered.X])
    combined_obs = pd.concat([nuclei_grouped_filtered.obs, nuclei_expanded_filtered.obs], axis=0)
    
    adata_final = ad.AnnData(X=combined_matrix, obs=combined_obs, var = nuclei_expanded_filtered.var)
    adata_final.uns["spatial"] = nuclei_grouped_filtered.uns["spatial"]
    spatial_1 = np.array(nuclei_grouped_filtered.obsm["spatial"])
    spatial_2 = np.array(nuclei_expanded_filtered.obsm["spatial"])
    spatial_concat = np.concatenate([spatial_1, spatial_2], axis=0)
    adata_final.obsm["spatial"] = spatial_concat
    adata_final.raw = adata_final

df_class_RCTD.index = df_class_RCTD.index.astype(str)
adata_final.obs["cell_types"] = df_class_RCTD.loc[adata_final.obs.index, "class"]

# --- save the adata ---

adata_final.write_h5ad(merged_adata_file)