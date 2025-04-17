# ======================================================================
# == load libraries ==
# ======================================================================

from stardist import random_label_cmap, _draw_polygons
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse
import bin2cell as b2c
from stardist.models import StarDist2D
from csbdeep.utils import normalize
import cv2
import os
import pickle
import anndata
import geopandas as gpd
from tifffile import imread, imwrite
from stardist.models import StarDist2D
from shapely.geometry import Polygon, Point
from shapely.affinity import scale
from scipy.spatial.distance import pdist
import scrublet as scr
import itertools
import anndata as ad
from pandas.api.types import CategoricalDtype
import random
import sys

# ======================================================================
# == custom function definition ==
# ======================================================================

from custom_functions import remove_destripe_artifacts, nuclei_detection, add_obs_variables, cell_geometry, max_diameter

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# Redirect stderr/stdout to Snakemake log file
log_file = snakemake.log[0]
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

print(f"--- Running bin2cell for sample {snakemake.wildcards.sample_id} ---")

# 1. Get input from Snakemake object
input_adata = snakemake.input.adata
print(f"Input adata: {input_adata}")

input_parquet = snakemake.input.parquet
print(f"Input parquet: {input_parquet}")

input_stardist = snakemake.input.stardist
print(f"Input stardist: {input_stardist}")

species = snakemake.params.species
print(f"species: {species}")

# stardist parameters. Consider whether to add them
# snakemake.params.prob_thresh
# snakemake.params.nms_thresh

# 2. Get output from Snakemake object
output_nuclei_grouped_pkl = snakemake.output.nuclei_grouped_pkl
print(f"Output nuclei_grouped_pkl path: {output_nuclei_grouped_pkl}")

output_nuclei_expanded_pkl = snakemake.output.nuclei_expanded_pkl
print(f"Output nuclei_expanded_pkl path: {output_nuclei_expanded_pkl}")

output_nuclei_grouped_adata = snakemake.output.nuclei_grouped_adata
print(f"Output nuclei_grouped_adata path: {output_nuclei_grouped_adata}")

output_nuclei_expanded_adata = snakemake.output.nuclei_expanded_adata
print(f"Output nuclei_expanded_pkl path: {output_nuclei_expanded_adata}")

# ======================================================================
# == load input ==
# ======================================================================

# load parquet
with open(input_stardist, 'rb') as f:  # open a text file
    polys = pickle.load(f)

# adata
adata = sc.read_h5ad(input_adata)

# parquet
tissue_position_file = input_parquet

# ======================================================================
# == processing ==
# ======================================================================

# --- preprocessing adata ---

sc.pp.filter_genes(adata, min_cells=3)  # ok filtro sui geni
sc.pp.filter_cells(adata, min_counts=0) # non filtro i pixel perchè andranno aggregati
b2c.destripe(adata,adjust_counts=True)
adata = remove_destripe_artifacts(adata)

# --- add segmentation labels to h5ad object, map stardist modeled nuclei on the adata ---

adata, gdf = nuclei_detection(
    polys = polys, adata=adata,
    tissue_position_file = tissue_position_file
    )
adata.obs.id = adata.obs.id.astype(int)

# --- group nuclei based on labels "id" ---

# species = "Mm"
nuclei_grouped = b2c.bin_to_cell(adata, labels_key="id", spatial_keys=["spatial"])

# annotate
nuclei_grouped = add_obs_variables(nuclei_grouped, "nucleus", species=species)
nuclei_grouped.obs["id"] = np.array(nuclei_grouped.obs.index)
nuclei_grouped.obs['geometry'] = gdf.loc[nuclei_grouped.obs.id]['geometry']
nuclei_grouped.obs['max_diameter'] = [max_diameter(x) for x in nuclei_grouped.obs['geometry']]
nuclei_grouped.obs['zero_mt'] = nuclei_grouped.obs.pct_counts_mt == 0
print('Using ' + str(np.round(nuclei_grouped.X.sum() / adata.X.sum() * 100)) + '% of the total read counts')

# --- remove empty nuclei and point geometries ---

idx = nuclei_grouped.obs.features_per_nucleus >= 2
print('Number of nuclei with at least 2 features: ' + str(np.sum(idx)))
nuclei_grouped = nuclei_grouped[idx]
idx = [geom.geom_type  == 'Polygon' for geom in nuclei_grouped.obs['geometry']]
print('Number of nuclei with polygon geometry: ' + str(np.sum(idx)))
nuclei_grouped = nuclei_grouped[idx]

# --- expand nuclei into cells ---

# reset to "0" (i.e. not assigned) pixels assigned to non-valid nulclei
adata.obs.loc[~adata.obs.id.isin(nuclei_grouped.obs.index.astype(int)),'id'] = 0

# expand nuclei
#adata.obs.id = adata.obs.id.astype(int)   # to expand needs to be numeric
expand = 2      # 2 bins = 4 microns, forse conviene metterlo come parametro in input per fare diversi tentativi di espansione
expansion_label = "id_exp_"+str(expand)
b2c.expand_labels(adata,
                  labels_key='id',
                  expanded_labels_key=expansion_label,
                  max_bin_distance = expand
                 )
adata.obs[expansion_label] = adata.obs[expansion_label].astype(int)

# --- group cells based on labels ---

expanded_nuclei = b2c.bin_to_cell(adata, labels_key=expansion_label, spatial_keys=["spatial"])

# --- define new geometries and filter out cells (and corresponding nuclei) with non-Polygon geometry ---

cell_labels = set(adata.obs[expansion_label])
cell_labels.remove(0)
cell_geoms = {x : cell_geometry(adata[adata.obs[expansion_label] == x].obs['geometry']) for x in cell_labels}
idx = [geom.geom_type  == 'Polygon' for geom in cell_geoms.values()]
print('Number of cells with polygon geometry: ' + str(np.sum(idx)))
nuclei_grouped = nuclei_grouped[idx]
expanded_nuclei = expanded_nuclei[idx]
cell_geoms = {int(x) : cell_geoms[int(x)] for x in expanded_nuclei.obs_names}

# --- annotate cells ---

expanded_nuclei = add_obs_variables(expanded_nuclei, "cell", species=species)
expanded_nuclei.obs["id"] = np.array(expanded_nuclei.obs.index)
expanded_nuclei.obs['geometry'] = [cell_geoms[int(x)] for x in expanded_nuclei.obs.index]
expanded_nuclei.obs['max_diameter'] = [max_diameter(x) for x in expanded_nuclei.obs['geometry']]
expanded_nuclei.obs['zero_mt'] = expanded_nuclei.obs.pct_counts_mt == 0
print('Using ' + str(np.round(expanded_nuclei.X.sum() / adata.X.sum() * 100)) + '% of the total read counts')

# --- save geometry on disk ---

nuclei_grouped_geometry = nuclei_grouped.obs['geometry']
expanded_nuclei_geometry = expanded_nuclei.obs['geometry']
with open(output_nuclei_grouped_pkl, "wb") as f:  # open a text file,
    pickle.dump(nuclei_grouped_geometry, f),
with open(output_nuclei_expanded_pkl, "wb") as f:  # open a text file,
    pickle.dump(expanded_nuclei_geometry, f)

# --- remove geometries and save h5ad on disk of nuclei and cells ---

del nuclei_grouped.obs['geometry']
del expanded_nuclei.obs['geometry']
nuclei_grouped.write_h5ad(output_nuclei_grouped_adata)
expanded_nuclei.write_h5ad(output_nuclei_expanded_adata)
