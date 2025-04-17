# ======================================================================
# == load the libraries ==
# ======================================================================

from stardist.data import test_image_nuclei_2d
from stardist.plot import render_label
from stardist.models import StarDist2D
from shapely.geometry import Polygon, Point
import matplotlib.image
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as patches
import pandas as pd
import numpy as np
import seaborn as sns
import PIL
import tifffile
import sys
import pickle

# ======================================================================
# == definition of custom functions ==
# ======================================================================

# import the function file
from custom_functions import custom_normalize

# ======================================================================
# == Snakemake integation ==
# ======================================================================

# Redirect stderr/stdout to Snakemake log file
log_file = snakemake.log[0]
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

print(f"--- Running Stardist for sample {snakemake.wildcards.sample_id} ---")

# 1. Get input from Snakemake object
img_input_path = snakemake.input.raw_image
print(f"Input source image path: {img_input_path}")

hist_img = snakemake.params.hist_img
print(f"Processing an {hist_img} image")

# stardist parameters. Consider whether to add them
# snakemake.params.prob_thresh
# snakemake.params.nms_thresh

# 2. Get output from Snakemake object
output_pkl_path = snakemake.output.pkl
print(f"Output pkl path: {output_pkl_path}")

# ======================================================================
# == processing ==
# ======================================================================

# --- read in the input image ---

# img hires
print(f"Reading raw image from {img_input_path}...")
PIL.Image.MAX_IMAGE_PIXELS = 1718032108
img = tifffile.imread(img_input_path)
img = img[:,:,0:3]
print("Reading raw image completed")

# --- run custom normalization ---

print("Custom image normalization start")
normalized_img = custom_normalize(img)
print("Custom image normalization completed")

# --- select approprote stardist model ---

if hist_img == "H&E":
    model_stardist = "2D_versatile_he"
elif hist_img == "IF":
    model_stardist = "2D_versatile_fluo"

model = StarDist2D.from_pretrained(model_stardist)

# --- run the model ---

# here is running with default parameters for prob_thresh and nms_thresh
print("Stardist run start")
labels, polys = model.predict_instances_big(
    normalized_img, axes='YXC', block_size=4096, min_overlap=128, context=128,
    normalizer=None, # n_tiles=(4,4,1),
    prob_thresh = .05)
print("Stardist run completed")

# ======================================================================
# == save the output ==
# ======================================================================

print(f"Saving pkl object to {output_pkl_path}...")
with open(output_pkl_path, "wb") as f:  # open a text file
    pickle.dump(polys, f)
print("pkl Saved")