import sys
import os
import bin2cell as b2c


# Redirect stderr/stdout to Snakemake log file
log_file = snakemake.log[0]
sys.stderr = open(log_file, "w")
sys.stdout = open(log_file, "w")

print(f"--- Starting AnnData generation for sample {snakemake.wildcards.sample_id} ---")

# 1. Get input paths from Snakemake object
# These correspond to the 'input:' section of the rule
# path_seq_data = snakemake.input.seq_data + "/square_002um/"
path_seq_data = os.path.join(snakemake.input.seq_data, "square_002um")
path_source_image = snakemake.input.raw_image

print(f"Input Visium data directory: {path_seq_data}")
print(f"Input source image path: {path_source_image}")

# 2. Get output path from Snakemake object
output_adata_path = snakemake.output.adata

print(f"Output AnnData path: {output_adata_path}")

# 3. generate the AnnData object
print("Reading Visium data using bin2cell.read_visium...")
adata = b2c.read_visium(path_seq_data, source_image_path=path_source_image)
print("Data loaded successfully.")

print("Making variable names unique...")
adata.var_names_make_unique()
print("Variable names updated.")

# Save the AnnData object to the path defined in the rule's output
print(f"Saving AnnData object to {output_adata_path}...")


adata.write_h5ad(output_adata_path, compression="gzip")
print("AnnData object saved successfully.")