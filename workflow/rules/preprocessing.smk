rule downloadData:
    input:
        sample_file = config["samples_csv"],
    output:
        counts_tar = temp(config["out_location"] + "rawData/{sample_id}/binned_outputs.tar.gz"),
        counts_dir= temp(directory(config["out_location"] + "rawData/{sample_id}/binned_outputs")),
        rawImage = temp(config["out_location"] + "rawData/{sample_id}/rawImage.img"),
        parquet_keep = config["out_location"] + "rawData/{sample_id}/tissue_positions.parquet"
        # counts_tar = temp(f"{config["out_location"]}rawData/{sample_id}/binned_outputs.tar.gz"),
        # counts_dir= temp(directory(f"{config["out_location"]}rawData/{sample_id}/binned_outputs")),
        # rawImage = temp(f"{config["out_location"]}rawData/{sample_id}/rawImage.img")
    log:
        'logs/{sample_id}/downloadData.log'
    benchmark:
        'benchmarks/{sample_id}/downloadData.txt'
    threads:
        config["CPU_download"]
    resources:
        mem_gb = config["RAM_download"]
    params:
        # Get links from the SAMPLES dictionary (ensure SAMPLES is loaded)
        link_counts = lambda wildcards: SAMPLES[wildcards.sample_id]['link_counts'],
        link_rawImage = lambda wildcards: SAMPLES[wildcards.sample_id]['link_rawImage'],
        # data_count = config["out_location"] + "rawData/{sample_id}/binned_outputs.tar.gz",
        extract_parent_dir = directory(config["out_location"] + "rawData/{sample_id}/"),
        # extract_parent_dir = directory(f"{config["out_location"]}rawData/{sample_id}/")
        parquet_input = config["out_location"] + "rawData/{sample_id}/binned_outputs/square_002um/spatial/tissue_positions.parquet"
    shell:
        '''
        # Exit on error, undefined variable, or pipe failure
        set -euo pipefail

        echo "Starting download and extraction for {wildcards.sample_id}"

        # Download counts data, referencing output path
        echo "Downloading counts data..."
        wget --quiet --show-progress -O {output.counts_tar} {params.link_counts}

        # Download raw image, referencing output path
        echo "Downloading raw image..."
        wget --quiet --show-progress -O {output.rawImage} {params.link_rawImage}

        # Extract counts data into the target PARENT directory
        # tar -C extracts *into* the specified directory.
        # Assuming the tarball contains a 'binned_outputs' folder,
        # extracting to its parent will create the desired output directory.
        echo "Extracting counts data..."
        tar -xzf {output.counts_tar} -C {params.extract_parent_dir}

        # keep the tissue_positions.parquet file to the output directory
        cp {params.parquet_input} {output.parquet_keep}

        echo "Download and extraction for {wildcards.sample_id} finished."
        '''

rule generateAdata:
    input:
        seq_data = rules.downloadData.output.counts_dir,
        raw_image = rules.downloadData.output.rawImage
    output:
        adata = config["out_location"] + "adata/{sample_id}_adata.h5ad"
    conda:
        config["env_bin2cell"]
    log:
        'logs/{sample_id}/generateAdata.log'
    benchmark:
        'benchmarks/{sample_id}/generateAdata.txt'
    threads:
        config["CPU_adata"]
    resources:
        mem_gb = config["RAM_adata"]
    script:
        "../scripts/00_generate_adata.py"

rule runStardist:
    input:
        raw_image = rules.downloadData.output.rawImage
    output:
        pkl = config["out_location"] + "stardist/{sample_id}_nuclei_polys.pkl"
    conda:
        config["env_stardist"]
    log:
        'logs/{sample_id}/runStardist.log'
    benchmark:
        'benchmarks/{sample_id}/runStardist.txt'
    threads:
        config["CPU_stardist"]
    resources:
        mem_gb = config["RAM_stardist"]
    params:
        hist_img = lambda wildcards: SAMPLES[wildcards.sample_id]['hist_img']
    script:
        "../scripts/01_run_stardist.py"

rule runBin2cell:
    input:
        adata = rules.generateAdata.output.adata,
        parquet = rules.downloadData.output.parquet_keep,
        stardist = rules.runStardist.output.pkl
    output:
        nuclei_grouped_pkl = config["out_location"] + "bin2cell/{sample_id}_nuclei_grouped_geometry.pkl",
        nuclei_expanded_pkl = config["out_location"] + "bin2cell/{sample_id}_nuclei_expanded_geometry.pkl",
        nuclei_grouped_adata = config["out_location"] + "bin2cell/{sample_id}_nuclei_grouped_geometry.h5ad",
        nuclei_expanded_adata = config["out_location"] + "bin2cell/{sample_id}_nuclei_expanded_geometry.h5ad"
    conda:
        config["env_bin2cell"]
    log:
        'logs/{sample_id}/runBin2cell.log'
    benchmark:
        'benchmarks/{sample_id}/runBin2cell.txt'
    threads:
        config["CPU_bin2cell"]
    resources:
        mem_gb = config["RAM_bin2cell"]
    params:
        species = lambda wildcards: SAMPLES[wildcards.sample_id]['species']
    script:
        "../scripts/02_run_bin2cell.py"

rule runRCTD:
    input:
        ref = lambda wildcards: SAMPLES[wildcards.sample_id]['RCTD_ref'],
        # nuclei_grouped_adata = rules.runBin2cell.output.nuclei_grouped_adata,
        # nuclei_expanded_adata = rules.runBin2cell.output.nuclei_grouped_adata
        # keep this only for testing
        nuclei_grouped_adata = "/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/Mouse_Embryo/Mouse_Embryo_nuclei_grouped_tiny.rds",
        nuclei_expanded_adata = "/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/Mouse_Embryo/Mouse_Embryo_expanded_nuclei_tiny.rds"
    output:
        csv_filters = config["out_location"] + "RCTD/{sample_id}_RCTD_filters.csv"
    conda:
        config["env_RCTD"]
    log:
        'logs/{sample_id}/runRCTD.log'
    benchmark:
        'benchmarks/{sample_id}/runRCTD.txt'
    threads:
        config["CPU_RCTD"]
    resources:
        mem_gb = config["RAM_RCTD"]
    params:
        # ref = lambda wildcards: SAMPLES[wildcards.sample_id]['RCTD_ref']
    script:
        "../scripts/03_run_RCTD.R"

rule runAdataMerge:
    input:
        # csv_filters = rules.runRCTD.output.csv_filters,
        # nuclei_grouped_adata = rules.runBin2cell.output.nuclei_grouped_adata,
        # nuclei_expanded_adata = rules.runBin2cell.output.nuclei_grouped_adata
        csv_filters = "/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/Mouse_Embryo/Mouse_Embryo_RCTD_filters.csv",
        nuclei_grouped_adata = "/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/Mouse_Embryo/Mouse_Embryo_nuclei_grouped.h5ad",
        nuclei_expanded_adata = "/media/edo/ExtremeSSD/training/test_snakemake/test_cosr_spatial/data/Mouse_Embryo/Mouse_Embryo_expanded_nuclei.h5ad"

    output:
        merged_adata = config["out_location"] + "adata/{sample_id}_adata_final.h5ad",
    conda:
        config["env_bin2cell"]
    log:
        'logs/{sample_id}/runAdataMerge.log'
    benchmark:
        'benchmarks/{sample_id}/runAdataMerge.txt'
    threads:
        config["CPU_AdataMerge"]
    resources:
        mem_gb = config["RAM_AdataMerge"]
    params:
        # ref = lambda wildcards: SAMPLES[wildcards.sample_id]['RCTD_ref']
    script:
        "../scripts/04_run_mergeAdata.py"