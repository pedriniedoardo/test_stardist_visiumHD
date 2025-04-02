rule downloadData:
    input:
        sample_file = config["samples_csv"],
    output:
        counts_tar = temp(config["out_location"] + "rawData/{sample_id}/binned_outputs.tar.gz"),
        counts_dir= temp(directory(config["out_location"] + "rawData/{sample_id}/binned_outputs")),
        rawImage = temp(config["out_location"] + "rawData/{sample_id}/rawImage.img"),
    log:
        'logs/{sample_id}/downloadData.log'
    benchmark:
        'benchmarks/{sample_id}/downloadData.txt'
    params:
        # Get links from the SAMPLES dictionary (ensure SAMPLES is loaded)
        link_counts = lambda wildcards: SAMPLES[wildcards.sample_id]['link_counts'],
        link_rawImage = lambda wildcards: SAMPLES[wildcards.sample_id]['link_rawImage'],
        # data_count = config["out_location"] + "rawData/{sample_id}/binned_outputs.tar.gz",
        extract_parent_dir = directory(config["out_location"] + "rawData/{sample_id}/"),
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

        echo "Download and extraction for {wildcards.sample_id} finished."
        '''

rule generateAdata:
    input:
        seq_data = rules.downloadData.output.counts_dir,
        raw_image = rules.downloadData.output.rawImage
    output:
        adata=config["out_location"] + "adata/{sample_id}_adata.h5ad"
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
