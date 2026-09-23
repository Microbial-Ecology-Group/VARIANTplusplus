Getting started with VARIANT++
-----------------

To get started, view the [Installation document](installation.md) to determine the best way to install VARIANT++ on your computing cluster.

Next, we will run a small sample dataset that comes with the pipeline source code, so we won't need to specify any input paths — they're already set as defaults in `params.config`. In this example, we'll assume conda is available in your computing environment.

During execution, the required tool dependencies come from the `VARIANT++_env` conda environment (see [installation.md](installation.md)).

```bash
# If you followed the instructions in the installation document, you must now navigate to the VARIANTplusplus directory
cd VARIANTplusplus

# Activate the conda environment created during installation
conda activate VARIANT++_env

# Run FastQC + MultiQC on the bundled test reads to check their quality first.
nextflow run main_VARIANT++.nf -profile local --pipeline eval_qc

# Explore the "test_results/" directory to view pipeline outputs
ls test_results/
```

The GSV pipeline runs as a series of steps (GSV_1 through GSV_5), each consuming the output of the previous one, or as a single `full_GSV_pipeline` run. View the [step-by-step guide](VARIANT++_step_by_step.md) for the full walkthrough and the [configuration doc](configuration.md) for details on every parameter.

```bash
# Step 1: QC trim + merge the bundled test reads with Trimmomatic and FLASH.
nextflow run main_VARIANT++.nf -profile local --pipeline GSV_1 --output test_results

# View the merged reads
ls test_results/Flash_reads/

# Continue on to deduplication (Step 2), host removal (Step 3), Kraken2
# extraction (Step 4), and GSV classification (Step 5) — see
# docs/VARIANT++_step_by_step.md for the exact --merged_reads glob to use
# at each step.
```

Alternatively, once you have a Kraken2 database and a host genome configured, you can run the entire pipeline in one command:

```bash
# You can use "--pipeline full_GSV_pipeline" to run all GSV steps end-to-end.
nextflow run main_VARIANT++.nf -profile local --pipeline full_GSV_pipeline --output test_GSV_output -w work_GSV
```
