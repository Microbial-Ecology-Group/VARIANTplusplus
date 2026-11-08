[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)


# VARIANT++ bioinformatic pipeline

VARIANT++ is a Nextflow pipeline for characterizing the profile of genomic sequence variants (GSVs) of a bacterial species of interest directly from raw metagenomic sequencing reads. GSVs are clusters of closely related genomes (grouped by average nucleotide identity) within a single species.

The pipeline performs read QC and merging, deduplication, host DNA removal, and Kraken2 filtration to extract reads from the target taxa, then classifies those reads into GSVs using Themisto and mSWEEP. By default, VARIANT++ is set up for *Mannheimia haemolytica*, but this repository also includes the code and a tutorial to benchmark and build GSVs for any species of interest.

With VARIANT++, you will obtain GSV abundance results for each sample that are combined into a count matrix that can be analyzed using any statistical and mathematical techniques that can operate on a matrix of observations.

More Information
----------------

- [Getting started](docs/GettingStarted.md)
- [Installation](docs/installation.md)
- [Usage](docs/usage.md)
- [VARIANT++ step by step guide](docs/VARIANT++_step_by_step.md)
- [Configuration](docs/configuration.md)
- [Output](docs/output.md)
- [Dependencies](docs/dependencies.md)
- [Software Requirements](docs/requirements.md)
- [GSV benchmarking tutorial (building GSVs for a new species)](docs/GSV_benchmarking_tutorial.md)
- [FAQs](docs/FAQs.md)
- [Changelog](docs/CHANGELOG.md)
- [Contact](docs/contact.md)



## VARIANT++ demonstration

If anaconda is already installed and nextflow is working, we'll just need to download the VARIANT++ github repository. Please review the [installation document](docs/installation.md) for alternative methods to install VARIANT++ in your computing environment.

```bash
# Install mamba for faster installation
conda install mamba -n base -c conda-forge
```

Clone the VARIANT++ repository.

```bash
git clone https://github.com/Microbial-Ecology-Group/VARIANTplusplus.git
```

Navigate into the VARIANT++ repository and run the test command.
```bash
cd VARIANTplusplus

# Print the help message with the available pipelines
nextflow run main_VARIANT++.nf --pipeline help

# Run FastQC + MultiQC on the bundled test reads using the conda profile.
nextflow run main_VARIANT++.nf -profile conda --pipeline eval_qc

# The first time this can take 5-10 mins (or more) depending on your internet speed because it is installing a conda environment. Subsequent runs will skip this step automatically.
```
Now, you can check out the results in the newly created "test_results" directory. See [Getting started](docs/GettingStarted.md) to run the remaining GSV steps on the test data.

# Using VARIANT++ to analyze your data

VARIANT++ is customizable to suit your computing needs and analyze your data. Primarily, the ```-profile``` parameter allows you to choose between running VARIANT++ using a singularity container, docker container, anaconda packages, or a local installation of your software, and the ```--pipeline``` parameter selects which steps to run (e.g. ```full_GSV_pipeline``` or the individual ```GSV_1``` to ```GSV_5``` steps).
All parameters used to control how VARIANT++ analyzes your data can also be changed as needed in the [params.config](params.config) file or on the command line. For full information, review this [configuration document.](docs/configuration.md)


Below is a brief example, running the full pipeline with the default parameters:

```nextflow run main_VARIANT++.nf -profile conda --pipeline full_GSV_pipeline --kraken_db /path/to/kraken_db/```

To change the reads that were analyzed, you should specify the ```--reads``` parameter. Here, we can use regular expressions to point to your samples in a different directory.
```bash
nextflow run main_VARIANT++.nf -profile conda --pipeline full_GSV_pipeline --kraken_db /path/to/kraken_db/ --reads "path/to/your/reads/*_R{1,2}.fastq.gz"
```

To analyze a species other than *M. haemolytica*, follow the [GSV benchmarking tutorial](docs/GSV_benchmarking_tutorial.md) to build and benchmark your GSVs, then update the Themisto index, mSWEEP clustering file, and Kraken2 taxid in [params.config](params.config).
