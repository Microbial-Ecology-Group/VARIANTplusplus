Overview
--------
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A50.25.1-brightgreen.svg)](https://www.nextflow.io/)

# *Under construction*

# VARIANT++ bioinformatic pipeline
(https://megares.meglab.org/)

VARIANT++ is a bioinformatic pipeline meant to aid in the analysis of raw sequencing reads to characterize the profile of phylogenetic sequence cluster variants (PSVs) for a certain bacterial species. In this repository you'll find code and tutorials to replicate our PSV annotation scheme using kraken2 with any species of interest.

With VARIANT++, you will obtain alignment count files for each sample that are combined into a count matrix that can be analyzed using any statistical and mathematical techniques that can operate on a matrix of observations.

More Information
----------------

- [Installation](https://github.com/Microbial-Ecology-Group/VARIANTplusplus/blob/master/docs/installation.md)
- [Usage](https://github.com/Microbial-Ecology-Group/VARIANTplusplus/blob/master/docs/usage.md)
- [Configuration](https://github.com/Microbial-Ecology-Group/VARIANTplusplus/blob/master/docs/configuration.md)
- [Output](https://github.com/Microbial-Ecology-Group/VARIANTplusplus/blob/master/docs/output.md)
- [Dependencies](https://github.com/Microbial-Ecology-Group/VARIANTplusplus/blob/master/docs/dependencies.md)
- [Software Requirements](https://github.com/Microbial-Ecology-Group/VARIANTplusplus/blob/master/docs/requirements.md)
- [FAQs](https://github.com/Microbial-Ecology-Group/VARIANTplusplus/blob/master/docs/FAQs.md)
- [Details on VARIANT++ updates](https://github.com/Microbial-Ecology-Group/VARIANTplusplus/blob/master/docs/update_details.md)
- [Contact](https://github.com/Microbial-Ecology-Group/VARIANTplusplus/blob/master/docs/contact.md)



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

# Run command to perform the demonstration pipeline using the conda profile.
nextflow run main_VARIANT++.nf -profile conda

# The first time this can take 5-10 mins (or more) depending on your internet speed because it is installing a conda environment. Subsequent runs will skip this step automatically.
```
Now, you can check out the results in the newly created "test_results" directory.

# Using VARIANT++ to analyze your data

VARIANT++ is customizable to suit your computing needs and analyze your data. Primarily, the ```-profile``` paramater allows you to choose between running VARIANT++ using a singularity container, docker container, anaconda packages, or a local installation of your software. 
All parameters used to control how VARIANT++ analyzes your data can also be changed as needed in a variety of ways. For full information, review this [configuration document.](docs/configuration.md)


Below is a brief example, the default parameters were run using this command:

```nextflow run main_VARIANT++.nf -profile conda```

To change the reads that were analyzed, you should specify the ```--reads`` parameters. Here, we can use regular expressions to point to your samples in a different directory.
```bash
nextflow run main_VARIANT++.nf -profile conda --reads "path/to/your/reads/*_R{1,2}.fastq.gz" 
```
