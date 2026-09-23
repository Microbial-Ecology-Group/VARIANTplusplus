# Installation overview
This section will help you get started with running the VARIANT++ pipeline. This tutorial assumes you will be running the pipeline from a POSIX compatible system such as Linux, Solaris, or OS X.

There are several ways to install and run VARIANT++, based on what is easiest for your computing cluster. Based on our experience, we recommend the **conda** installation below — it's the actively-tested path and the one used throughout the [step-by-step guide](VARIANT++_step_by_step.md) and the [GSV benchmarking tutorial](GSV_benchmarking_tutorial.md).

Usually, you'll have to install nextflow, unless it's available to be loaded as a module in your HPC.
* [Install Nextflow](#installing-nextflow)

To make all of the bioinformatic tool dependencies available for use with VARIANT++, we have a few options:
* [Run VARIANT++ with Anaconda](#run-variant-with-anaconda) (recommended)
    * [Install miniconda without "sudo" permissions](#installing-miniconda-without-sudo-permissions)
* [Run VARIANT++ with Singularity](#run-variant-using-singularity)
* [Run VARIANT++ with Docker](#run-variant-using-docker)
* [Run VARIANT++ with locally installed tools](#local-installation-of-tools)

## Installing nextflow

```bash
# username and host address
$ ssh [USER]@[HOST]

# Check if you have nextflow installed,
$ nextflow -h

# If not available, install Nextflow
$ curl -s https://get.nextflow.io | bash
# If you do not have curl installed, try wget
# $ wget -qO- https://get.nextflow.io | bash

# give write permissions to user
chmod u+x nextflow

# move nextflow executable to a folder in your $PATH environment variable. For example:
mv nextflow $HOME/bin
```

## Run VARIANT++ with Anaconda
Requirements:
* Nextflow
* Anaconda or Miniconda

```bash
# Download the VARIANT++ repository
git clone https://github.com/Microbial-Ecology-Group/VARIANTplusplus.git

# Navigate into the directory
cd VARIANTplusplus

# Install mamba for faster installation
conda install mamba -n base -c conda-forge

# Create the VARIANT++ conda environment
conda env create -f envs/VARIANT++_env.yaml
conda activate VARIANT++_env

# Run VARIANT++, specifying the "local" profile since dependencies are now on your $PATH
nextflow run main_VARIANT++.nf -profile local --pipeline eval_qc

# If your computing cluster uses the slurm scheduler, modify "run_VARIANT++_slurm.sbatch" to
# accurately request computing resources, then run it using:
sbatch run_VARIANT++_slurm.sbatch
```

Alternatively, Nextflow can build the conda environment for you automatically via the `conda`/`conda_slurm` profiles (no manual `conda env create` needed):

```bash
nextflow run main_VARIANT++.nf -profile conda --pipeline eval_qc
```

### Installing miniconda without "sudo" permissions

We will go over a typical pipeline setup scenario in which you connect to a remote server, install miniconda (or use a local installation of anaconda), and download the pipeline source code. In cases where the Anaconda installation on your computing cluster is not updated or you are experiencing errors while installing packages, we recommend miniconda. Use [this site](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html) for further information on installing miniconda in a user-writable directory.

```bash
# Download miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.1.0-1-Linux-x86_64.sh
# Run installation, follow default options. Depending on your computing cluster,
# consider changing the install location to somewhere other than your home directory,
# which can have storage limits.
bash Miniconda3-py310_23.1.0-1-Linux-x86_64.sh

# Download the VARIANT++ repository
git clone https://github.com/Microbial-Ecology-Group/VARIANTplusplus.git

# Navigate into the directory
cd VARIANTplusplus

# Create and activate the conda environment
conda env create -f envs/VARIANT++_env.yaml
conda activate VARIANT++_env

# Run VARIANT++ using the "local" profile. If your computing cluster uses the
# slurm scheduler, use the "local_slurm" profile instead.
nextflow run main_VARIANT++.nf -profile local --pipeline eval_qc
```

There's one tool, Themisto, that ships as a compiled binary in the repo (`bin/themisto`) rather than through conda. To make it available on your `$PATH`:

```bash
cd bin/
pwd
echo 'export PATH="/path/to/your/VARIANTplusplus/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

## Run VARIANT++ using Singularity

Requirements:
* Nextflow
* Singularity

```bash
# Download the VARIANT++ repository
git clone https://github.com/Microbial-Ecology-Group/VARIANTplusplus.git
cd VARIANTplusplus

# Run command with singularity profile
nextflow run main_VARIANT++.nf -profile singularity --pipeline eval_qc
```

**Known limitation:** the `singularity`/`singularity_slurm`/`docker` profiles currently point at `enriquedoster/amrplusplus:latest`, the container image inherited from the AMR++ fork. It has not been verified to include Themisto, mSWEEP, or mGEMS, which the GSV_5/GSV_5_mGEMS steps require. Until an updated VARIANT++ image is published, the conda installation above is the recommended path.

## Run VARIANT++ using Docker

Requirements:
* Nextflow
* Docker

```bash
git clone https://github.com/Microbial-Ecology-Group/VARIANTplusplus.git
cd VARIANTplusplus

nextflow run main_VARIANT++.nf -profile docker --pipeline eval_qc
```

See the known limitation noted above — verify the container includes Themisto/mSWEEP/mGEMS before relying on it for GSV_5.

## Local installation of tools

Requirements:
* All [software requirements](requirements.md)
* Nextflow

If none of the above options work for your computing cluster, configure `config/local.config` to specify the absolute path to each required bioinformatic tool, or add them to your `$PATH` (e.g. by loading the appropriate modules). Then run with the "local" profile:

```bash
nextflow run main_VARIANT++.nf -profile local --pipeline eval_qc
```
