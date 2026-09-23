Software Requirements
---------------------
To run VARIANT++, you will need the following tools installed on your server or local machine.

  - Anaconda or Miniconda (Required)
    - Visit this website for further information: https://docs.anaconda.com/anaconda/install/
  - Java 8+ (Required)
  - Nextflow (Required)

NOTE: If you choose not to install anaconda, you will need to download each of the required dependencies and add their executable paths to your `.bashrc` file to run the pipeline. A list of these dependencies can be found in the [Dependencies](dependencies.md) section of this document. Themisto specifically ships as a compiled binary in `bin/themisto` — see [installation.md](installation.md) for how to add it to your `$PATH`.

Once Java, Nextflow, and Anaconda/Miniconda are available, see [installation.md](installation.md) to create the `VARIANT++_env` conda environment and run a first pipeline step:

```bash
nextflow run main_VARIANT++.nf -profile conda --pipeline eval_qc
```
