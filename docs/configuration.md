# Contents

* [Configuration](#configuration)
* [Customize environmental variables using profiles](#customize-environment-variables-using-profiles)
* [Customize VARIANT++ pipeline parameters](#customize-variant-pipeline-parameters)
  * [Modifying the params.config file](#modifying-the-paramsconfig-file)
  * [Modifying parameters using the command-line](#modifying-parameters-using-the-command-line)
    * [Analyzing your samples](#analyzing-your-samples)
    * [Running with Kraken2](#running-with-kraken2)
* [Selecting the right pipeline](#selecting-the-right-pipeline)

## Configuration
-------------

The pipeline source code comes with two configuration files that can be used to set environment variables and default command-line options. These configuration files can be found in the root source code directory and are called **nextflow.config** and **params.config**.

The **nextflow.config** file mainly contains parameters regarding how VARIANT++ will run on your computing cluster using the ```-profile``` parameter.

The **params.config** contains parameters that control which files are being analyzed and parameters for the software in the pipeline. Setting the variables in the **params.config** beforehand may be useful in situations when you do not want to specify a long list of options from the command line or want to have a separate file for each project. You can modify these files, save the changes, and run the pipeline directly. More details below.


## Customize Environment Variables using profiles
----------------------------------------------

The **nextflow.config** contains a section that allows the use of environment "profiles" when running VARIANT++. Further information for each profile can be found within the `/config` directory. In brief, profiles allow control over how the pipeline is run on different computing clusters. We recommend the "conda" profile, which builds a single conda environment (`envs/VARIANT++_env.yaml`) containing all the required bioinformatic tools.

We make the following profiles available to suit your computing needs: "local", "local_slurm", "conda", "conda_slurm", "singularity", "singularity_slurm", and "docker". You specify which profile to use with the `-profile` flag.


```bash
profiles {
  local {
    includeConfig "config/local.config"
  }
  local_slurm {
    includeConfig "config/local_slurm.config"
  }
  conda {
    includeConfig "config/conda.config"
    conda.enabled = true
    conda.cacheDir = "$baseDir/envs/"
    conda.useMamba = true
    conda.createTimeout = '30 min'
  }
  docker {
    includeConfig "config/local.config"
    docker.enabled = true
    process.container = 'enriquedoster/amrplusplus:latest'
  }
  singularity {
    includeConfig "config/singularity.config"
    singularity.enabled = true
    singularity.autoMounts = true
    singularity.cacheDir = "$baseDir/envs/"
  }
  conda_slurm {
    includeConfig "config/conda_slurm.config"
    conda.cacheDir = "$baseDir/envs/"
    conda.enabled = true
    conda.useMamba = true
    conda.createTimeout = '30 min'
  }
   singularity_slurm {
    includeConfig "config/singularity_slurm.config"
    singularity.enabled = true
    singularity.autoMounts = true
    singularity.cacheDir = "$baseDir/envs/"
  }
}
```

Note: the `docker`/`singularity`/`singularity_slurm` profiles currently point at `enriquedoster/amrplusplus:latest`, a container image inherited from the AMR++ fork. It has not been verified to include Themisto, mSWEEP, or mGEMS — the `conda` profile is the actively-tested path for VARIANT++ today (see [installation.md](installation.md)).

## Customize VARIANT++ pipeline parameters
------------------------------

The params section allows you to set the different command-line options that can be used within the pipeline. Here, you can specify input/output options, trimming options, and classification options.

### Modifying the params.config file
Below is a list of the parameters VARIANT++ uses by default. They can be found in the `params.config` file in the main directory. These parameters can be modified by changing this file or specifying any of these parameters on the command line using a double dash, like this: `--reads "path/to/your/reads/*_R{1,2}.fastq.gz"`. Otherwise, change the parameters in the `params.config` file prior to running the VARIANT++ pipeline.

These are all of the parameters used by VARIANT++ (see `params.config` for the authoritative, up-to-date list):
```bash
params {
    /* Display help message */
    help = false

    /* Location of forward and reverse read pairs */
    reads = "${baseDir}/data/raw/*_R{1,2}.fastq.gz"

    /* Location of reference/host genome */
    host = "${baseDir}/data/host/chr21.fasta.gz"

    /* Optionally, the location of bwa host index files (path + wildcard *) */
    host_index = ""

    /* Output directory */
    output = "test_results"

    /* Default memory to run clumpify (GSV_2 deduplication) */
    clumpify_mem_gb = 8

    /* Kraken2 confidence score */
    kraken_confidence = 0.1

    /* Kraken2 database location */
    kraken_db = ""
    krakendb_inter = ""
    kraken_options = ""

    /* Kraken2 db for genus/species confirmation */
    confirmation_db = ""

    /* Optional flags for extract_kraken_reads.py */
    extract_reads_taxid = "75985"
    extract_reads_options_single = "--include-children"
    extract_reads_options_double = "--include-children"

    /* Location of reference genome directory */
    genome_ref_dir = ""

    coverage_threshold = 0.0001
    dedup_sam = "Y"

    /* Number of threads */
    threads = 4

    /* Trimmomatic trimming parameters */
    adapters = "${baseDir}/data/adapters/nextera.fa"
    leading = 3
    trailing = 3
    slidingwindow = "4:15"
    minlen = 36

    /* multiQC config directory */
    multiqc = "$baseDir/data/multiqc"

    /* Optional flag for running GSV subworkflows starting with merged/unmerged reads (GSV_2 - GSV_5) */
    merged_reads = "${params.output}/Flash_reads/*.{extendedFrags,notCombined}.fastq.gz"

    QC_dir = "${params.output}/QC_trimming/Paired/"
    QC_prefix = "QC_trimmed"

    /* Themisto pseudoalignment index (used by GSV_5) */
    themisto_index = "$baseDir/data/themisto"
    themisto_index_prefix = "2025_themisto_index_no"
    clustering_file = "$baseDir/data/themisto/2025_Mh_msweep_annotations_k8.tsv"

    /* mGEMS binning (GSV_5_mGEMS) */
    run_mgems = true
    msweep_write_probs = true
    mgems_write_assignment_table = true
    mgems_min_abundance = 0.01
}
```
### Modifying parameters using the command-line

#### Analyzing your samples
------
If you intend to run multiple samples in parallel, you must specify a glob pattern for your sequence data as shown for the **reads** parameter. For more information on globs, please see this related [article](https://en.wikipedia.org/wiki/Glob_(programming)).

For example, the default parameters can be used to run the pipeline with this command:

```bash
nextflow run main_VARIANT++.nf -profile conda --pipeline GSV_1
```

This will run the default samples through the pipeline, as seen under the `--reads` parameter. To change the reads being analyzed, specify the `--reads` parameter on the command line:

```bash
nextflow run main_VARIANT++.nf -profile conda --pipeline GSV_1 --reads "path/to/your/reads/*_R{1,2}.fastq.gz"
```

#### Running with Kraken2
-----
GSV_4 uses Kraken2 to extract reads belonging to your target species before classification. You need to point `--kraken_db` at a Kraken2 database on your system — see [docs/VARIANT++_step_by_step.md](VARIANT++_step_by_step.md) for how to download the coreNT database used for GSV benchmarking.

```bash
nextflow run main_VARIANT++.nf -profile conda --pipeline GSV_4 --kraken_db /path/to/your/kraken2_db
```

## Selecting the right pipeline

VARIANT++ lets you run different components of the GSV pipeline at a time by specifying the `--pipeline` flag.

Main pipeline option
  * Run all GSV steps end-to-end (QC trim + merge > dedup > host removal > Kraken2 extraction > Themisto/mSWEEP classification)
    ```bash
    --pipeline full_GSV_pipeline
    ```

Pipeline components
  * Evaluate raw-read QC with FastQC/MultiQC
    ```bash
    --pipeline eval_qc
    ```
  * Merge paired-end reads with FLASH only
    ```bash
    --pipeline merge
    ```
  * GSV_1: QC trimming (Trimmomatic) + merge reads (FLASH)
    ```bash
    --pipeline GSV_1
    ```
  * GSV_2: deduplicate merged/unmerged reads (BBMap clumpify)
    ```bash
    --pipeline GSV_2
    ```
  * GSV_3: remove host DNA (BWA)
    ```bash
    --pipeline GSV_3
    ```
  * GSV_4: extract target-species reads with Kraken2
    ```bash
    --pipeline GSV_4
    ```
  * GSV_5: classify GSVs with Themisto + mSWEEP
    ```bash
    --pipeline GSV_5
    ```
  * GSV_5_mGEMS: classify and bin reads per GSV with Themisto + mSWEEP + mGEMS
    ```bash
    --pipeline GSV_5_mGEMS
    ```
