Usage
-----

### Display Help Message

The `help` parameter (or `--pipeline help`) displays the available pipeline options and commands.
```
nextflow run main_VARIANT++.nf --help
```

# Parameter selection
VARIANT++ comes with a default selection of parameters to perform a demonstration using example data provided in the "data/" directory. The example command below uses two types of parameters.

```
nextflow run main_VARIANT++.nf -profile conda --pipeline GSV_1
```

The `-profile` parameter only has a single dash (`-`), meaning it corresponds to a nextflow-specific parameter, unlike `--pipeline`. Examples of parameters with one dash include `-profile`, `-resume`, and `-config`. Further details on using the profile parameter, which determines how VARIANT++ runs on a computing cluster, can be found in the [configuration document](configuration.md). Below, we'll see how to change the pipeline-specific parameters, which are denoted using two dashes (`--`), such as `--pipeline` and `--reads`.

The VARIANT++ pipeline pulls information from various sources to determine the correct parameters for running the pipeline. This is the order in which Nextflow prioritizes parameters it receives:

1. Parameters specified on the command line (--something value)

2. Parameters provided using the -params-file option (params.config by default)

3. Config file specified using the -c my_config option (e.g. config/local.config)

4. The config file named nextflow.config in the current directory

5. The config file named nextflow.config in the workflow project directory

6. The config file $HOME/.nextflow/config

7. Values defined within the pipeline script itself (e.g. main_VARIANT++.nf)



## File Inputs

### Set custom sequence data

The `reads` parameter accepts paired-end sequence files in standard fastq/fastq.gz format (used by GSV_1, `eval_qc`, and `merge`).
```
$ nextflow run main_VARIANT++.nf --pipeline GSV_1 --reads "data/raw/*_R{1,2}.fastq.gz"
```

### Set merged/unmerged reads for downstream GSV steps

GSV_2 through GSV_5 pick up from the merged/unmerged FASTQ pairs produced by the previous step, via the `merged_reads` parameter (note the required single quotes around the glob):
```
$ nextflow run main_VARIANT++.nf --pipeline GSV_2 --merged_reads 'test_results/Flash_reads/*.{extendedFrags,notCombined}.fastq.gz'
```

### Set host genome

The `host` parameter accepts a fasta-formatted host genome, used by GSV_3 to remove host DNA.
```
$ nextflow run main_VARIANT++.nf --pipeline GSV_3 --host "data/host/chr21.fasta.gz"
```

### Set Kraken2 database

The `kraken_db` parameter points GSV_4 at a Kraken2 database used to extract reads belonging to your target species.
```
$ nextflow run main_VARIANT++.nf --pipeline GSV_4 --kraken_db "/path/to/kraken2_db"
```

### Set Themisto index and clustering file

GSV_5 and GSV_5_mGEMS classify reads against a Themisto pseudoalignment index and a mSWEEP clustering file.
```
$ nextflow run main_VARIANT++.nf --pipeline GSV_5 \
    --themisto_index "data/themisto" \
    --themisto_index_prefix "2025_themisto_index_no" \
    --clustering_file "data/themisto/2025_Mh_msweep_annotations_k8.tsv"
```

### Set adapter file

The `adapters` parameter accepts a fasta-formatted adapter file, used by Trimmomatic in GSV_1.
```
$ nextflow run main_VARIANT++.nf --pipeline GSV_1 --adapters "data/adapters/nextera.fa"
```

## File Outputs

### Set output and work directories

The `--output` parameter writes the results to the specified directory. As a Nextflow variable, the `-w` parameter only requires one dash and determines where the temporary files will be directed. Upon completing the run, you can delete the temporary work directory.
```
$ nextflow run main_VARIANT++.nf --pipeline GSV_1 --output "test_results/" -w "work_dir/"
```

## Resume a pipeline run

If the pipeline run is cancelled or stopped for whatever reason, using the same command with the addition of the `-resume` flag will attempt to pick up where the pipeline stopped. The "work" directory can take a lot of storage space and we recommend deleting it after completion of the pipeline.

```
$ nextflow run main_VARIANT++.nf --pipeline GSV_1 --output "test_results/" -w "work_dir/" -resume
```

## Trimming Options

### Set custom trimming parameters for Trimmomatic (GSV_1)

```
$ nextflow run main_VARIANT++.nf --pipeline GSV_1 \
    --reads "data/raw/*_R{1,2}.fastq.gz" \
    --leading 3 \
    --trailing 3 \
    --minlen 36 \
    --slidingwindow "4:15" \
    --adapters "data/adapters/nextera.fa" \
    --output "test_results/"
```

## Kraken2 Options

### Set the Kraken2 confidence score used during extraction (GSV_4)

```
$ nextflow run main_VARIANT++.nf --pipeline GSV_4 \
    --kraken_db "/path/to/kraken2_db" \
    --kraken_confidence 0.1 \
    --output "test_results/"
```

## Set number of threads to use for each process (when possible)

```
$ nextflow run main_VARIANT++.nf --pipeline GSV_1 --threads 8
```

## Run VARIANT++ with SLURM

Using the `local_slurm` (or `conda_slurm`/`singularity_slurm`) profile submits each individual process in the pipeline as its own SLURM job. Adding the `-bg` flag runs the Nextflow head process itself in the background:

```
nextflow run main_VARIANT++.nf -profile local_slurm --pipeline full_GSV_pipeline -bg > bg-log.out
```
