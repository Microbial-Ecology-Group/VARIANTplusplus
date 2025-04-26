# VARIANT++ step by step guide

Make sure that you have:
1. An updated VARIANT++_env and github repository
** Install it with conda from the VARIANT++ directory with: 
```
git clone https://github.com/Microbial-Ecology-Group/VARIANTplusplus.git

cd VARIANTplusplus/

conda env create -f envs/VARIANT++_env.yaml 

conda activate VARIANT++_env
# if that doesn't work, use "source" instead of "conda" for the command above.

```
2. The coreNT kraken database. It's about 240 GB and can be downloaded like this:
```
# Downlaod database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20241228.tar.gz

# Make directory for db contents
mkdir -p k2_core_nt_20241228

# unzip it 
tar -xzvf k2_core_nt_20241228.tar.gz -C k2_core_nt_20241228/ 
```


## Step 1: QC trimming, merge reads, and deduplicate

Parameters that have to change:
* `--pipeline` ==> `--pipeline GSV_1`
* `--reads`  ==> `--reads "/path/to/your/reads/*R{1,2}.fastq.gz"`

Defaults for Trimmomatic
* `--leading` = 3
* `--trailing` = 3
* `--slidingwindow` = "4:15"
* `--minlen` = 36

Optional
* `--threads` = 4


Example command:

```
 nextflow run main_VARIANT++.nf --pipeline GSV_1 --output GSV_analysis --reads "data/raw/*_R{1,2}.fastq.gz" 
```


## Step 2: Remove host DNA

Parameters that have to change:
* `--pipeline` ==> `--pipeline GSV_2`
* `--merged_reads`  ==> `--merged_reads 'GSV_analysis/Deduped_reads/*_{merged,unmerged}.dedup.fastq.gz'`
    * If you named you used "--output GSV_analysis", then the command below should work with your data, otherwise just change it to match your output directory name.
    * Also note, that this parameter requires the use of single quotes `'`, anything else will not work. 
* `host` ==> `--host "/path/to/your/host/chr21.fasta.gz"`

Example command:
```
nextflow run main_VARIANT++.nf --pipeline GSV_2 --output GSV_analysis --merged_reads 'GSV_analysis/Deduped_reads/*_{merged,unmerged}.dedup.fastq.gz'
```




## Step 3: Filter reads with kraken

Parameters that have to change:
* `--pipeline` ==> `--pipeline GSV_3`
* `--merged_reads`  ==> `--merged_reads 'GSV_analysis/Deduped_reads/*{merged,unmerged}.dedup.fastq.gz'`
* `--kraken_db` ==> `--kraken_db /path/to/k2_core_nt_20241228`


Example command:

```
nextflow run main_VARIANT++.nf --pipeline GSV_3 --output GSV_analysis --merged_reads 'GSV_analysis/HostRemoval/NonHostFastq/*_{merged,unmerged}.non.host.fastq.gz'
```

## Step 4: Perform classification with themisto and mSweep

Parameters that have to change:
* `--pipeline` ==> `--pipeline GSV_4`
* `--merged_reads`  ==> `--merged_reads 'GSV_analysis/MicrobiomeAnalysis/Kraken/extracted_reads/*_{merged,unmerged}.dedup.fastq.gz'`


Example command:
```
nextflow run main_VARIANT++.nf --pipeline GSV_4 --output GSV_analysis --merged_reads 'GSV_analysis/MicrobiomeAnalysis/Kraken/extracted_reads/*_Mh_extracted_{merged,unmerged}.fastq.gz'
```

## Explore the results

Check the "Results" folder for the kraken analytic matrix and the mSweep results. The "mSweep_results_summary.tsv" file contains all results for each the merged and unmerged reads, but the "mSweep_results_count_matrix.tsv" file has a count matrix with the combined results.