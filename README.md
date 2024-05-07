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

# Optional flags

## SNP verification

VARIANT++ now works in conjuction with a [custom SNP verification software](https://github.com/Isabella136/VARIANTPlusPlus_SNP) to evaluate alignments to gene accessions requiring SNP confirmation to confer resistance. To include this workflow, include the ```--snp Y``` flag in your command like this:

```bash
nextflow run main_VARIANT++.nf -profile conda --snp Y
```
This will create with the standard count table (VARIANT_analytic_matrix.csv) in addition to a count matrix with SNP confirmed counts (SNPconfirmed_VARIANT_analytic_matrix.csv).

## Deduplicated counts

Another option is to include results for deduplicated counts by using the ```--deduped Y``` flag in your command.

```bash
nextflow run main_VARIANT++.nf -profile conda --snp Y --deduped Y
```

With this flag, VARIANT++ will extract the deduplicated alignments to MEGARes also output a count matrix with deduplicated counts. Since also we included the ```--snp Y``` flag, we will end up with 4 total output count matrices.

# Choosing the right pipeline

VARIANT++ analyzes data by combining workflows that takes a set of sequencing reads through various bioinformatic software. We recommend our standard VARIANT++ pipeline as a comprehensive way to start from raw sequencing reads, QC assessment, host DNA removal, and resistome analysis with MEGARes. However, users might only want to replicate portions of the pipeline and have more control over their computing needs. Using the ```--pipeline``` parameter, users can now change how VARIANT++ runs.



## Pipeline workflows
*  omitting the ```--pipeline``` flag or using ```--pipeline demo```    
    * Simple demonstration on test data

* ```--pipeline standard_VARIANT```   
    * Steps: QC trimming > Host DNA removal > Resistome alignment > Resistome results

* ```--pipeline fast_VARIANT```
    * This workflow simply skips host removal to speed up analysis.
    * Steps: QC trimming > Resistome alignment > Resistome results

* ```--pipeline standard_VARIANT_wKraken```
    * This workflow adds microbiome analysis with kraken. It requires having a local kraken database. The minikraken_8GB_202003 will be downloaded automatically and requires ~8GB of space. Otherwise, you can specify the location to your own database with the flag, ```--kraken_db "/Path/to/KrakenDb/"```
    * Steps:
        * QC trimming > Host DNA removal > Resistome alignment > Resistome results 
        * Non-host reads > Microbiome analysis

## Pipeline subworkflows
* ```--pipeline eval_qc```  
    * Evaluate sample QC 
* ```--pipeline trim_qc```  
    * QC trimming using trimmomatic 
* ```--pipeline rm_host```  
    * Align reads to host DNA using bwa and remove contaminants 
* ```--pipeline resistome```  
    * Align reads to MEGARes using bwa, perform rarefaction and resistome analysis
* ```--pipeline kraken```  
    * Classify reads taxonomically using kraken.
* ```--pipeline bam_resistome```
    * This will run the resistome pipeline starting with bam files from a previous alignment to MEGARes.
    * Need to include ```--bam_files "Path/to/BAM/*.bam"``` in the command line.

## Example command
In the following example, we'll choose to run the standard VARIANT++ workflow, which includes QC trimming, host removal, and Resistome analysis. Since we included the ```--snp Y --deduped Y``` flags, we'll also get ouput for deduped counts and SNP confirmed counts.

Alternatively, you can modify all of these variables and more in the "params.config" file which will be loaded automatically. Just make sure to include the "-profile" and "--pipeline" flags. More information [in this document](docs/configuration.md)

```bash
# Remember to update the --reads flag to match your read location
nextflow run main_VARIANT++.nf -profile conda --pipeline standard_VARIANT --reads "path/to/your/reads/*_R{1,2}.fastq.gz" --snp Y --deduped Y
```

# Repeating all steps

## Step 1; Kraken Database for BRD (morax, pasteur, mycoplasma) was created. These Genomes were downloaded using NCBI genome download tool. Refer to https://anaconda.org/bioconda/ncbi-genome-download 

module load Kraken2/2.0.8-beta-foss-2018b-Perl-5.28.0
module load WebProxy
kraken2-build --download-taxonomy --db BRD_bacteria_kraken_db --use-ftp

#Adding Bacteria Genomes to build the database (genomes downloaded from NCBI <refer to step 1> was stored inside refseq_bacteria/refseq/bacteria)
find refseq_bacteria/refseq/bacteria/*/ -name '*.fna.gz' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db BRD_bacteria_kraken_db

#Adding BRD family genomes into kraken database (genomes downloaded from NCBI was stored inside the directory family_genome_database/)
find family_genome_database/*/ -name '*.fna.gz' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db BRD_bacteria_kraken_db

#Building Kraken Database (Note; Once you have added all the necessary taxonamy and genomes you can run the command below)
kraken2-build --build --db BRD_bacteria_kraken_db --threads 24


Step 2; Kraken Database with Mannheimia genomes 

module load Kraken2/2.0.8-beta-foss-2018b-Perl-5.28.0
module load WebProxy
kraken2-build --download-taxonomy --db kraken_db_Mh --use-ftp

#Adding Mh fasta files to the database (All Mh genomes downloaded from the NCBI was stored inside Mh_all_strains_2023_02 directory)
find /scratch/user/sal/Genome_Databases/Mh_all_strains_2023_02 -name '*.fna' -print0 | xargs -0 -I{} -n1 kraken2-build --add-to-library {} --db kraken_db_Mh
kraken2-build --build --db kraken_db_Mh --threads 24

Step 3; Transferring Files 

Using Globus file transfer fastq files out of sequencer was transferred into the server/hpc.

Step 4; Classification

Using Kraken2 database created in step 2 to classify fastq files from step 3

kraken2 --db /scratch/user/sal/kraken_projects/kraken_db_Mh --paired <read1> <read2> --report s1.report --output s1.output

#Replicate the step above with all the paired reads

Step 5; Once the output was analyzed next step is to "re-extract" the "new-fastq-files" using the reports and output files we got from step 4. We used Kraken tools for this step

module load WebProxy 
module load Biopython/
python /scratch/group/morley_group/bin/KrakenTools/extract_kraken_reads.py -k s1.output --report s1.report --taxid 75984 --include-children --include-parents --fastq-output -s1 <read1> -s2 <read2> -o extracted-r1.fastq -o2 extracted-r2.fastq 

Step 6; These fastq files off of the kraken-extraction tools were run against the kraken_Mh_Database made on step 2. This was done once we were sure from the analysis against BRD database that most of the sequence were pointing towards the Mannheimia. 

module load Kraken2/2.0.8-beta-foss-2018b-Perl-5.28.0 
kraken2 --db /scratch/user/sal/kraken_projects/kraken_db_Mh --paired <read1> <read2> --report extracted_s3.report --output extracted_s3.output

Step 7; We repeat step 5 but this time excluding --include-parents   

Step 8; Using Humann to see the functional profile for the fastq files generated using Kraken Tools from step 5 or 7. 

Note; Make sure to load all the necessary module listed below. The database necessary to run Humann has already been pre-downloaded in group scratch space in /scratch/group/vero_research/chocophlan_databases/chocophlan/ & --protein-database /scratch/group/vero_research/humann2_databases/uniref/

module load GCC/11.3.0 OpenMPI/4.1.4 MetaPhlAn/4.0.6-with_dbs
module load humann/3.6-mpa_vOct22
humann3 --input /scratch/group/vero_research/Mh_test_reads/s3_Mh_classified_r1.fastq --input <extracted fastq from step 5 or 7> --output s3_output --nucleotide-database /scratch/group/vero_research/chocophlan_databases/chocophlan/ --protein-database /scratch/group/vero_research/humann2_databases/uniref/

