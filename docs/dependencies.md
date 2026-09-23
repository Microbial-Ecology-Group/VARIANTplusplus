Dependencies
------------

VARIANT++ uses a variety of open-source tools. The tools used and version specifics (where pinned) are provided below. These dependencies can be installed manually on your local computing cluster, or by changing the `-profile` parameter, most of them can be installed automatically via conda (`envs/VARIANT++_env.yaml`) — see [installation.md](installation.md).

### Trimmomatic
  - Description: removes low quality base pairs and adapter sequences from raw sequence data (GSV_1).
  - Version: 0.39
  - DOI: https://doi.org/10.1093/bioinformatics/btu170

### FastQC
  - Description: a quality control tool for high throughput sequence data (`--pipeline eval_qc`).
  - Version: 0.11.8

### multiQC
  - Description: aggregates FastQC reports across samples into a single report (`--pipeline eval_qc`).
  - DOI: https://doi.org/10.1093/bioinformatics/btw354

### FLASH
  - Description: merges overlapping paired-end reads into a single extended fragment (GSV_1).
  - DOI: https://doi.org/10.1093/bioinformatics/btr507

### BBMap (clumpify.sh)
  - Description: clumpify.sh deduplicates merged/unmerged reads (GSV_2).

### BWA
  - Description: aligns reads to a reference genome; used here to remove host DNA (GSV_3).
  - Version: 0.7.17
  - DOI: https://doi.org/10.1093/bioinformatics/btp324

### Samtools
  - Description: manipulates and extracts information from SAM/BAM alignment files (GSV_3).
  - Version: 1.15.1
  - DOI: https://doi.org/10.1093/bioinformatics/btp352

### Kraken2
  - Description: a fast taxonomic sequence classifier used to extract reads belonging to the target species (GSV_4).
  - Version: 2.1.2
  - DOI: https://doi.org/10.1186/gb-2014-15-3-r46

### seqkit
  - Description: used for read-count statistics at several GSV steps.
  - DOI: https://doi.org/10.1371/journal.pone.0163962

### Themisto
  - Description: pseudoaligns reads against a colored de Bruijn graph index of reference genomes; the compiled binary ships in `bin/themisto` (GSV_5).
  - DOI: https://doi.org/10.1093/bioinformatics/btad233

### mSWEEP
  - Description: estimates relative abundances of GSV clusters from Themisto pseudoalignments, and can bin reads per cluster (GSV_5, GSV_5_mGEMS).
  - DOI: https://doi.org/10.1371/journal.pcbi.1009041

### mGEMS
  - Description: extracts per-cluster read bins from mSWEEP's `--bin-reads` output (GSV_5_mGEMS).
  - DOI: https://doi.org/10.1093/bioinformatics/btab003

### pigz
  - Description: parallel gzip, used to compress host-removed FASTQ output (GSV_3).
