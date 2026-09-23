Output
------

All intermediate outputs produced from each module of this pipeline are provided as flat files that can be viewed in a text editor. These files are copied from the root **work/** directory created by Nextflow, so if disk space is a concern, this directory should be deleted once you have your results, as it can get quite large.

Directory Structure
-------------------

The output directories created by the pipeline are named after the GSV step that produced them, and are populated incrementally as you run `--pipeline GSV_1` through `--pipeline GSV_5` (or `GSV_5_mGEMS`), or all at once with `--pipeline full_GSV_pipeline`. Each file output is prefixed with the sample name.

Below is an example of the results created by a full run against the bundled test samples (`S1_test`, `S2_test`, `S3_test`), including the optional `--pipeline eval_qc` step.

```bash
test_results
├── QC_analysis                                  # --pipeline eval_qc
│   ├── FastQC
│   │   ├── S1_test_fastqc_logs
│   │   │   ├── S1_test_R1_fastqc.html
│   │   │   └── S1_test_R1_fastqc.zip
│   │   ├── S2_test_fastqc_logs
│   │   └── S3_test_fastqc_logs
│   └── MultiQC_stats
│       ├── multiqc_report.html
│       ├── multiqc_general_stats.txt
│       └── multiqc_data
├── QC_trimming                                  # GSV_1 (Trimmomatic)
│   ├── Paired
│   │   ├── S1_test.1P.fastq.gz
│   │   └── S1_test.2P.fastq.gz
│   └── Unpaired
│       ├── S1_test.1U.fastq.gz
│       └── S1_test.2U.fastq.gz
├── Flash_reads                                  # GSV_1 (FLASH merge)
│   ├── S1_test.extendedFrags.fastq.gz
│   ├── S1_test.notCombined.fastq.gz
│   ├── S1_test.hist
│   └── S1_test.log
├── Deduped_reads                                # GSV_2 (BBMap clumpify)
│   ├── S1_test_merged.dedup.fastq.gz
│   ├── S1_test_unmerged.dedup.fastq.gz
│   └── S1_test.dedupe_clumpify.stats.log
├── HostRemoval                                  # GSV_3 (BWA)
│   └── NonHostFastq
│       ├── S1_test.merged.non.host.fastq.gz
│       └── S1_test.unmerged.non.host.fastq.gz
├── MicrobiomeAnalysis                           # GSV_4 (Kraken2 extraction)
│   └── Kraken
│       ├── standard
│       │   ├── S1_test.merged.kraken.raw
│       │   └── S1_test.unmerged.kraken.raw
│       ├── standard_report
│       │   ├── S1_test.merged.kraken.report
│       │   └── S1_test.unmerged.kraken.report
│       └── extracted_reads
│           ├── S1_test_extracted_merged.fastq.gz
│           └── S1_test_extracted_unmerged.fastq.gz
├── Filtered_pseudoaligned_reads                 # GSV_5 (Themisto)
│   ├── S1_test_pseudoaligned_merged.fastq.gz
│   └── S1_test_pseudoaligned_unmerged.fastq.gz
├── mSWEEP_results                               # GSV_5 (mSWEEP)
│   ├── S1_test.merged.msweep_abundances.txt
│   └── S1_test.unmerged.msweep_abundances.txt
├── GSV_binned_reads                             # GSV_5_mGEMS only
│   └── S1_test_GSV_1.fastq.gz
└── Results
    ├── Stats
    │   ├── trimmomatic.stats
    │   └── Raw_reads.txt
    ├── kraken_analytic_matrix.csv               # GSV_4 (Kraken2 per-sample counts)
    ├── mSweep_results_summary.tsv               # GSV_5 (per-sample, per-read-type GSV calls)
    └── mSweep_results_count_matrix.tsv          # GSV_5 (combined GSV count matrix)
```

Check the `Results` folder for the Kraken2 analytic matrix and the mSWEEP results. `mSweep_results_summary.tsv` contains results for both the merged and unmerged reads separately, while `mSweep_results_count_matrix.tsv` has the combined count matrix you can load into R for downstream analysis (alpha/beta diversity, etc.) — see the [step-by-step guide](VARIANT++_step_by_step.md#explore-the-results).
