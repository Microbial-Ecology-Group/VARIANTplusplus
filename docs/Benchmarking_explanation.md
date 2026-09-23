# `generate_simulated_samples_byConf.py` – Detailed Guide
_Automated generator for GSV-benchmarking bash scripts_

---

## 1. Purpose

Genomic Sequence Variant (**GSV**) benchmarking evaluates how well different classifier settings distinguish closely related genomes.  
This script clusters reference genomes from a target species into GSVs, simulates metagenomic samples from those references, and tests multiple Kraken2 and mSWEEP configurations **to maximise recall and precision**.  
Off-target genomes from other species in the same family are processed in parallel to measure **false-positive rates**, allowing the pipeline to pinpoint the best combination of:

- Number of GSV clusters (annotation/classification resolution)
- Kraken2 classifier confidence thresholds  
- Count filtration to account for false positive classification

### Terminology

| Term | Meaning |
|------|---------|
| **classification_cluster** | The `k_col` annotation used by mSWEEP to classify reads (determines how many groups mSWEEP partitions into) |
| **On-target** | Simulated reads from genomes belonging to the target species (`numGSV > 0`). For these samples, reads are both simulated from AND classified against the same `k_col`. |
| **Off-target** | Simulated reads from genomes of non-target species in the same family (`numGSV == 0`). One sample per iteration, classified against ALL annotation clusters. |

---

## 2. Quick Usage

Follow the GSV_benchmarking_tutorial.md document for detailed instructions.
In short, you prepare a bunch of files and databases. Then you run the benchmarking script to take that information and create a series of 3 scripts which will make simulated samples and take them through our VARIANT++ classification workflow.

```bash
python generate_simulated_samples_byConf.py --benchmarking_params params.txt
```

The random seed is fixed (`random.seed(2)`) for reproducibility.

---

## 3. What the Script Generates

| Bash Script                           | Stage                             | Key Tasks                                                  |
|---------------------------------------|-----------------------------------|------------------------------------------------------------|
| **script_build_reads_\<PREFIX\>.sh**    | Sample generation                 | • Concatenate target / off-target genomes<br>• Simulate reads with ISS (NovaSeq error model)<br>• Merge R1/R2 with FLASH<br>• Standardise filenames & headers with `rename_sample_files.py` and `rename_headers.py`<br>• Concatenate partial samples per iteration into one merged FASTQ |
| **script_kraken_extract_split_\<PREFIX\>.sh** | Taxonomic classification           | • Classify merged FASTQs at each `--confidence`<br>• Extract target-taxon reads with `extract_kraken_reads.py`<br>• Split back into individual samples and separate merged vs unmerged reads with `split_extracted_mergedreads_and_unpaired.py` |
| **script_themisto_msweep_\<PREFIX\>.sh** | GSV quantification                | • Pseudo-align extracted reads with Themisto (per individual sample)<br>• Run mSWEEP separately for merged & unmerged reads to estimate GSV abundances |

---

## 4. Workflow Breakdown

### Step 1: Parse configuration
Reads `params.txt` (key=value format) with paths, iteration counts, GSV list, Kraken confidences, etc. This file can be found in `VARIANT++/bin/database_creation/Benchmarking_code`. List values are parsed using `ast.literal_eval()`.

### Step 2: Create directories
```
<PREFIX>_cat_reads/     # concatenated iteration-level FASTQs for Kraken2
<PREFIX>_split_reads/   # individual sample unmerged reads (post-Kraken extraction)
<PREFIX>_merged_reads/  # individual sample merged reads (post-Kraken extraction)
<PREFIX>_results/       # mSWEEP annotation files, Themisto outputs, mSWEEP results
```

### Step 3: Generate mSWEEP annotation files
Creates one `k_#_msweep.txt` file per annotation column from the metadata TSV. Each file has one line per genome in the Themisto index, listing its GSV cluster ID for that resolution.

### Step 4: Write `script_build_reads` — Sample simulation

Loops over `num_GSV_list × num_iters`:

**Off-target samples (numGSV == 0):**
- ONE sample per iteration (no k_col loop — non-target genomes don't have GSV annotations)
- Randomly select `num_genomes_per_GSV` genomes from the non-target genome directory
- Concatenate into a single reference FASTA
- Simulate paired-end reads with ISS (NovaSeq error model, fragment length 311 ± 50 bp)
- Merge R1/R2 with FLASH (max overlap 120 bp), then concatenate merged and unmerged reads
- Run `rename_sample_files.py` to record sample metadata, then `rename_headers.py` to embed the sample identity into read headers
- The sample will be classified against ALL k_col annotations in the mSWEEP step

**On-target samples (numGSV > 0):**
- One sample per k_col per iteration
- For each `k_col`, randomly select `numGSV` distinct GSV groups from that annotation
- For each selected GSV, randomly pick `num_genomes_per_GSV` reference genomes
- Concatenate all selected genomes into a single reference FASTA
- Genome headers are renamed to include their GSV ID (using `rename_headers.py`) so ISS distributes reads across all GSVs
- Same ISS → FLASH → rename pipeline as off-target

**Merging per iteration:**
After all partial samples are created for a given (numGSV, iteration), they are concatenated into a single iteration-level FASTQ: `GSV_<N>_iter_<I>_ALL_renamed.fastq.gz`. This reduces the number of Kraken2 runs needed since all partial samples are classified together.

Read counts per partial sample are randomly chosen from `num_reads_options`.

The shared simulation pipeline (ISS → FLASH → rename) is handled by the `write_simulate_and_rename()` helper function to eliminate code duplication between on-target and off-target branches.

### Step 5: Write `script_kraken_extract_split` — Classification & extraction

For each (numGSV, iteration, confidence):
1. **Kraken2** classifies the concatenated iteration-level FASTQ at the given `--confidence` threshold
2. **extract_kraken_reads.py** extracts reads matching `extract_reads_taxid` (and optionally children taxa via `--include-children`)
3. **split_extracted_mergedreads_and_unpaired.py** de-multiplexes the extracted reads back into individual samples using the embedded header information, and separates FLASH-merged reads from unmerged reads

This produces per-sample, per-confidence FASTQ files in the `_merged_reads/` and `_split_reads/` directories.

### Step 6: Write `script_themisto_msweep` — Pseudo-alignment & abundance estimation

For each individual sample × confidence combination:
1. **Themisto** pseudo-aligns reads against the pre-built index (separately for merged and unmerged reads)
2. **mSWEEP** estimates GSV abundances using the Themisto output and the appropriate annotation file

**Key distinction:**
- **On-target samples**: classified only against their matching `k_col` annotation (the one used to select the GSVs during simulation)
- **Off-target samples**: classified against ALL `k_col` annotations, to test FP rates across every classification resolution

Both on-target and off-target samples use the same Themisto and mSWEEP commands (no `--themisto-mode` differences).

### Step 7: Finish
Print a summary of the three generated bash scripts.

---

## 5. Key Parameters (`params.txt`)

| Key                   | Type         | Example               | Description                                        |
|-----------------------|--------------|-----------------------|----------------------------------------------------|
| `target_genome_dir`   | str          | `/data/targets/`      | Directory of `.fna[.gz]` files for the target species |
| `nontarget_genome_dir`| str          | `/data/offtarget/`    | Directory of off-target family member genomes      |
| `input_file`          | str          | `ani_clusters.tsv`    | TSV with genome filenames + GSV annotation per `k_col` |
| `bin_dir`             | str          | `/usr/local/bin/`     | Path to helper scripts (`rename_sample_files.py`, `rename_headers.py`, `split_extracted_mergedreads_and_unpaired.py`) |
| `themisto_index`      | str          | `/indices/themisto/`  | Prefix of existing Themisto index                  |
| `krakendb`            | str          | `/databases/kraken2/` | Kraken2 database directory                         |
| `kraken_confidence`   | list / int   | `[0, 0.1, 0.5]`      | Confidence thresholds to test                      |
| `extract_reads_taxid` | str / int    | `75985`               | NCBI taxonomy ID for target taxon extraction       |
| `extract_reads_options` | str        | `--include-children`  | Additional flags for `extract_kraken_reads.py`     |
| `num_iters`           | int          | `100`                 | Number of replicate iterations per `numGSV` setting |
| `num_GSV_list`        | list[int]    | `[0, 1, 3, 5, 7, 9]` | GSV counts to test (0 = off-target)                |
| `num_reads_options`   | list[int]    | `[10000, 20000]`      | Possible read counts per partial sample (randomly chosen) |
| `num_genomes_per_GSV` | int          | `5`                   | Genomes sampled per GSV group (or from non-target pool) |
| `threads`             | int          | `48`                  | CPU threads for all tools                          |
| `output_name_prefix`  | str          | `Salmonella_bench`    | Prefix for all output directories and scripts      |
| `tmp_build`           | str          | `$TMPDIR`             | Scratch directory (can use environment variable)   |

---

## 6. Expected Outputs

```
<PREFIX>_cat_reads/      # concatenated iteration-level FASTQs for Kraken2
<PREFIX>_split_reads/    # per-sample unmerged FASTQs (by confidence)
<PREFIX>_merged_reads/   # per-sample merged FASTQs (by confidence)
<PREFIX>_results/        # mSWEEP annotation files, Themisto outputs, mSWEEP abundance tables
script_build_reads_<PREFIX>.sh
script_kraken_extract_split_<PREFIX>.sh
script_themisto_msweep_<PREFIX>.sh
```

---

## 7. Design Notes

### Off-target sample handling
Off-target samples (`numGSV == 0`) use non-target genomes from the same taxonomic family. Since these reads should not match any GSV, they serve as a negative control. One off-target sample is created per iteration (no k_col loop — non-target genomes don't have GSV annotations), and it is classified by mSWEEP against every annotation cluster to measure false positive rates at each resolution.

### On-target source = classification
For on-target samples, the k_col used to select the simulated GSVs is the same k_col used for mSWEEP classification. There is no separate "source_cluster" vs "classification_cluster" distinction for on-target data — the classification_cluster in the filename tells you both what was simulated and what was used for classification.

### Filename conventions
Partial sample names encode simulation parameters:
- On-target: `{iteration}XX{numGSV}_GSVsXX{k_col}XX{GSV_IDs}`
- Off-target: `{iteration}XX0_GSVsXXoff_target`

The `XX` delimiter separates fields for downstream parsing.

---

## 8. Dependencies

- Python ≥ 3.6  
- **ISS** (InSilicoSeq — read simulator)  
- **FLASH** (read merger)  
- **Kraken2** ≥ 2.1  
- **extract_kraken_reads.py** (KrakenTools)  
- **Themisto** ≥ 2.5  
- **mSWEEP** ≥ 1.9  
- **pigz** (parallel gzip)
- Helper scripts in `bin_dir`:
  - `rename_sample_files.py` — records sample metadata for downstream parsing
  - `rename_headers.py` — embeds sample identity into FASTQ read headers
  - `split_extracted_mergedreads_and_unpaired.py` — de-multiplexes extracted reads back into individual samples and separates merged from unmerged

---

## 9. Troubleshooting Tips

- **Empty output?** Verify `num_GSV_list` values exist in your TSV annotations.  
- **File-not-found errors**: Check all paths in `params.txt`; the script checks for both `.fna` and `.fna.gz` extensions.  
- **Kraken2 memory issues**: Add `--memory-mapping` to reduce RAM, but significantly increases runtime.
- **mSWEEP annotation errors**: Ensure each `k_#_msweep.txt` has exactly one entry per genome in the Themisto index, in the same order.
- **Reproducibility**: The random seed is hardcoded as `random.seed(2)`. Change this value or make it a parameter if you want different random samples.

---