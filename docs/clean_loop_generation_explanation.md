# Sample Generation Loop Structure

## What's being iterated over

```
for numGSV in num_GSV_list:            # e.g., [0, 1, 3, 5, 7, 9]
    for iteration in range(1, num_iters+1):  # e.g., 1..100
        
        if numGSV == 0:
            → 1 off-target partial sample (no k_col loop)
        
        else:
            for k_col in k_columns:    # e.g., [k_3, k_5, k_7, k_9]
                → 1 on-target partial sample per k_col
        
        → all partials concatenated into 1 iteration-level merged FASTQ
```

## Concrete example

Using these parameters:

| Parameter | Value |
|-----------|-------|
| `num_GSV_list` | `[0, 1, 3, 5, 7, 9]` (6 levels) |
| `k_columns` | `[k_3, k_5, k_7, k_9]` (4 columns) |
| `num_iters` | `100` |
| `kraken_confidence` | `[0, 0.1, 0.2, 0.3, 0.4, 0.5]` (6 levels) |

### Partial samples per iteration (iter = 1)

| numGSV | k_col loop? | Partial samples created | Partial sample names |
|--------|-------------|------------------------|----------------------|
| 0 | No | **1** | `1XX0_GSVsXXoff_target` |
| 1 | Yes (4 k_cols) | **4** | `1XX1_GSVsXXk_3XX{gsv}`, `...k_5...`, `...k_7...`, `...k_9...` |
| 3 | Yes (4 k_cols) | **4** | `1XX3_GSVsXXk_3XX{gsv_gsv_gsv}`, etc. |
| 5 | Yes (4 k_cols) | **4** | `1XX5_GSVsXXk_3XX{...}`, etc. |
| 7 | Yes (4 k_cols) | **4** | `1XX7_GSVsXXk_3XX{...}`, etc. |
| 9 | Yes (4 k_cols) | **4** | `1XX9_GSVsXXk_3XX{...}`, etc. |
| **Total** | | **21** | |

**Per iteration: 1 + (5 on-target numGSV levels × 4 k_cols) = 21 unique partial samples**

> Note: if numGSV > number of GSVs available in a k_col (e.g., numGSV=9 but k_3 only has 3 groups), the script still creates the sample — it just caps at `min(numGSV, available)`. The sample is created with fewer GSVs than requested.

### Iteration-level merged FASTQs per iteration

Each (numGSV, iteration) pair produces one concatenated FASTQ containing all its partial samples:

| numGSV | Partials merged | Merged FASTQ name |
|--------|----------------|-------------------|
| 0 | 1 partial | `GSV_0_iter_1_ALL_renamed.fastq.gz` |
| 1 | 4 partials | `GSV_1_iter_1_ALL_renamed.fastq.gz` |
| 3 | 4 partials | `GSV_3_iter_1_ALL_renamed.fastq.gz` |
| 5 | 4 partials | `GSV_5_iter_1_ALL_renamed.fastq.gz` |
| 7 | 4 partials | `GSV_7_iter_1_ALL_renamed.fastq.gz` |
| 9 | 4 partials | `GSV_9_iter_1_ALL_renamed.fastq.gz` |
| **Total** | | **6 merged FASTQs per iteration** |

These merged FASTQs are what Kraken2 classifies. After extraction, the split script de-multiplexes them back into the 21 individual partial samples using the embedded read headers.

## Grand totals across all iterations

| Count | Formula | Value |
|-------|---------|-------|
| Unique partial samples | `num_iters × [1 + (on-target numGSV levels × k_cols)]` | 100 × 21 = **2,100** |
| Iteration-level merged FASTQs | `num_iters × len(num_GSV_list)` | 100 × 6 = **600** |
| Kraken2 runs | `merged FASTQs × len(kraken_confidence)` | 600 × 6 = **3,600** |
| mSWEEP runs (merged + unmerged) | see below | **28,800** |

### mSWEEP run breakdown

Each partial sample is run through Themisto once per confidence, then mSWEEP twice (merged + unmerged reads):

| Sample type | Partials per iter | Annotations classified against | mSWEEP calls per partial per conf | Per iter per conf |
|-------------|-------------------|-------------------------------|----------------------------------|-------------------|
| Off-target | 1 | ALL 4 k_cols | 4 × 2 = 8 | 8 |
| On-target | 20 | Matching k_col only | 1 × 2 = 2 | 40 |
| **Total per iter per conf** | | | | **48** |

- Themisto runs: 2,100 partials × 6 confs × 2 (merged + unmerged) = **25,200**
- mSWEEP runs: (20 on-target × 2 + 1 off-target × 4 × 2) × 100 iters × 6 confs = (40 + 8) × 600 = **28,800**

## General formula

```
N_k       = len(k_columns)
N_gsv     = len(num_GSV_list)          # includes 0
N_on      = N_gsv - 1                  # on-target numGSV levels (excluding 0)
N_iter    = num_iters
N_conf    = len(kraken_confidence)

partials_per_iter    = 1 + (N_on × N_k)
merged_fastqs        = N_gsv × N_iter
kraken_runs          = merged_fastqs × N_conf
themisto_runs        = partials_per_iter × N_iter × N_conf × 2
msweep_runs          = N_iter × N_conf × [(N_on × N_k × 2) + (1 × N_k × 2)]
                     = N_iter × N_conf × 2 × N_k × (N_on + 1)
                     = N_iter × N_conf × 2 × N_k × N_gsv
```