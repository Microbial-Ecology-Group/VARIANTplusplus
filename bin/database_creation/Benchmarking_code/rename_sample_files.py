#!/usr/bin/env python3
import os
import sys
import gzip
import argparse

"""
rename_sample_files.py

Usage Example:
  python rename_sample_files.py --input /path/to/1XX3_PSVsXXk_4XX1_2_3.fastq.gz --results-dir /path/to/results

This script:
  1) Parses a FASTQ or FASTQ.GZ file line by line.
  2) For every header line ('@...'), extracts the PSV from between the 1st and 2nd '|' chars (if present).
  3) Counts how many reads belong to each unique PSV.
  4) If no PSVs are found, label the sample "off_target" with the total read count.
  5) Builds a suffix like "PSV1.5000_PSV2.3000" (or "off_target.12345") sorted by PSV name.
  6) Writes a line to results-dir/counts_{basename}.txt:
       <originalFileName>\t<newFileName>
     where <newFileName> is the original file name but with everything after the last 'XX'
     replaced by that new suffix. The extension ".fastq.gz" is re-added.

No files are actually renamed. It's just a reference file mapping old to new names.
"""

def parse_args():
    p = argparse.ArgumentParser(description="Parse FASTQ headers to count PSVs, store old->new name in a text file.")
    p.add_argument("--input", required=True, help="Path to .fastq/.fastq.gz file to parse. We'll not rename; just record.")
    p.add_argument("--results-dir", required=True, help="Directory where we write counts_<basename>.txt")
    return p.parse_args()

def main():
    args = parse_args()
    in_file = args.input
    results_dir = args.results_dir

    # Ensure results_dir exists
    if not os.path.isdir(results_dir):
        print(f"[ERROR] The results-dir '{results_dir}' does not exist or is not a directory.")
        sys.exit(1)

    # Open FASTQ
    if in_file.endswith(".gz"):
        fin = gzip.open(in_file, "rt")
    else:
        fin = open(in_file, "r")

    psv_counts = {}
    total_reads = 0  # We'll track how many '@' (header) lines appear

    # For each header line starting with "@", parse the PSV from splitted[1].
    for line in fin:
        line=line.rstrip("\n")
        if line.startswith("@"):
            total_reads += 1
            splitted=line.split("|")
            if len(splitted) >= 3:
                psv_id = splitted[1]  # e.g. "PSV2"
                psv_counts[psv_id] = psv_counts.get(psv_id,0) + 1

    fin.close()

    # If we found no PSVs => treat as off_target w/ total_reads
    if not psv_counts:
        print(f"[INFO] No PSVs found in {in_file} => marking this as off_target with total read count = {total_reads}")
        psv_counts["off_target"] = total_reads

    # Build suffix: e.g. "PSV1.5000_PSV2.3000" or "off_target.12345"
    items = sorted(psv_counts.items(), key=lambda x: x[0]) 
    suffix_pieces = [f"{psv}.{count}" for (psv, count) in items]
    new_suffix = "_".join(suffix_pieces)

    # Original file name (without path)
    old_basename = os.path.basename(in_file)

    # Remove extension ".fastq.gz" or ".fastq"
    if old_basename.endswith(".fastq.gz"):
        core = old_basename[:-9]  # remove .fastq.gz
        ext  = ".fastq.gz"
    elif old_basename.endswith(".fastq"):
        core = old_basename[:-6]
        ext  = ".fastq"
    else:
        core = old_basename
        ext  = ""

    # Insert new suffix after last 'XX'
    idx = core.rfind("XX")
    if idx == -1:
        # no 'XX', just append
        new_core = f"{core}.{new_suffix}"
    else:
        prefix = core[:idx]  # up to 'XX'
        new_core = f"{prefix}XX{new_suffix}"

    new_name = new_core + ext  # final new filename

    # Write the row to results-dir/counts_<basename>.txt
    # We'll store: oldBasename <tab> newName
    outfilename = os.path.join(results_dir, f"counts_{old_basename}.txt")

    with open(outfilename, "w") as outf:
        outf.write(f"{old_basename}\t{new_name}\n")

    print(f"[INFO] Wrote PSV counts-based new name to: {outfilename}")
    print(f"[INFO] Original => {old_basename}")
    print(f"[INFO] Proposed => {new_name}")

if __name__=="__main__":
    main()
