#!/usr/bin/env python3
import sys
import gzip

def rename_headers(input_fasta, sample_id, output_fasta):
    """
    For each header:
      >oldHeader
    becomes
      >[uniqueInt]|[sample_id]|oldHeader
    """
    if input_fasta.endswith(".gz"):
        fin = gzip.open(input_fasta, "rt")
    else:
        fin = open(input_fasta, "r")

    if output_fasta.endswith(".gz"):
        fout = gzip.open(output_fasta, "wt")
    else:
        fout = open(output_fasta, "w")

    unique_counter = 1
    for line in fin:
        if line.startswith(">"):
            old = line[1:].strip()
            fout.write(f">{unique_counter}|{sample_id}|{old}\n")
            unique_counter += 1
        else:
            fout.write(line)

    fin.close()
    fout.close()

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rename_headers(args.input, args.sample_id, args.output)

if __name__ == "__main__":
    main()
