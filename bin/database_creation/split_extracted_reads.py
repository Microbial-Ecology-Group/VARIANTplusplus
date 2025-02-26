#!/usr/bin/env python3
import sys
import os
import gzip

def split_extracted_reads(input_fasta, output_dir):
    """
    Split each read by the 'sample_id' portion in headers:
      >uniqueInt|sample_id|oldHeader
    Writes to output_dir/{sample_id}.fasta.gz
    """
    if input_fasta.endswith(".gz"):
        fin = gzip.open(input_fasta, "rt")
    else:
        fin = open(input_fasta, "r")

    file_handles = {}  # sample_id -> gz handle

    def get_handle(sample_id):
        if sample_id not in file_handles:
            out_path = os.path.join(output_dir, f"{sample_id}")
            fhandle = gzip.open(out_path, "wt")
            file_handles[sample_id] = fhandle
        return file_handles[sample_id]

    current_sample_id = None
    for line in fin:
        if line.startswith(">"):
            header = line[1:].strip()  # e.g. "45|100XX0_PSVsXXoff_targetXX20000XX.fasta.gz|M01835:42:..."
            parts = header.split("|", 2)  # up to 3 parts
            if len(parts) != 3:
                print(f"WARNING: malformed header: {header}")
                current_sample_id = None
                continue
            unique_id, sample_id, old_header = parts
            current_sample_id = sample_id
            out_f = get_handle(sample_id)
            out_f.write(f">{unique_id}|{old_header}\n")
        else:
            if current_sample_id is not None:
                out_f = get_handle(current_sample_id)
                out_f.write(line)

    fin.close()
    for f in file_handles.values():
        f.close()

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_fasta", required=True)
    parser.add_argument("--output_dir", required=True)
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    split_extracted_reads(args.input_fasta, args.output_dir)

if __name__ == "__main__":
    main()
