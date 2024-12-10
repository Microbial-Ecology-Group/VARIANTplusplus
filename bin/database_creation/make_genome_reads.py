#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import gzip

# Function to extract segments from the genome
def extract_segments(genome_file, segment_length, output_file):
    # Determine if the file is gzipped
    open_func = gzip.open if genome_file.endswith('.gz') else open
    print(f"Subsetting {genome_file}")

    with open_func(genome_file, 'rt') as handle:
        with open(output_file, 'w') as out_handle:
            for record in SeqIO.parse(handle, 'fasta'):
                genome_length = len(record.seq)
                # Skip fragments smaller than the segment length
                if genome_length < segment_length:
                    print(f"Skipping {record.id} due to insufficient length: {genome_length} bp")
                    continue
                
                num_segments = genome_length // segment_length
                # Calculate intervals for even distribution
                interval = genome_length / num_segments

                for i in range(num_segments):
                    start = int(i * interval)
                    end = start + segment_length
                    segment = record.seq[start:end]
                    
                    # Create a unique segment identifier based on sequence name and position
                    segment_id = f"{record.id}_Segment_{i+1}"
                    out_handle.write(f">{segment_id}\n{segment}\n")

    print(f"Done with {genome_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract evenly distributed segments from a genome.")
    parser.add_argument("genome_file", help="Path to the genome FASTA file (can be gzipped).")
    parser.add_argument("segment_length", type=int, help="Length of the segments to be extracted.")
    parser.add_argument("output_file", help="Path to the output file where segments will be written.")
    args = parser.parse_args()

    extract_segments(args.genome_file, args.segment_length, args.output_file)
