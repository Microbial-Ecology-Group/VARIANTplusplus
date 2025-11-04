#!/usr/bin/env python3
import sys
import gzip
import os
import argparse

"""
rename_headers.py

- Autodetect FASTA vs. FASTQ.
- If FASTA: same as before (each header line => new unique ID).
- If FASTQ: lines where line_index % 4 == 0 => read header lines.
  We'll parse out the base name ignoring trailing /1 or /2 so that
  both forward and reverse read get the SAME unique integer ID.

Example:
  Original:
    @ABC|mySample|ReadXYZ/1
  We parse parted = ["ABC","mySample","ReadXYZ/1"] -> baseName="ReadXYZ"
  dictionary[ "ReadXYZ" ] => assigned integer, e.g. 42
  => rename => @42|mySample|ReadXYZ/1

Thus, both ReadXYZ/1 and ReadXYZ/2 => "42|mySample|...".
"""

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input",     required=True, help="Input FASTA/FASTQ, possibly gzipped.")
    p.add_argument("--sample_id", required=True, help="Sample ID to embed in headers.")
    p.add_argument("--output",    required=True, help="Output file name, possibly .gz as well.")
    return p.parse_args()

def detect_format(file_path):
    """
    Reads first non-blank line (in text mode).
    If it starts with '>', => FASTA
    If it starts with '@', => FASTQ
    else => None
    """
    if file_path.endswith(".gz"):
        fin = gzip.open(file_path, "rt")
    else:
        fin = open(file_path, "r")

    first_line=None
    for line in fin:
        line=line.strip()
        if line:
            first_line=line
            break
    fin.close()

    if not first_line:
        return None
    if first_line.startswith(">"):
        return "FASTA"
    elif first_line.startswith("@"):
        return "FASTQ"
    else:
        return None

def rename_fasta(in_file, sample_id, out_file):
    """
    FASTA rename:
      >oldHeader => >uniqueCounter|sample_id|oldHeader
    No paired logic here, just each new header => increment uniqueCounter.
    """
    if in_file.endswith(".gz"):
        fin = gzip.open(in_file,"rt")
    else:
        fin = open(in_file,"r")

    if out_file.endswith(".gz"):
        fout = gzip.open(out_file,"wt")
    else:
        fout = open(out_file,"w")

    unique_counter = 1
    for line in fin:
        if line.startswith(">"):
            old_header = line[1:].strip()
            new_header = f">{unique_counter}|{sample_id}|{old_header}"
            fout.write(new_header+"\n")
            unique_counter += 1
        else:
            fout.write(line)
    fin.close()
    fout.close()

def rename_fastq_paired(in_file, sample_id, out_file):
    """
    FASTQ rename: line_index % 4 == 0 => '@someHeader'.
    We'll parse parted => [uniqueInt, sampleID, baseNameAndSuffix],
    then extract the baseName ignoring trailing /1 or /2,
    so that /1 and /2 share the same integer from a dictionary.
    """
    if in_file.endswith(".gz"):
        fin = gzip.open(in_file,"rt")
    else:
        fin = open(in_file,"r")

    if out_file.endswith(".gz"):
        fout = gzip.open(out_file,"wt")
    else:
        fout = open(out_file,"w")

    # baseName -> assigned integer
    assigned_ids = {}
    next_id = 1

    line_index=0
    for line in fin:
        line=line.rstrip("\n")

        if line_index % 4 == 0 and line.startswith("@"):
            # parse parted => "uniqueInt|sampleID|stuff"
            parted = line[1:].split("|",2)
            if len(parted)<2:
                # can't parse => do default
                new_header=f"@{next_id}|{sample_id}|{line[1:]}"
                fout.write(new_header+"\n")
                next_id+=1
            else:
                # parted => [someUnique, thisSample, theRest]
                # theRest might end with /1 or /2 => let's strip that out to find baseName
                uniquePart, theirSampleID = parted[0], parted[1]
                theRest = parted[2] if len(parted)==3 else ""

                # If theirSampleID != sample_id, you might keep or override. We'll ignore that for now.
                # Next: parse baseName ignoring trailing /1 or /2
                baseName=theRest
                readSuffix=""
                if baseName.endswith("/1"):
                    baseName= baseName[:-2]  # remove /1
                    readSuffix="/1"
                elif baseName.endswith("/2"):
                    baseName= baseName[:-2]
                    readSuffix="/2"

                # check dictionary
                if baseName not in assigned_ids:
                    assigned_ids[baseName] = next_id
                    next_id +=1
                assignedInt = assigned_ids[baseName]

                # build new header => e.g. "@42|sample_id|BaseName/1"
                new_header= f"@{assignedInt}|{sample_id}|{baseName}{readSuffix}"
                fout.write(new_header+"\n")

        else:
            fout.write(line+"\n")

        line_index+=1
    fin.close()
    fout.close()

def main():
    args= parse_args()
    in_file = args.input
    out_file= args.output
    sample_id = args.sample_id

    fmt = detect_format(in_file)
    if not fmt:
        print(f"[ERROR] Could not detect FASTA or FASTQ from '{in_file}'. Empty or unknown.")
        sys.exit(1)

    if fmt=="FASTA":
        rename_fasta(in_file, sample_id, out_file)
        print(f"[INFO] Renamed FASTA headers => {out_file}")
    elif fmt=="FASTQ":
        rename_fastq_paired(in_file, sample_id, out_file)
        print(f"[INFO] Renamed FASTQ headers => {out_file}, pairing /1 and /2 with same ID.")
    else:
        print(f"[ERROR] Unexpected format => {fmt}")
        sys.exit(1)

if __name__=="__main__":
    main()