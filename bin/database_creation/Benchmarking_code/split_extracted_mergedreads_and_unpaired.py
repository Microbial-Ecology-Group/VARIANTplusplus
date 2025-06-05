#!/usr/bin/env python3
import os
import sys
import gzip
import argparse

"""
split_extracted_mergedreads_and_unpaired.py

Similar to split_extracted_reads_and_pair.py, but:
- We separate by sample_id (from second field in the header).
- We define "merged" vs "unmerged":
  * If parted[2] ends with "/1" or "/2" => unmerged
  * else => merged
- We write 2 files per sample:
  * {output_merged_dir}/{sample_id}_merged.fastq.gz
  * {output_unmerged_dir}/{sample_id}_unmerged.fastq.gz
- Summaries => <prefix>_mergedreads_summary.tsv in output_merged_dir,
  one row per sample with columns:
    sample_id, total_reads, merged_count, forward_count, reverse_count
- If input is FASTA => we store them in the "unmerged" directory one file per sample
  (keeping the old style of "FASTA => each sample => unpaired/unmerged"). 
  Summaries treat all of them as "unmerged" or you can treat them all as "merged"
  depending on your preference. The script below treats them as "unmerged."
"""

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input_fasta", required=True,
                   help="Path to input FASTA or FASTQ (gz or not). We'll autodetect format.")
    p.add_argument("--output_merged_dir", required=True,
                   help="Directory where merged reads are stored per sample if FASTQ.")
    p.add_argument("--output_unmerged_dir", required=True,
                   help="Directory where unmerged reads (or all FASTA) are stored.")
    return p.parse_args()

def open_text_read(fname):
    """Opens a file (possibly .gz) in text mode."""
    if fname.endswith(".gz"):
        return gzip.open(fname, "rt")
    else:
        return open(fname, "r")

def open_text_write(fname):
    """Opens a file (possibly .gz) in write text mode."""
    if fname.endswith(".gz"):
        return gzip.open(fname, "wt")
    else:
        return open(fname, "w")

def detect_format(file_path):
    """
    Reads the first non-blank line => '>' => FASTA, '@' => FASTQ, or None if uncertain.
    """
    with open_text_read(file_path) as fin:
        for line in fin:
            line=line.strip()
            if line:
                if line.startswith(">"):
                    return "FASTA"
                elif line.startswith("@"):
                    return "FASTQ"
                else:
                    return None
    return None

def strip_file_prefix(filename):
    """Removes .gz, .fastq, .fasta from the end to create a prefix."""
    prefix = os.path.basename(filename)
    if prefix.endswith(".gz"):
        prefix = prefix[:-3]
    if prefix.endswith(".fastq"):
        prefix = prefix[:-6]
    elif prefix.endswith(".fasta"):
        prefix = prefix[:-6]
    return prefix

def split_fasta(in_file, unmerged_dir):
    """
    For FASTA => old style: each sample => {unmerged_dir}/{sample_id}.fasta.gz
    We'll treat all as "unmerged." Summaries => forward=0, reverse=0, merged=0, total=someCount, etc.
    """
    os.makedirs(unmerged_dir, exist_ok=True)

    handles = {}  # sample_id => file handle
    # to track stats => { sample_id : [countOfSequences] }
    sample_counts = {}

    def get_handle(sample_id):
        if sample_id not in handles:
            outpath = os.path.join(unmerged_dir, f"{sample_id}.fasta.gz")
            fh = gzip.open(outpath, "wt")
            handles[sample_id] = fh
            sample_counts[sample_id] = 0
        return handles[sample_id]

    current_sample_id = None
    with open_text_read(in_file) as fin:
        for line in fin:
            if line.startswith(">"):
                # parse
                header = line[1:].strip()
                parts = header.split("|", 2)
                if len(parts) != 3:
                    print(f"[WARN:FASTA] malformed header => {header}", file=sys.stderr)
                    current_sample_id = None
                    continue
                _, sample_id, old_header = parts
                current_sample_id = sample_id
                out_f = get_handle(sample_id)
                out_f.write(f">{_}|{old_header}\n")
                sample_counts[sample_id] += 1
            else:
                if current_sample_id is not None:
                    out_f = get_handle(current_sample_id)
                    out_f.write(line)

    for fh in handles.values():
        fh.close()

    return sample_counts  # sample_id => count

def split_fastq_merged_unmerged(in_file, merged_dir, unmerged_dir):
    """
    FASTQ => We parse each read's 4-line chunk. parted => parted[0], parted[1], parted[2].
      parted[1] => sample_id
      parted[2] => if ends with /1 => forward => unmerged, /2 => reverse => unmerged, else => merged.
    We'll store them in memory for each sample => merged_blocks, unmerged_blocks.
    Then we'll write 2 files per sample => sample_id_merged.fastq.gz, sample_id_unmerged.fastq.gz
    We'll also track how many are forward or reverse or merged per sample for the summary.
    """
    os.makedirs(merged_dir, exist_ok=True)
    os.makedirs(unmerged_dir, exist_ok=True)

    # sample_data[sample_id] => {
    #   "merged_blocks": [],
    #   "unmerged_blocks": [],
    #   "merged_count":0,
    #   "forward_count":0,
    #   "reverse_count":0,
    #   "total":0
    # }
    sample_data = {}

    def ensure_sample_struct(sid):
        if sid not in sample_data:
            sample_data[sid] = {
                "merged_blocks": [],
                "unmerged_blocks": [],
                "merged_count": 0,
                "forward_count":0,
                "reverse_count":0,
                "total": 0
            }
        return sample_data[sid]

    with open_text_read(in_file) as fin:
        block = []
        for line in fin:
            block.append(line)
            if len(block)==4:
                header_line = block[0].rstrip("\n")
                if not header_line.startswith("@"):
                    # we treat it as unmerged forward, but no sample_id? We'll store in "unknown" sample?
                    # or skip it. We'll do "unknown" sample to avoid crashing
                    sample_id = "unknown"
                    parted = []
                    read_type = "forward"
                else:
                    parted = header_line[1:].split("|",2)
                    if len(parted)<2:
                        sample_id = "unknown"
                        read_type = "forward"
                    else:
                        sample_id = parted[1]
                    # parted[2] => if endswith /1 => forward, /2 => reverse => unmerged, else => merged
                    rest = parted[2] if len(parted)==3 else ""
                    if rest.endswith("/1"):
                        read_type="forward"
                    elif rest.endswith("/2"):
                        read_type="reverse"
                    else:
                        read_type="merged"
                st = ensure_sample_struct(sample_id)
                st["total"] +=1

                if read_type=="merged":
                    st["merged_blocks"].append(block[:])
                    st["merged_count"]+=1
                else:
                    # unmerged
                    st["unmerged_blocks"].append(block[:])
                    if read_type=="forward":
                        st["forward_count"]+=1
                    else:
                        st["reverse_count"]+=1
                block=[]

    # Now we write out 2 files per sample
    # sample_id_merged.fastq.gz => merged_dir
    # sample_id_unmerged.fastq.gz => unmerged_dir
    for sid, data in sample_data.items():
        if not sid:
            sid = "unknown"
        sid = sid.replace('.fastq.gz','')
        print(sid)
        merged_out  = os.path.join(merged_dir,   f"{sid}_merged.fastq.gz")
        unmerged_out= os.path.join(unmerged_dir, f"{sid}_unmerged.fastq.gz")
        fh_m = gzip.open(merged_out,   "wt")
        fh_u = gzip.open(unmerged_out, "wt")

        for blk in data["merged_blocks"]:
            fh_m.write("".join(blk))
        for blk in data["unmerged_blocks"]:
            fh_u.write("".join(blk))

        fh_m.close()
        fh_u.close()

    return sample_data

def main():
    parser = argparse.ArgumentParser(
        description="Split extracted reads into merged vs unmerged, per sample, then produce summary."
    )
    parser.add_argument("--input_fasta", required=True,
                        help="Path to input (FASTA or FASTQ).")
    parser.add_argument("--output_merged_dir", required=True,
                        help="Directory for storing merged reads per sample.")
    parser.add_argument("--output_unmerged_dir", required=True,
                        help="Directory for storing unmerged reads per sample.")
    args = parser.parse_args()

    # figure out prefix from input file
    prefix = strip_file_prefix(args.input_fasta)
    # summary => {prefix}_mergedreads_summary.tsv in output_merged_dir
    summary_path = os.path.join(args.output_merged_dir, f"{prefix}_mergedreads_summary.tsv")

    # detect format
    ftype = detect_format(args.input_fasta)
    if not ftype:
        print(f"[ERROR] can't detect format from", args.input_fasta, file=sys.stderr)
        sys.exit(1)

    summary_rows = []  # one row per sample

    if ftype=="FASTA":
        print(f"[INFO] Detected FASTA => each sample => unmerged (like original).", file=sys.stderr)
        sample_counts = split_fasta(args.input_fasta, args.output_unmerged_dir)
        # We produce one row per sample => sample_id, total_reads, merged_count=0, forward_count=0, reverse_count=0 or do we treat them all as unmerged?
        # We'll treat them as all unmerged => merged=0
        for sid, countVal in sample_counts.items():
            row = [
                sid,                # sample_id
                str(countVal),      # total_reads
                "0",                # merged_count
                "0",                # forward_count
                "0"                 # reverse_count
            ]
            summary_rows.append(row)

    else:
        print(f"[INFO] Detected FASTQ => merged vs unmerged logic, by parted[2] endswith /1 or /2", file=sys.stderr)
        sample_data = split_fastq_merged_unmerged(args.input_fasta, args.output_merged_dir, args.output_unmerged_dir)
        # sample_data[sample_id] => total, merged_count, forward_count, reverse_count
        for sid, data in sample_data.items():
            row = [
                sid,
                str(data["total"]),
                str(data["merged_count"]),
                str(data["forward_count"]),
                str(data["reverse_count"])
            ]
            summary_rows.append(row)

    # Write summary => columns: sample_id, total_reads, merged_count, forward_count, reverse_count
    # place in output_merged_dir
    with open(summary_path,"w") as outf:
        outf.write("\t".join(["sample_id","total_reads","merged_count","forward_count","reverse_count"])+"\n")
        for row in summary_rows:
            outf.write("\t".join(row)+"\n")

    print(f"[INFO] Wrote summary => {summary_path}", file=sys.stderr)

if __name__=="__main__":
    main()
