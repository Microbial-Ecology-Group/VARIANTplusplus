#!/usr/bin/env python3
"""
split_script_iteration_markers.py

Split a long bash script into N chunks at marker lines and wrap each
chunk in an SBATCH launcher.

Example
-------
python split_script_iteration_markers.py \
        --input_file bigScript.sh \
        --num_output_files 20 \
        --name_prefix myBatch \
        --mem 250G \
        --time 48:00:00 \
        --cpus 64
"""
import os
import argparse

# ──────────────── argument parsing ─────────────────────────────────────────
parser = argparse.ArgumentParser(
    description="Split a long bash script by marker lines & generate SBATCH wrappers."
)
parser.add_argument("-i", "--input_file", required=True,  help="Original long script")
parser.add_argument("-n", "--num_output_files", required=True, type=int,
                    help="Number of chunks/scripts to create")
parser.add_argument("-p", "--name_prefix", required=True,
                    help="Prefix for the generated files")

# optional SBATCH tuning
parser.add_argument("--mem",  default="40G",      help="SBATCH --mem value")
parser.add_argument("--time", default="24:00:00", help="SBATCH -t wall time")
parser.add_argument("--cpus", default=48, type=int, help="SBATCH --cpus-per-task value")
args = parser.parse_args()

# ──────────────── SBATCH template ─────────────────────────────────────────
SBATCH_TEMPLATE = (
    "#!/bin/bash\n"
    "#SBATCH -J {jobname}\n"
    "#SBATCH -o scripts/log_{jobname}.out\n"
    "#SBATCH -t {time}\n"
    "#SBATCH --mem={mem}\n"
    "#SBATCH --nodes=1\n"
    "#SBATCH --ntasks=1\n"
    "#SBATCH --cpus-per-task={cpus}\n"
    "#SBATCH --ntasks-per-node=1\n\n"
)

# directories
SCRIPTS_DIR = "scripts"
os.makedirs(SCRIPTS_DIR, exist_ok=True)

MEGASAMPLE_MARKER = "# - Starting megasample"
ITERATION_MARKER  = "# --- Starting iteration"

# ──────────────── main logic ──────────────────────────────────────────────
def split_by_iteration_markers(input_file: str, n_out: int, prefix: str):
    with open(input_file) as fh:
        lines = fh.readlines()

    total_lines = len(lines)
    total_megasample = sum(1 for l in lines if l.startswith(MEGASAMPLE_MARKER))
    total_iteration  = sum(1 for l in lines if l.startswith(ITERATION_MARKER))

    marker = (MEGASAMPLE_MARKER if total_megasample else ITERATION_MARKER)
    print(f"Using marker '{marker}'")

    # group between markers
    groups, curr, m_total, s_total = [], [], 0, 0
    for l in lines:
        if l.startswith(marker):
            m_total += 1
            if curr:
                groups.append(curr)
            curr = [l]
        else:
            curr.append(l)
        if l.startswith("# Unique Sample"):
            s_total += 1
    if curr:
        groups.append(curr)

    # round-robin distribution
    chunks = [[] for _ in range(n_out)]
    for idx, grp in enumerate(groups):
        chunks[idx % n_out].extend(grp)

    # write output scripts & SBATCH wrappers
    for idx, chunk in enumerate(chunks, start=1):
        cmd = f"{SCRIPTS_DIR}/{prefix}_part_{idx}.sh"
        sb  = f"{SCRIPTS_DIR}/{prefix}_part_{idx}.sbatch"
        job = f"{prefix}_{idx}"

        m_cnt = sum(1 for l in chunk if l.startswith(marker))
        s_cnt = sum(1 for l in chunk if l.startswith("# Unique Sample"))

        with open(cmd, "w") as cf:
            cf.write(
                f"echo \"[{prefix}] Script {idx}/{n_out}: "
                f"{m_cnt} markers of {m_total}; "
                f"{s_cnt} samples of {s_total}.\"\n\n"
            )
            s_idx = 0
            for l in chunk:
                if l.startswith("# Unique Sample"):
                    s_idx += 1
                    cf.write(f"echo \"Sample {s_idx}/{s_cnt} in part {idx}.\"\n")
                cf.write(l)

        with open(sb, "w") as sf:
            sf.write(SBATCH_TEMPLATE.format(
                jobname=job, mem=args.mem, time=args.time, cpus=args.cpus
            ))
            sf.write(f"echo 'Launching {cmd} (part {idx})'\n")
            sf.write(f"bash {cmd}\n")

        print(f"Created {cmd} & {sb} ({len(chunk)} lines)")

    # detailed summary
    print("\n=== Global Summary ===")
    print(f"Total lines in '{input_file}': {total_lines}")
    print(f"Megasample marker lines:       {total_megasample}")
    print(f"Iteration marker lines:        {total_iteration}")
    print(f"Processed marker lines ('{marker}'): {m_total}")
    print(f"Total '# Unique Sample' lines: {s_total}\n")

# ──────────────── entry point ─────────────────────────────────────────────
if __name__ == "__main__":
    split_by_iteration_markers(args.input_file, args.num_output_files, args.name_prefix)
