#!/usr/bin/env python3
import sys

def split_by_iteration_markers(input_file, num_output_files, name_prefix):
    """
    Reads 'input_file' line by line, grouping lines between occurrences of
    either:
      1) "# - Starting megasample creation" (if present), OR
      2) "# --- Starting iteration" (fallback).

    Distributes these groups across 'num_output_files' round-robin.
    For each output index i:
      - Writes a "command file": name_prefix_part_{i}.sh
      - Writes a minimal SBATCH file: name_prefix_part_{i}.sbatch

    Also tracks the chosen marker lines and "# Unique Sample" lines, inserting
    echo statements for progress.

    Finally, prints how many marker lines (either iteration or megasample) and
    sample lines were found in total, plus distribution across the output files.
    """

    # 0) First pass: detect if there's a "megasample" marker
    with open(input_file, 'r') as fcheck:
        all_lines = fcheck.readlines()

    megasample_marker = "# - Starting megasample"
    iteration_marker  = "# --- Starting iteration"

    # Determine which marker to use
    use_megasample_marker = any(line.startswith(megasample_marker) for line in all_lines)
    if use_megasample_marker:
        marker = megasample_marker
        print(f"Using marker '{marker}' because lines start with '# - Starting megasample'")
    else:
        marker = iteration_marker
        print(f"No lines start with '# - Starting megasample', using '{marker}' instead.")

    # 1) Split lines into groups based on the chosen marker
    groups = []
    current_group = []

    marker_count_global = 0
    sample_count_global = 0

    for line in all_lines:
        if line.startswith(marker):
            marker_count_global += 1
            if current_group:
                groups.append(current_group)
            current_group = [line]
        else:
            current_group.append(line)

        if line.startswith("# Unique Sample"):
            sample_count_global += 1

    # Don't forget the last group if it's non-empty
    if current_group:
        groups.append(current_group)

    # 2) Distribute groups among the output files (round-robin)
    file_contents = [[] for _ in range(num_output_files)]
    for idx, grp in enumerate(groups):
        file_index = idx % num_output_files
        file_contents[file_index].extend(grp)

    # 3) For each file, create command + SBATCH scripts
    for i in range(num_output_files):
        part_index = i + 1

        cmd_filename = f"{name_prefix}_part_{part_index}.sh"
        sbatch_filename = f"{name_prefix}_part_{part_index}.sbatch"

        # Count marker & sample lines in this part
        marker_count_file = 0
        sample_count_file = 0
        for line in file_contents[i]:
            if line.startswith(marker):
                marker_count_file += 1
            if line.startswith("# Unique Sample"):
                sample_count_file += 1

        # Prepare the command file
        with open(cmd_filename, 'w') as cmd_f:
            # Optional header for summary in this part
            cmd_f.write(
                f"echo \"[{name_prefix}] Script {part_index}/{num_output_files}: "
                f"{marker_count_file} marker(s) in this file of {marker_count_global} total; "
                f"{sample_count_file} sample(s) in this file of {sample_count_global} total.\"\n\n"
            )

            # Insert an echo line before each "# Unique Sample"
            sample_index_in_file = 0
            for line in file_contents[i]:
                if line.startswith("# Unique Sample"):
                    sample_index_in_file += 1
                    cmd_f.write(
                        f"echo \"Now running sample {sample_index_in_file} / {sample_count_file} "
                        f"in script {part_index}.\"\n"
                    )

                cmd_f.write(line)

        # Prepare the minimal SBATCH file
        sbatch_header = (
            "#!/bin/bash\n"
            f"#SBATCH -J {name_prefix}_{part_index}\n"
            f"#SBATCH -o log_{name_prefix}_{part_index}.out\n"
            "#SBATCH -t 24:00:00\n"
            "#SBATCH --mem=40G\n"
            "#SBATCH --nodes=1\n"
            "#SBATCH --ntasks=1\n"
            "#SBATCH --cpus-per-task=48\n"
            "#SBATCH --ntasks-per-node=1\n\n"
        )

        with open(sbatch_filename, 'w') as sb_f:
            sb_f.write(sbatch_header)
            sb_f.write(f"echo 'Launching {cmd_filename} in script {part_index}'\n")
            sb_f.write(f"bash {cmd_filename}\n")

        print(
            f"Created {cmd_filename} (commands) and {sbatch_filename} (SBATCH). "
            f"Lines in part={len(file_contents[i])} "
            f"({marker_count_file} marker lines, {sample_count_file} samples)."
        )

    # 4) Print global summary
    print("\n=== Global Summary ===")
    print(f"Total marker lines ('{marker}') found in '{input_file}': {marker_count_global}")
    print(f"Total sample lines ('# Unique Sample') found in '{input_file}': {sample_count_global}")

def main():
    # Expect 3 arguments: input_file, num_output_files, name_prefix
    if len(sys.argv) != 4:
        print("Usage: python split_script_iteration_markers.py <input_file> <num_output_files> <name_prefix>")
        sys.exit(1)

    input_file = sys.argv[1]
    num_output_files = int(sys.argv[2])
    name_prefix = sys.argv[3]

    split_by_iteration_markers(input_file, num_output_files, name_prefix)

if __name__ == "__main__":
    main()
