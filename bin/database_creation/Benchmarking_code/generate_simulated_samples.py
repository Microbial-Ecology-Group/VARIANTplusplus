#!/usr/bin/env python3
import os
import random
import argparse

"""
Usage:
    python generate_simulated_samples.py --benchmarking_params params.txt

This script generates 3 shell scripts:
  1) script_build_reads_{PREFIX}.sh
     - merges selected reference genomes into one FASTA
     - uses inSilicoSeq (ISS) to simulate FASTQ reads (single or paired-end)
     - renames and merges them for each iteration & condition
  2) script_kraken_extract_split_{PREFIX}.sh
     - runs Kraken2 to classify the final FASTQ
     - optionally extracts reads belonging to a given taxid
     - splits the extracted reads if needed
  3) script_themisto_msweep_{PREFIX}.sh
     - runs themisto pseudoalignment
     - runs mSWEEP for abundance estimation

All key config variables are read from a 'benchmarking_params' file, specifying:
- genome directories, iteration counts, read lengths, etc.
"""

def parse_benchmarking_params(params_file):
    """
    Parse a key=value text file containing these variables:
      - target_genome_dir, nontarget_genome_dir, input_file
      - themisto_index, krakendb
      - kraken_confidence, extract_reads_taxid, extract_reads_options
      - num_iters, num_PSV_list, num_reads_options
      - threads, segment_length, num_genomes_per_PSV
      - output_name_prefix, tmp_build

    Lines starting with '#' or blank lines are ignored.
    For lists, we expect Python list syntax, e.g.:  num_PSV_list=[0,1,3].
    Returns a dict of {key -> parsed_value}.
    """
    config = {}
    with open(params_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Ignore blank or commented lines
            if not line or line.startswith('#'):
                continue
            # Must have key=value format
            if '=' not in line:
                continue
            key, val = line.split('=', 1)
            key = key.strip()
            val = val.strip()

            # Attempt to parse as list or integer
            if val.startswith('[') or val.startswith('('):
                try:
                    config[key] = eval(val)  # parse list/tuple
                except:
                    config[key] = val  # fallback as string
            else:
                # Maybe it's an integer
                if val.isdigit():
                    config[key] = int(val)
                else:
                    # else keep as string (paths, floats, etc.)
                    config[key] = val
    return config


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build scripts for mSWEEP+Themisto pipeline, using inSilicoSeq "
            "to generate FASTQ reads from concatenated references."
        )
    )
    parser.add_argument(
        '--benchmarking_params',
        type=str,
        required=True,
        help='Path to the text file containing all configurable variables.'
    )
    args = parser.parse_args()

    # -----------------------------------------------------
    # 1) Parse the external benchmarking params
    # -----------------------------------------------------
    params = parse_benchmarking_params(args.benchmarking_params)

    # Retrieve the key variables from the dict
    target_genome_dir = params.get("target_genome_dir")
    nontarget_genome_dir = params.get("nontarget_genome_dir")
    input_file = params.get("input_file")

    themisto_index = params.get("themisto_index")
    krakendb = params.get("krakendb")

    # For Kraken2 + read extraction
    kraken_confidence = str(params.get("kraken_confidence", 0))
    extract_reads_taxid = str(params.get("extract_reads_taxid", "75985"))
    extract_reads_options = params.get("extract_reads_options", "--include-children")

    # For iteration logic
    num_iters = params.get("num_iters", 100)
    num_PSV_list = params.get("num_PSV_list", [0,1,3,5,7,9])
    num_reads_options = params.get("num_reads_options", [10000,20000])
    threads = params.get("threads", 48)
    segment_length = params.get("segment_length", 150)  # Not directly used by iss
    num_genomes_per_PSV = params.get("num_genomes_per_PSV", 5)

    # Output prefixes / directories
    output_name_prefix = params.get("output_name_prefix", "set1_100i_NTcore_0conf")
    tmp_build = params.get("tmp_build", "$TMPDIR")

    # -----------------------------------------------------
    # 2) Derivative paths and output script names
    # -----------------------------------------------------
    script_build_reads = f"script_build_reads_{output_name_prefix}.sh"
    script_kraken_extract_split = f"script_kraken_extract_split_{output_name_prefix}.sh"
    script_run_themisto_msweep = f"script_themisto_msweep_{output_name_prefix}.sh"

    # Folders
    merged_reads_dir = f"{output_name_prefix}_merged_reads"
    split_reads_dir = f"{output_name_prefix}_split_reads"
    results_dir = f"{output_name_prefix}_results"

    # Ensure output directories exist
    os.makedirs(merged_reads_dir, exist_ok=True)
    os.makedirs(split_reads_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    # If tmp_build is $TMPDIR, we skip creation (on HPC, $TMPDIR is set)
    if tmp_build.startswith("$"):
        print(f"Skipping creation of tmp_build directory, using literal {tmp_build}")
    else:
        os.makedirs(tmp_build, exist_ok=True)
        print(f"Created actual tmp_build directory: {tmp_build}")

    # -----------------------------------------------------
    # Read the input metadata file, skipping first column
    # The first column is the genome filename, the subsequent columns are k_# ...
    # e.g.:
    #   GenomeFilename  k_2  k_3  k_4 ...
    #   genomeX.fna.gz  PSV1 ...
    # We'll read them into (genome_file_name, genome_name, [k2, k3, k4, ...])
    # -----------------------------------------------------
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # The second line onward are actual data lines
    k_columns = lines[0].strip().split("\t")[1:]  # skip the first col name

    genome_data = []
    for line in lines[1:]:
        parts = line.strip().split("\t")
        genome_file_name = parts[0]
        # e.g. "genomeXYZ.fna.gz" => "genomeXYZ.fna"
        genome_name = ".".join(genome_file_name.split('.')[:2])
        genome_name = os.path.basename(genome_name)

        ani_annotations = parts[1:]  # the k_col labels
        genome_data.append((genome_file_name, genome_name, ani_annotations))

    # -----------------------------------------------------
    # Create annotation files for mSWEEP (one per k_col)
    # e.g. "k_2_msweep.txt" => each line is that row's annotation
    # -----------------------------------------------------
    for ani_idx, k_col in enumerate(k_columns):
        msweep_annot_file = os.path.join(results_dir, f"{k_col}_msweep.txt")
        with open(msweep_annot_file, 'w') as out_annot:
            for (_gf, _gn, annotations) in genome_data:
                out_annot.write(f"{annotations[ani_idx]}\n")

    # Helper function to locate a .fna or .fna.gz file
    def get_full_path(base_dir, filename):
        path_uncompressed = os.path.join(base_dir, filename)
        path_gz = path_uncompressed + '.gz'
        if filename.endswith('.gz'):
            # e.g. "genomeXYZ.fna.gz"
            if os.path.exists(path_uncompressed):
                return path_uncompressed
            alt_uncompressed = path_uncompressed.replace('.gz', '')
            if os.path.exists(alt_uncompressed):
                return alt_uncompressed
            return None
        else:
            # e.g. "genomeXYZ.fna"
            if os.path.exists(path_uncompressed):
                return path_uncompressed
            if os.path.exists(path_gz):
                return path_gz
            return None

    # We'll produce three scripts
    psv_to_renamed = {psv: [] for psv in num_PSV_list}
    sample_map = {}

    # -----------------------------------------------------
    # SCRIPT 1: Build iteration reads
    # Concatenate references => run inSilicoSeq => produce FASTQ
    # rename headers => merges
    # -----------------------------------------------------
    print("Building script:", script_build_reads)
    with open(script_build_reads, 'w') as build_out:

        sample_counter = 0

        for num_PSV in num_PSV_list:
            build_out.write(f"\n# - Starting megasample creation for num_PSV={num_PSV}\n")

            for iter_num in range(1, num_iters + 1):
                build_out.write(f"\n# --- Starting iteration {iter_num} for num_PSV={num_PSV} ---\n")

                # We'll build a single merged references file
                temp_merged_fna = f"{tmp_build}/iter_{iter_num}_PSV_{num_PSV}_merged.fna"
                build_out.write(f"rm -f {temp_merged_fna}\n")

                if num_PSV == 0:
                    continue # remove this line to include off-target genomes
                    # Off-target logic: pick non-target genomes
                    all_nontarget_genomes = [
                        f for f in os.listdir(nontarget_genome_dir)
                        if f.endswith(".fna") or f.endswith(".fna.gz")
                    ]
                    chosen = random.sample(
                        all_nontarget_genomes,
                        min(num_genomes_per_PSV, len(all_nontarget_genomes))
                    )

                    # Cat or zcat each chosen genome
                    for gfName in chosen:
                        real_path = get_full_path(nontarget_genome_dir, gfName)
                        if gfName.endswith(".gz"):
                            build_out.write(f"zcat {real_path} >> {temp_merged_fna}\n")
                        else:
                            build_out.write(f"cat {real_path} >> {temp_merged_fna}\n")

                    # Then pick a single random read count
                    final_reads_count = random.choice(num_reads_options)
                    # The final name => e.g. "1XX0_PSVsXXoff_targetXX10000.fastq.gz"
                    final_name = f"{iter_num}XX0_PSVsXXoff_targetXX{final_reads_count}.fastq.gz"

                    # Run inSilicoSeq => produce FASTQ
                    build_out.write(
                        f"iss generate --genomes {temp_merged_fna} "
                        f"--n_reads {final_reads_count} "
                        f"--model Miseq-20 "
                        f"--cpus {threads} "
                        f"--output {tmp_build}/iter_{iter_num}_0PSV_{k_col}_iss "
                        f"--compress --quiet\n"
                    )
                    # Combine single-end or paired-end reads into final FASTQ
                    build_out.write(
                        f"cat {tmp_build}/iter_{iter_num}_0PSV_{k_col}_iss_R* > {tmp_build}/{final_name}\n"
                    )
                    # Clean up
                    build_out.write(f"rm {temp_merged_fna}\n")
                    build_out.write(f"rm {tmp_build}/iter_{iter_num}_0PSV_{k_col}_iss_R*.fastq.gz\n")

                    # rename headers
                    renamed_name = final_name.replace(".fastq.gz", "_renamed.fastq.gz")
                    build_out.write(
                        f"python3 rename_headers.py --input {tmp_build}/{final_name} "
                        f'--sample_id "{final_name}" '
                        f"--output {tmp_build}/{renamed_name}\n"
                    )
                    build_out.write(f"rm {tmp_build}/{final_name}\n")

                    # Store references for merging & classification
                    psv_to_renamed[0].append(renamed_name)
                    sample_map[final_name] = "off_target"
                    sample_counter += 1
                    build_out.write(f"# Unique Sample {sample_counter}: num_PSV=0\n")

                else:
                    # num_PSV > 0 => we do it per k_col
                    for ani_idx, k_col in enumerate(k_columns):
                        unique_psvs = list(set(g[2][ani_idx] for g in genome_data))
                        # pick which psvs are in this iteration
                        selected = random.sample(unique_psvs, min(num_PSV, len(unique_psvs)))

                        # single read count for entire iteration
                        some_read_count = random.choice(num_reads_options)
                        psv_sums = {}

                        # For each chosen psv, store that read count & cat references
                        for psv_val in selected:
                            gwp = [x for x in genome_data if x[2][ani_idx] == psv_val]
                            if not gwp:
                                build_out.write(f"# WARNING: no genomes for PSV={psv_val} in {k_col}\n")
                                continue

                            # Assign read count
                            psv_sums[psv_val] = some_read_count

                            # pick up to num_genomes_per_PSV from that label
                            chosen_genomes = random.sample(gwp, min(num_genomes_per_PSV, len(gwp)))
                            for (gFile, gName, _anns) in chosen_genomes:
                                real_path = get_full_path(target_genome_dir, gFile)
                                if real_path is None:
                                    build_out.write("# WARNING: Could not find genome file:\n")
                                    build_out.write(f"# {target_genome_dir}/{gFile}\n")
                                    continue

                                if gFile.endswith(".gz"):
                                    build_out.write(f"zcat {real_path} >> {temp_merged_fna}\n")
                                else:
                                    build_out.write(f"cat {real_path} >> {temp_merged_fna}\n")

                        # Build naming string => "PSV1.10000_PSV2.10000"
                        psv_strings = [f"{p}.{psv_sums[p]}" for p in sorted(psv_sums.keys())]
                        psv_info_str = "_".join(psv_strings)

                        final_reads_count = sum(psv_sums.values())
                        # e.g. "1XX2_PSVsXXk_3XXPSV1.10000_PSV5.10000.fastq.gz"
                        final_name = f"{iter_num}XX{num_PSV}_PSVsXX{k_col}XX{psv_info_str}.fastq.gz"

                        # inSilicoSeq
                        build_out.write(
                            f"iss generate --genomes {temp_merged_fna} "
                            f"--n_reads {final_reads_count} "
                            f"--model Miseq-20 "
                            f"--cpus {threads} "
                            f"--output {tmp_build}/iter_{iter_num}_{num_PSV}_PSVs_{k_col}_iss "
                            f"--compress --quiet\n"
                        )
                        build_out.write(
                            f"cat {tmp_build}/iter_{iter_num}_{num_PSV}_PSVs_{k_col}_iss_R* > {tmp_build}/{final_name}\n"
                        )
                        build_out.write(f"rm {temp_merged_fna}\n")
                        build_out.write(f"rm {tmp_build}/iter_{iter_num}_{num_PSV}_PSVs_{k_col}_iss_R*.fastq.gz\n")

                        # rename headers
                        renamed_name = final_name.replace(".fastq.gz", "_renamed.fastq.gz")
                        build_out.write(
                            f"python3 rename_headers.py --input {tmp_build}/{final_name} "
                            f'--sample_id "{final_name}" '
                            f"--output {tmp_build}/{renamed_name}\n"
                        )
                        build_out.write(f"rm {tmp_build}/{final_name}\n")

                        psv_to_renamed[num_PSV].append(renamed_name)
                        sample_map[final_name] = k_col
                        sample_counter += 1
                        build_out.write(f"# Unique Sample {sample_counter}: num_PSV={num_PSV}, k_col={k_col}\n")

            # After finishing all iterations for this num_PSV,
            # merge them into a single file => e.g. "PSV_1_ALL_renamed.fastq.gz"
            build_out.write(f"\n# Merge all iteration-level renamed files for psv={num_PSV} into merged_reads_dir\n")
            merged_name = f"PSV_{num_PSV}_ALL_renamed.fastq.gz"
            if psv_to_renamed[num_PSV]:
                cat_line = (
                    "cat " + " ".join(f"{tmp_build}/{nm}" for nm in psv_to_renamed[num_PSV]) +
                    f" > {merged_reads_dir}/{merged_name}\n"
                )
                build_out.write(cat_line)
            else:
                build_out.write(f"# No iteration files for psv={num_PSV}, skipping merge.\n")

    print("\nDone. The script_build_reads is complete.")

    # -----------------------------------------------------
    # SCRIPT 2: Classification, extraction, splitting
    # -----------------------------------------------------
    print("Building script:", script_kraken_extract_split)
    with open(script_kraken_extract_split, 'w') as cl_out:
        for psv in num_PSV_list:
            # We'll use the merged final from script 1
            merged_name = f"PSV_{psv}_ALL_renamed.fastq.gz"
            merged_path = f"{merged_reads_dir}/{merged_name}"

            kr_report = f"{tmp_build}/PSV_{psv}.kraken.report"
            kr_raw = f"{tmp_build}/PSV_{psv}.kraken.raw"
            # We'll store the extracted reads in a .fastq => then compress
            extracted_name = f"PSV_{psv}_extracted.fastq"
            extracted_path = f"{tmp_build}/{extracted_name}"
            splitted_subdir = split_reads_dir

            cl_out.write(f"\n# - Starting megasample classification for num_PSV={psv}\n")
            cl_out.write(f"# --- Classification for psv={psv} ---\n")

            # Run kraken2
            kr_cmd = (
                f"kraken2 --db {krakendb} --threads {threads} "
                f"--confidence {kraken_confidence} "
                f"--report {kr_report} "
                f"--output {kr_raw} {merged_path}\n"
            )
            cl_out.write(kr_cmd)

            # Extract reads matching 'extract_reads_taxid'
            extr_cmd = (
                f"extract_kraken_reads.py -k {kr_raw} --max 1000000000 "
                f"--report {kr_report} --taxid {extract_reads_taxid} {extract_reads_options} "
                f"-s {merged_path} -o {extracted_path}\n"
            )
            cl_out.write(extr_cmd)

            # Compress the extracted fastq
            cl_out.write(f"pigz --processes {threads} {extracted_path}\n")

            # Split the extracted .fastq.gz if needed
            # If your script expects the input as "--input_fasta", adapt accordingly
            split_cmd = (
                f"python3 split_extracted_reads.py "
                f"--input_fasta {extracted_path}.gz "
                f"--output_dir {splitted_subdir}\n"
            )
            cl_out.write(split_cmd)

    # -----------------------------------------------------
    # SCRIPT 3: themisto + mSWEEP
    # -----------------------------------------------------
    print("Building script:", script_run_themisto_msweep)
    all_final_names = list(sample_map.keys())

    with open(script_run_themisto_msweep, 'w') as run_out:
        for final_name, col_info in sample_map.items():
            # e.g. final_name = "1XX1_PSVsXXk_2XXPSV1.10000.fastq.gz"
            base_name = final_name.replace(".fastq.gz", "")
            parts = base_name.split("XX")
            if len(parts) != 4:
                print(f"WARNING: final_name {final_name} doesn't have 4 parts with 'XX'")
                continue

            iteration_num = parts[0]
            second = parts[1]  # e.g. "0_PSVs" or "1_PSVs"
            if "_" in second:
                psv_str = second.split("_")[0]
            else:
                print(f"WARNING: cannot parse psv from {second}")
                continue

            try:
                psv_int = int(psv_str)
            except ValueError:
                print(f"WARNING: cannot parse psv as int from {psv_str}")
                continue

            # The splitted read set => e.g. "1XX1_PSVsXXk_2XXPSV1.10000.fastq.gz" => in split_reads_dir
            splitted_path = os.path.join(split_reads_dir, f"{final_name}")
            themisto_out = os.path.join(results_dir, f"{base_name}.themisto_output")

            run_out.write(f"\n# --- Starting iteration {iteration_num}, running themisto + mSWEEP for sample: {base_name}, psv={psv_int} ---\n")

            # themisto pseudoalign
            themisto_cmd = (
                f"themisto pseudoalign -q {splitted_path} "
                f"-i {themisto_index} "
                f"--temp-dir tmp -t {threads} --gzip-output --sort-output-lines "
                f"-o {themisto_out} {splitted_path}\n"
            )
            run_out.write(themisto_cmd)

            # If it's off_target, we run mSWEEP for every k_col
            if col_info == "off_target":
                for kc in k_columns:
                    msweep_annot = os.path.join(results_dir, f"{kc}_msweep.txt")
                    msweep_out = os.path.join(results_dir, f"{base_name}.{kc}.msweep_output")
                    ms_cmd = (
                        f"mSWEEP --themisto {themisto_out}.gz "
                        f"-i {msweep_annot} -t {threads} "
                        f"-o {msweep_out} --verbose\n"
                    )
                    run_out.write(ms_cmd)
            else:
                # Otherwise, we run it for the specific col_info that produced this sample
                msweep_annot = os.path.join(results_dir, f"{col_info}_msweep.txt")
                msweep_out = os.path.join(results_dir, f"{base_name}.msweep_output")
                ms_cmd = (
                    f"mSWEEP --themisto {themisto_out}.gz "
                    f"-i {msweep_annot} -t {threads} "
                    f"-o {msweep_out} --verbose\n"
                )
                run_out.write(ms_cmd)

    print("\nDone generating 3 scripts:\n"
          f"1) {script_build_reads} (handles iteration-level merges -> inSilicoSeq -> .fastq.gz)\n"
          f"2) {script_kraken_extract_split} (Kraken2 classification + read extraction)\n"
          f"3) {script_run_themisto_msweep} (Themisto + mSWEEP)\n"
          f"Annotation files are created in {results_dir}.\n")


if __name__ == "__main__":
    main()
