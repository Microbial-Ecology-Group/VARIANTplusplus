#!/usr/bin/env python3
import os
import random
import argparse
random.seed(2)

"""
Usage:
    python generate_simulated_samples_byConf.py --benchmarking_params params.txt

Description:
- This script reads user config (benchmarking_params),
- Creates 3 bash scripts:
  1) script_build_reads_{PREFIX}.sh
       * For each (numGSV, iteration, k_col):
         - If numGSV==0 => pick non-target => cat references => run ISS => run flash => final partial_name => 
           rename_sample_files.py => rename_headers.py ...
         - If numGSV>0 => pick multiple GSVs => cat references => run ISS => run flash => final partial_name => 
           rename_sample_files.py => rename_headers.py ...
       * merges partial_name files across k_col => final iteration-level file for classification
  2) script_kraken_extract_split_{PREFIX}.sh
       * Classifies each iteration-level final merged file with Kraken2 across multiple kraken_confidence values
       * Extracts reads using extract_kraken_reads.py
       * Splits merged/unmerged reads
  3) script_themisto_msweep_{PREFIX}.sh
       * Runs themisto pseudoalign on each partial iteration-level sample per confidence,
         then runs mSWEEP separately for merged & unmerged reads.
"""

def parse_benchmarking_params(params_file):
    """
    Reads a key=value file with config variables.
    Returns a dict.
    """
    config = {}
    with open(params_file, 'r') as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' not in line:
                continue
            key, val = line.split('=', 1)
            key = key.strip()
            val = val.strip()

            if val.startswith('[') or val.startswith('('):
                try:
                    config[key] = eval(val)
                except:
                    config[key] = val
            elif val.isdigit():
                config[key] = int(val)
            else:
                config[key] = val
    return config

def get_full_path(base_dir, filename):
    """
    Locates a .fna or .fna.gz file.
    Returns the appropriate existing path or None if not found.
    """
    path_uncompressed = os.path.join(base_dir, filename)
    path_gz = path_uncompressed + '.gz'
    if filename.endswith('.gz'):
        if os.path.exists(path_uncompressed):
            return path_uncompressed
        alt_uncompressed = path_uncompressed.replace('.gz','')
        if os.path.exists(alt_uncompressed):
            return alt_uncompressed
        return None
    else:
        if os.path.exists(path_uncompressed):
            return path_uncompressed
        if os.path.exists(path_gz):
            return path_gz
        return None

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build scripts for pipeline that merges references for multiple GSVs into a single .fna, "
            "runs ISS once (and flash to merge R1/R2), merges partial samples across k_col for classification, "
            "and calls rename_sample_files.py => rename_headers.py. The themisto + mSWEEP sections remain intact."
        )
    )
    parser.add_argument(
        '--benchmarking_params',
        type=str,
        required=True,
        help='Path to the text file containing all configurable variables.'
    )
    args = parser.parse_args()

    # 1) Parse external benchmarking params
    params = parse_benchmarking_params(args.benchmarking_params)

    # Retrieve user config
    target_genome_dir = params.get("target_genome_dir")
    nontarget_genome_dir = params.get("nontarget_genome_dir")
    input_file = params.get("input_file")
    bin_dir = params.get("bin_dir")
    themisto_index = params.get("themisto_index")
    krakendb = params.get("krakendb")

    kraken_confidence = params.get("kraken_confidence", [0])
    if not isinstance(kraken_confidence, list):
        kraken_confidence = [kraken_confidence]

    extract_reads_taxid = str(params.get("extract_reads_taxid", "75985"))
    extract_reads_options = params.get("extract_reads_options", "--include-children")
    num_iters = params.get("num_iters", 100)
    num_GSV_list = params.get("num_GSV_list", [0,1,3,5,7,9])
    num_reads_options = params.get("num_reads_options", [10000,20000])
    threads = params.get("threads", 48)
    segment_length = params.get("segment_length", 150)
    num_genomes_per_GSV = params.get("num_genomes_per_GSV", 5)
    output_name_prefix = params.get("output_name_prefix", "default_output_prefix")
    tmp_build = params.get("tmp_build", "$TMPDIR")

    # Output scripts
    script_build_reads = f"script_build_reads_{output_name_prefix}.sh"
    script_kraken_extract_split = f"script_kraken_extract_split_{output_name_prefix}.sh"
    script_run_themisto_msweep = f"script_themisto_msweep_{output_name_prefix}.sh"

    # Directories
    cat_reads_dir = f"{output_name_prefix}_cat_reads"
    split_reads_dir = f"{output_name_prefix}_split_reads"
    results_dir = f"{output_name_prefix}_results"
    merged_extracted_reads_dir = f"{output_name_prefix}_merged_reads"

    os.makedirs(cat_reads_dir, exist_ok=True)
    os.makedirs(split_reads_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    if not tmp_build.startswith("$"):
        os.makedirs(tmp_build, exist_ok=True)

    # Read metadata file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    k_columns = lines[0].strip().split('\t')[1:]
    genome_data = []
    for line in lines[1:]:
        parts = line.strip().split('\t')
        gf = parts[0]
        gn = ".".join(gf.split('.')[:2])
        gn = os.path.basename(gn)
        ann = parts[1:]
        genome_data.append((gf, gn, ann))

    # Create annotation files
    for idx, k_col in enumerate(k_columns):
        msweep_annot_file = os.path.join(results_dir, f"{k_col}_msweep.txt")
        with open(msweep_annot_file, 'w') as out_annot:
            for (_gf, _gn, annotations) in genome_data:
                out_annot.write(f"{annotations[idx]}\n")

    sample_map = {}
    iteration_to_renamed = {}

    # --------- 1) script_build_reads -----------
    print("Building script:", script_build_reads)
    with open(script_build_reads, 'w') as build_out:

        for numGSV in num_GSV_list:
            build_out.write(f"\n# ====== Building samples for numGSV={numGSV} ======\n")

            for iteration in range(1, num_iters+1):
                build_out.write(f"\n# - Starting megasample creation for numGSV={numGSV}\n")
                build_out.write(f"# --- Iteration {iteration} ---\n")

                iter_key = (numGSV, iteration)
                iteration_to_renamed[iter_key] = []

                for k_idx, k_col in enumerate(k_columns):
                    build_out.write(f"\n# (k_col={k_col})\n")

                    if numGSV == 0:
                        # Off-target => pick references => run ISS => run flash => partial_name
                        all_nt = [
                            f for f in os.listdir(nontarget_genome_dir)
                            if f.endswith(".fna") or f.endswith(".fna.gz")
                        ]
                        chosen_nt = random.sample(all_nt, min(num_genomes_per_GSV, len(all_nt)))

                        tmp_cat = f"{tmp_build}/iter_{iteration}_0GSVs_{k_col}.fna"
                        build_out.write(f"rm -f {tmp_cat}\n")
                        for gfName in chosen_nt:
                            real_path = get_full_path(nontarget_genome_dir, gfName)
                            if gfName.endswith('.gz'):
                                build_out.write(f"zcat {real_path} >> {tmp_cat}\n")
                            else:
                                build_out.write(f"cat {real_path} >> {tmp_cat}\n")

                        read_count = random.choice(num_reads_options)
                        # Use the same naming approach as the >0 GSVs branch
                        iss_prefix = f"{tmp_build}/iter_{iteration}XX{numGSV}_GSVsXX{k_col}_iss"
                        nodir_prefix = f"iter_{iteration}XX{numGSV}_GSVsXX{k_col}_iss"

                        build_out.write(
                            f"iss generate --genomes {tmp_cat} --n_reads {read_count} --model NovaSeq "
                            f"--cpus {threads} --output {iss_prefix} --compress --quiet --fragment-length 311 --fragment-length-sd 50\n"
                        )

                        partial_name = f"{iteration}XX{numGSV}_GSVsXX{k_col}XXoff_target"

                        # Use flash to merge R1 & R2
                        build_out.write(
                            f"flash -M 120 -o {partial_name} -d {tmp_build} --interleaved-output -z -t {threads} "
                            f"{tmp_build}/{nodir_prefix}_R1.fastq.gz {tmp_build}/{nodir_prefix}_R2.fastq.gz\n"
                        )
                        build_out.write(
                            f"cat {tmp_build}/{partial_name}.extendedFrags.fastq.gz "
                            f"{tmp_build}/{partial_name}.notCombined.fastq.gz "
                            f"> {tmp_build}/{partial_name}.fastq.gz\n"
                        )

                        # rename_sample_files.py
                        build_out.write(
                            f"python3 {bin_dir}/rename_sample_files.py --input {tmp_build}/{partial_name}.fastq.gz "
                            f"--results-dir {results_dir}\n"
                        )

                        renamed_part = f"{partial_name}_renamed.fastq.gz"
                        build_out.write(
                            f"python3 {bin_dir}/rename_headers.py --input {tmp_build}/{partial_name}.fastq.gz "
                            f"--sample_id \"{partial_name}.fastq.gz\" "
                            f"--output {tmp_build}/{renamed_part}\n"
                        )
                        build_out.write(f"rm {tmp_build}/{partial_name}.fastq.gz\n")
                        build_out.write(f"rm {iss_prefix}*tmp*\n")

                        iteration_to_renamed[iter_key].append(renamed_part)
                        sample_map[partial_name] = "off_target"

                    else:
                        # numGSV>0 => combine multiple GSVs => single ISS => flash => partial_name
                        all_col_GSVs = list(set(x[2][k_idx] for x in genome_data))
                        if not all_col_GSVs:
                            build_out.write(f"# WARNING: no GSV in k_col={k_col}\n")
                            continue

                        selected_GSVs = random.sample(all_col_GSVs, min(numGSV, len(all_col_GSVs)))
                        build_out.write(f"# selected GSVs => {selected_GSVs}\n")

                        tmp_cat = f"{tmp_build}/iter_{iteration}_{k_col}_GSVs_combined.fna"
                        build_out.write(f"rm -f {tmp_cat}\n")

                        GSV_labels = []
                        for GSV_val in selected_GSVs:
                            matched = [xx for xx in genome_data if xx[2][k_idx] == GSV_val]
                            chosen_g = random.sample(matched, min(num_genomes_per_GSV, len(matched)))

                            GSV_fna = f"{tmp_build}/iter_{iteration}_{GSV_val}_unrenamed.fna"
                            build_out.write(f"rm -f {GSV_fna}\n")
                            for (gFile, gName, _anns) in chosen_g:
                                real_path = get_full_path(target_genome_dir, gFile)
                                if gFile.endswith('.gz'):
                                    build_out.write(f"zcat {real_path} >> {GSV_fna}\n")
                                else:
                                    build_out.write(f"cat {real_path} >> {GSV_fna}\n")

                            GSV_fna_ren = f"{tmp_build}/iter_{iteration}_{GSV_val}_renamed.fna"
                            build_out.write(
                                f"python3 {bin_dir}/rename_headers.py --input {GSV_fna} --sample_id \"{GSV_val}\" "
                                f"--output {GSV_fna_ren}\n"
                            )
                            build_out.write(f"cat {GSV_fna_ren} >> {tmp_cat}\n")
                            build_out.write(f"rm {GSV_fna} {GSV_fna_ren}\n")
                            GSV_labels.append(GSV_val)

                        read_count = random.choice(num_reads_options)
                        iss_prefix = f"{tmp_build}/iter_{iteration}XX{numGSV}_GSVsXX{k_col}_iss"
                        nodir_prefix= f"iter_{iteration}XX{numGSV}_GSVsXX{k_col}_iss"
                        build_out.write(
                            f"iss generate --genomes {tmp_cat} --n_reads {read_count} --model NovaSeq "
                            f"--cpus {threads} --output {iss_prefix} --compress --quiet --fragment-length 311 --fragment-length-sd 50\n"
                        )

                        GSV_str = "_".join(GSV_labels)
                        partial_name= f"{iteration}XX{numGSV}_GSVsXX{k_col}XX{GSV_str}"
                        build_out.write(
                            f"flash -M 120 -o {partial_name} -d {tmp_build} --interleaved-output -z -t {threads} "
                            f"{tmp_build}/{nodir_prefix}_R1.fastq.gz {tmp_build}/{nodir_prefix}_R2.fastq.gz\n"
                        )
                        build_out.write(
                            f"cat {tmp_build}/{partial_name}.extendedFrags.fastq.gz "
                            f"{tmp_build}/{partial_name}.notCombined.fastq.gz "
                            f"> {tmp_build}/{partial_name}.fastq.gz\n"
                        )

                        # rename_sample_files.py
                        build_out.write(
                            f"python3 {bin_dir}/rename_sample_files.py --input {tmp_build}/{partial_name}.fastq.gz "
                            f"--results-dir {results_dir}\n"
                        )

                        renamed_part= f"{partial_name}_renamed.fastq.gz"
                        build_out.write(
                            f"python3 {bin_dir}/rename_headers.py --input {tmp_build}/{partial_name}.fastq.gz "
                            f"--sample_id \"{partial_name}.fastq.gz\" "
                            f"--output {tmp_build}/{renamed_part}\n"
                        )
                        build_out.write(f"rm {tmp_build}/{partial_name}.fastq.gz\n")
                        build_out.write(f"rm {iss_prefix}*tmp*\n")

                        iteration_to_renamed[iter_key].append(renamed_part)
                        sample_map[partial_name] = k_col

                # merges partial_name across k_col => iteration-level final
                final_merged_iter= f"GSV_{numGSV}_iter_{iteration}_ALL_renamed.fastq.gz"
                partial_list= iteration_to_renamed[iter_key]
                if partial_list:
                    build_out.write(
                      f"cat {' '.join(f'{tmp_build}/{nm}' for nm in partial_list)} "
                      f"> {cat_reads_dir}/{final_merged_iter}\n"
                    )
                else:
                    build_out.write("# no iteration-level partial => skip merges.\n")


    print("Building:", script_kraken_extract_split)
    with open(script_kraken_extract_split, 'w') as cl_out:
        for numGSV in num_GSV_list:
            for iteration in range(1, num_iters+1):
                final_merged = f"GSV_{numGSV}_iter_{iteration}_ALL_renamed.fastq.gz"
                final_merged_path = f"{cat_reads_dir}/{final_merged}"

                for conf in kraken_confidence:
                    #conf_str = str(conf).replace('.', 'p')
                    conf_str = str(conf)
                    kr_report = f"{tmp_build}/{final_merged}.conf{conf_str}.kraken.report"
                    kr_raw = f"{tmp_build}/{final_merged}.conf{conf_str}.kraken.raw"
                    extracted_name = final_merged.replace('.fastq.gz', f"_conf{conf_str}_extracted.fastq")
                    extracted_path = f"{tmp_build}/{extracted_name}"
                    cl_out.write(f"\n# - Starting megasample classification for num_GSV={numGSV}, iteration={iteration}\n")
                    cl_out.write(f"\n# Kraken2 classification at confidence {conf}\n")
                    cl_out.write(
                        f"kraken2 --db {krakendb} --threads {threads} --confidence {conf} "
                        f"--report {kr_report} --output {kr_raw} {final_merged_path}\n"
                    )

                    cl_out.write(
                        f"extract_kraken_reads.py -k {kr_raw} --max 1000000000 --report {kr_report} "
                        f"--taxid {extract_reads_taxid} {extract_reads_options} "
                        f"-s {final_merged_path} -o {extracted_path} --fastq-output\n"
                    )

                    cl_out.write(f"pigz --processes {threads} {extracted_path}\n")

                    cl_out.write(
                        f"python3 {bin_dir}/split_extracted_reads_by_conf.py --input_fasta {extracted_path}.gz "
                        f"--output_merged_dir {merged_extracted_reads_dir} "
                        f"--output_unmerged_dir {split_reads_dir} --label conf{conf_str} \n"
                    )

    print("Building:", script_run_themisto_msweep)
    with open(script_run_themisto_msweep, 'w') as run_out:
        for partial_name, col_info in sample_map.items():
            for conf in kraken_confidence:
                #conf_str = str(conf).replace('.', 'p')
                conf_str = str(conf)
                merged_paired_path = os.path.join(merged_extracted_reads_dir, partial_name + f"_conf{conf_str}_merged.fastq.gz")
                split_paired_path  = os.path.join(split_reads_dir, partial_name + f"_conf{conf_str}_unmerged.fastq.gz")
                themisto_out = os.path.join(results_dir, partial_name + f"_conf{conf_str}.themisto_output")

                run_out.write(f"\n# Themisto + mSWEEP for {partial_name} at conf={conf}\n")
                run_out.write(
                    f"themisto pseudoalign -q {merged_paired_path} -i {themisto_index} --temp-dir {tmp_build} -t {threads} "
                    f"--gzip-output --sort-output-lines -o {themisto_out}.merged\n"
                )
                run_out.write(
                    f"themisto pseudoalign -q {split_paired_path} -i {themisto_index} --temp-dir {tmp_build} -t {threads} "
                    f"--gzip-output --sort-output-lines -o {themisto_out}.unmerged\n"
                )

                if col_info == "off_target":
                    for col_name in k_columns:
                        msweep_annot = os.path.join(results_dir, f"{col_name}_msweep.txt")
                        msweep_out   = os.path.join(results_dir, f"{partial_name}_conf{conf_str}.{col_name}.msweep_output")
                        run_out.write(
                            f"mSWEEP --themisto {themisto_out}.merged.gz -i {msweep_annot} -t {threads} -o {msweep_out}.merged --verbose --themisto-mode intersection\n"
                        )
                        run_out.write(
                            f"mSWEEP --themisto {themisto_out}.unmerged.gz -i {msweep_annot} -t {threads} -o {msweep_out}.unmerged --verbose --themisto-mode intersection\n"
                        )
                else:
                    msweep_annot = os.path.join(results_dir, f"{col_info}_msweep.txt")
                    msweep_out   = os.path.join(results_dir, f"{partial_name}_conf{conf_str}.{col_info}.msweep_output")
                    run_out.write(
                        f"mSWEEP --themisto {themisto_out}.merged.gz -i {msweep_annot} -t {threads} -o {msweep_out}.merged --verbose\n"
                    )
                    run_out.write(
                        f"mSWEEP --themisto {themisto_out}.unmerged.gz -i {msweep_annot} -t {threads} -o {msweep_out}.unmerged --verbose\n"
                    )

    print(f"\nDONE generating:\n- {script_build_reads}\n- {script_kraken_extract_split}\n- {script_run_themisto_msweep}\n")

if __name__ == "__main__":
    main()
