import random
import os

# usage: python make_benchmarking_samples_unique_PSVs_wIters
# Working directory needs to include the "make_genome_reads.py"
# Output is two bash scripts, these need to be run from the same working directory. 
# 5 output directories will be made in your working dir. 

# Script fuction:
# Samples are generated from randomly selected genomes, short sample reads are made, those are subset randomly based on a "num_reads_options" variable, then classified with kraken, filtered 
# based on the "extract_reads_taxid" variable, then classified using themisto, finally the corresponding annotation set will be used to parse the results with mSWEEP. 

# Define variables
# Input directory for target genomes
target_genome_dir = "/scratch/user/enriquedoster/Mh_validation_paper_results/Mh_genomes/Mh_genbank_contig/"
# Input directory for non-target genomes
nontarget_genome_dir = "/scratch/user/enriquedoster/Mh_validation_paper_results/NonMh_Past_genomes/"
# This metadata file comes from running cutree on R on ANI distances.
# The input file containing genome names and ANI annotations, first column is full genome file name (no path) with no header, the other column headers are k_2 - k_#.
# make sure the genome names are in the same order as the genomes in your themisto index
input_file = "ANI_groupings_2-30.tsv"  
# Input themisto index
themisto_index = "/scratch/group/vero_research/databases/themisto_index/themisto_index/themisto_index"
# Input kraken database used for classification and filtering
krakendb = "/scratch/group/vero_research/databases/kraken2/PlusPF"  # Update with the actual Kraken2 database path
# Variables for kraken
kraken_confidence = 0  # Confidence threshold for Kraken2 (0 is least stringent)
extract_reads_taxid = "75984"  # TaxID to extract reads, change this to your taxaID from NCBI
extract_reads_options_single = ""  # Optional flag space for extract_kraken_reads.py

# Benchmarking characteristics 
num_iters = 100  # Number of iterations for each num_PSV
num_PSV_list = [0, 3, 6, 9] # How many unique PSVs to be represented in your benchmarking samples. 
num_reads_options = [5000, 20000]  # Possible number of reads to subset
threads = 48 # I used 48 because I would hog a node for each set of commands in multiple sbatch scripts
segment_length = 150  # Reads segment length for make_genome_reads.py

# Name for first output script with all commands to make, filter, and classify benchmarking samples
output_script = "run_seqtk_commands_uniquePSVs_mh.sh"
# Second output script with commands to run themisto and mSWEEP
run_them_msweep_script = "run_them_msweep_mh.sh"  # New script for Themisto and mSWEEP commands
# Variable for adding a useful naming suffix
output_name_suffix = "75984" # You can name this whatever you want and it'll go at the end of your output directories

#
##
# Nothing else to input below here ##############
##
#

# Makes output folder names
output_dir = f"output_genome_reads_{output_name_suffix}"  # Directory for generated genome reads
output_random_reads_dir = f"output_random_reads_{output_name_suffix}"
output_msweep_annotations = f"output_msweep_annotations_{output_name_suffix}" #  Msweep annotations from your main metadata file in the format to work with mSWEEP
output_themisto_alignments = f"output_themisto_alignments_{output_name_suffix}"
output_msweep_results = f"output_msweep_results_{output_name_suffix}"


# Create output directories if they don't exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(output_random_reads_dir, exist_ok=True)
os.makedirs(output_msweep_annotations, exist_ok=True)
os.makedirs(output_themisto_alignments, exist_ok=True)
os.makedirs(output_msweep_results, exist_ok=True)

# Read input metadata file
with open(input_file, 'r') as file:
    lines = file.readlines()

# Parse the header to identify k_ columns
k_columns = lines[0].strip().split("\t")
print(k_columns)

# Parse the metadata to get genome names
genome_data = []
for line in lines[1:]:
    parts = line.strip().split("\t")
    genome_file_name = parts[0]
    genome_name = ".".join(genome_file_name.split('.')[:2])  # Extract everything to the left of the second "."
    ani_annotations = parts[1:]
    genome_data.append((genome_file_name, genome_name, ani_annotations))

# Write bash commands for running "make_genome_reads.py" to make reads from genomes, seqtk randomly picking reads, Kraken2 classification and extraction, and concatenation to output file
with open(output_script, 'w') as out_file:
    out_file.write("#!/bin/bash\n\n")

    # List to store final output read names for subsequent processing
    final_output_read_names = []

    # Iterate over each PSV number
    for num_PSV in num_PSV_list:
        for iter_num in range(1, num_iters + 1):  # Adding iteration layer
            if num_PSV == 0: # To handle "non-target" genomes 
                # For PSV 0, select 10 random genomes from all files in the non target genomes directory
                all_nontarget_genomes = [f for f in os.listdir(nontarget_genome_dir) if f.endswith('.fna') or f.endswith('.fna.gz')]
                selected_nontarget_genomes = random.sample(all_nontarget_genomes, 10)

                num_nontarget_reads = 0
                for genome_file_name in selected_nontarget_genomes:
                    # Generate genome reads using make_genome_reads.py
                    genome_name = ".".join(genome_file_name.split('.')[:2])
                    genome_reads_file = f"{genome_name}_reads.fasta"
                    make_reads_cmd = f"python3 make_genome_reads.py {nontarget_genome_dir}{genome_file_name} {segment_length} {output_dir}/{genome_reads_file}\n"
                    out_file.write(make_reads_cmd)

                    # Kraken2 command
                    kraken_report_file = f"{output_dir}/{genome_name}.kraken.report"
                    kraken_raw_file = f"{output_dir}/{genome_name}.kraken.raw"
                    kraken_cmd = f"kraken2 --db {krakendb} --confidence {kraken_confidence} --threads {threads} --report {kraken_report_file} --output {kraken_raw_file} {output_dir}/{genome_reads_file}\n"
                    out_file.write(kraken_cmd)

                    # Extract Kraken reads command with updated format
                    extracted_reads_file = f"{output_dir}/{genome_name}_extracted.fasta"
                    extract_cmd = (
                        f"extract_kraken_reads.py -k {kraken_raw_file} --max 1000000000 --report {kraken_report_file} "
                        f"--taxid {extract_reads_taxid} -s {output_dir}/{genome_reads_file} -o {extracted_reads_file}\n"
                    )
                    out_file.write(extract_cmd)

                    # Compress extracted reads using pigz
                    pigz_cmd = f"pigz --processes {threads} {extracted_reads_file}\n"
                    out_file.write(pigz_cmd)

                    # Randomly select number of reads for seqtk
                    num_reads = random.choice(num_reads_options)
                    num_nontarget_reads += num_reads

                    # Generate seqtk commands for the reads with new naming convention
                    seqtk_output_name = f"{iter_num}_{num_PSV}_{num_reads}_{genome_name}.fasta.gz"  # Modified to include iteration number
                    seqtk_cmd = f"seqtk sample -s100 {extracted_reads_file}.gz {num_reads} > {output_dir}/{seqtk_output_name}\n"
                    out_file.write(seqtk_cmd)

                output_reads = f"off_target.{num_nontarget_reads}"

                # Concatenate the reads
                cat_cmd = f"cat {' '.join([f'{output_dir}/{iter_num}_{num_PSV}_{num_reads}_{genome_name}.fasta.gz' for genome_file_name in selected_nontarget_genomes])} > {output_random_reads_dir}/{iter_num}XX0_PSVsXXoff_targetXX{output_reads}XX.fasta.gz\n"
                out_file.write(cat_cmd)

                # Add to final output read names for subsequent processing
                final_output_read_names.append(f"{iter_num}XX0_PSVsXXoff_targetXX{output_reads}XX.fasta.gz")

            else:
                # For other PSV numbers, select random genomes and associated k_column values
                for ani_set_index, k_col in enumerate(k_columns):
                    # Identify all unique PSV values in the current k_# column
                    unique_psvs = list(set(genome[2][ani_set_index] for genome in genome_data))  # Convert to list to avoid deprecation warning

                    # Randomly select a subset of unique PSVs according to num_PSV
                    selected_psvs = random.sample(unique_psvs, min(num_PSV, len(unique_psvs)))

                    print(f"starting with {k_col}")
                    print(f"unique_psvs {unique_psvs}")
                    selected_genomes = []
                    seqtk_info = []  # To store PSV and num_reads for final concatenated file name
                    for psv in selected_psvs:
                        # Select genomes that have the current PSV value
                        genomes_with_psv = [genome for genome in genome_data if genome[2][ani_set_index] == psv]
                        if genomes_with_psv:
                            selected_genome = random.choice(genomes_with_psv)
                            selected_genomes.append((selected_genome, psv))

                    # Generate the output_reads name based on selected PSVs 
                    print(selected_psvs)

                    # Create msweep annotation file for current k_col
                    # No row names or headers, must match order of genomes in themisto database
                    msweep_annotation_file = os.path.join(output_msweep_annotations, f"{k_col}_msweep.txt")
                    with open(msweep_annotation_file, 'w') as msweep_file:
                        for genome_file_name, genome_name, annotations in genome_data:
                            msweep_file.write(f"{annotations[ani_set_index]}\n")

                    cat_cmd_input_files = []

                    for (genome_file_name, genome_name, _), psv in selected_genomes:
                        # Generate genome reads using make_genome_reads.py
                        genome_reads_file = f"{genome_name}_reads.fasta"
                        make_reads_cmd = f"python3 make_genome_reads.py {target_genome_dir}{genome_file_name} {segment_length} {output_dir}/{genome_reads_file}\n"
                        out_file.write(make_reads_cmd)

                        # Kraken2 command
                        kraken_report_file = f"{output_dir}/{genome_name}.kraken.report"
                        kraken_raw_file = f"{output_dir}/{genome_name}.kraken.raw"
                        kraken_cmd = f"kraken2 --db {krakendb} --confidence {kraken_confidence} --threads {threads} --report {kraken_report_file} --output {kraken_raw_file} {output_dir}/{genome_reads_file}\n"
                        out_file.write(kraken_cmd)

                        # Extract Kraken reads command with updated format
                        extracted_reads_file = f"{output_dir}/{genome_name}_extracted.fasta"
                        extract_cmd = (
                            f"extract_kraken_reads.py -k {kraken_raw_file} --max 1000000000 --report {kraken_report_file} "
                            f"--taxid {extract_reads_taxid} -s {output_dir}/{genome_reads_file} -o {extracted_reads_file}\n"
                        )
                        out_file.write(extract_cmd)

                        # Compress extracted reads using pigz
                        pigz_cmd = f"pigz --processes {threads} {extracted_reads_file}\n"
                        out_file.write(pigz_cmd)

                        # Randomly select number of reads for seqtk
                        num_reads = random.choice(num_reads_options)
                        # Generate seqtk commands for the reads with new naming convention
                        seqtk_output_name = f"{iter_num}_{psv}_{num_reads}.{k_col}.{genome_name}.fasta.gz"  # Include iteration number
                        seqtk_cmd = f"seqtk sample -s100 {extracted_reads_file}.gz {num_reads} > {output_dir}/{seqtk_output_name}\n"
                        out_file.write(seqtk_cmd)

                        # Store PSV and num_reads for final concatenated file name
                        seqtk_info.append(f"{psv}.{num_reads}")

                        # Prepare input files for cat command
                        cat_cmd_input_files.append(f"{output_dir}/{seqtk_output_name}")

                    # Concatenate the reads, using PSV and num_reads for each genome
                    final_output_reads = f"{iter_num}XX{num_PSV}_PSVsXX{k_col}XX{'_'.join(seqtk_info)}XX.fasta.gz"
                    cat_cmd = f"cat {' '.join(cat_cmd_input_files)} > {output_random_reads_dir}/{final_output_reads}\n"
                    out_file.write(cat_cmd)
                    print(f"done with {k_col}")

                    # Add to final output read names for subsequent processing
                    final_output_read_names.append(final_output_reads)

# Write Themisto alignment and mSWEEP commands to separate script
with open(run_them_msweep_script, 'w') as them_msweep_out_file:
    them_msweep_out_file.write("#!/bin/bash\n\n")

    # Generate Themisto alignment and mSWEEP commands for each concatenated file
    for final_output in final_output_read_names:
        print(final_output)
        base_name = final_output.replace("XX.fasta.gz", "")
        print(base_name)

        # Parse the base name into iter_num, numPSVs, k_col, and genome_info
        iter_num, numPSVs, k_col, genome_info = base_name.split("XX")
        print(numPSVs)

        # Themisto alignment command
        themisto_output_file = os.path.join(output_themisto_alignments, f"{base_name}.themisto_output")
        themisto_cmd = f"themisto pseudoalign -q {output_random_reads_dir}/{final_output} -i {themisto_index} --temp-dir tmp -t {threads} --gzip-output --sort-output-lines -o {themisto_output_file} {output_random_reads_dir}/{final_output}\n"
        them_msweep_out_file.write(themisto_cmd)

        if k_col == "off_target":
            # For 'off_target', loop through all k_#_msweep.txt annotations
            for k_column in k_columns:
                msweep_output_file = os.path.join(output_msweep_results, f"{base_name}.{k_column}.msweep_output")
                msweep_annotation_file = os.path.join(output_msweep_annotations, f"{k_column}_msweep.txt")
                msweep_cmd = f"mSWEEP --themisto {themisto_output_file}.gz -i {msweep_annotation_file} -t {threads} -o {msweep_output_file} --verbose\n"
                them_msweep_out_file.write(msweep_cmd)
        else:
            # Handle other k_col values
            msweep_output_file = os.path.join(output_msweep_results, f"{base_name}.msweep_output")
            msweep_annotation_file = os.path.join(output_msweep_annotations, f"{k_col}_msweep.txt")
            msweep_cmd = f"mSWEEP --themisto {themisto_output_file}.gz -i {msweep_annotation_file} -t {threads} -o {msweep_output_file} --verbose\n"
            them_msweep_out_file.write(msweep_cmd)

print(f"Commands have been written to {output_script}")
print(f"Themisto and mSWEEP commands have been written to {run_them_msweep_script}")
