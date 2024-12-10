import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def parse_msweep_results(input_dir, output_file):
    # Lists to store all data
    long_form_data = []
    summarized_data = []

    # Identify all .txt files in the directory
    files = glob.glob(os.path.join(input_dir, '*.txt'))
    print(f"Found {len(files)} files in the directory: {input_dir}")

    if not files:
        print("No .txt files found. Please check the directory path.")
        return

    # Process each file
    for file in files:
        # Parse the file name
        base_name = os.path.basename(file).replace(".msweep_output_abundances.txt", "")  # Remove file name suffix
        print(f"Processing file: {file}, Base name: {base_name}")

        # Parse the base name into iter_num, numPSVs, k_col, and genome_info
        try:
            iter_num, numPSVs, k_col, genome_info = base_name.split("XX")
        except ValueError:
            print(f"File name {base_name} does not conform to the expected format 'XX' delimiter.")
            continue
        
        # Initialize variables for the current file
        num_reads = None
        num_aligned = None
        expected_PSVs = {}

        # Parse genome_info to extract expected_PSVs and expected_reads
        if numPSVs == "0_PSVs":
            expected_PSVs = "off_target"
            expected_psv_keys = set()  # Set expected_psv_keys to an empty set for off_target
        else:
            genome_parts = genome_info.split("_")
            print(genome_parts)
            for part in genome_parts:
                try:
                    psv, reads = part.split(".")
                    expected_PSVs[psv] = int(reads)
                except ValueError:
                    print(f"Unexpected format in genome_info '{genome_info}': '{part}' does not split correctly.")
                    continue
            expected_psv_keys = set(expected_PSVs.keys())  # Set expected_psv_keys for PSVs

        # Calculate the total expected reads for calculating expected_relative_abund
        total_expected_reads = sum(expected_PSVs.values()) if expected_PSVs != "off_target" else 0

        try:
            # Open and read the file
            with open(file, 'r') as f:
                for line in f:
                    line = line.strip()

                    if line.startswith("#num_reads"):
                        num_reads = int(line.split('\t')[1].strip())
                    elif line.startswith("#num_aligned"):
                        num_aligned = int(line.split('\t')[1].strip())
                    elif not line.startswith("#") and line.strip():
                        parts = line.split('\t')
                        if len(parts) == 2:
                            psv_group = parts[0].strip()
                            rel_abund = float(parts[1].strip())
                            count = num_aligned * rel_abund

                            # Calculate the expected reads and differences
                            expected_reads = expected_PSVs.get(psv_group, 0) if expected_PSVs != "off_target" else 0
                            difference = expected_reads - count
                            expected_relative_abund = expected_reads / total_expected_reads if total_expected_reads > 0 else 0
                            rel_abund_difference = rel_abund - expected_relative_abund

                            # Add data to the list
                            long_form_data.append({
                                'file': base_name,
                                'iter_num': iter_num,
                                'numPSVs': numPSVs,
                                'k_col': k_col,
                                'genome_info': genome_info,
                                'num_reads': num_reads,
                                'num_aligned': num_aligned,
                                'psv_group': psv_group,
                                'rel_abund': rel_abund,
                                'count': count,  # Add the count column
                                'expected_PSVs': expected_PSVs if expected_PSVs != "off_target" else {},  # Include expected_PSVs
                                'expected_reads': expected_reads,  # Include expected reads for the psv_group
                                'difference': difference,  # Include difference between expected reads and count
                                'expected_relative_abund': expected_relative_abund,  # Include expected relative abundance
                                'rel_abund_difference': rel_abund_difference  # Include difference between rel_abund and expected_relative_abund
                            })
        except Exception as e:
            print(f"Error processing file {file}: {e}")
            continue

        print(f"Finished processing file: {file}")

    # Create a DataFrame from the collected data
    df_long_form = pd.DataFrame(long_form_data)
    
    # Save the long-form data table to a file
    df_long_form.to_csv(f"long_form_{output_file}", index=False, sep='\t')
    print(f"Long-form data table saved to long_form_{output_file}")

    # Process summarized data
    grouped = df_long_form.groupby(['file', 'iter_num'])
    
    for (file_name, iter_num), group in grouped:
        # Initialize variables for summary
        expected_PSVs = group['expected_PSVs'].iloc[0]
        psv_groups = set(group.loc[group['rel_abund'] > 0.01, 'psv_group'])
        perfect_match = "N"
        all_found = "N"
        false_positive_psv = "N"
        false_negative = "N"
        false_positive_mh = "N"
        true_negative = "N"
        precision = recall = f1_score = 0.0
        false_positive_count = 0

        # Logic for match checking
        if expected_PSVs == "off_target":
            false_positive_count = len(psv_groups)  # Number of false positives
            if len(psv_groups) > 0:
                false_positive_mh = "Y"
            elif len(psv_groups) == 0:
                true_negative = "Y"
        else:
            tp = len(expected_psv_keys & psv_groups)  # True positives: Correctly identified PSVs
            fp = len(psv_groups - expected_psv_keys)  # False positives: Incorrectly identified PSVs
            fn = len(expected_psv_keys - psv_groups)  # False negatives: Missed PSVs

            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

            if expected_psv_keys == psv_groups:
                perfect_match = "Y"
            if expected_psv_keys <= psv_groups:
                all_found = "Y"
            if fp > 0:
                false_positive_psv = "Y"
            if fn > 0:
                false_negative = "Y"
        
        summarized_data.append({
            'file': file_name,
            'iter_num': iter_num,
            'numPSVs': group['numPSVs'].iloc[0],
            'k_col': group['k_col'].iloc[0],
            'expected_PSVs': str(expected_PSVs),  # Add expected_PSVs as a string
            'psv_groups': str(psv_groups),  # Add psv_groups as a string
            'Perfect_match': perfect_match,
            'All_found': all_found,
            'False_positive_PSV': false_positive_psv,
            'False_negative': false_negative,
            'False_positive_mh': false_positive_mh,
            'True_negative': true_negative,
            'False_positive_count': false_positive_count,
            'Precision': precision,
            'Recall': recall,
            'F1_score': f1_score
        })

    # Create a DataFrame from the summarized data
    df_summarized = pd.DataFrame(summarized_data)
    
    # Save the summarized data to a file
    df_summarized.to_csv(f"summarized_{output_file}", index=False, sep='\t')
    print(f"Summarized data table saved to summarized_{output_file}")

    # Calculate precision, recall, F1 score, and false positives for numPSVs == 0
    sn_sp_data = []

    # Calculate metrics for each k_col
    for k_col in df_summarized['k_col'].unique():
        subset = df_summarized[df_summarized['k_col'] == k_col]
        avg_precision = subset['Precision'].mean()
        avg_recall = subset['Recall'].mean()
        avg_f1_score = subset['F1_score'].mean()
        avg_false_positives = subset[subset['numPSVs'] == '0_PSVs']['False_positive_count'].mean()

        sn_sp_data.append({
            'k_col': k_col,
            'Precision': avg_precision,
            'Recall': avg_recall,
            'F1_score': avg_f1_score,
            'False_positive_rate': avg_false_positives
        })

    # Create a DataFrame for precision, recall, F1 score, and false positive rate
    df_sn_sp = pd.DataFrame(sn_sp_data)
    
    # Save the precision, recall, F1 score, and false positive rate data to a file
    df_sn_sp.to_csv(f"sn_sp_{output_file}", index=False, sep='\t')
    print(f"Precision, Recall, F1 Score, and False Positive Rate data saved to sn_sp_{output_file}")
    
    # Exclude 'off_target' rows from the plot
    df_sn_sp_filtered = df_sn_sp[df_sn_sp['k_col'] != 'off_target'].copy()
    
    # Extract the numerical part from 'k_col' and sort by it
    df_sn_sp_filtered.loc[:, 'k_col_number'] = df_sn_sp_filtered['k_col'].str.extract('(\d+)').astype(int)
    df_sn_sp_filtered = df_sn_sp_filtered.sort_values('k_col_number')
    
    # Plot Precision, Recall, F1 Score, and False Positive Rate by k_col, excluding 'off_target'
    plt.figure(figsize=(10, 6))
    plt.plot(df_sn_sp_filtered['k_col'], df_sn_sp_filtered['Precision'], marker='o', label='Precision')
    plt.plot(df_sn_sp_filtered['k_col'], df_sn_sp_filtered['Recall'], marker='o', label='Recall')
    plt.plot(df_sn_sp_filtered['k_col'], df_sn_sp_filtered['F1_score'], marker='o', label='F1 Score')
    plt.plot(df_sn_sp_filtered['k_col'], df_sn_sp_filtered['False_positive_rate'], marker='o', label='False Positive Rate')
    plt.xlabel('k_col')
    plt.ylabel('Value')
    plt.title('Precision, Recall, F1 Score, and False Positive Rate by k_col')
    plt.legend()
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'sn_sp_plot_{output_file.replace(".txt", ".png")}')
    plt.show()
    print(f"Plot saved as sn_sp_plot_{output_file.replace('.txt', '.png')}")

if __name__ == "__main__":
    import argparse

    # Setup argument parser
    parser = argparse.ArgumentParser(description='Parse mSWEEP result files into a long-form data table and calculate precision, recall, F1 score, and false positive rate.')
    parser.add_argument('input_dir', type=str, help='Directory containing the mSWEEP .txt files')
    parser.add_argument('output_file', type=str, help='Base name for the output files')

    # Parse arguments
    args = parser.parse_args()

    # Run the parsing function
    parse_msweep_results(args.input_dir, args.output_file)
