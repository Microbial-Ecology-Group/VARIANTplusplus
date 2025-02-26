#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def parse_msweep_results(input_dirs, output_file):
    """
    Parses mSWEEP result files from multiple directories, filtering by rel_abund > 0.001. 
    Summarizes precision, recall, F1, stores psv_abundance_diff, and includes num_reads/num_aligned 
    in summary. Then plots mean & std dev as line plots with error bars.
    Adds logic to parse total expected PSVs and expected reads from the filename.
    """

    long_form_data = []
    summarized_data = []

    # Gather all .txt files from all specified directories
    all_files = []
    for d in input_dirs:
        matched = glob.glob(os.path.join(d, '*.txt'))
        all_files.extend(matched)
        print(f"In directory '{d}', found {len(matched)} .txt files.")

    total_found = len(all_files)
    print(f"\nTotal .txt files found across all directories: {total_found}")
    if not all_files:
        print("No .txt files found. Please check your directory paths.")
        return

    # Thresholds
    REL_ABUND_THRESHOLD = 0.001
    DETECTION_THRESHOLD = 0  # (count >= 0 means detected, not used here explicitly but noted)

    # Process each file
    for file in all_files:
        base_name = os.path.basename(file).replace(".msweep_output_abundances.txt", "")
        #print(f"\nProcessing file: {file}, base_name: {base_name}")

        # Attempt to parse
        try:
            iter_num, numPSVs, k_col, genome_info = base_name.split("XX")
        except ValueError:
            print(f"File name {base_name} does not conform to 'XX' delimiter format.")
            continue

        # Keep track of expected PSVs
        expected_PSVs = {}
        expected_psv_keys = set()
        expected_pvs_count = 0
        expected_total_reads = 0

        if numPSVs == "0_PSVs":
            # off_target
            expected_PSVs = "off_target"
        else:
            # parse e.g. "PSV1.10000_PSV2.20000" ...
            parts = genome_info.split("_")
            for part in parts:
                try:
                    psv, reads_str = part.split(".")
                    reads_int = int(reads_str)
                    expected_PSVs[psv] = reads_int
                except ValueError:
                    print(f"Unexpected format in genome_info '{genome_info}' for chunk '{part}'")
                    continue
            if isinstance(expected_PSVs, dict):
                expected_psv_keys = set(expected_PSVs.keys())
                expected_pvs_count = len(expected_psv_keys)
                expected_total_reads = sum(expected_PSVs.values())

        # Prepare tracking
        num_reads = None
        num_aligned = None

        # psv_dict => psv->count for rel_abund>0.001
        psv_dict = {}

        try:
            with open(file, 'r') as f_in:
                for line in f_in:
                    line = line.strip()
                    if line.startswith("#num_reads"):
                        num_reads = int(line.split('\t')[1].strip())
                    elif line.startswith("#num_aligned"):
                        num_aligned = int(line.split('\t')[1].strip())
                    elif not line.startswith("#") and line:
                        parts2 = line.split('\t')
                        if len(parts2) == 2:
                            psv_group = parts2[0].strip()
                            rel_abund = float(parts2[1].strip())
                            if rel_abund > REL_ABUND_THRESHOLD:
                                count = (num_aligned or 0) * rel_abund
                                psv_dict[psv_group] = count

                                expected_reads = 0
                                if isinstance(expected_PSVs, dict):
                                    expected_reads = expected_PSVs.get(psv_group, 0)
                                diff_reads = expected_reads - count

                                exp_rel_ab = (expected_reads / expected_total_reads 
                                              if expected_total_reads > 0 else 0)
                                rel_ab_diff = rel_abund - exp_rel_ab

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
                                    'count': count,
                                    'expected_PSVs': (
                                        expected_PSVs if isinstance(expected_PSVs, dict) else {}
                                    ),
                                    'expected_reads': expected_reads,
                                    'difference': diff_reads,
                                    'expected_relative_abund': exp_rel_ab,
                                    'rel_abund_difference': rel_ab_diff
                                })

        except Exception as e:
            print(f"Error reading file {file}: {e}")
            continue

        #print(f"Finished reading file: {file}")

    # 1) Long-form DF
    df_long_form = pd.DataFrame(long_form_data)
    df_long_form.to_csv(f"{output_file}_long_form.txt", sep='\t', index=False)
    print(f"Saved {output_file}_long_form.txt")

    # 2) Summarize
    grouped = df_long_form.groupby(['file','iter_num'], as_index=False)

    for (fname, iteration), group in grouped:
        # from the group, pick the first row's metadata
        n_reads = group['num_reads'].iloc[0]
        n_aligned = group['num_aligned'].iloc[0]
        num_psvs_str = group['numPSVs'].iloc[0]
        k_col_val = group['k_col'].iloc[0]
        genome_info_val = group['genome_info'].iloc[0]

        # Reconstruct expected dict from the group
        possible_exp = group['expected_PSVs'].iloc[0]
        if isinstance(possible_exp, dict):
            expected_dict = possible_exp
            expected_keys = set(expected_dict.keys())
        else:
            expected_dict = {}
            expected_keys = set()

        # Build psv_counts => sum of "count"
        detected_psv_counts = group.groupby('psv_group')['count'].sum().to_dict()

        # Build psv_abundance_diff => for each psv in expected_dict
        psv_abundance_diff = {}
        for psv_k, exp_val in expected_dict.items():
            detected = detected_psv_counts.get(psv_k, 0)
            psv_abundance_diff[psv_k] = exp_val - detected

        # Precision/Recall
        if num_psvs_str == "0_PSVs":
            # off_target
            fp = len(detected_psv_counts)  # all are false positives
            tp = 0
            fn = 0
            precision = 0.0 if fp > 0 else 1.0
            recall = 0.0
            f1 = 0.0
            false_positive_count = fp
            perfect_match = "Y" if fp == 0 else "N"
            all_found = "N"
            false_positive_psv = "Y" if fp > 0 else "N"
            false_negative = "N"
            false_positive_mh = "Y" if fp > 0 else "N"
            true_negative = "Y" if fp == 0 else "N"
        else:
            # Non-zero
            detected_keys = set(detected_psv_counts.keys())
            tp = len(expected_keys & detected_keys)
            fp = len(detected_keys - expected_keys)
            fn = len(expected_keys - detected_keys)
            precision = tp/(tp+fp) if (tp+fp) > 0 else 0
            recall = tp/(tp+fn) if (tp+fn) > 0 else 0
            f1 = 2*(precision*recall)/(precision+recall) if (precision+recall) > 0 else 0
            false_positive_count = fp
            perfect_match = "Y" if expected_keys == detected_keys else "N"
            all_found = "Y" if expected_keys <= detected_keys else "N"
            false_positive_psv = "Y" if fp > 0 else "N"
            false_negative = "Y" if fn > 0 else "N"
            false_positive_mh = "N"
            true_negative = "N"

        summarized_data.append({
            'file': fname,
            'iter_num': iteration,
            'num_reads': n_reads,
            'num_aligned': n_aligned,
            'numPSVs': num_psvs_str,
            'k_col': k_col_val,
            'expected_pvs_count': len(expected_dict),
            'expected_total_reads': sum(expected_dict.values()) if expected_dict else 0,
            'expected_PSVs': str(expected_dict),
            'psv_groups': str(detected_psv_counts),
            'psv_abundance_diff': str(psv_abundance_diff),
            'Perfect_match': perfect_match,
            'All_found': all_found,
            'False_positive_PSV': false_positive_psv,
            'False_negative': false_negative,
            'False_positive_mh': false_positive_mh,
            'True_negative': true_negative,
            'False_positive_count': false_positive_count,
            'Precision': precision,
            'Recall': recall,
            'F1_score': f1
        })

    df_summarized = pd.DataFrame(summarized_data)

    # Round numeric columns
    for col in ['Precision','Recall','F1_score']:
        df_summarized[col] = df_summarized[col].round(1)

    # Save summary
    sum_name = f"{output_file}_summarized.txt"
    df_summarized.to_csv(sum_name, sep='\t', index=False)
    print(f"Saved {sum_name}")

    # -------------------------------------------------
    # Compute means & std per k_col => sn_sp_{output_file}
    # -------------------------------------------------
    sn_sp_data = []
    unique_k_cols = df_summarized['k_col'].unique()
    for kc in unique_k_cols:
        sub = df_summarized[df_summarized['k_col'] == kc]
        if sub.empty:
            continue
        # Means
        m_prec = sub['Precision'].mean()
        m_rec = sub['Recall'].mean()
        m_f1 = sub['F1_score'].mean()
        # Std
        s_prec = sub['Precision'].std()
        s_rec = sub['Recall'].std()
        s_f1 = sub['F1_score'].std()

        # Off-target => false_positive_count
        off_sub = sub[sub['numPSVs'] == '0_PSVs']
        avg_fp = off_sub['False_positive_count'].mean() if not off_sub.empty else 0

        sn_sp_data.append({
            'k_col': kc,
            'Precision_mean': round(m_prec, 1) if pd.notnull(m_prec) else None,
            'Precision_std': round(s_prec, 1) if pd.notnull(s_prec) else None,
            'Recall_mean': round(m_rec, 1) if pd.notnull(m_rec) else None,
            'Recall_std': round(s_rec, 1) if pd.notnull(s_rec) else None,
            'F1_score_mean': round(m_f1, 1) if pd.notnull(m_f1) else None,
            'F1_score_std': round(s_f1, 1) if pd.notnull(s_f1) else None,
            'False_positive_rate': avg_fp
        })

    df_sn_sp = pd.DataFrame(sn_sp_data)
    sp_name = f"{output_file}_sn_sp.txt"
    df_sn_sp.to_csv(sp_name, sep='\t', index=False)
    print(f"Saved {sp_name}")

    # Filter out "off_target"
    df_sn_sp_filt = df_sn_sp[df_sn_sp['k_col'] != 'off_target'].copy()
    # Convert k_col to numeric for sorting if possible
    df_sn_sp_filt['k_col_number'] = df_sn_sp_filt['k_col'].str.extract('(\\d+)').astype(float)
    df_sn_sp_filt.sort_values('k_col_number', inplace=True)

    # line plot
    x_labels = df_sn_sp_filt['k_col'].tolist()
    x = range(len(x_labels))

    prec_mean = df_sn_sp_filt['Precision_mean']
    prec_std = df_sn_sp_filt['Precision_std']
    rec_mean = df_sn_sp_filt['Recall_mean']
    rec_std = df_sn_sp_filt['Recall_std']
    f1_mean = df_sn_sp_filt['F1_score_mean']
    f1_std = df_sn_sp_filt['F1_score_std']
    fp_rate = df_sn_sp_filt['False_positive_rate']

    plt.figure(figsize=(8, 6))
    plt.errorbar(x, prec_mean, yerr=prec_std, fmt='-o', label='Precision')
    plt.errorbar(x, rec_mean, yerr=rec_std, fmt='-o', label='Recall')
    plt.errorbar(x, f1_mean, yerr=f1_std, fmt='-o', label='F1 Score')
    plt.plot(x, fp_rate, '-o', label='False Positive Rate', color='red')

    plt.xticks(x, x_labels, rotation=45)
    plt.xlabel('k_col')
    plt.ylabel('Value')
    plt.title('Precision, Recall, F1, and False Positive Rate by k_col')
    plt.legend()
    plt.tight_layout()

    fig_name = f"{output_file}_sn_sp_plot.png"
    plt.savefig(fig_name)
    plt.show()
    print(f"Saved line plot with error bars as {fig_name}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Parse mSWEEP result files from multiple directories, filter rel_abund>0.001, "
            "compute precision/recall/f1, store psv_abundance_diff, parse total expected PSVs & reads, "
            "and plot line with error bars."
        )
    )
    # Allow multiple directories
    parser.add_argument(
        'input_dirs', 
        type=str, 
        nargs='+', 
        help='One or more directories containing .txt mSWEEP outputs'
    )
    parser.add_argument('output_file', type=str, help='Base name for output files (e.g., results)')
    args = parser.parse_args()

    parse_msweep_results(args.input_dirs, args.output_file)


if __name__ == "__main__":
    main()
