#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def parse_msweep_results(input_dirs, output_file, filter_mode,
                         rel_abund_threshold=0.01, detection_threshold=10):
    """
    Parses mSWEEP result files from multiple directories.

    Filtering:
      - If filter_mode == 'rel_abund':
           Keep only PSVs with rel_abund > rel_abund_threshold.
      - If filter_mode == 'count':
           Keep only PSVs with (rel_abund * num_aligned) > detection_threshold.
      - If filter_mode == 'count_rel_abund':
           Keep only PSVs where (rel_abund * num_aligned) > (num_aligned * rel_abund_threshold).

    Summaries:
      - Once a PSV passes the filter, we compute 'count' = rel_abund * num_aligned
        for it, store that in the long_form data, and eventually compute precision, recall, etc.
      - 'found_psvs' column shows all psv_group->raw rel_abund pairs from the file (before filtering).
      - 'Num_false_negative' column shows the integer count of false negatives.
    """

    long_form_data = []
    summarized_data = []

    # Dictionary to store all PSVs (with raw rel_abund) for each file => found_psvs_no_filter
    # keyed by the base_name we extract from the file
    found_psvs_all = {}

    # Collect all .txt files from all specified directories
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

    # ------------------------------------------------------
    # Main loop: Parse each file's output
    # ------------------------------------------------------
    for file in all_files:
        base_name = os.path.basename(file).replace(".msweep_output_abundances.txt", "")

        # Attempt to parse expected format => iter_num + "XX" + numPSVs + "XX" + k_col + "XX" + genome_info
        try:
            iter_num, numPSVs, k_col, genome_info = base_name.split("XX")
        except ValueError:
            print(f"File name {base_name} does not conform to 'XX' delimiter format.")
            continue

        # Parse expected PSVs (unless '0_PSVs' => off_target)
        expected_PSVs = {}
        if numPSVs == "0_PSVs":
            expected_PSVs = "off_target"
            expected_total_reads = 0
        else:
            parts = genome_info.split("_")
            for part in parts:
                try:
                    psv, reads_str = part.split(".")
                    reads_int = int(reads_str)
                    expected_PSVs[psv] = reads_int
                except ValueError:
                    print(f"Unexpected format in genome_info '{genome_info}' for chunk '{part}'")
                    continue
            expected_total_reads = sum(expected_PSVs.values()) if expected_PSVs else 0

        num_reads = None
        num_aligned = None

        # psv_dict => "detected" PSVs (i.e., pass the chosen filter)
        psv_dict = {}

        # found_psvs_no_filter => dictionary of ALL PSVs & their *raw* rel_abund from the file
        found_psvs_no_filter = {}

        # ------------------------------------------------------
        # Read the .txt file
        # ------------------------------------------------------
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

                            # Store raw rel_abund in found_psvs_no_filter
                            found_psvs_no_filter[psv_group] = rel_abund

                            # Calculate raw count
                            count_val = (num_aligned or 0) * rel_abund

                            # -------------------------
                            #  Choose which filter to apply
                            # -------------------------
                            if filter_mode == 'rel_abund':
                                # PSV is "kept" if rel_abund > rel_abund_threshold
                                if rel_abund > rel_abund_threshold:
                                    psv_dict[psv_group] = count_val

                            elif filter_mode == 'count':
                                # PSV is "kept" if count_val > detection_threshold
                                if count_val > detection_threshold:
                                    psv_dict[psv_group] = count_val

                            elif filter_mode == 'count_rel_abund':
                                # PSV is "kept" if (rel_abund * num_aligned) > (num_aligned * rel_abund_threshold)
                                sample_read_abund_threshold = (num_aligned or 0) * rel_abund_threshold
                                if count_val > sample_read_abund_threshold:
                                    psv_dict[psv_group] = count_val

                            else:
                                raise ValueError(f"Unknown filter_mode: {filter_mode}")

            # After reading the file, store found_psvs_no_filter so we can retrieve it in summary
            found_psvs_all[base_name] = found_psvs_no_filter

            # For all PSVs that passed the filter, record them in long_form_data
            for psv_group, cval in psv_dict.items():
                if isinstance(expected_PSVs, dict):
                    expected_reads = expected_PSVs.get(psv_group, 0)
                else:
                    expected_reads = 0

                diff_reads = expected_reads - cval
                exp_rel_ab = (expected_reads / expected_total_reads) if expected_total_reads > 0 else 0
                # Recompute the "actual" rel_ab based on the retained count
                new_rel_ab = cval / (num_aligned or 1)
                rel_ab_diff = new_rel_ab - exp_rel_ab

                long_form_data.append({
                    'file': base_name,
                    'iter_num': iter_num,
                    'numPSVs': numPSVs,
                    'k_col': k_col,
                    'genome_info': genome_info,
                    'num_reads': num_reads,
                    'num_aligned': num_aligned,
                    'psv_group': psv_group,
                    'rel_abund': new_rel_ab,
                    'count': cval,
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

    # ------------------------------------------------
    # 1) Long-form DF
    # ------------------------------------------------
    df_long_form = pd.DataFrame(long_form_data)
    df_long_form_filename = f"{output_file}_long_form.txt"
    df_long_form.to_csv(df_long_form_filename, sep='\t', index=False)
    print(f"Saved {df_long_form_filename}")

    # ------------------------------------------------
    # 2) Summarize
    # ------------------------------------------------
    grouped = df_long_form.groupby(['file','iter_num'], as_index=False)

    for (fname, iteration), group in grouped:
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

        # psv_abundance_diff => for each psv in expected_dict
        psv_abundance_diff = {}
        for psv_k, exp_val in expected_dict.items():
            detected = detected_psv_counts.get(psv_k, 0)
            psv_abundance_diff[psv_k] = exp_val - detected

        # Precision/Recall, false negative count
        if num_psvs_str == "0_PSVs":
            # off_target => everything is a false positive
            fp = len(detected_psv_counts)
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
            detected_keys = set(detected_psv_counts.keys())
            tp = len(expected_keys & detected_keys)
            fp = len(detected_keys - expected_keys)
            fn = len(expected_keys - detected_keys)  # <--- This is our false negatives
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

        # Retrieve the dict of ALL PSVs from found_psvs_all
        # (raw rel_abund as reported by mSWEEP)
        found_psvs_dict = found_psvs_all.get(fname, {})

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
            'found_psvs': str(found_psvs_dict),
            'psv_abundance_diff': str(psv_abundance_diff),
            'Perfect_match': perfect_match,
            'All_found': all_found,
            'False_positive_PSV': false_positive_psv,
            'False_negative': false_negative,
            'False_positive_mh': false_positive_mh,
            'True_negative': true_negative,
            'False_positive_count': false_positive_count,
            'False_negative_count': fn,
            'Precision': round(precision, 1),
            'Recall': round(recall, 1),
            'F1_score': round(f1, 1)
        })

    df_summarized = pd.DataFrame(summarized_data)
    sum_name = f"{output_file}_summarized.txt"
    df_summarized.to_csv(sum_name, sep='\t', index=False)
    print(f"Saved {sum_name}")

    # -------------------------------------------------
    # 3) Compute means & std per k_col => _sn_sp.txt
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

    # -------------------------------------------------
    # 4) Create a line plot of Precision/Recall/F1 + False Positive Rate
    # -------------------------------------------------
    df_sn_sp_filt = df_sn_sp[df_sn_sp['k_col'] != 'off_target'].copy()
    # Convert k_col to numeric for sorting if possible
    df_sn_sp_filt['k_col_number'] = df_sn_sp_filt['k_col'].str.extract(r'(\d+)').astype(float)
    df_sn_sp_filt.sort_values('k_col_number', inplace=True)

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
            "Parse mSWEEP result files from multiple directories. "
            "Choose either relative-abundance or count-based filtering."
        )
    )

    parser.add_argument(
        '--input_dirs',
        type=str, 
        nargs='+', 
        required=True,
        help='One or more directories containing .txt mSWEEP outputs'
    )
    parser.add_argument(
        '--output_file',
        type=str,
        required=True,
        help='Base name for output files (e.g., results)'
    )
    parser.add_argument(
        '--filter-mode',
        type=str,
        choices=['rel_abund', 'count', 'count_rel_abund'],
        default='rel_abund',
        help=(
            "Choose filtering mode:\n"
            "  'rel_abund' => only keep PSVs with rel_abund > rel_abund_threshold.\n"
            "  'count' => only keep PSVs with (rel_abund * num_aligned) > count_threshold.\n"
            "  'count_rel_abund' => only keep PSVs where (rel_abund * num_aligned) > (num_aligned * rel_abund_threshold).\n"
            "Default is 'rel_abund'."
        )
    )
    parser.add_argument(
        '--rel-abund-threshold',
        type=float,
        default=0.01,
        help='Relative abundance threshold (default=0.01)'
    )
    parser.add_argument(
        '--count-threshold',
        type=float,
        default=10,
        help='Count threshold (default=10)'
    )

    args = parser.parse_args()

    parse_msweep_results(
        input_dirs=args.input_dirs,
        output_file=args.output_file,
        filter_mode=args.filter_mode,
        rel_abund_threshold=args.rel_abund_threshold,
        detection_threshold=args.count_threshold
    )


if __name__ == "__main__":
    main()
