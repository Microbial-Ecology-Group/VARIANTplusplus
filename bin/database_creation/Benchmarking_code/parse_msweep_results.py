#!/usr/bin/env python3
"""
Parse & summarise mSWEEP abundance tables, supporting several filtering modes.

2025-05 • adds `rel_abund_by_GSV`
"""

import os, glob, argparse, re, pandas as pd
from collections import defaultdict

# ---------------------------------------------------------------------
# editable per-GSV relative-abundance thresholds
# ---------------------------------------------------------------------
#GSV_rel_map = {1: 0.000639 , 2: 0.0000215,3:0.000522 ,4:0.000588 ,5:0.000204 ,6:0.000294 ,7:0.0000616 ,8:0.0007}          # 99 quantile COnf 0.1
#GSV_rel_map = {1: 0.00063 , 2: 0.00002,3:0.00052 ,4:0.00059 ,5:0.0002 ,6:0.00029 ,7:0.00006 ,8:0.0007}          # ← tweak as required
#GSV_rel_map = {1: 0.00107 , 2:0.000777,3:0.000953 ,4:0.00293 ,5:0.00145 ,6:0.000600 ,7:0.000598 ,8:0.00140}          # ← conf 0.1, max values
GSV_rel_map = {1: 0.00138 , 2: 0.0000653,3:0.000743 ,4:0.000739 ,5:0.000427 ,6:0.00104 ,7:0.000283 ,8:0.000765} # Conf 0.0, 99 quantile




# ───────── helpers ────────────────────────────────────────────────────
def apply_filter(df, mode, rel_thr, cnt_thr, denom_reads, GSV_rel_map=None):
    """
    Return a filtered (and possibly modified) copy of *df*.

    df must contain columns: rel_abund, abs_count  (and taxon if mode needs it)
    """
    if mode == "rel_abund":
        keep = df.rel_abund > rel_thr

    elif mode == "count":
        keep = df.abs_count > cnt_thr

    elif mode == "count_rel_abund":
        keep = df.abs_count > denom_reads * rel_thr

    elif mode == "sub_count":
        df["abs_count"] -= cnt_thr
        keep = df.abs_count > 0

    elif mode == "sub_count_rel_abund":
        df["abs_count"] -= denom_reads * rel_thr
        keep = df.abs_count > 0

    elif mode == "rel_abund_combined":
        keep = df.rel_abund > rel_thr

    # --- NEW MODE -----------------------------------------------------
    elif mode == "rel_abund_by_GSV":
        # for each row subtract GSV-specific threshold
        def adjust(row):
            thr_rel = GSV_rel_map.get(int(row.taxon), rel_thr)
            return max(0, row.abs_count - denom_reads * thr_rel)

        df["abs_count"] = df.apply(adjust, axis=1)
        keep = df.abs_count > 0
    # -----------------------------------------------------------------
    else:                        # safety net
        keep = df.rel_abund > rel_thr

    return df[keep]


# ───────── main routine ───────────────────────────────────────────────
def parse_msweep_results(input_dirs, output_file, filter_mode,
                         rel_abund_threshold=0.01, detection_threshold=10):

    summarized_rows = []
    found_GSVs_all  = defaultdict(dict)
    file_bucket     = defaultdict(lambda: {"merged": None, "unmerged": None})

    # 1 · collect all *_abundances.txt
    all_files = []
    for d in input_dirs:
        hits = glob.glob(os.path.join(d, "*.txt"))
        all_files.extend((f, d) for f in hits)
        print(f"[INFO] {len(hits):4} txt files in {d}")
    if not all_files:
        raise SystemExit("[ERROR] no .txt files found.")

    # 2 · per-file parsing ------------------------------------------------
    for path, srcdir in all_files:
        if path.endswith(".merged_abundances.txt"):
            rtype = "merged"
            base  = os.path.basename(path).replace(
                        ".msweep_output.merged_abundances.txt", "")
        elif path.endswith(".unmerged_abundances.txt"):
            rtype = "unmerged"
            base  = os.path.basename(path).replace(
                        ".msweep_output.unmerged_abundances.txt", "")
        else:
            continue

        # meta …
        temp_name  = base.split(".", 1)[0]
        clean_name = re.sub(r"_conf\d+(\.\d+)?", "", temp_name)
        counts_fn  = os.path.join(os.path.dirname(path),
                                  f"counts_{clean_name}.fastq.gz.txt")

        counts_name = ""
        if os.path.exists(counts_fn):
            with open(counts_fn) as fh:
                for ln in fh:
                    cols = ln.rstrip().split("\t")
                    if len(cols) >= 2:
                        counts_name = cols[1].replace(".fastq.gz", "")
                        break
        try:
            iter_num, numGSVs, k_col, genome_info = counts_name.split("XX")
        except Exception:                                   # fallback
            iter_num, numGSVs, k_col, genome_info = base.split("XX")

        # expected GSVs / total input reads
        if numGSVs == "0_GSVs":
            exp_GSVs     = "off_target"
            total_input  = sum(int(p.split(".")[1]) for p in genome_info.split("_") if "." in p)
        else:
            tmp          = {p.split(".")[0]: int(p.split(".")[1])
                            for p in genome_info.split("_") if "." in p}
            exp_GSVs     = tmp
            total_input  = sum(tmp.values())

        # read abundance table
        n_reads = n_aligned = None
        rows    = []
        with open(path) as fh:
            for ln in fh:
                ln = ln.rstrip()
                if ln.startswith("#num_reads"):
                    n_reads = int(ln.split("\t")[1])
                elif ln.startswith("#num_aligned"):
                    n_aligned = int(ln.split("\t")[1])
                elif ln and not ln.startswith("#"):
                    GSV, rel = ln.split("\t")
                    rows.append((GSV, float(rel)))

        df = pd.DataFrame(rows, columns=["taxon", "rel_abund"])
        df["abs_count"] = df.rel_abund * (n_aligned or 0)
        found_GSVs_dict = dict(zip(df.taxon, df.rel_abund.round(6)))

        # per-file filtering only for classic rel_abund mode
        if filter_mode == "rel_abund":
            df_kept  = apply_filter(df.copy(), filter_mode,
                                    rel_abund_threshold, detection_threshold,
                                    n_aligned or 1)
            kept_dic = dict(zip(df_kept.taxon, df_kept.abs_count.round(2)))
        else:
            df_kept, kept_dic = df.copy(), {}

        status_flag = "passed" if kept_dic else "filtered"

        found_GSVs_all[(base, rtype, srcdir)] = {
            "found_GSV_groups": found_GSVs_dict,
            "counts_name"     : counts_name
        }
        file_bucket[(base, iter_num, srcdir)][rtype] = {
            "df_raw"        : df,
            "df_kept"       : df_kept,
            "kept_dict"     : kept_dic,
            "status"        : status_flag,
            "num_reads"     : n_reads,
            "num_aligned"   : n_aligned,
            "total_input"   : total_input,
            "k_col"         : k_col,
            "expected_GSVs" : exp_GSVs
        }

        summarized_rows.append({
            "SourceDir"          : srcdir,
            "file"               : base,
            "read_type"          : rtype,
            "iter_num"           : iter_num,
            "counts_name"        : counts_name,
            "num_reads"          : n_reads,
            "num_aligned"        : n_aligned,
            "total_input_reads"  : total_input,
            "numGSVs"            : numGSVs,
            "k_col"              : k_col,
            "expected_GSVs"      : exp_GSVs,
            "filtered_GSV_groups": kept_dic,
            "found_GSV_groups"   : found_GSVs_dict,
            "filter_mode"        : filter_mode,
            "rel_abund_threshold": rel_abund_threshold,
            "count_threshold"    : detection_threshold,
            "status"             : status_flag
        })

    # 3 · build combined rows --------------------------------------------
    for (base, iter_num, srcdir), sub in file_bucket.items():
        m, u = sub["merged"], sub["unmerged"]
        if m is None and u is None:
            continue

        # choose tables to merge
        dfs_merge = ([x["df_kept"] for x in (m, u) if x] if filter_mode == "rel_abund"
                     else [x["df_raw"]  for x in (m, u) if x])
        if not dfs_merge:
            continue

        df_sum = (pd.concat(dfs_merge, ignore_index=True)
                    .groupby("taxon", as_index=False)["abs_count"].sum())

        total_aligned = (m["num_aligned"] if m else 0) + (u["num_aligned"] if u else 0)
        total_input   =  m["total_input"] if m else (u["total_input"] if u else 0)

        denom_ra = total_aligned or 1
        df_sum["rel_abund"] = df_sum.abs_count / denom_ra

        # POST-MERGE FILTERING -------------------------------------------
        if   filter_mode == "rel_abund_combined":
            df_sum = apply_filter(df_sum, filter_mode,
                                  rel_abund_threshold, detection_threshold,
                                  denom_ra)
        elif filter_mode in {"count", "count_rel_abund",
                             "sub_count", "sub_count_rel_abund",
                             "rel_abund_by_GSV"}:
            if   filter_mode == "count_rel_abund":
                denom = total_aligned
            elif filter_mode in {"sub_count_rel_abund", "rel_abund_by_GSV"}:
                denom = total_input
            else:
                denom = 1
            df_sum = apply_filter(df_sum, filter_mode,
                                  rel_abund_threshold, detection_threshold,
                                  denom, GSV_rel_map)

        comb_dic   = dict(zip(df_sum.taxon, df_sum.abs_count.round(2)))
        status_flag = "passed" if comb_dic else "filtered"

        summarized_rows.append({
            "SourceDir"          : srcdir,
            "file"               : base,
            "read_type"          : "combined",
            "iter_num"           : iter_num,
            "counts_name"        : ";".join(filter(None, [
                                        found_GSVs_all.get((base,"merged"  ,srcdir),{}).get("counts_name",""),
                                        found_GSVs_all.get((base,"unmerged",srcdir),{}).get("counts_name","")])),
            "num_reads"          : (m["num_reads"]  if m else 0) + (u["num_reads"]  if u else 0),
            "num_aligned"        : total_aligned,
            "total_input_reads"  : total_input,
            "numGSVs"            : sub["merged"]["k_col"] if m else u["k_col"],
            "k_col"              : m["k_col"] if m else u["k_col"],
            "expected_GSVs"      : m["expected_GSVs"] if m else u["expected_GSVs"],
            "filtered_GSV_groups": comb_dic,
            "found_GSV_groups"   : {},
            "filter_mode"        : filter_mode,
            "rel_abund_threshold": rel_abund_threshold,
            "count_threshold"    : (
                total_input * rel_abund_threshold
                if filter_mode in {"sub_count_rel_abund", "rel_abund_by_GSV"}
                else detection_threshold),
            "denominator_choice" : (
                "NA" if filter_mode in {"rel_abund", "rel_abund_combined"}
                else "num_aligned" if filter_mode == "count_rel_abund"
                else "total_input_reads"),
            "status"             : status_flag
        })

    # 4 · save ------------------------------------------------------------
    out_df = pd.DataFrame(summarized_rows)
    out_df.to_csv(f"{output_file}_summarized.txt", sep="\t", index=False)
    print(f"[INFO] Saved {output_file}_summarized.txt")


# ───────── CLI ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Summarise mSWEEP results with multiple filtering modes."
    )
    ap.add_argument("--input_dirs", nargs="+", required=True,
                    help="One or more directories containing *_abundances.txt")
    ap.add_argument("--output_file", required=True,
                    help="Output file prefix (no extension)")
    ap.add_argument("--filter-mode",
                    choices=["rel_abund", "rel_abund_combined",
                             "count", "count_rel_abund",
                             "sub_count_rel_abund", "sub_count",
                             "rel_abund_by_GSV"],
                    default="rel_abund_by_GSV")
    ap.add_argument("--rel-abund-threshold", type=float, default=0.01,
                    help="Global relative-abundance threshold")
    ap.add_argument("--count-threshold",     type=float, default=10,
                    help="Absolute-count threshold for count-based modes")
    args = ap.parse_args()

    parse_msweep_results(args.input_dirs, args.output_file,
                         args.filter_mode,
                         args.rel_abund_threshold,
                         args.count_threshold)
