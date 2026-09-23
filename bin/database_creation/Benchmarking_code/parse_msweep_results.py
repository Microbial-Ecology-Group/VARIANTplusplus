#!/usr/bin/env python3
"""
Parse & summarise mSWEEP abundance tables, supporting several filtering modes.

The set of samples is defined by the counts_<sample>.fastq.gz.txt files that
rename_sample_files.py writes for EVERY simulated sample, not by the
*_abundances.txt files that only exist for samples that survived Kraken2
extraction.  Samples with no Mannheimia reads (no extraction -> no themisto ->
no mSWEEP output) are emitted with status = "not_analyzed" so they still count
towards the denominator.
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


# ───────── name parsing helpers ───────────────────────────────────────
CONF_RE   = re.compile(r"_conf(\d+(?:\.\d+)?)")
COUNTS_RE = re.compile(r"^counts_(.+)\.fastq\.gz\.txt$")


def split_base(base):
    """
    '1XX0_GSVsXXk_10XXoff_target_conf0.1.k_10'
        -> ('1XX0_GSVsXXk_10XXoff_target', '0.1', 'k_10')

    Returns (clean_sample_name, kraken_conf_string, msweep_column).
    Splitting on the '_conf' token (instead of the first '.') keeps
    fractional confidences such as 0.1 intact.
    """
    m = CONF_RE.search(base)
    if not m:
        return base.split(".", 1)[0], "", ""
    clean = base[:m.start()]
    rest  = base[m.end():].lstrip(".")      # mSWEEP annotation column, if any
    return clean, m.group(1), rest


def read_counts_file(path):
    """Return the 2nd column of a counts_*.txt with '.fastq.gz' stripped."""
    try:
        with open(path) as fh:
            for ln in fh:
                cols = ln.rstrip("\n").split("\t")
                if len(cols) >= 2:
                    return cols[1].replace(".fastq.gz", "")
    except OSError:
        pass
    return ""


def parse_sample_metadata(counts_name, fallback):
    """
    counts_name looks like '1XX0_GSVsXXk_10XXoff_target.815631' or
    '3XX2_GSVsXXk_10XXGSV1.4021_GSV5.3980'.
    Returns (iter_num, numGSVs, k_col, expected_GSVs, total_input_reads).
    """
    src = counts_name if counts_name else fallback
    try:
        iter_num, numGSVs, k_col, genome_info = src.split("XX")
    except ValueError:
        try:
            iter_num, numGSVs, k_col, genome_info = fallback.split("XX")
        except ValueError:
            return "", "", "", {}, 0

    if numGSVs == "0_GSVs":
        exp_GSVs    = "off_target"
        total_input = sum(int(p.split(".")[1])
                          for p in genome_info.split("_") if "." in p)
    else:
        tmp         = {p.split(".")[0]: int(p.split(".")[1])
                       for p in genome_info.split("_") if "." in p}
        exp_GSVs    = tmp
        total_input = sum(tmp.values())
    return iter_num, numGSVs, k_col, exp_GSVs, total_input


def discover_counts(counts_dirs):
    """
    Every simulated sample, keyed by (dir, clean_sample_name).
    This is the true denominator: one record per sample that ISS actually made.
    """
    registry = {}
    for d in counts_dirs:
        hits = glob.glob(os.path.join(d, "counts_*.fastq.gz.txt"))
        for path in hits:
            m = COUNTS_RE.match(os.path.basename(path))
            if not m:
                continue
            clean       = m.group(1)
            counts_name = read_counts_file(path)
            it, ng, kc, exp, tot = parse_sample_metadata(counts_name, clean)
            registry[(d, clean)] = {
                "counts_path"  : path,
                "counts_name"  : counts_name,
                "iter_num"     : it,
                "numGSVs"      : ng,
                "k_col"        : kc,
                "expected_GSVs": exp,
                "total_input"  : tot,
            }
        print(f"[INFO] {len(hits):4} counts_*.fastq.gz.txt files in {d}")
    return registry


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
                         rel_abund_threshold=0.01, detection_threshold=10,
                         counts_dirs=None, confidences=None,
                         missing_read_types=("combined",),
                         report_missing=True):

    summarized_rows = []
    found_GSVs_all  = defaultdict(dict)
    file_bucket     = defaultdict(lambda: {"merged": None, "unmerged": None})

    # 0 · full sample universe from the counts files ----------------------
    strict_dir      = counts_dirs is None          # counts live beside results
    counts_registry = discover_counts(counts_dirs or input_dirs)
    analyzed_keys   = set()                        # (dir?, sample, conf)
    confs_seen      = set()

    def akey(srcdir, clean, conf):
        return (srcdir, clean, conf) if strict_dir else (clean, conf)

    # 1 · collect all *_abundances.txt
    all_files = []
    for d in input_dirs:
        hits = glob.glob(os.path.join(d, "*.txt"))
        all_files.extend((f, d) for f in hits)
        print(f"[INFO] {len(hits):4} txt files in {d}")
    if not all_files and not counts_registry:
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
        clean_name, conf_str, msweep_col = split_base(base)
        counts_fn = os.path.join(os.path.dirname(path),
                                 f"counts_{clean_name}.fastq.gz.txt")

        counts_name = read_counts_file(counts_fn) if os.path.exists(counts_fn) else ""
        iter_num, numGSVs, k_col, exp_GSVs, total_input = \
            parse_sample_metadata(counts_name, clean_name)

        # this sample DID reach mSWEEP at this confidence
        analyzed_keys.add(akey(srcdir, clean_name, conf_str))
        confs_seen.add(conf_str)

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
            "numGSVs"       : numGSVs,
            "expected_GSVs" : exp_GSVs
        }

        summarized_rows.append({
            "SourceDir"          : srcdir,
            "file"               : base,
            "sample"             : clean_name,
            "read_type"          : rtype,
            "iter_num"           : iter_num,
            "kraken_confidence"  : conf_str,
            "msweep_col"         : msweep_col,
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

        clean_name, conf_str, msweep_col = split_base(base)

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
            "sample"             : clean_name,
            "read_type"          : "combined",
            "iter_num"           : iter_num,
            "kraken_confidence"  : conf_str,
            "msweep_col"         : msweep_col,
            "counts_name"        : ";".join(filter(None, [
                                        found_GSVs_all.get((base,"merged"  ,srcdir),{}).get("counts_name",""),
                                        found_GSVs_all.get((base,"unmerged",srcdir),{}).get("counts_name","")])),
            "num_reads"          : (m["num_reads"]  if m else 0) + (u["num_reads"]  if u else 0),
            "num_aligned"        : total_aligned,
            "total_input_reads"  : total_input,
            "numGSVs"            : m["numGSVs"] if m else u["numGSVs"],
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

    # 3b · samples that never produced any mSWEEP output ------------------
    n_missing = 0
    if report_missing and counts_registry:
        if confidences:
            conf_list = [str(c) for c in confidences]
        else:
            conf_list = sorted(confs_seen) or [""]
            print(f"[INFO] kraken confidences inferred from filenames: {conf_list}")

        for (srcdir, clean), meta in sorted(counts_registry.items()):
            for conf in conf_list:
                if akey(srcdir, clean, conf) in analyzed_keys:
                    continue
                n_missing += 1
                for rtype in missing_read_types:
                    summarized_rows.append({
                        "SourceDir"          : srcdir,
                        "file"               : f"{clean}_conf{conf}" if conf else clean,
                        "sample"             : clean,
                        "read_type"          : rtype,
                        "iter_num"           : meta["iter_num"],
                        "kraken_confidence"  : conf,
                        "msweep_col"         : "",
                        "counts_name"        : meta["counts_name"],
                        "num_reads"          : 0,
                        "num_aligned"        : 0,
                        "total_input_reads"  : meta["total_input"],
                        "numGSVs"            : meta["numGSVs"],
                        "k_col"              : meta["k_col"],
                        "expected_GSVs"      : meta["expected_GSVs"],
                        "filtered_GSV_groups": {},
                        "found_GSV_groups"   : {},
                        "filter_mode"        : filter_mode,
                        "rel_abund_threshold": rel_abund_threshold,
                        "count_threshold"    : (
                            meta["total_input"] * rel_abund_threshold
                            if filter_mode in {"sub_count_rel_abund", "rel_abund_by_GSV"}
                            else detection_threshold),
                        "denominator_choice" : (
                            "NA" if filter_mode in {"rel_abund", "rel_abund_combined"}
                            else "num_aligned" if filter_mode == "count_rel_abund"
                            else "total_input_reads"),
                        "status"             : "not_analyzed"
                    })

    # 4 · save ------------------------------------------------------------
    out_df = pd.DataFrame(summarized_rows)
    out_df.to_csv(f"{output_file}_summarized.txt", sep="\t", index=False)

    n_expected = len(counts_registry)
    n_analyzed = len({k[-2] if strict_dir else k[0] for k in analyzed_keys})
    print(f"[INFO] samples in counts files      : {n_expected}")
    print(f"[INFO] samples with mSWEEP output   : {n_analyzed}")
    print(f"[INFO] sample x confidence dropouts : {n_missing} (status=not_analyzed)")
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
    ap.add_argument("--counts-dirs", nargs="+", default=None,
                    help="Directories holding counts_*.fastq.gz.txt "
                         "(default: same as --input_dirs)")
    ap.add_argument("--confidences", nargs="+", default=None,
                    help="Kraken2 confidence values every sample was run at, "
                         "written exactly as they appear in filenames "
                         "(e.g. 0 0.1). Default: inferred from the "
                         "*_abundances.txt that do exist.")
    ap.add_argument("--missing-read-types", nargs="+",
                    choices=["merged", "unmerged", "combined"],
                    default=["combined"],
                    help="Which read_type rows to emit for un-analysed samples")
    ap.add_argument("--no-missing", action="store_true",
                    help="Do not emit not_analyzed rows (old behaviour)")
    args = ap.parse_args()

    parse_msweep_results(args.input_dirs, args.output_file,
                         args.filter_mode,
                         args.rel_abund_threshold,
                         args.count_threshold,
                         counts_dirs=args.counts_dirs,
                         confidences=args.confidences,
                         missing_read_types=tuple(args.missing_read_types),
                         report_missing=not args.no_missing)