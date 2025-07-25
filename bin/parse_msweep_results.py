#!/usr/bin/env python3
import argparse, gzip, glob, os, sys, re
from io import StringIO
from collections import defaultdict
import pandas as pd

"""
Slim parser for mSWEEP abundance files.

* Accepts
    --msweep_dir   directory with  <sample>.<merged|unmerged>.msweep_abundances.txt
    --reads_dir    directory with  <sample>_<merged|unmerged>.non.host.fastq(.gz)
* Counts reads in each FASTQ and uses that number as the denominator when
  converting relative-abundance to absolute counts.
* Applies one of **six** filter modes
     rel_abund | rel_abund_combined | count | count_rel_abund | sub_count | sub_count_rel_abund | rel_abund_by_GSV
* Produces
  ─ <out>_summary.tsv       – one row per file **plus** one “combined” row per sample
      • found_GSV_groups      – {taxon: rel_abund, …} before filtering
      • filtered_GSV_groups   – {taxon: abs_count, …} after filtering
  ─ <out>_count_matrix.tsv  – taxa × samples (combined only) absolute-count matrix
"""


# ────────────────────────────────────────────────────────────────────────
# GSV-specific relative-abundance thresholds
# (edit these numbers whenever you need new cut-offs)
# ────────────────────────────────────────────────────────────────────────
#GSV_rel_map = {1: 0.00107, 2: 0.000777, 3: 0.000953, 4: 0.00293,
#               5: 0.00145, 6: 0.000600, 7: 0.000598, 8: 0.00140} # Conf 0.1, 99quantile
#GSV_rel_map = {1: 0.000639 , 2: 0.0000215,3:0.000522 ,4:0.000588 ,5:0.000204 ,6:0.000294 ,7:0.0000616 ,8:0.0007} # Conf 0.1, 99 quantile

GSV_rel_map = {1: 0.00138 , 2: 0.0000653,3:0.000743 ,4:0.000739 ,5:0.000427 ,6:0.00104 ,7:0.000283 ,8:0.000765} # Conf 0.0, 99 quantile



# ───────── helpers ──────────────────────────────────────────────────────
def count_fastq_reads(fq):
    op = gzip.open if fq.endswith(".gz") else open
    with op(fq, "rt") as fh:
        return sum(1 for _ in fh) // 4


def load_msweep(path):
    """
    Read *_abundances.txt  →  (num_reads, num_aligned, DataFrame[taxon, rel_abund])
    """
    n_reads = n_aligned = None
    table   = []
    with open(path) as fh:
        for ln in fh:
            if   ln.startswith("#num_reads"):
                n_reads = int(ln.split(":")[1].strip())
            elif ln.startswith("#num_aligned"):
                n_aligned = int(ln.split(":")[1].strip())
            elif ln.startswith("#"):
                continue
            elif ln.strip():
                table.append(ln)

    if n_aligned is None:
        raise ValueError(f"{path}: missing #num_aligned")

    df = pd.read_csv(StringIO("".join(table)), sep=r"\s+", header=None,
                     names=["taxon", "rel_abund"])
    return n_reads, n_aligned, df


def apply_filter(df, mode, rel_thr, cnt_thr, denom_reads, GSV_rel_map=None):
    """
    Return a filtered (possibly modified) copy of *df*.
    `df` must already contain columns: rel_abund, abs_count
    """
    if mode in {"rel_abund", "rel_abund_combined"}:
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

    # ── NEW: rel_abund_by_GSV ──────────────────────────────────────────
    elif mode == "rel_abund_by_GSV":
        # subtract GSV-specific floor for every row, keep if > 0
        def adjusted(row):
            floor_rel = GSV_rel_map.get(int(row.taxon), rel_thr)   # fallback → global
            return max(0, row.abs_count - denom_reads * floor_rel)

        df["abs_count"] = df.apply(adjusted, axis=1)
        keep = df.abs_count > 0
    # ------------------------------------------------------------------
    else:                       # safety net
        keep = df.rel_abund > rel_thr

    return df[keep]


# ───────── core ─────────────────────────────────────────────────────────
def parse(msweep_dir, reads_dir, out_prefix,
          mode, rel_thr=0.01, count_thr=10):

    # bookkeeping for the summary table
    denominator_choice = (
        "NA"                 if mode in {"rel_abund", "rel_abund_combined"} else
        "num_aligned"        if mode == "count_rel_abund" else
        "nonhost_reads"      if mode in {"sub_count_rel_abund", "rel_abund_by_GSV"} else
        "constant"
    )

    ms_files = sorted(glob.glob(os.path.join(msweep_dir, "*.msweep_abundances.txt")))
    if not ms_files:
        sys.exit(f"No mSWEEP files in {msweep_dir}")

    # per-file FASTQ read counts – only required for *_rel_abund* count modes
    need_reads     = mode in {"count_rel_abund", "sub_count_rel_abund", "rel_abund_by_GSV"}
    fastq_lookup   = {}
    if need_reads:
        for fq in glob.glob(os.path.join(reads_dir, "*.non.host.fastq*")):
            base = os.path.basename(fq)
            #print(base)
            try:
                sample, rtype = re.match(r"(.+?)_(merged|unmerged)$",
                                     base.replace(".non.host.fastq.gz", "")
                                    ).groups()
                #print(sample, rtype)
                fastq_lookup[(sample, rtype)] = fq
            except ValueError:
                print(f"[WARN] unexpected FASTQ skipped: {base}")

    # print(fastq_lookup) show dictionary of reads
    summary_rows, count_matrix = [], defaultdict(dict)
    combined_bucket            = {}

    # ------------------------------------------------------------------ #
    # 1 · per-file pass (immediate filtering only for classic rel_abund) #
    # ------------------------------------------------------------------ #
    for mfile in ms_files:
        base = os.path.basename(mfile)
        try:
            sample_id, read_type, *_ = base.split(".")
        except ValueError:
            print(f"[WARN] bad mSWEEP name {base}; skipped")
            continue

        # non-host read count if needed
        if need_reads:
            fq = fastq_lookup.get((sample_id, read_type))
            if fq is None:
                print(f"[WARN] FASTQ missing for {sample_id} {read_type}")
                continue
            nonhost_reads = count_fastq_reads(fq)
        else:
            nonhost_reads = 0

        n_reads_hdr, n_aligned_hdr, df = load_msweep(mfile)
        df["abs_count"] = df.rel_abund * n_aligned_hdr
        found_GSVs_full = dict(zip(df.taxon, df.rel_abund.round(6)))

        # per-file filtering (only for classic rel_abund)
        if mode == "rel_abund":
            df_kept = apply_filter(df.copy(), mode, rel_thr, count_thr,
                                   n_aligned_hdr or 1)
            kept_dict = dict(zip(df_kept.taxon, df_kept.abs_count.round(2)))
            kept_n    = len(df_kept)
        else:
            df_kept, kept_dict, kept_n = df.copy(), {}, 0

        # store summary
        summary_rows.append({
            "sample"              : sample_id,
            "read_type"           : read_type,
            "num_reads_hdr"       : n_reads_hdr,
            "num_aligned_hdr"     : n_aligned_hdr,
            "nonhost_reads"       : nonhost_reads,
            "found_taxa_n"        : len(df),
            "kept_taxa_n"         : kept_n,
            "filter_mode"         : mode,
            "found_GSV_groups"    : found_GSVs_full,
            "filtered_GSV_groups" : kept_dict,
            "denominator_choice"  : denominator_choice
        })

        bucket = combined_bucket.setdefault(
            sample_id,
            {"reads": 0, "hdr_reads": 0, "hdr_aligned": 0,
             "tables_filtered": [], "tables_raw": []}
        )
        bucket["hdr_reads"]       += n_reads_hdr or 0
        bucket["hdr_aligned"]     += n_aligned_hdr or 0
        bucket["reads"]           += nonhost_reads
        bucket["tables_filtered"].append(df_kept)
        bucket["tables_raw"].append(df)

    # ------------------------------------------------------------------ #
    # 2 · build combined rows                                            #
    # ------------------------------------------------------------------ #
    for sample_id, info in combined_bucket.items():
        denom_reads     = info["reads"]          # non-host reads
        hdr_aligned_sum = info["hdr_aligned"]    # for rel_abund calc

        tables_to_merge = (info["tables_filtered"] if mode == "rel_abund"
                           else info["tables_raw"])
        if not tables_to_merge:
            continue

        df_sum = (pd.concat(tables_to_merge, ignore_index=True)
                    .groupby("taxon", as_index=False)["abs_count"].sum())

        df_sum["rel_abund"] = df_sum.abs_count / (hdr_aligned_sum or 1)

        # post-merge filtering
        if   mode == "rel_abund_combined":
            df_sum = apply_filter(df_sum, mode, rel_thr, count_thr,
                                  hdr_aligned_sum)
        elif mode in {"count", "count_rel_abund",
                      "sub_count", "sub_count_rel_abund",
                      "rel_abund_by_GSV"}:
            if   mode == "count_rel_abund":
                denom = hdr_aligned_sum
            else:                               # sub*  OR  rel_abund_by_GSV
                denom = denom_reads
            df_sum = apply_filter(df_sum, mode, rel_thr, count_thr,
                                  denom, GSV_rel_map)

        comb_filtered = dict(zip(df_sum.taxon, df_sum.abs_count.round(2)))

        summary_rows.append({
            "sample"              : sample_id,
            "read_type"           : "combined",
            "num_reads_hdr"       : info["hdr_reads"],
            "num_aligned_hdr"     : hdr_aligned_sum,
            "nonhost_reads"       : denom_reads,
            "found_taxa_n"        : len(df_sum),
            "kept_taxa_n"         : len(comb_filtered),
            "filter_mode"         : mode,
            "found_GSV_groups"    : {},
            "filtered_GSV_groups" : comb_filtered,
            "denominator_choice"  : denominator_choice
        })

        # matrix (combined rows only)
        for taxon, cnt in comb_filtered.items():
            count_matrix[taxon][f"{sample_id}_combined"] = cnt

    # ------------------------------------------------------------------ #
    # 3 · output                                                         #
    # ------------------------------------------------------------------ #
    pd.DataFrame(summary_rows).to_csv(f"{out_prefix}_summary.tsv",
                                      sep="\t", index=False)
    pd.DataFrame(count_matrix).fillna(0).T.sort_index().sort_index(axis=1) \
        .to_csv(f"{out_prefix}_count_matrix.tsv", sep="\t")

    print("Created:")
    print(f"  {out_prefix}_summary.tsv")
    print(f"  {out_prefix}_count_matrix.tsv")


# ───────── CLI ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Parse & filter mSWEEP abundances; build combined matrix."
    )
    ap.add_argument("--msweep_dir", required=True)
    ap.add_argument("--reads_dir",  required=True)
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument("--filter-mode", default="rel_abund_by_GSV",
        choices=["rel_abund", "rel_abund_combined",
                 "count", "count_rel_abund",
                 "sub_count", "sub_count_rel_abund",
                 "rel_abund_by_GSV"])
    ap.add_argument("--rel-thr",   type=float, default=0.01)
    ap.add_argument("--count-thr", type=float, default=10)
    args = ap.parse_args()

    parse(args.msweep_dir, args.reads_dir, args.out_prefix,
          args.filter_mode, args.rel_thr, args.count_thr)