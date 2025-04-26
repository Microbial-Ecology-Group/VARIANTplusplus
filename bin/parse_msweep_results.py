#!/usr/bin/env python3
import argparse, gzip, glob, os, sys
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
* Applies one of five filter modes
     rel_abund | count | count_rel_abund | sub_count | sub_count_rel_abund
* Produces
  ─ <out>_summary.tsv       – one row per file **plus** one “combined” row per sample
  ─ <out>_count_matrix.tsv  – taxa × samples (combined only) absolute-count matrix
"""



# usage: python bin/parse_msweep_results.py --msweep_dir test_results/mSWEEP_results/ --reads_dir test_results/HostRemoval/NonHostFastq/ -o test_matrix --filter-mode sub_count_rel_abund --rel-thr 0.01


# ───────── helpers ───────────────────────────────────────────────────────
def count_fastq_reads(fq):
    op = gzip.open if fq.endswith(".gz") else open
    with op(fq, "rt") as fh:
        return sum(1 for _ in fh) // 4

def load_msweep(path):
    """
    Parse an mSWEEP *_abundances.txt file that looks like:

      #mSWEEP_version: …
      #num_reads:     49
      #num_aligned:   49
      #c_id   mean_theta
      1       0.92585
      …

    Returns
    -------
    n_reads       : int
    n_aligned     : int
    df_abund      : pandas.DataFrame with columns (taxon, rel_abund)
    """
    n_reads = n_aligned = None; table = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#num_reads"):
                n_reads = int(ln.split(":")[1].strip())
            elif ln.startswith("#num_aligned"):
                n_aligned = int(ln.split(":")[1].strip())
            elif ln.startswith("#"): continue
            elif ln.strip(): table.append(ln)
    if n_aligned is None:
        raise ValueError(f"{path}: missing #num_aligned")
    df = pd.read_csv(StringIO("".join(table)), sep=r"\s+", header=None,
                     names=["taxon", "rel_abund"])
    return n_reads, n_aligned, df
# ───────── core ──────────────────────────────────────────────────────────
def parse(msweep_dir, reads_dir, out_prefix,
          mode, rel_thr=0.01, count_thr=10):

    ms_files = sorted(glob.glob(os.path.join(msweep_dir,
                                             "*.msweep_abundances.txt")))
    if not ms_files:
        sys.exit(f"No mSWEEP files in {msweep_dir}")

    # map (sample, merged|unmerged) ▶ FASTQ
    fastq_lookup = {}
    for fq in glob.glob(os.path.join(reads_dir, "*.non.host.fastq*")):
        base = os.path.basename(fq)
        try:
            sample, rtype = base.replace(".non.host.fastq.gz","").rsplit("_",1)
            fastq_lookup[(sample, rtype)] = fq
        except ValueError:
            print(f"[WARN] unexpected FASTQ skipped: {base}")

    summary_rows, count_matrix = [], defaultdict(dict)
    combined_bucket = {}     # sample ▶ {"reads":int, "table":[DF,…] }

    # ── per-file pass ────────────────────────────────────────────────────
    for mfile in ms_files:
        base = os.path.basename(mfile)
        try:
            sample_id, read_type, *_ = base.split(".")
        except ValueError:
            print(f"[WARN] bad mSWEEP name {base}; skip"); continue

        fq = fastq_lookup.get((sample_id, read_type))
        if fq is None:
            print(f"[WARN] FASTQ missing for {sample_id} {read_type}"); continue

        n_reads_hdr, n_aligned_hdr, df = load_msweep(mfile)
        nonhost_reads = count_fastq_reads(fq)

        # ----- absolute counts based on **num_aligned** ------------------
        df["abs_count"] = df.rel_abund * n_aligned_hdr

        # ----- filtering -------------------------------------------------
        if   mode == "rel_abund":
            keep = df.rel_abund > rel_thr
        elif mode == "count":
            keep = df.abs_count  > count_thr
        elif mode == "count_rel_abund":
            keep = df.abs_count  > nonhost_reads*rel_thr
        elif mode == "sub_count":
            df["abs_count"] -= count_thr;              keep = df.abs_count>0
        elif mode == "sub_count_rel_abund":
            df["abs_count"] -= nonhost_reads*rel_thr;   keep = df.abs_count>0
        else:
            sys.exit(f"unknown filter mode {mode}")

        df_kept = df[keep]

        summary_rows.append({
            "sample"        : sample_id,
            "read_type"     : read_type,
            "num_reads_hdr" : n_reads_hdr,
            "num_aligned_hdr": n_aligned_hdr,
            "nonhost_reads" : nonhost_reads,
            "found_taxa"    : len(df),
            "kept_taxa"     : len(df_kept),
            "filter_mode"   : mode
        })

        # bucket for combined
        b = combined_bucket.setdefault(sample_id, {"reads":0,"table":[]})
        b["reads"] += nonhost_reads
        b["table"].append(df)          # keeps abs_count column

    # ── combined rows per sample ─────────────────────────────────────────
    for sample_id, info in combined_bucket.items():
        comb_reads = info["reads"]
        if comb_reads == 0: continue

        df_all = pd.concat(info["table"], ignore_index=True)
        df_comb = df_all.groupby("taxon", as_index=False)["abs_count"].sum()
        df_comb["rel_abund"] = df_comb.abs_count / comb_reads

        # apply same filter on combined
        if   mode == "rel_abund":
            keep = df_comb.rel_abund > rel_thr
        elif mode == "count":
            keep = df_comb.abs_count  > count_thr
        elif mode == "count_rel_abund":
            keep = df_comb.abs_count  > comb_reads*rel_thr
        elif mode == "sub_count":
            df_comb["abs_count"] -= count_thr;           keep = df_comb.abs_count>0
        elif mode == "sub_count_rel_abund":
            df_comb["abs_count"] -= comb_reads*rel_thr;  keep = df_comb.abs_count>0

        df_kept = df_comb[keep]

        summary_rows.append({
            "sample"        : sample_id,
            "read_type"     : "combined",
            "num_reads_hdr" : "",
            "num_aligned_hdr": "",
            "nonhost_reads" : comb_reads,
            "found_taxa"    : len(df_comb),
            "kept_taxa"     : len(df_kept),
            "filter_mode"   : mode
        })

        for _, r in df_kept.iterrows():
            count_matrix[r.taxon][f"{sample_id}_combined"] = round(r.abs_count,2)

    # ── outputs ──────────────────────────────────────────────────────────
    pd.DataFrame(summary_rows).to_csv(f"{out_prefix}_summary.tsv",
                                      sep="\t", index=False)

    pd.DataFrame(count_matrix).fillna(0).T\
        .sort_index().sort_index(axis=1) \
        .to_csv(f"{out_prefix}_count_matrix.tsv", sep="\t")

    print("Created:")
    print(f"  {out_prefix}_summary.tsv")
    print(f"  {out_prefix}_count_matrix.tsv")

# ───────── CLI ───────────────────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Parse & filter mSWEEP abundances; build combined matrix.")
    ap.add_argument("--msweep_dir", required=True)
    ap.add_argument("--reads_dir",  required=True)
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument("--filter-mode", default="rel_abund",
        choices=["rel_abund","count","count_rel_abund",
                 "sub_count","sub_count_rel_abund"])
    ap.add_argument("--rel-thr",   type=float, default=0.01)
    ap.add_argument("--count-thr", type=float, default=10)
    a = ap.parse_args()

    parse(a.msweep_dir, a.reads_dir, a.out_prefix,
          a.filter_mode, a.rel_thr, a.count_thr)
