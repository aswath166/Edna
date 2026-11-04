"""
Quality Control (QC) and sequence profiling module for eDNA reads.
Computes Phred score metrics (for FASTQ), length distribution,
GC% distributions, homopolymer run detection, and quality filtering diagnostics.
"""

import numpy as np


def compute_qc_metrics(records):
    """
    Analyzes a list of parsed record dicts from fasta_parser.
    Returns a structured dictionary of QC metrics and distribution histograms.
    """
    if not records:
        return {
            "total_reads": 0,
            "total_bases": 0,
            "mean_length": 0.0,
            "min_length": 0,
            "max_length": 0,
            "mean_gc": 0.0,
            "has_quality_scores": False,
            "length_distribution": {"bins": [], "counts": []},
            "gc_distribution": {"bins": [], "counts": []},
            "quality_summary": None,
        }

    lengths = np.array([r["length"] for r in records], dtype=int)
    gc_values = np.array([r["gc_pct"] for r in records], dtype=float)
    total_reads = len(records)
    total_bases = int(np.sum(lengths))

    # GC distribution (0% to 100% in 10 bins)
    gc_hist, gc_bin_edges = np.histogram(gc_values, bins=10, range=(0.0, 100.0))
    gc_bins = [f"{int(gc_bin_edges[i])}-{int(gc_bin_edges[i+1])}%" for i in range(len(gc_hist))]

    # Length distribution (up to 8 bins)
    min_len = int(np.min(lengths))
    max_len = int(np.max(lengths))
    if min_len == max_len:
        len_bins = [f"{min_len} bp"]
        len_counts = [total_reads]
    else:
        num_len_bins = min(8, max(2, max_len - min_len))
        len_hist, len_bin_edges = np.histogram(lengths, bins=num_len_bins)
        len_bins = [
            f"{int(round(len_bin_edges[i]))}-{int(round(len_bin_edges[i+1]))} bp"
            for i in range(len(len_hist))
        ]
        len_counts = [int(c) for c in len_hist]

    # Quality scores (if available from FASTQ)
    all_quals = []
    has_quals = False
    for r in records:
        if r.get("quality_scores"):
            has_quals = True
            all_quals.extend(r["quality_scores"])

    quality_summary = None
    if has_quals and all_quals:
        qual_arr = np.array(all_quals, dtype=float)
        mean_q = float(np.mean(qual_arr))
        q20_pct = float(np.mean(qual_arr >= 20.0) * 100.0)
        q30_pct = float(np.mean(qual_arr >= 30.0) * 100.0)
        
        # Q-score histogram (0 to 45 in bins of 5)
        q_hist, q_bin_edges = np.histogram(qual_arr, bins=9, range=(0.0, 45.0))
        q_bins = [f"Q{int(q_bin_edges[i])}-{int(q_bin_edges[i+1])}" for i in range(len(q_hist))]
        
        quality_summary = {
            "mean_phred_score": round(mean_q, 2),
            "q20_percentage": round(q20_pct, 2),
            "q30_percentage": round(q30_pct, 2),
            "distribution": {
                "bins": q_bins,
                "counts": [int(c) for c in q_hist],
            },
        }

    # Sequence complexity check: count reads with long homopolymers (>= 8 consecutive identical bases)
    homopolymer_count = 0
    for r in records:
        seq = r["sequence"]
        if any(base * 8 in seq for base in "ACGT"):
            homopolymer_count += 1

    return {
        "total_reads": total_reads,
        "total_bases": total_bases,
        "mean_length": round(float(np.mean(lengths)), 1),
        "median_length": int(np.median(lengths)),
        "min_length": min_len,
        "max_length": max_len,
        "mean_gc": round(float(np.mean(gc_values)), 2),
        "std_gc": round(float(np.std(gc_values)), 2),
        "homopolymer_flagged_reads": homopolymer_count,
        "has_quality_scores": has_quals,
        "length_distribution": {
            "bins": len_bins,
            "counts": len_counts,
        },
        "gc_distribution": {
            "bins": gc_bins,
            "counts": [int(c) for c in gc_hist],
        },
        "quality_summary": quality_summary,
    }


def filter_records(records, min_length=30, max_n_pct=5.0):
    """
    Filters out reads that are too short or have high ambiguous 'N' base content.
    Returns (filtered_records, removed_count).
    """
    valid = []
    removed = 0
    for r in records:
        seq = r["sequence"]
        if len(seq) < min_length:
            removed += 1
            continue
        n_count = seq.count("N")
        if (n_count / len(seq)) * 100.0 > max_n_pct:
            removed += 1
            continue
        valid.append(r)
    return valid, removed
