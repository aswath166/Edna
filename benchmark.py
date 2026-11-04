"""
Empirical Benchmark Module: Alignment (BLAST-like) vs Embedding-based Classification.

Provides quantitative scientific evidence for final year thesis evaluation:
1. Simulates standard sequence identity thresholding (e.g. 97% species, 90% genus).
2. Quantifies how many environmental reads are discarded as 'Unassigned/Dark Taxa' by alignment.
3. Proves how many of those discarded reads are successfully clustered into candidate novel taxa by BioEmbed.
4. Compares biodiversity underestimation (BLAST vs Embedding).
"""

from collections import Counter
from ml.biodiversity import shannon_diversity, simpson_diversity


def _approx_sequence_identity(seq1, seq2):
    """
    Computes an approximate sequence identity fraction using k-mer jaccard overlap
    and length matching. Runs instantaneously in pure Python without external BLAST binaries.
    """
    if not seq1 or not seq2:
        return 0.0
    l1, l2 = len(seq1), len(seq2)
    len_ratio = min(l1, l2) / max(l1, l2)
    if len_ratio < 0.6:
        return 0.0

    k = 4
    kmers1 = set(seq1[i:i + k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i + k] for i in range(len(seq2) - k + 1))
    if not kmers1 or not kmers2:
        return 0.0

    overlap = len(kmers1.intersection(kmers2))
    union = len(kmers1.union(kmers2))
    jaccard = overlap / union if union > 0 else 0.0
    # Jaccard to estimated identity transform: Id ~= Jaccard^(1/k) * len_ratio
    est_identity = min(1.0, (jaccard ** (1.0 / k)) * (0.85 + 0.15 * len_ratio))
    return round(float(est_identity), 4)


def run_alignment_benchmark(sequences, reference_seqs_by_taxon, embedding_results):
    """
    Compares traditional alignment classification against BioEmbed embedding results.

    reference_seqs_by_taxon: dict {taxon_name: [list of reference sequence strings]}
    embedding_results: list of dicts from analyze_novelty
    """
    n_reads = len(sequences)
    if n_reads == 0:
        return {}

    blast_results = []
    blast_assigned_taxa = []
    blast_unassigned_count = 0

    for i, seq in enumerate(sequences):
        best_taxon = None
        best_id = 0.0

        for taxon, ref_list in reference_seqs_by_taxon.items():
            for ref_seq in ref_list:
                ident = _approx_sequence_identity(seq, ref_seq)
                if ident > best_id:
                    best_id = ident
                    best_taxon = taxon

        # Standard barcoding criteria:
        # >= 97% identity: Species match
        # 90% - 97%: Genus match
        # < 90%: Discarded / Unassigned
        if best_id >= 0.90 and best_taxon:
            blast_call = best_taxon
            blast_status = "assigned_species" if best_id >= 0.97 else "assigned_genus"
            blast_assigned_taxa.append(best_taxon)
        else:
            blast_call = "Unassigned"
            blast_status = "discarded"
            blast_unassigned_count += 1

        blast_results.append({
            "index": i,
            "blast_call": blast_call,
            "blast_identity": round(best_id * 100.0, 1),
            "blast_status": blast_status,
        })

    # Alignment ecological metrics (traditionally calculated only on assigned reads)
    blast_abundance = Counter(blast_assigned_taxa)
    blast_richness = len(blast_abundance)
    blast_shannon = shannon_diversity(blast_abundance)
    blast_simpson = simpson_diversity(blast_abundance)

    # Embedding results summary
    bioembed_known = sum(1 for r in embedding_results if r["category"] == "known")
    bioembed_sister = sum(1 for r in embedding_results if r["category"] == "sister_taxon")
    bioembed_novel = sum(1 for r in embedding_results if r["category"] == "novel_taxon")
    bioembed_artifacts = sum(1 for r in embedding_results if r["category"] == "artifact_outlier")

    # Dark taxa rescue count: reads discarded by BLAST that BioEmbed structured into candidate novel clusters
    rescued_reads = 0
    for b_res, emb_res in zip(blast_results, embedding_results):
        if b_res["blast_status"] == "discarded" and emb_res["category"] in ("novel_taxon", "sister_taxon"):
            rescued_reads += 1

    return {
        "total_reads": n_reads,
        "blast_assigned_count": n_reads - blast_unassigned_count,
        "blast_assigned_pct": round(((n_reads - blast_unassigned_count) / n_reads) * 100.0, 1),
        "blast_discarded_count": blast_unassigned_count,
        "blast_discarded_pct": round((blast_unassigned_count / n_reads) * 100.0, 1),
        "blast_richness": blast_richness,
        "blast_shannon": blast_shannon,
        "blast_simpson": blast_simpson,
        "bioembed_classification": {
            "known_taxa": bioembed_known,
            "sister_taxa": bioembed_sister,
            "novel_taxa": bioembed_novel,
            "artifacts": bioembed_artifacts,
        },
        "rescued_dark_reads_count": rescued_reads,
        "rescued_dark_reads_pct": round((rescued_reads / n_reads) * 100.0, 1),
        "per_read_comparison": [
            {
                "index": i,
                "blast_call": blast_results[i]["blast_call"],
                "blast_identity": blast_results[i]["blast_identity"],
                "blast_status": blast_results[i]["blast_status"],
                "bioembed_category": embedding_results[i]["category"],
                "bioembed_nearest": embedding_results[i]["nearest_taxon"],
                "rescued": (
                    blast_results[i]["blast_status"] == "discarded"
                    and embedding_results[i]["category"] in ("novel_taxon", "sister_taxon")
                ),
            }
            for i in range(min(50, n_reads))  # Sample up to 50 for the table
        ],
    }
