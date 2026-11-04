"""
Ecological Biodiversity Suite (Alpha & Beta Diversity) for eDNA metabarcoding.

Computes:
- Alpha Diversity:
  * Observed Richness (S)
  * Chao1 Richness Estimator (accounts for rare singletons f1 & doubletons f2)
  * Shannon-Wiener Diversity (H')
  * Simpson Diversity (1 - D)
  * Pielou's Evenness (J')
  * Berger-Parker Dominance Index (d = N_max / N)
  * Whittaker Rank-Abundance curve
- Beta Diversity:
  * Bray-Curtis Dissimilarity between community samples
  * Jaccard Community Distance & Taxonomic Overlap
"""

import math
from collections import Counter


def compute_abundance(labels):
    """Counts sequence abundance for each taxonomic/cluster label."""
    return Counter(labels)


def species_richness(abundance_counts):
    """Observed species richness (number of distinct taxonomic units)."""
    return len(abundance_counts)


def chao1_richness(abundance_counts):
    """
    Chao1 non-parametric asymptotic richness estimator:
    Estimates true species richness including undetected rare species based on
    singletons (f1, taxa seen once) and doubletons (f2, taxa seen twice).
    """
    s_obs = len(abundance_counts)
    f1 = sum(1 for c in abundance_counts.values() if c == 1)
    f2 = sum(1 for c in abundance_counts.values() if c == 2)

    if f2 > 0:
        chao1 = s_obs + (f1 * f1) / (2.0 * f2)
    else:
        chao1 = s_obs + (f1 * max(0, f1 - 1)) / 2.0
    return round(float(chao1), 2)


def shannon_diversity(abundance_counts):
    """Shannon-Wiener entropy index (H = -sum(p * ln(p)))."""
    total = sum(abundance_counts.values())
    if total == 0:
        return 0.0
    h = 0.0
    for count in abundance_counts.values():
        if count <= 0:
            continue
        p = count / total
        h -= p * math.log(p)
    return round(float(h), 4)


def simpson_diversity(abundance_counts):
    """Gini-Simpson index (1 - D): probability that two random reads belong to different taxa."""
    total = sum(abundance_counts.values())
    if total <= 1:
        return 0.0
    numerator = sum(n * (n - 1) for n in abundance_counts.values())
    denominator = total * (total - 1)
    d = numerator / denominator
    return round(float(1.0 - d), 4)


def pielou_evenness(abundance_counts):
    """Pielou's evenness (J = H' / ln(S)). Ranges from 0 (uneven) to 1 (perfectly even)."""
    richness = len(abundance_counts)
    if richness <= 1:
        return 0.0
    h = shannon_diversity(abundance_counts)
    h_max = math.log(richness)
    return round(float(h / h_max), 4) if h_max > 0 else 0.0


def berger_parker_dominance(abundance_counts):
    """Berger-Parker dominance index: relative abundance of the single most dominant taxon."""
    total = sum(abundance_counts.values())
    if total == 0:
        return 0.0
    max_count = max(abundance_counts.values()) if abundance_counts else 0
    return round(float(max_count / total), 4)


def compute_rank_abundance(abundance_counts, taxon_names=None):
    """
    Generates Whittaker rank-abundance distribution for visualization.
    Returns sorted list of {rank, label, count, relative_abundance_pct}.
    """
    total = sum(abundance_counts.values()) or 1
    sorted_items = sorted(abundance_counts.items(), key=lambda x: x[1], reverse=True)
    whittaker = []
    for rank, (cid, count) in enumerate(sorted_items, start=1):
        name = None
        if taxon_names and cid in taxon_names:
            name = taxon_names[cid]
        elif cid < 0:
            name = f"Outlier lineage ({cid})"
        else:
            name = f"Cluster {cid}"

        whittaker.append({
            "rank": rank,
            "cluster_id": int(cid),
            "label": name,
            "count": int(count),
            "relative_abundance_pct": round((count / total) * 100.0, 2),
        })
    return whittaker


def compute_all_indices(labels, taxon_names=None):
    """Computes full suite of Alpha diversity metrics and Whittaker rank-abundance."""
    abundance = compute_abundance(labels)
    rank_abundance = compute_rank_abundance(abundance, taxon_names)
    return {
        "richness": species_richness(abundance),
        "chao1_richness": chao1_richness(abundance),
        "shannon_diversity": shannon_diversity(abundance),
        "simpson_diversity": simpson_diversity(abundance),
        "pielou_evenness": pielou_evenness(abundance),
        "berger_parker_dominance": berger_parker_dominance(abundance),
        "abundance_per_cluster": {int(k): int(v) for k, v in abundance.items()},
        "rank_abundance": rank_abundance,
    }


def compute_beta_diversity(sample_a_labels, sample_b_labels, name_a="Sample A", name_b="Sample B"):
    """
    Computes Beta diversity between two eDNA communities:
    - Bray-Curtis Dissimilarity (abundance-weighted turnover)
    - Jaccard Community Distance (presence/absence turnover)
    - Shared vs private taxa breakdown
    """
    counts_a = compute_abundance(sample_a_labels)
    counts_b = compute_abundance(sample_b_labels)

    all_taxa = set(counts_a.keys()).union(set(counts_b.keys()))
    if not all_taxa:
        return {
            "bray_curtis_dissimilarity": 0.0,
            "jaccard_distance": 0.0,
            "shared_taxa_count": 0,
            "unique_to_a_count": 0,
            "unique_to_b_count": 0,
            "sample_a_name": name_a,
            "sample_b_name": name_b,
        }

    # Bray-Curtis
    sum_min = 0
    sum_a = sum(counts_a.values())
    sum_b = sum(counts_b.values())
    for taxon in all_taxa:
        sum_min += min(counts_a.get(taxon, 0), counts_b.get(taxon, 0))

    denom = sum_a + sum_b
    bray_curtis = 1.0 - (2.0 * sum_min / denom) if denom > 0 else 0.0

    # Jaccard
    shared = set(counts_a.keys()).intersection(set(counts_b.keys()))
    jaccard_sim = len(shared) / len(all_taxa) if len(all_taxa) > 0 else 1.0
    jaccard_dist = 1.0 - jaccard_sim

    return {
        "bray_curtis_dissimilarity": round(float(bray_curtis), 4),
        "jaccard_distance": round(float(jaccard_dist), 4),
        "shared_taxa_count": len(shared),
        "unique_to_a_count": len(set(counts_a.keys()) - shared),
        "unique_to_b_count": len(set(counts_b.keys()) - shared),
        "sample_a_name": name_a,
        "sample_b_name": name_b,
        "taxa_overlap_summary": [
            {
                "taxon_id": int(t),
                "count_a": counts_a.get(t, 0),
                "count_b": counts_b.get(t, 0),
            }
            for t in sorted(list(all_taxa), key=lambda x: counts_a.get(x, 0) + counts_b.get(x, 0), reverse=True)
        ],
    }
