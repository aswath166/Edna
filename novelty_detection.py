"""
Multi-Tiered Novelty & Taxon Divergence Detection Engine.
Clusters sequence embeddings with HDBSCAN and analyzes divergence against
known reference centroids.

Categorizes reads into 4 biological classes:
1. 'known': High-confidence match to a reference taxon (d <= KNOWN_TAXON_THRESHOLD).
2. 'sister_taxon': Putative novel sister lineage / divergent species in the same family/genus
   (KNOWN_TAXON_THRESHOLD < d <= SISTER_TAXON_THRESHOLD).
3. 'novel_taxon': Candidate uncharacterized 'dark taxon' forming a coherent cluster far
   from all cataloged references (d > SISTER_TAXON_THRESHOLD).
4. 'artifact_outlier': Isolated read with high GLOSH outlier score or noise label (-1),
   representing potential chimeric reads, sequence artifacts, or rare stray reads.
"""

import numpy as np
from sklearn.cluster import HDBSCAN

from config import (
    MIN_CLUSTER_SIZE,
    MIN_SAMPLES,
    KNOWN_TAXON_THRESHOLD,
    SISTER_TAXON_THRESHOLD,
)

OUTLIER_SCORE_THRESHOLD = 0.75


def cluster_sequences(embeddings, min_cluster_size=MIN_CLUSTER_SIZE, min_samples=MIN_SAMPLES):
    """
    Fits HDBSCAN on the embedding matrix.
    Returns:
        labels: np.ndarray of cluster IDs (-1 = noise)
        outlier_scores: np.ndarray of GLOSH outlier probabilities [0.0, 1.0]
    """
    n = len(embeddings)
    if n < min_cluster_size:
        return np.array([-1] * n), np.zeros(n, dtype=float)

    # In scikit-learn >= 1.3, HDBSCAN is built-in
    clusterer = HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        metric="euclidean",
        copy=True,
    )
    labels = clusterer.fit_predict(embeddings)
    outlier_scores = getattr(clusterer, "outlier_scores_", np.zeros(n, dtype=float))
    return labels, outlier_scores


def _distance_to_cluster_mates(embeddings, labels):
    """Computes distance of each clustered point to leave-one-out cluster centroid."""
    n = len(embeddings)
    dists = np.full(n, np.nan)
    for i in range(n):
        if labels[i] == -1:
            continue
        mates = [j for j in range(n) if labels[j] == labels[i] and j != i]
        if not mates:
            continue
        centroid = embeddings[mates].mean(axis=0)
        dists[i] = np.linalg.norm(embeddings[i] - centroid)
    return dists


def analyze_novelty(
    embeddings,
    labels,
    outlier_scores=None,
    reference_centroids=None,
    known_threshold=KNOWN_TAXON_THRESHOLD,
    sister_threshold=SISTER_TAXON_THRESHOLD,
):
    """
    Evaluates every read against cluster structure and reference centroids.
    Returns:
        results: list of dicts with detailed classification metadata.
    """
    embeddings = np.asarray(embeddings, dtype=float)
    n = len(embeddings)
    if outlier_scores is None:
        outlier_scores = np.zeros(n, dtype=float)

    dists = _distance_to_cluster_mates(embeddings, labels)
    valid_dists = dists[~np.isnan(dists)]
    typical_dist = float(np.median(valid_dists)) if len(valid_dists) > 0 else 0.15

    results = []
    for i in range(n):
        vec = embeddings[i]
        label = int(labels[i])
        score = float(outlier_scores[i]) if i < len(outlier_scores) else 0.0

        # Outlier flag from HDBSCAN noise label or high GLOSH score
        is_noise = (label == -1)
        is_statistical_outlier = (score >= OUTLIER_SCORE_THRESHOLD)
        far_from_mates = (
            not np.isnan(dists[i])
            and typical_dist > 0
            and dists[i] > typical_dist * 2.8
        )
        is_artifact = bool(is_noise or is_statistical_outlier or far_from_mates)

        nearest_taxon = None
        nearest_dist = None
        ref_matches = []

        if reference_centroids and len(reference_centroids) > 0:
            for name, centroid in reference_centroids.items():
                d = float(np.linalg.norm(vec - centroid))
                ref_matches.append((name, d))
            ref_matches.sort(key=lambda x: x[1])
            nearest_taxon, nearest_dist = ref_matches[0]

        # Determine 4-tier classification category
        if is_artifact:
            category = "artifact_outlier"
            novel = True
            confidence = round(max(0.2, 1.0 - score), 3)
            status_desc = "Stray read or sequencing artifact (isolated/noise)"
        elif nearest_dist is not None:
            if nearest_dist <= known_threshold:
                category = "known"
                novel = False
                confidence = round(max(0.0, 1.0 - (nearest_dist / known_threshold)), 3)
                status_desc = f"Identified as {nearest_taxon} (high confidence)"
            elif nearest_dist <= sister_threshold:
                category = "sister_taxon"
                novel = True
                confidence = round(
                    max(0.0, 1.0 - ((nearest_dist - known_threshold) / (sister_threshold - known_threshold))),
                    3,
                )
                status_desc = f"Putative novel sister lineage of {nearest_taxon}"
            else:
                category = "novel_taxon"
                novel = True
                confidence = round(min(1.0, nearest_dist / sister_threshold), 3)
                status_desc = f"Candidate dark taxon cluster (closest: {nearest_taxon}, d={nearest_dist:.2f})"
        else:
            # No reference centroids loaded
            if is_artifact:
                category = "artifact_outlier"
                novel = True
                confidence = 0.5
                status_desc = "Isolated outlier read"
            else:
                category = "novel_taxon"
                novel = True
                confidence = 0.85
                status_desc = f"Uncharacterized cluster {label} (no reference DB)"

        results.append({
            "index": i,
            "cluster": label,
            "category": category,
            "novel": novel,
            "possibly_novel": novel,  # backward compatibility
            "is_outlier": is_artifact,
            "outlier_score": round(score, 3),
            "nearest_taxon": nearest_taxon,
            "nearest_dist": round(nearest_dist, 4) if nearest_dist is not None else None,
            "confidence": confidence,
            "status_desc": status_desc,
        })

    return results


def assign_final_labels(labels, novelty):
    """
    Generates unified operational taxonomic unit (OTU) labels:
    - Coherent clusters (known or novel taxa) keep non-negative IDs (0, 1, 2, ...).
    - Outliers / artifacts are assigned unique negative IDs (-1, -2, -3, ...)
      so they count as distinct single-read groups in diversity estimators.
    """
    final = []
    next_outlier = -1
    for label, info in zip(labels, novelty):
        if info["is_outlier"] or info["category"] == "artifact_outlier":
            final.append(next_outlier)
            next_outlier -= 1
        else:
            final.append(int(label))
    return final


def compute_cluster_summaries(embeddings, final_labels, novelty):
    """
    Computes centroid, size, and consensus taxonomic call for each cluster.
    """
    embeddings = np.asarray(embeddings, dtype=float)
    clusters = {}
    for i, cid in enumerate(final_labels):
        clusters.setdefault(cid, []).append(i)

    summaries = []
    for cid, indices in clusters.items():
        sub_vecs = embeddings[indices]
        centroid = sub_vecs.mean(axis=0)
        # Dispersion: average distance to cluster centroid
        dispersion = float(np.mean([np.linalg.norm(v - centroid) for v in sub_vecs])) if len(indices) > 1 else 0.0

        # Consensus taxon call
        calls = [novelty[idx]["nearest_taxon"] for idx in indices if novelty[idx]["nearest_taxon"]]
        dominant_taxon = max(set(calls), key=calls.count) if calls else "Uncharacterized"
        dominant_cat = max(
            set(novelty[idx]["category"] for idx in indices),
            key=lambda c: sum(1 for idx in indices if novelty[idx]["category"] == c),
        )

        summaries.append({
            "cluster_id": int(cid),
            "size": len(indices),
            "dispersion": round(dispersion, 4),
            "dominant_taxon": dominant_taxon,
            "dominant_category": dominant_cat,
            "is_novel": dominant_cat in ("novel_taxon", "sister_taxon", "artifact_outlier"),
        })

    summaries.sort(key=lambda x: x["size"], reverse=True)
    return summaries
