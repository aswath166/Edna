"""
Squashes high-dimensional sequence embeddings down to 2D coordinates (x, y)
for interactive visualization.

Supports:
- UMAP projection (non-linear manifold preserving local & global structure)
- PCA fallback (deterministic, instantaneous linear projection)
- Joint projection of sequence reads and reference taxon centroids as reference anchors!
"""

import threading
import numpy as np
from sklearn.decomposition import PCA

PROJECTION_METHOD = "umap"  # "umap" or "pca"

_umap_ready = threading.Event()


def _pca_projection(matrix):
    n = len(matrix)
    if n < 2:
        return np.zeros((n, 2))
    components = min(2, matrix.shape[1], n)
    xy = PCA(n_components=components, random_state=42).fit_transform(matrix)
    if xy.shape[1] == 1:
        xy = np.hstack([xy, np.zeros((n, 1))])
    return xy


def _run_projection(matrix, method=PROJECTION_METHOD):
    matrix = np.asarray(matrix, dtype=float)
    n = len(matrix)
    if n < 2:
        return np.zeros((n, 2))

    xy = None
    if method == "umap" and n >= 5 and _umap_ready.is_set():
        try:
            import umap
            reducer = umap.UMAP(
                n_components=2,
                n_neighbors=max(2, min(15, n - 1)),
                min_dist=0.3,
                random_state=42,
            )
            xy = reducer.fit_transform(matrix)
        except Exception:
            xy = None

    if xy is None:
        xy = _pca_projection(matrix)

    return xy


def project_2d(embeddings, method=PROJECTION_METHOD):
    """Returns a list of [x, y] coordinates, one per embedding vector."""
    if len(embeddings) == 0:
        return []
    xy = _run_projection(embeddings, method=method)
    return [[round(float(x), 4), round(float(y), 4)] for x, y in xy]


def project_with_references(embeddings, reference_centroids=None, method=PROJECTION_METHOD):
    """
    Jointly projects sequence reads and known reference centroids into the same 2D coordinate space.
    Returns:
        read_coords: list of [x, y] for each sequence read
        ref_coords: dict {taxon_name: [x, y]} for reference anchors
    """
    n_reads = len(embeddings)
    if n_reads == 0:
        return [], {}

    ref_names = []
    ref_vecs = []
    if reference_centroids:
        for name, vec in reference_centroids.items():
            ref_names.append(name)
            ref_vecs.append(vec)

    if ref_vecs:
        combined_matrix = np.vstack([embeddings, np.array(ref_vecs)])
        all_xy = _run_projection(combined_matrix, method=method)
        read_xy = all_xy[:n_reads]
        ref_xy = all_xy[n_reads:]
        
        read_coords = [[round(float(x), 4), round(float(y), 4)] for x, y in read_xy]
        ref_coords = {
            ref_names[i]: [round(float(ref_xy[i][0]), 4), round(float(ref_xy[i][1]), 4)]
            for i in range(len(ref_names))
        }
        return read_coords, ref_coords
    else:
        read_coords = project_2d(embeddings, method=method)
        return read_coords, {}


def warm_up():
    """Compiles UMAP in a background thread at server start."""
    def _run():
        try:
            rng = np.random.default_rng(0)
            import umap
            umap.UMAP(n_components=2, n_neighbors=5, min_dist=0.3, random_state=0).fit_transform(
                rng.random((10, 16))
            )
            _umap_ready.set()
        except Exception:
            pass

    threading.Thread(target=_run, daemon=True).start()