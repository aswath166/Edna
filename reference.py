"""
Builds and caches the reference taxa library for eDNA classification.

Reads reference FASTA files from data/reference_embeddings/.
Extracts taxon names, optional taxonomic lineages, sequence strings, and
computes reference centroid vectors in embedding space.
"""

import os
import re
import numpy as np
from Bio import SeqIO

from ml.embedding import embed_batch

_ACCESSION = re.compile(r"^[A-Za-z]{1,2}_?\d+(\.\d+)?$")
_FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".txt")

_cache = {
    "signature": None,
    "centroids": {},
    "sequences": {},
    "lineages": {},
}


def taxon_from_header(record):
    """Extracts a readable taxon name from a FASTA record header."""
    raw = record.description.strip()
    # Strip bracketed metadata such as [lineage=...] or [phylum=...]
    clean_raw = re.sub(r"\[.*?\]", "", raw).strip()

    # Check for SILVA/Greengenes style hierarchy separated by semicolons
    if ";" in clean_raw:
        parts = [p.strip() for p in clean_raw.split(";") if p.strip()]
        last_part = parts[-1].replace("_", " ")
        if len(last_part.split()) >= 1:
            return last_part

    # One-word header such as ">Salmo_salar"
    if len(clean_raw.split()) == 1:
        return clean_raw.replace("_", " ")

    # NCBI-style header: drop accession if present, take first two words
    words = clean_raw.split()
    if _ACCESSION.match(words[0]):
        words = words[1:]
    if not words:
        return record.id
    return " ".join(words[:2]).replace("_", " ")


def lineage_from_header(record):
    """Extracts taxonomic hierarchy if available in header description."""
    raw = record.description.strip()
    if ";" in raw:
        return " > ".join([p.strip().replace("_", " ") for p in raw.split(";") if p.strip()])
    # Check for bracketed taxonomy [lineage=...] or [phylum=...]
    match = re.search(r"\[lineage=([^\]]+)\]", raw, re.IGNORECASE)
    if match:
        return match.group(1).replace(";", " > ")
    return None


def _signature(reference_dir):
    if not os.path.isdir(reference_dir):
        return ()
    files = sorted(
        f for f in os.listdir(reference_dir) if f.lower().endswith(_FASTA_EXTENSIONS)
    )
    return tuple((f, os.path.getmtime(os.path.join(reference_dir, f))) for f in files)


def load_reference_data(reference_dir):
    """
    Returns (centroids, sequences_by_taxon, lineages_by_taxon).
    Cached and rebuilt only when files in reference_dir are modified.
    """
    if not os.path.isdir(reference_dir):
        return {}, {}, {}

    signature = _signature(reference_dir)
    if signature == _cache["signature"] and _cache["centroids"]:
        return _cache["centroids"], _cache["sequences"], _cache["lineages"]

    by_taxon = {}
    lineages = {}

    for filename, _ in signature:
        path = os.path.join(reference_dir, filename)
        with open(path, "r", encoding="utf-8-sig") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                seq = str(record.seq).upper().strip()
                if len(seq) < 20:
                    continue
                taxon = taxon_from_header(record)
                lineage = lineage_from_header(record)
                if lineage and taxon not in lineages:
                    lineages[taxon] = lineage
                by_taxon.setdefault(taxon, []).append(seq)

    centroids = {}
    for taxon, seqs in by_taxon.items():
        embeddings = embed_batch(seqs)
        # Average into single centroid vector and L2-normalize
        mean_vec = embeddings.mean(axis=0)
        norm = np.linalg.norm(mean_vec)
        if norm > 1e-8:
            mean_vec /= norm
        centroids[taxon] = mean_vec

    _cache["signature"] = signature
    _cache["centroids"] = centroids
    _cache["sequences"] = by_taxon
    _cache["lineages"] = lineages

    return centroids, by_taxon, lineages


def load_reference_centroids(reference_dir):
    """Backward-compatible helper returning only the centroids dict."""
    centroids, _, _ = load_reference_data(reference_dir)
    return centroids
