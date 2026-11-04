import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
UPLOAD_DIR = os.path.join(BASE_DIR, "data", "uploads")
REFERENCE_DIR = os.path.join(BASE_DIR, "data", "reference_embeddings")
PRESETS_DIR = os.path.join(BASE_DIR, "data", "sample_presets")

# "kmer" runs immediately with no extra downloads.
# "transformer" uses DNABERT-2 (requires torch + transformers installed).
EMBEDDING_MODE = "kmer"

# k-mer settings
KMER_SIZE = 4
USE_CANONICAL_KMERS = True

# HDBSCAN clustering parameters
MIN_CLUSTER_SIZE = 3
MIN_SAMPLES = 2

# Distance thresholds in normalized embedding space for multi-tier novelty:
# 1. d <= KNOWN_TAXON_THRESHOLD -> Confirmed Known Taxon
# 2. KNOWN_TAXON_THRESHOLD < d <= SISTER_TAXON_THRESHOLD -> Putative Novel Sister Lineage
# 3. d > SISTER_TAXON_THRESHOLD in compact cluster -> Candidate Novel / Dark Taxon
# 4. Outlier score > OUTLIER_FACTOR -> Artifact / Stray Read
KNOWN_TAXON_THRESHOLD = 0.35
SISTER_TAXON_THRESHOLD = 0.65
NOVELTY_DISTANCE_THRESHOLD = 0.60

ALLOWED_EXTENSIONS = {"fasta", "fa", "fna", "fastq", "fq", "txt"}
MAX_UPLOAD_SIZE_MB = 50
