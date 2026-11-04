"""
Turns raw DNA sequences into numeric embedding vectors so we can cluster
and compare them — this is the core difference from a BLAST-alignment
approach.

Modes:
- kmer: Frequency vector of canonical (strand-invariant) or standard k-mers.
  Fast, deterministic, normalized, and runs completely offline with no dependencies.
- transformer: DNABERT-2 embeddings (rich contextual nucleotide representations).
  Lazy-loaded only when torch + transformers are installed.
"""

import itertools
import numpy as np

from config import EMBEDDING_MODE, KMER_SIZE, USE_CANONICAL_KMERS

_BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def _reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def _build_kmer_vocab(k=KMER_SIZE, canonical=USE_CANONICAL_KMERS):
    """Generates the unique k-mer vocabulary. Canonical merges reverse complements."""
    all_raw = ["".join(p) for p in itertools.product(_BASES, repeat=k)]
    if not canonical:
        vocab = sorted(all_raw)
    else:
        unique_canonical = set()
        for kmer in all_raw:
            rc = _reverse_complement(kmer)
            unique_canonical.add(min(kmer, rc))
        vocab = sorted(list(unique_canonical))
    index_map = {kmer: i for i, kmer in enumerate(vocab)}
    return vocab, index_map


_KMER_VOCAB, _KMER_INDEX = _build_kmer_vocab(KMER_SIZE, USE_CANONICAL_KMERS)


def _kmer_embedding(sequence, k=KMER_SIZE, canonical=USE_CANONICAL_KMERS):
    """
    Normalized k-mer frequency vector with L2 normalization for robust distance metrics.
    Strand-invariant when canonical is True.
    """
    global _KMER_VOCAB, _KMER_INDEX
    # If parameters match global cache, use cached index
    if k == KMER_SIZE and canonical == USE_CANONICAL_KMERS:
        index_map = _KMER_INDEX
        vocab_len = len(_KMER_VOCAB)
    else:
        vocab, index_map = _build_kmer_vocab(k, canonical)
        vocab_len = len(vocab)

    vec = np.zeros(vocab_len, dtype=np.float32)
    n_kmers = 0
    seq_len = len(sequence)

    for i in range(seq_len - k + 1):
        kmer = sequence[i:i + k]
        if "N" in kmer:
            continue
        if canonical:
            rc = _reverse_complement(kmer)
            kmer = min(kmer, rc)
        idx = index_map.get(kmer)
        if idx is not None:
            vec[idx] += 1.0
            n_kmers += 1

    if n_kmers > 0:
        vec /= n_kmers
        # L2-normalize vector to unit sphere for stable Euclidean & cosine clustering
        norm = np.linalg.norm(vec)
        if norm > 1e-8:
            vec /= norm

    return vec


_transformer_model = None
_transformer_tokenizer = None


def _load_transformer():
    """Lazy-loads DNABERT-2 only when transformer mode is actually invoked."""
    global _transformer_model, _transformer_tokenizer
    if _transformer_model is not None:
        return
    try:
        import torch  # noqa: F401
        from transformers import AutoModel, AutoTokenizer
    except ImportError as e:
        raise RuntimeError(
            "transformer mode requires torch and transformers to be installed. "
            "Please switch EMBEDDING_MODE to 'kmer' in config.py or install torch."
        ) from e

    model_name = "zhihan1996/DNABERT-2-117M"
    _transformer_tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    _transformer_model = AutoModel.from_pretrained(model_name, trust_remote_code=True)
    _transformer_model.eval()


def _transformer_embedding(sequence):
    import torch
    _load_transformer()
    with torch.no_grad():
        inputs = _transformer_tokenizer(sequence, return_tensors="pt")["input_ids"]
        hidden_states = _transformer_model(inputs)[0]
        # Mean-pool across sequence tokens
        embedding = hidden_states.mean(dim=1).squeeze().cpu().numpy()
        norm = np.linalg.norm(embedding)
        if norm > 1e-8:
            embedding /= norm
    return embedding


def embed_sequence(sequence, mode=None, k=KMER_SIZE, canonical=USE_CANONICAL_KMERS):
    """Embeds a single DNA sequence string into a numeric vector."""
    mode = mode or EMBEDDING_MODE
    if mode == "kmer":
        return _kmer_embedding(sequence, k=k, canonical=canonical)
    elif mode == "transformer":
        return _transformer_embedding(sequence)
    else:
        raise ValueError(f"Unknown EMBEDDING_MODE: {mode}")


def embed_batch(sequences, mode=None, k=KMER_SIZE, canonical=USE_CANONICAL_KMERS):
    """Embeds an iterable of DNA sequence strings, returns a 2D numpy array (N x D)."""
    return np.array([embed_sequence(seq, mode=mode, k=k, canonical=canonical) for seq in sequences])
