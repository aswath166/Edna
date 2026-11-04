# BioEmbed-eDNA: Environmental DNA Taxonomy & Dark Taxa Discovery Platform

**Final-Year Project / Capstone Research Platform in Bioinformatics & Machine Learning**

BioEmbed-eDNA identifies known taxa, detects putative sister lineages, flags **candidate novel / dark taxa**, filters sequencing artifacts, and computes true ecological biodiversity indices (Alpha & Beta diversity) across environmental DNA (eDNA) metabarcoding reads using continuous sequence embeddings and density clustering.

---

## 🌟 Key Research Innovations & Defendable Features

1. **Embedding-Based Sequence Representation (Canonical $k$-mers & DNABERT-2)**:
   - Projects sequences into a 136-dimensional canonical $k$-mer simplex invariant to reverse-complement strand flips, with $L_2$-norm sphere normalization.
   - Optional deep contextual nucleotide representations via **DNABERT-2**.
2. **Quality Control (QC) & Sequence Profiling Engine**:
   - Automated profiling of Phred quality scores ($Q20$, $Q30$ rates), GC content distribution, amplicon length dispersion, and homopolymer complexity flags.
3. **Multi-Tiered Biological Divergence Categorization**:
   - Classifies reads into 4 categories:
     - **Confirmed Known Taxa** ($d \le \theta_{known}$)
     - **Putative Novel Sister Lineages** ($\theta_{known} < d \le \theta_{sister}$)
     - **Candidate Dark Taxon Clusters** ($d > \theta_{sister}$)
     - **Chimeric Reads / Sequencing Artifacts** (HDBSCAN GLOSH outlier scores $> 0.75$)
4. **Comprehensive Ecological Diversity Suite**:
   - **Alpha Diversity**: Observed Richness ($S$), Chao1 Asymptotic Estimator ($S_{Chao1}$), Shannon-Wiener Entropy ($H'$), Gini-Simpson Index ($1 - D$), Pielou's Evenness ($J'$), Berger-Parker Dominance ($d$), and Whittaker Rank-Abundance curves.
   - **Beta Diversity**: Multi-sample community turnover via **Bray-Curtis Dissimilarity** and **Jaccard Community Distance**.
5. **Empirical Benchmark vs. Traditional Alignment (BLAST)**:
   - Quantifies the percentage of environmental reads discarded as "Unassigned Dark Taxa" by traditional alignment ($<90\%$ identity threshold).
   - Proves the exact number of discarded reads rescued by BioEmbed into structured candidate novel taxa.
   - Measures biodiversity underestimation ($\Delta S$ and $\Delta H'$).
6. **Interactive 2D Sonar Sequence Manifold**:
   - Joint UMAP/PCA projection plotting sequence reads alongside **Reference Anchor Beacons (★)** with hover tooltips and table synchronization.
7. **1-Click Preset Demos & Academic Thesis Dossier**:
   - Built-in demonstration datasets (Marine Coastal eDNA, River Estuary Hydro-eDNA).
   - Export full results as CSV or print publication-ready thesis summaries.

---

## 🚀 Quickstart & Setup

### Requirements
- **Python 3.10 or 3.11** (Recommended: Python 3.11)
- Windows / macOS / Linux

### Installation & Launch
```bash
# Clone or open repository folder
cd edna-project

# Run application
py -3.11 app.py
```
Then navigate to: **`http://127.0.0.1:5000`** in your web browser.

---

## 🧪 Automated Testing

Execute the comprehensive test suite verifying the QC engine, canonical embedding simplex, HDBSCAN clustering, Alpha/Beta diversity formulations, and Flask REST endpoints:

```bash
py -3.11 -m unittest tests/test_pipeline.py
```

---

## 📂 Project Architecture

```
edna-project/
├── app.py                      # Flask REST API + routes
├── config.py                   # Central hyperparameters & distance thresholds
├── requirements.txt            # Dependency specification
├── README.md                   # Project documentation & thesis guide
├── ml/
│   ├── embedding.py            # Canonical strand-invariant k-mers & DNABERT-2
│   ├── qc.py                   # FASTQ Q-score, GC%, length distribution engine
│   ├── novelty_detection.py    # HDBSCAN clustering, GLOSH outlier scoring, 4-tier taxonomy
│   ├── biodiversity.py         # Chao1, Shannon, Simpson, Pielou, Whittaker, Bray-Curtis
│   ├── benchmark.py            # Quantitative BLAST alignment vs BioEmbed evaluation
│   ├── projection.py           # UMAP & PCA 2D projection with reference anchors
│   └── reference.py            # Biological reference database loader & centroid engine
├── utils/
│   └── fasta_parser.py         # Resilient FASTA/FASTQ parser with Phred extraction
├── templates/
│   └── index.html              # Modern 7-tab scientific analytical dashboard
├── static/
│   ├── css/style.css           # Biotech theme stylesheet
│   └── js/dashboard.js         # Interactive charts, Sonar map, CSV export, preset runner
├── data/
│   ├── reference_embeddings/   # Reference marker gene library (18S, 16S, COI)
│   ├── sample_presets/         # 1-click demo datasets (Marine coastal, River estuary)
│   └── uploads/                # User-uploaded sequence files
└── tests/
    └── test_pipeline.py        # Automated test suite
```

---

## 🎓 Academic Defense Notes for Viva / Final Review

1. **Why k-mers instead of BLAST?**
   - Traditional alignment algorithms reject sequences that do not match known reference databases above stringent identity thresholds. In real-world eDNA, up to 80% of reads from understudied ecosystems fall into the "dark taxa" category. BioEmbed embeds every read into a continuous geometric manifold, grouping uncharacterized reads by biological kinship rather than discarding them.
2. **Why Canonical k-mers?**
   - DNA is double-stranded and high-throughput sequencing can read either strand randomly. Standard k-mers treat $5'\to 3'$ and $3'\to 5'$ reverse complements as different features. Canonical k-mers ($\min(k, \text{revcomp}(k))$) enforce strand invariance, reducing dimensionality from 256 to 136 ($k=4$) while eliminating strand artifact noise.
3. **What is Chao1 and why is it essential for eDNA?**
   - Observed species richness ($S_{obs}$) in environmental sequencing is always an underestimate due to limited sequencing depth. Chao1 ($S_{Chao1} = S_{obs} + \frac{f_1^2}{2f_2}$) uses singletons ($f_1$) and doubletons ($f_2$) to asymptotically estimate unobserved rare species.
