import os
import uuid
import time
import io
import csv
from flask import Flask, render_template, request, jsonify, redirect, Response
from werkzeug.utils import secure_filename

import config
from utils.fasta_parser import parse_sequence_file, validate_file_extension
from ml.embedding import embed_batch
from ml.qc import compute_qc_metrics
from ml.reference import load_reference_data
from ml.novelty_detection import (
    cluster_sequences,
    analyze_novelty,
    assign_final_labels,
    compute_cluster_summaries,
)
from ml.biodiversity import compute_all_indices, compute_beta_diversity
from ml.benchmark import run_alignment_benchmark
from ml.projection import project_with_references, warm_up

app = Flask(__name__)
app.config["TEMPLATES_AUTO_RELOAD"] = True
app.config["SEND_FILE_MAX_AGE_DEFAULT"] = 0
app.config["MAX_CONTENT_LENGTH"] = config.MAX_UPLOAD_SIZE_MB * 1024 * 1024

os.makedirs(config.UPLOAD_DIR, exist_ok=True)
os.makedirs(config.REFERENCE_DIR, exist_ok=True)
os.makedirs(config.PRESETS_DIR, exist_ok=True)

# Start UMAP background JIT compiler warm up
warm_up()


@app.route("/")
def index():
    return render_template("index.html", cache_bust=int(time.time()))


@app.route("/dashboard")
def dashboard():
    return redirect("/")


@app.route("/api/presets", methods=["GET"])
def get_presets():
    """Lists available demo sample datasets."""
    presets = [
        {
            "id": "marine_coastal_edna",
            "filename": "marine_coastal_edna.fastq",
            "name": "Marine Coastal eDNA (Fish + Novel Deep Protists)",
            "description": "Metabarcoding sample with Atlantic salmon, cod, and uncharacterized dark protists.",
        },
        {
            "id": "river_estuary_edna",
            "filename": "river_estuary_edna.fastq",
            "name": "River Estuary Hydro-eDNA (Freshwater Biota)",
            "description": "Freshwater river sample with Daphnia, algae, aquatic bacteria, and novel ciliates.",
        },
    ]
    return jsonify({"presets": presets})


def _process_records(records, reference_centroids, reference_seqs, reference_lineages):
    """Core analytical pipeline for eDNA records."""
    qc_metrics = compute_qc_metrics(records)
    sequences = [r["sequence"] for r in records]
    seq_ids = [r["id"] for r in records]

    # Feature representation
    embeddings = embed_batch(sequences)

    # Density clustering
    labels, outlier_scores = cluster_sequences(
        embeddings,
        min_cluster_size=config.MIN_CLUSTER_SIZE,
        min_samples=config.MIN_SAMPLES,
    )

    # Multi-tier novelty & dark taxa divergence
    novelty = analyze_novelty(
        embeddings,
        labels,
        outlier_scores=outlier_scores,
        reference_centroids=reference_centroids,
        known_threshold=config.KNOWN_TAXON_THRESHOLD,
        sister_threshold=config.SISTER_TAXON_THRESHOLD,
    )

    final_labels = assign_final_labels(labels, novelty)
    cluster_summaries = compute_cluster_summaries(embeddings, final_labels, novelty)

    # Map cluster id to taxon name for Whittaker rank abundance
    taxon_names = {}
    for cs in cluster_summaries:
        taxon_names[cs["cluster_id"]] = cs["dominant_taxon"]

    # Ecological biodiversity indices
    alpha_indices = compute_all_indices(final_labels, taxon_names=taxon_names)

    # 2D projection with reference beacons
    read_coords, ref_coords = project_with_references(embeddings, reference_centroids)

    # Benchmark vs traditional BLAST-like alignment
    benchmark = run_alignment_benchmark(sequences, reference_seqs, novelty)

    results = []
    for i in range(len(sequences)):
        info = novelty[i]
        match = info["nearest_taxon"] if info["category"] == "known" else None
        sister_match = info["nearest_taxon"] if info["category"] == "sister_taxon" else None
        lineage = reference_lineages.get(info["nearest_taxon"], "Unspecified")

        results.append({
            "id": seq_ids[i],
            "cluster": int(final_labels[i]),
            "raw_cluster": int(labels[i]),
            "category": info["category"],
            "possibly_novel": info["novel"],
            "closest_known_taxon": match or sister_match or info["nearest_taxon"],
            "taxonomic_lineage": lineage,
            "nearest_dist": info["nearest_dist"],
            "confidence": info["confidence"],
            "status_desc": info["status_desc"],
            "outlier_score": info["outlier_score"],
            "length": records[i]["length"],
            "gc_pct": records[i]["gc_pct"],
            "x": read_coords[i][0] if i < len(read_coords) else 0.0,
            "y": read_coords[i][1] if i < len(read_coords) else 0.0,
        })

    return {
        "qc_metrics": qc_metrics,
        "biodiversity_indices": alpha_indices,
        "cluster_summaries": cluster_summaries,
        "reference_anchors": ref_coords,
        "benchmark": benchmark,
        "sequence_results": results,
    }


@app.route("/api/analyze", methods=["POST"])
def analyze():
    """
    Accepts an uploaded FASTA/FASTQ or a preset ID and executes the complete
    BioEmbed pipeline: QC -> Canonical Embeddings -> HDBSCAN -> Multi-Tier Novelty ->
    Alpha Diversity -> 2D UMAP Map -> BLAST Benchmark.
    """
    preset_id = request.form.get("preset_id")
    filepath = None

    if preset_id:
        preset_filename = f"{preset_id}.fastq"
        preset_path = os.path.join(config.PRESETS_DIR, preset_filename)
        if not os.path.isfile(preset_path):
            preset_filename = f"{preset_id}.fasta"
            preset_path = os.path.join(config.PRESETS_DIR, preset_filename)

        if not os.path.isfile(preset_path):
            return jsonify({"error": f"Preset '{preset_id}' not found"}), 404
        filepath = preset_path
    else:
        if "file" not in request.files:
            return jsonify({"error": "No file uploaded"}), 400

        uploaded_file = request.files["file"]
        if uploaded_file.filename == "":
            return jsonify({"error": "Empty filename"}), 400

        if not validate_file_extension(uploaded_file.filename, config.ALLOWED_EXTENSIONS):
            return jsonify({"error": "Unsupported file format. Please upload FASTA or FASTQ."}), 400

        filename = secure_filename(uploaded_file.filename)
        unique_name = f"{uuid.uuid4().hex}_{filename}"
        filepath = os.path.join(config.UPLOAD_DIR, unique_name)
        uploaded_file.save(filepath)

    try:
        records = parse_sequence_file(filepath)
    except Exception as e:
        return jsonify({
            "error": "Could not read this file as FASTA/FASTQ. Ensure records start with '@' or '>'.",
            "details": str(e),
        }), 400

    if not records:
        return jsonify({"error": "No valid DNA sequences found in file."}), 400

    try:
        ref_centroids, ref_seqs, ref_lineages = load_reference_data(config.REFERENCE_DIR)
        pipeline_data = _process_records(records, ref_centroids, ref_seqs, ref_lineages)

        response_payload = {
            "total_sequences": len(records),
            "reference_taxa_loaded": len(ref_centroids),
            "reference_taxa_names": list(ref_centroids.keys()),
            "qc_metrics": pipeline_data["qc_metrics"],
            "biodiversity_indices": pipeline_data["biodiversity_indices"],
            "cluster_summaries": pipeline_data["cluster_summaries"],
            "reference_anchors": pipeline_data["reference_anchors"],
            "benchmark": pipeline_data["benchmark"],
            "sequence_results": pipeline_data["sequence_results"],
        }
        return jsonify(response_payload)
    except Exception as e:
        return jsonify({"error": f"Pipeline analysis failed: {str(e)}"}), 500


@app.route("/api/compare_samples", methods=["POST"])
def compare_samples():
    """
    Computes comparative Beta diversity (Bray-Curtis dissimilarity, Jaccard distance,
    community turnover) between two eDNA samples.
    """
    try:
        file_a = request.files.get("file_a")
        file_b = request.files.get("file_b")
        preset_a = request.form.get("preset_a", "marine_coastal_edna")
        preset_b = request.form.get("preset_b", "river_estuary_edna")

        path_a = os.path.join(config.PRESETS_DIR, f"{preset_a}.fastq") if not file_a else None
        path_b = os.path.join(config.PRESETS_DIR, f"{preset_b}.fastq") if not file_b else None

        if file_a:
            path_a = os.path.join(config.UPLOAD_DIR, f"comp_a_{uuid.uuid4().hex}_{secure_filename(file_a.filename)}")
            file_a.save(path_a)
        if file_b:
            path_b = os.path.join(config.UPLOAD_DIR, f"comp_b_{uuid.uuid4().hex}_{secure_filename(file_b.filename)}")
            file_b.save(path_b)

        records_a = parse_sequence_file(path_a)
        records_b = parse_sequence_file(path_b)

        ref_centroids, ref_seqs, ref_lineages = load_reference_data(config.REFERENCE_DIR)
        data_a = _process_records(records_a, ref_centroids, ref_seqs, ref_lineages)
        data_b = _process_records(records_b, ref_centroids, ref_seqs, ref_lineages)

        labels_a = [r["cluster"] for r in data_a["sequence_results"]]
        labels_b = [r["cluster"] for r in data_b["sequence_results"]]

        name_a = file_a.filename if file_a else preset_a.replace("_", " ").title()
        name_b = file_b.filename if file_b else preset_b.replace("_", " ").title()

        beta = compute_beta_diversity(labels_a, labels_b, name_a=name_a, name_b=name_b)

        return jsonify({
            "sample_a": {
                "name": name_a,
                "reads": len(records_a),
                "richness": data_a["biodiversity_indices"]["richness"],
                "shannon": data_a["biodiversity_indices"]["shannon_diversity"],
                "simpson": data_a["biodiversity_indices"]["simpson_diversity"],
            },
            "sample_b": {
                "name": name_b,
                "reads": len(records_b),
                "richness": data_b["biodiversity_indices"]["richness"],
                "shannon": data_b["biodiversity_indices"]["shannon_diversity"],
                "simpson": data_b["biodiversity_indices"]["simpson_diversity"],
            },
            "beta_diversity": beta,
        })
    except Exception as e:
        return jsonify({"error": f"Comparative analysis failed: {str(e)}"}), 500


@app.route("/api/export_csv", methods=["POST"])
def export_csv():
    """Generates and downloads a CSV export of all analyzed read records."""
    data = request.get_json(silent=True) or {}
    results = data.get("sequence_results", [])
    if not results:
        return jsonify({"error": "No results available to export"}), 400

    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow([
        "Sequence_ID",
        "Assigned_Cluster",
        "Classification_Category",
        "Closest_Reference_Taxon",
        "Taxonomic_Lineage",
        "Divergence_Distance",
        "Confidence_Score",
        "Status_Description",
        "Read_Length_bp",
        "GC_Percentage",
        "Outlier_Score",
        "UMAP_X",
        "UMAP_Y",
    ])

    for r in results:
        writer.writerow([
            r.get("id"),
            r.get("cluster"),
            r.get("category"),
            r.get("closest_known_taxon"),
            r.get("taxonomic_lineage"),
            r.get("nearest_dist"),
            r.get("confidence"),
            r.get("status_desc"),
            r.get("length"),
            r.get("gc_pct"),
            r.get("outlier_score"),
            r.get("x"),
            r.get("y"),
        ])

    csv_data = output.getvalue()
    return Response(
        csv_data,
        mimetype="text/csv",
        headers={"Content-Disposition": "attachment;filename=edna_taxonomy_results.csv"},
    )


if __name__ == "__main__":
    app.run(debug=True, host="127.0.0.1", port=5000)
