from __future__ import annotations

import math
import os
import json
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Optional
from urllib.error import URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

import anndata as ad
import numpy as np
import pandas as pd
from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform
from scipy import sparse


ROOT = Path(__file__).parent
DATASETS_DIR = Path(os.environ.get("DATASETS_DIR", ROOT / "datasets"))
WORMBASE_CACHE_DIR = Path(os.environ.get("WORMBASE_CACHE_DIR", ROOT / "cache" / "wormbase"))
MAX_BROWSER_CELLS = int(os.environ.get("MAX_BROWSER_CELLS", "3000"))
MAX_BROWSER_GENES = int(os.environ.get("MAX_BROWSER_GENES", "250"))
GENE_SYMBOL_COLUMNS = (
    "gene_symbols",
    "gene_symbol",
    "gene_name",
    "gene_names",
    "gene_short_name",
    "feature_name",
    "feature_names",
    "symbol",
    "name",
)
CLUSTER_COLUMNS = ("leiden", "louvain", "cluster", "clusters", "cell_type", "celltype", "annotation")
HVG_COLUMNS = ("highly_variable", "highly_variable_features")
COMMON_MARKER_GENES = (
    "aap-1",
    "act-1",
    "ama-1",
    "ced-1",
    "daf-2",
    "daf-16",
    "elt-2",
    "hlh-1",
    "myo-2",
    "myo-3",
    "rab-3",
    "rgef-1",
    "unc-17",
    "unc-25",
    "unc-47",
    "vit-2",
    "CD3D",
    "CD3E",
    "IL7R",
    "CCR7",
    "MS4A1",
    "CD79A",
    "CD74",
    "LYZ",
    "S100A8",
    "S100A9",
    "FCGR3A",
    "NKG7",
    "GNLY",
    "GZMB",
    "PPBP",
    "MKI67",
    "TOP2A",
    "EPCAM",
    "KRT8",
    "COL1A1",
    "PECAM1",
    "VWF",
)

app = FastAPI(title="Cell Portal")
GLOBAL_CLUSTER_CACHE: dict[str, dict[str, Any]] = {}


@app.get("/api/health")
async def health() -> dict[str, Any]:
    return {
        "status": "ok",
        "datasetsDir": str(DATASETS_DIR),
        "datasetCount": len(list(DATASETS_DIR.glob("*.h5ad"))) if DATASETS_DIR.exists() else 0,
    }


@app.post("/api/upload-h5ad")
async def upload_h5ad(file: UploadFile = File(...)) -> dict[str, Any]:
    if not file.filename or not file.filename.lower().endswith(".h5ad"):
        raise HTTPException(status_code=400, detail="Upload a .h5ad AnnData file.")

    with tempfile.NamedTemporaryFile(suffix=".h5ad") as tmp:
        content = await file.read()
        tmp.write(content)
        tmp.flush()

        try:
            return h5ad_to_portal_payload(Path(tmp.name), file.filename)
        except Exception as exc:  # noqa: BLE001 - preserve parser details for the UI
            raise HTTPException(status_code=400, detail=f"Could not read H5AD file: {exc}") from exc


@app.get("/api/datasets")
async def list_datasets() -> list[dict[str, Any]]:
    DATASETS_DIR.mkdir(exist_ok=True)
    datasets = []
    for path in sorted(DATASETS_DIR.glob("*.h5ad")):
        datasets.append(
            {
                "id": path.name,
                "label": dataset_label(path.name),
                "sizeBytes": path.stat().st_size,
            }
        )
    return datasets


@app.get("/api/datasets/{dataset_id}")
async def load_dataset(dataset_id: str) -> dict[str, Any]:
    path = DATASETS_DIR / dataset_id
    if path.parent != DATASETS_DIR or not path.exists() or path.suffix.lower() != ".h5ad":
        raise HTTPException(status_code=404, detail="Downloaded dataset not found.")
    try:
        return h5ad_to_portal_payload(path, dataset_label(path.name))
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(status_code=400, detail=f"Could not read H5AD file: {exc}") from exc


@app.get("/api/wormbase/gene/{gene_name}")
async def wormbase_gene(gene_name: str, gene_id: Optional[str] = None) -> dict[str, Any]:
    resolved_id = gene_id if gene_id and gene_id.startswith("WBGene") else None
    if not resolved_id:
        raise HTTPException(
            status_code=404,
            detail="WormBase lookup needs a WBGene identifier from the loaded H5AD dataset.",
        )
    try:
        return wormbase_summary(gene_name, resolved_id)
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(status_code=502, detail=f"WormBase lookup failed: {exc}") from exc


@app.post("/api/global-clustering")
async def global_clustering(payload: dict[str, Any]) -> dict[str, Any]:
    """Compute a pairwise neuron similarity heatmap and dendrogram from the loaded expression view."""
    rows = payload.get("rows") or []
    genes = payload.get("genes") or []
    labels = payload.get("labels") or []
    classes = payload.get("classes") or []
    metric = str(payload.get("metric") or "pearson").lower()
    method = str(payload.get("method") or "average").lower()
    cache_key = str(payload.get("cacheKey") or "")

    if cache_key and cache_key in GLOBAL_CLUSTER_CACHE:
        return GLOBAL_CLUSTER_CACHE[cache_key]
    if not rows or not genes:
        raise HTTPException(status_code=400, detail="Global clustering needs loaded cells and numeric gene columns.")
    if len(rows) > 600:
        raise HTTPException(status_code=400, detail="Global clustering is limited to 600 displayed cells for interactive use.")
    if metric not in {"pearson", "cosine"}:
        raise HTTPException(status_code=400, detail="Metric must be pearson or cosine.")
    if method not in {"average", "complete", "single", "ward"}:
        raise HTTPException(status_code=400, detail="Linkage method must be average, complete, single, or ward.")

    matrix = np.asarray([[float(row.get(gene) or 0.0) for gene in genes] for row in rows], dtype=np.float64)
    matrix = normalize_expression_matrix(matrix)
    labels = [clean_label(label) or f"cell_{index + 1}" for index, label in enumerate(labels)]
    classes = [clean_label(value) or "Unannotated" for value in classes]
    result = compute_global_clustering(matrix, labels, classes, metric, method)

    if cache_key:
        GLOBAL_CLUSTER_CACHE[cache_key] = result
        if len(GLOBAL_CLUSTER_CACHE) > 12:
            GLOBAL_CLUSTER_CACHE.pop(next(iter(GLOBAL_CLUSTER_CACHE)))
    return result


def h5ad_to_portal_payload(path: Path, label: str) -> dict[str, Any]:
    adata = ad.read_h5ad(path, backed="r")
    try:
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise HTTPException(status_code=400, detail="The AnnData object has no cells or genes.")
        cell_idx = select_cell_indices(adata.n_obs)
        all_gene_names = readable_gene_names(adata)
        gene_idx = select_gene_indices_fast(adata, all_gene_names)
        matrix = subset_matrix(adata, cell_idx, gene_idx)

        obs = adata.obs.iloc[cell_idx].reset_index(names="cell_id")
        obs = obs.replace({np.nan: ""})
        obs_columns = [str(column) for column in obs.columns]
        gene_names = [all_gene_names[index] for index in gene_idx]
        gene_metadata = gene_metadata_for_indices(adata, gene_names, gene_idx)

        rows: list[dict[str, Any]] = []
        for row_index in range(len(cell_idx)):
            item: dict[str, Any] = {}
            for column in obs_columns:
                item[column] = clean_value(obs.iloc[row_index][column])
            for gene_index, gene in enumerate(gene_names):
                item[gene] = round(float(matrix[row_index, gene_index]), 4)
            rows.append(item)

        response = {
            "label": label,
            "columns": obs_columns + gene_names,
            "metadataColumns": obs_columns,
            "geneColumns": gene_names,
            "geneMetadata": gene_metadata,
            "rows": rows,
            "projection": existing_projection(adata, cell_idx),
            "clusters": existing_clusters(adata, cell_idx),
            "source": {
                "cellsInFile": int(adata.n_obs),
                "genesInFile": int(adata.n_vars),
                "cellsReturned": int(len(cell_idx)),
                "genesReturned": int(len(gene_idx)),
                "geneNameSource": gene_name_source(adata),
                "mode": "fast_browser_view",
            },
        }
        return response
    finally:
        adata.file.close()


def normalize_expression_matrix(matrix: np.ndarray) -> np.ndarray:
    totals = matrix.clip(min=0).sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    normalized = np.log1p((matrix.clip(min=0) / totals) * 10000)
    normalized -= normalized.mean(axis=0, keepdims=True)
    std = normalized.std(axis=0, keepdims=True)
    std[std == 0] = 1
    return normalized / std


def compute_global_clustering(
    matrix: np.ndarray,
    labels: list[str],
    classes: list[str],
    metric: str,
    method: str,
) -> dict[str, Any]:
    similarity = similarity_matrix(matrix, metric)
    distance = np.clip(1 - similarity, 0, 2)
    np.fill_diagonal(distance, 0)
    condensed = squareform(distance, checks=False)
    linkage_method = "average" if method == "ward" and metric != "pearson" else method
    tree = linkage(condensed, method=linkage_method)
    tree_info = dendrogram(tree, no_plot=True)
    order = [int(index) for index in tree_info["leaves"]]
    ordered_similarity = similarity[np.ix_(order, order)]
    ordered_labels = [labels[index] for index in order]
    ordered_classes = [classes[index] for index in order]
    flat_clusters = [int(value) for value in fcluster(tree, t=min(12, len(labels)), criterion="maxclust")]
    ordered_cluster_ids = [flat_clusters[index] for index in order]

    return {
        "metric": metric,
        "method": linkage_method,
        "labels": ordered_labels,
        "classes": ordered_classes,
        "order": order,
        "similarity": round_matrix(ordered_similarity),
        "dendrogram": {
            "icoord": tree_info["icoord"],
            "dcoord": tree_info["dcoord"],
        },
        "clusters": ordered_cluster_ids,
        "blocks": class_blocks(ordered_classes, ordered_cluster_ids),
        "summary": {
            "cells": len(labels),
            "genes": int(matrix.shape[1]),
            "clusterCount": len(set(flat_clusters)),
        },
    }


def similarity_matrix(matrix: np.ndarray, metric: str) -> np.ndarray:
    if metric == "cosine":
        norms = np.linalg.norm(matrix, axis=1, keepdims=True)
        norms[norms == 0] = 1
        normalized = matrix / norms
        similarity = normalized @ normalized.T
    else:
        similarity = np.corrcoef(matrix)
    similarity = np.nan_to_num(similarity, nan=0.0, posinf=1.0, neginf=-1.0)
    similarity = np.clip(similarity, -1, 1)
    np.fill_diagonal(similarity, 1)
    return similarity


def round_matrix(matrix: np.ndarray) -> list[list[float]]:
    return [[round(float(value), 4) for value in row] for row in matrix]


def class_blocks(classes: list[str], cluster_ids: list[int]) -> list[dict[str, Any]]:
    blocks = []
    if not classes:
        return blocks
    start = 0
    for index in range(1, len(classes) + 1):
        if index == len(classes) or classes[index] != classes[start]:
            blocks.append(
                {
                    "start": start,
                    "end": index - 1,
                    "label": classes[start],
                    "cluster": cluster_ids[start],
                    "count": index - start,
                }
            )
            start = index
    return blocks


def dataset_label(filename: str) -> str:
    labels = {
        "cao2017.h5ad": "Cao et al. 2017 L2 sci-RNA-seq",
        "taylor2020_cengen.h5ad": "Taylor et al. 2021 CeNGEN neurons",
        "adult_wormseq_cds_baseline_post_sub.h5ad": "Adult C. elegans atlas / WormSeq",
    }
    return labels.get(filename, filename)


def gene_metadata_for_indices(adata: ad.AnnData, gene_names: list[str], gene_idx: np.ndarray) -> dict[str, dict[str, str]]:
    metadata: dict[str, dict[str, str]] = {}
    for name, index in zip(gene_names, gene_idx):
        var_row = adata.var.iloc[int(index)]
        gene_id = first_present(var_row, ("gene_id", "id")) or str(adata.var_names[int(index)])
        description = first_present(var_row, ("gene_description", "description", "brief_description"))
        metadata[name] = {
            "id": gene_id,
            "description": description or "",
        }
    return metadata


def first_present(row: pd.Series, columns: tuple[str, ...]) -> str:
    for column in columns:
        if column in row and not pd.isna(row[column]) and str(row[column]).strip():
            return str(row[column]).strip()
    return ""


def wormbase_summary(gene_name: str, gene_id: str) -> dict[str, Any]:
    WORMBASE_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = WORMBASE_CACHE_DIR / f"{gene_id}.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text())

    widget_names = ("overview", "gene_ontology", "expression", "interactions", "human_diseases", "location", "phenotypes")
    with ThreadPoolExecutor(max_workers=len(widget_names)) as executor:
        widgets = dict(zip(widget_names, executor.map(lambda name: fetch_wormbase_widget(gene_id, name), widget_names)))
    summary = summarize_wormbase_widgets(gene_name, gene_id, widgets)
    cache_path.write_text(json.dumps(summary, indent=2))
    return summary


def fetch_wormbase_widget(gene_id: str, widget: str) -> dict[str, Any]:
    url = f"https://rest.wormbase.org/rest/widget/gene/{quote(gene_id)}/{widget}"
    request = Request(url, headers={"Accept": "application/json", "User-Agent": "CellPortal/0.1"})
    try:
        with urlopen(request, timeout=12) as response:
            return json.loads(response.read().decode("utf-8"))
    except (TimeoutError, URLError, json.JSONDecodeError) as exc:
        return {"error": str(exc), "fields": {}}


def summarize_wormbase_widgets(gene_name: str, gene_id: str, widgets: dict[str, dict[str, Any]]) -> dict[str, Any]:
    overview_fields = widgets["overview"].get("fields", {})
    go_fields = widgets["gene_ontology"].get("fields", {})
    expression_fields = widgets["expression"].get("fields", {})
    interaction_fields = widgets["interactions"].get("fields", {})
    disease_fields = widgets["human_diseases"].get("fields", {})
    location_fields = widgets["location"].get("fields", {})
    phenotype_fields = widgets["phenotypes"].get("fields", {})

    go_summary = go_fields.get("gene_ontology_summary", {}).get("data", {}) or {}
    expression_labels = collect_labels(expression_fields.get("anatomy_expression", {}).get("data"))[:20]
    expression_notes = collect_texts(expression_fields.get("fourd_expression_movies", {}).get("data"))[:5]
    interactions = summarize_interactions(interaction_fields)
    diseases = collect_diseases(disease_fields)

    summary = {
        "gene": clean_label(overview_fields.get("name", {}).get("data", {}).get("label")) or gene_name,
        "id": gene_id,
        "description": field_text(overview_fields, "concise_description")
        or field_text(overview_fields, "legacy_manual_description"),
        "manualDescription": field_text(overview_fields, "legacy_manual_description"),
        "humanDiseaseRelevance": disease_relevance(overview_fields, disease_fields),
        "humanDiseases": diseases[:20],
        "homology": infer_homology(field_text(overview_fields, "concise_description")),
        "molecularFunction": go_terms(go_summary, "Molecular_function"),
        "pathways": go_terms(go_summary, "Biological_process"),
        "location": {
            "genomic": genomic_location(location_fields),
            "cellular": go_terms(go_summary, "Cellular_component"),
        },
        "expression": {
            "tissues": expression_labels,
            "notes": expression_notes,
        },
        "interactingGenes": interactions[:20],
        "phenotypes": collect_phenotypes(phenotype_fields)[:20],
        "source": f"https://wormbase.org/species/c_elegans/gene/{gene_id}",
    }
    return summary


def field_text(fields: dict[str, Any], name: str) -> str:
    data = fields.get(name, {}).get("data")
    if isinstance(data, dict):
        return str(data.get("text") or "").strip()
    if isinstance(data, list):
        return " ".join(str(item.get("text", "")) for item in data if isinstance(item, dict)).strip()
    return str(data or "").strip()


def disease_relevance(overview_fields: dict[str, Any], disease_fields: dict[str, Any]) -> str:
    return field_text(overview_fields, "human_disease_relevance") or field_text(disease_fields, "human_disease_relevance")


def go_terms(go_summary: dict[str, Any], key: str) -> list[str]:
    terms = []
    for item in go_summary.get(key, []) or []:
        descriptions = item.get("term_description") or []
        if isinstance(descriptions, dict):
            descriptions = [descriptions]
        for description in descriptions:
            label = clean_label(description.get("label") if isinstance(description, dict) else description)
            if label:
                terms.append(label)
    return dedupe_strings(terms)[:20]


def summarize_interactions(fields: dict[str, Any]) -> list[str]:
    values = []
    for item in fields.get("interactor_types", {}).get("data", []) or []:
        label = clean_label((item.get("interactor") or {}).get("label"))
        types = ", ".join(item.get("types", []) or [])
        if label:
            values.append(f"{label} ({types})" if types else label)
    edges = fields.get("interactions", {}).get("data", {}).get("edges", []) or []
    for edge in edges[:20]:
        affected = clean_label((edge.get("affected") or {}).get("label"))
        effector = clean_label((edge.get("effector") or {}).get("label"))
        other = affected if affected != effector else ""
        interaction_type = clean_label(edge.get("type"))
        if other:
            values.append(f"{other} ({interaction_type})" if interaction_type else other)
    return dedupe_strings(values)


def collect_diseases(fields: dict[str, Any]) -> list[str]:
    data = fields.get("human_diseases", {}).get("data", {}) or {}
    values = []
    for key in ("potential_model", "experimental_model", "disease_relevance"):
        values.extend(clean_label(item.get("label")) for item in data.get(key, []) or [] if isinstance(item, dict))
    return dedupe_strings([value for value in values if value])


def collect_phenotypes(fields: dict[str, Any]) -> list[str]:
    return dedupe_strings(collect_labels(fields)[:100])


def collect_labels(value: Any) -> list[str]:
    labels = []
    if isinstance(value, dict):
        label = clean_label(value.get("label"))
        item_class = value.get("class")
        if label and item_class in {"phenotype", "anatomy_term", "go_term", "do_term"}:
            labels.append(label)
        for child in value.values():
            labels.extend(collect_labels(child))
    elif isinstance(value, list):
        for item in value:
            labels.extend(collect_labels(item))
    return dedupe_strings(labels)


def collect_texts(value: Any) -> list[str]:
    texts = []
    if isinstance(value, dict):
        for key in ("details", "description"):
            text = clean_label(value.get(key))
            if text:
                texts.append(text)
        for child in value.values():
            texts.extend(collect_texts(child))
    elif isinstance(value, list):
        for item in value:
            texts.extend(collect_texts(item))
    return dedupe_strings(texts)


def genomic_location(fields: dict[str, Any]) -> str:
    data = fields.get("genomic_position", {}).get("data", []) or []
    if not data:
        return ""
    item = data[0]
    return clean_label(item.get("label") or item.get("pos_string"))


def infer_homology(description: str) -> list[str]:
    if not description:
        return []
    markers = ("Is an ortholog of", "Human ortholog")
    snippets = []
    for marker in markers:
        if marker in description:
            snippet = description[description.find(marker) :].split(".")[0]
            snippets.append(snippet)
    return dedupe_strings(snippets)


def clean_label(value: Any) -> str:
    if value is None:
        return ""
    return " ".join(str(value).split())


def dedupe_strings(values: list[str]) -> list[str]:
    seen = set()
    result = []
    for value in values:
        if not value or value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


@app.get("/")
async def index() -> FileResponse:
    return FileResponse(ROOT / "index.html")


@app.get("/app.js")
async def app_js() -> FileResponse:
    return FileResponse(ROOT / "app.js", media_type="application/javascript")


@app.get("/styles.css")
async def styles_css() -> FileResponse:
    return FileResponse(ROOT / "styles.css", media_type="text/css")


@app.get("/vendor/plotly-2.35.2.min.js")
async def plotly_js() -> FileResponse:
    return FileResponse(ROOT / "vendor" / "plotly-2.35.2.min.js", media_type="application/javascript")


def select_cell_indices(total: int) -> np.ndarray:
    if total <= MAX_BROWSER_CELLS:
        return np.arange(total)
    return np.arange(MAX_BROWSER_CELLS)


def select_gene_indices_fast(adata: ad.AnnData, gene_names: list[str]) -> np.ndarray:
    selected: list[int] = []
    gene_lookup = {name.upper(): index for index, name in enumerate(gene_names)}

    for marker in COMMON_MARKER_GENES:
        index = gene_lookup.get(marker.upper())
        if index is not None:
            selected.append(index)

    for column in HVG_COLUMNS:
        if column in adata.var.columns:
            flags = np.asarray(adata.var[column]).astype(bool)
            selected.extend(np.flatnonzero(flags).tolist())
            break

    selected.extend(range(adata.n_vars))
    return np.asarray(dedupe(selected)[:MAX_BROWSER_GENES], dtype=int)


def dedupe(values: list[int]) -> list[int]:
    seen: set[int] = set()
    result: list[int] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        result.append(value)
    return result


def existing_projection(adata: ad.AnnData, cell_idx: np.ndarray) -> list[dict[str, float]]:
    for key in ("X_umap", "UMAP", "X_tsne", "X_pca", "PCA", "Aligned"):
        if key not in adata.obsm:
            continue
        coords = np.asarray(adata.obsm[key][cell_idx])
        if coords.ndim == 2 and coords.shape[1] >= 2:
            return [{"x": round(float(row[0]), 4), "y": round(float(row[1]), 4)} for row in coords]
    return []


def existing_clusters(adata: ad.AnnData, cell_idx: np.ndarray) -> list[str]:
    for column in CLUSTER_COLUMNS:
        if column in adata.obs.columns:
            values = adata.obs.iloc[cell_idx][column]
            return [str(value) for value in values]
    return []


def subset_matrix(adata: ad.AnnData, cell_idx: np.ndarray, gene_idx: np.ndarray) -> np.ndarray:
    if len(cell_idx) and np.array_equal(cell_idx, np.arange(len(cell_idx))):
        return to_dense(adata.X[: len(cell_idx), gene_idx])
    cell_view = adata[cell_idx, :].to_memory()
    return to_dense(cell_view.X[:, gene_idx])


def readable_gene_names(adata: ad.AnnData) -> list[str]:
    source = gene_name_source(adata)
    if source == "var_names":
        raw_names = [str(gene) for gene in adata.var_names]
    else:
        raw_names = [clean_gene_name(value, fallback) for value, fallback in zip(adata.var[source], adata.var_names)]
    return make_unique(raw_names)


def gene_name_source(adata: ad.AnnData) -> str:
    for column in GENE_SYMBOL_COLUMNS:
        if column not in adata.var.columns:
            continue
        values = adata.var[column].dropna().astype(str).str.strip()
        if not values.empty and (values != "").any():
            return column
    return "var_names"


def clean_gene_name(value: Any, fallback: Any) -> str:
    if pd.isna(value) or str(value).strip() == "":
        return str(fallback)
    return str(value).strip()


def make_unique(names: list[str]) -> list[str]:
    counts: dict[str, int] = {}
    unique_names: list[str] = []
    for name in names:
        counts[name] = counts.get(name, 0) + 1
        if counts[name] == 1:
            unique_names.append(name)
        else:
            unique_names.append(f"{name}_{counts[name]}")
    return unique_names


def to_dense(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        return matrix.toarray()
    return np.asarray(matrix)


def clean_value(value: Any) -> Any:
    if pd.isna(value):
        return ""
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        number = float(value)
        if math.isfinite(number):
            return round(number, 4)
        return ""
    return str(value)
