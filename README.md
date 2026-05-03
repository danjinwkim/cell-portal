# Cell Portal

A local prototype web portal for single-cell RNA-seq exploration with natural-language dataset queries.

## What it does now

- Uploads CSV, TSV, TXT, or H5AD datasets.
- CSV/TSV files run entirely in the browser and use genes as rows with cell types or cells as columns.
- H5AD files are parsed by the Python backend into a fast browser-ready analysis view.
- H5AD gene labels prefer common readable symbols from `adata.var`, such as `gene_symbols`, `gene_name`, `feature_name`, or `symbol`, before falling back to `adata.var_names`.
- WormBase-backed H5AD genes can be queried for description, expression, homology, human disease, interactions, genomic/cellular location, pathways, and phenotypes.
- Detects numeric gene-expression columns automatically.
- Computes normalized expression, a PCA-style 2D projection, lightweight clustering, cluster markers, and group differential expression.
- Lets users color the interactive Plotly PCA projection by metadata, cluster, or gene expression.
- PCA points are hoverable and clickable; selecting a point shows cell ID, cluster, cell type, top expressed genes, PC1/PC2 coordinates, and a biological interpretation.
- Adds a **Global Clustering** tab with a pairwise transcriptomic similarity heatmap, top dendrogram, neuron-class block annotations, Pearson/cosine similarity options, and linked cluster highlighting back to the PCA view.
- Adds a multimodal neuron query panel that links transcriptomics to connectome edges, 2D/3D spatial coordinates, and developmental lineage records by `neuron_id`.
- Includes a natural-language query box for common analysis tasks:
  - `show marker genes for cluster 2`
  - `plot MS4A1 expression`
  - `which cells show highest CD3D expression?`
  - `top 10 cells expressing LYZ`
  - `describe daf-2`
  - `what diseases are associated with aap-1?`
  - `what genes interact with daf-2?`
  - `compare treated vs control`
  - `how many clusters are there?`

## Run locally, CSV/TSV only

```bash
python3 -m http.server 8000
```

Then open:

```text
http://localhost:8000
```

## Run locally with H5AD support

Install backend dependencies:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Start the backend:

```bash
.venv/bin/uvicorn server:app --reload --port 8000
```

Then open:

```text
http://localhost:8000
```

## Publish as a website

This app is ready to deploy as a Docker-backed FastAPI website. Because the bundled H5AD files are multi-GB, deploy it on a host that supports a Python web service plus persistent disk, such as Render, Railway, Fly.io, or a small VPS. Static-only hosts such as GitHub Pages, Netlify static, or basic Vercel static hosting will not support the H5AD backend.

### Option A: Docker on a VPS

Build the image from this folder:

```bash
docker build -t cell-portal .
```

Run it with the local `datasets/` folder mounted into the container:

```bash
docker run --rm -p 8000:8000 \
  -v "$PWD/datasets:/app/datasets" \
  -v "$PWD/cache:/app/cache" \
  cell-portal
```

Then point your domain or reverse proxy to:

```text
http://YOUR_SERVER_IP:8000
```

### Option B: Render

This repo includes `render.yaml` for a Docker web service with a persistent disk mounted at `/app/datasets`.

1. Push this project to a GitHub repository.
2. In Render, create a new Blueprint from that repository.
3. After the service is created, upload the `.h5ad` files into the mounted `/app/datasets` disk. The app will discover them through `GET /api/datasets`.
4. Open the Render service URL and load a downloaded H5AD from the sidebar.

The container reads these deployment environment variables:

```text
DATASETS_DIR=/app/datasets
WORMBASE_CACHE_DIR=/app/cache/wormbase
MAX_BROWSER_CELLS=3000
MAX_BROWSER_GENES=250
```

Health check endpoint:

```text
/api/health
```

The backend uses backed AnnData reading for H5AD uploads, then returns up to 3,000 cells and 250 selected genes by default. It prioritizes common interpretable marker genes and any `highly_variable` gene annotations. If the file already contains `X_umap`, `X_tsne`, `X_pca`, `leiden`, `louvain`, `cell_type`, or related annotations, the portal reuses them instead of recomputing everything in the browser.

The Global Clustering view posts the already-loaded browser expression view to `/api/global-clustering`; the backend normalizes expression, computes a Pearson or cosine similarity matrix, runs SciPy hierarchical clustering, caches the result by dataset/settings, and returns Plotly-ready dendrogram and heatmap data.

## Multimodal data integration

The backend module `multimodal_store.py` defines the unified schema:

```text
neuron_id
cluster
transcriptomic_neighbors
connectome_edges: [{target, chemical_weight, gap_weight}]
spatial_coordinates: {x, y, z}
lineage_path
parent_cell
```

Optional ingestion files can be added under `multimodal_data/`:

- `neurons.csv`: `neuron_id,neuron_class,cluster`
- `connectome_edges.csv`: `source,target,chemical_weight,gap_weight`
- `spatial_coordinates.csv`: `neuron_id,x,y,z`
- `lineage.json`: records with `neuron_id`, `parent_cell`, and `lineage_path`

If those files are absent, the app uses a small built-in adult neuron demo index so the API and UI remain usable. API endpoints:

- `GET /api/multimodal/neurons`
- `GET /api/multimodal/neuron/{neuron_id}`
- `GET /api/multimodal/connectivity/{neuron_id}`
- `GET /api/multimodal/spatial/nearest/{neuron_id}`
- `GET /api/multimodal/lineage/{neuron_id}`
- `POST /api/multimodal/query`

The natural-language query engine is intentionally deterministic in this prototype: it extracts neuron IDs and maps terms such as `synaptic`, `connected`, `near`, `lineage`, and `transcriptionally similar` to the structured endpoints above. This can later be replaced by LLM function calling while preserving the same API tools.

## Downloaded datasets

Place full `.h5ad` files in `datasets/` and restart or refresh the portal. The sidebar will list them under "Downloaded datasets" so they can be loaded directly from the backend without re-uploading through the browser.

Currently downloaded:

- `datasets/cao2017.h5ad`
- `datasets/taylor2020_cengen.h5ad`
- `datasets/adult_wormseq_cds_baseline_post_sub.h5ad`

## Dataset format

The fastest input format is a CSV/TSV expression matrix where genes are rows and cell types or cells are columns:

```csv
gene,T_cells,B_cells,Monocytes,NK_cells
CD3D,12,0,1,4
IL7R,8,0,0,2
MS4A1,0,14,1,0
LYZ,1,0,16,1
```

The first column should contain gene names. Every following column becomes a cell or cell-type profile in the portal. H5AD files are not affected by this CSV/TSV orientation because AnnData already stores cells and genes separately.

## Production upgrade path

For real research use, keep this portal UI and add a Python service:

- FastAPI API for uploads, dataset metadata, analysis jobs, and natural-language tool calls.
- Scanpy/AnnData support for `.h5ad`, 10x folders, Matrix Market, Loom, and large sparse matrices.
- Celery or RQ workers for long-running QC, normalization, integration, clustering, marker detection, and differential expression.
- PostgreSQL for metadata, S3/MinIO for uploaded data, and Redis for job state.
- LLM function calling where the model can only invoke curated tools such as `plot_umap`, `find_marker_genes`, `run_differential_expression`, and `query_gene_expression`.

The important boundary is that the LLM interprets user intent, while deterministic bioinformatics functions perform the scientific computation.
