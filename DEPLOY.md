# Deploy Cell Portal on Render

This project is configured as a Docker-backed FastAPI website. The frontend, API, H5AD parser, and WormBase integration run in one web service.

## 1. Push the code to GitHub

Do not commit the `.h5ad` files. They are ignored by `.gitignore` because the downloaded datasets are multi-GB and should live on Render's persistent disk.

```bash
git init
git add .
git commit -m "Deploy cell portal"
```

Create a GitHub repository, then push:

```bash
git remote add origin https://github.com/YOUR_USERNAME/cell-portal.git
git branch -M main
git push -u origin main
```

## 2. Create the Render service

1. Open Render and choose **New > Blueprint**.
2. Connect the GitHub repository.
3. Render will read `render.yaml` and create a Docker web service.
4. Use a paid instance size with enough RAM for AnnData backed reads and a persistent disk. The included `render.yaml` requests a 10 GB disk at `/app/datasets`.

## 3. Upload the H5AD files to Render's disk

Upload the datasets into:

```text
/app/datasets
```

Expected filenames:

```text
cao2017.h5ad
taylor2020_cengen.h5ad
adult_wormseq_cds_baseline_post_sub.h5ad
```

After upload, restart the Render service. The sidebar's "Downloaded datasets" dropdown should list all `.h5ad` files.

## 4. Verify

Open:

```text
https://YOUR_RENDER_SERVICE.onrender.com/api/health
```

Expected response:

```json
{
  "status": "ok",
  "datasetsDir": "/app/datasets",
  "datasetCount": 3
}
```

Then open the main website:

```text
https://YOUR_RENDER_SERVICE.onrender.com/
```

## Notes

- Static-only hosts are not enough because H5AD support requires Python, FastAPI, AnnData, NumPy, Pandas, and SciPy.
- Raw `.h5ad` files are not directly served by the website. Users load them through the backend API.
- The app accepts uploaded `.h5ad` files, but large uploads may be limited by the host's request size and timeout settings.
