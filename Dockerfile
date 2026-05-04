FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    DATASETS_DIR=/app/datasets \
    WORMBASE_CACHE_DIR=/app/cache/wormbase \
    MAX_BROWSER_CELLS=3000 \
    MAX_BROWSER_GENES=250

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends curl \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY server.py multimodal_store.py multimodal_ingest.py index.html app.js styles.css README.md ./
COPY vendor ./vendor

RUN mkdir -p /app/datasets /app/cache/wormbase

EXPOSE 8000

CMD ["sh", "-c", "uvicorn server:app --host 0.0.0.0 --port ${PORT:-8000}"]
