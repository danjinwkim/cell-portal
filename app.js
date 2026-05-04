const state = {
  fileName: "",
  rows: [],
  columns: [],
  genes: [],
  metadata: [],
  cellIdColumn: "",
  groupColumn: "",
  projection: [],
  clusters: [],
  serverProjection: [],
  serverClusters: [],
  geneMetadata: {},
  selectedCellIndex: null,
  activeView: "pca",
  globalClustering: null,
  globalHighlightIndices: new Set(),
  multimodalNeurons: [],
  multimodalHighlightIds: new Set(),
  markers: new Map(),
};

const palette = ["#1f7a8c", "#d88742", "#ba4a68", "#4d8d62", "#7353ba", "#247ba0", "#c44536", "#5b8e7d"];
const numericThreshold = 0.9;
const maxGlobalClusterCells = 350;
const maxGlobalClusterGenes = 120;
const globalClusterCache = new Map();
const deployedApiBase = "https://celegans-single-cell-rnaseq-explorer.vercel.app";
const apiBase = location.protocol === "file:" ? deployedApiBase : "";
const apiUrl = (path) => `${apiBase}${path}`;

const $ = (id) => document.getElementById(id);

const controls = {
  datasetFile: $("datasetFile"),
  serverDataset: $("serverDataset"),
  loadServerDataset: $("loadServerDataset"),
  loadExample: $("loadExample"),
  cellIdColumn: $("cellIdColumn"),
  groupColumn: $("groupColumn"),
  colorBy: $("colorBy"),
  geneOverlay: $("geneOverlay"),
  runAnalysis: $("runAnalysis"),
  pcaTab: $("pcaTab"),
  globalTab: $("globalTab"),
  globalMetric: $("globalMetric"),
  globalLinkage: $("globalLinkage"),
  runGlobalClustering: $("runGlobalClustering"),
  multimodalForm: $("multimodalForm"),
  multimodalInput: $("multimodalInput"),
  markerCluster: $("markerCluster"),
  groupA: $("groupA"),
  groupB: $("groupB"),
  chatForm: $("chatForm"),
  chatInput: $("chatInput"),
};

controls.datasetFile.addEventListener("change", async (event) => {
  const file = event.target.files[0];
  if (!file) return;
  if (file.name.toLowerCase().endsWith(".h5ad")) {
    await uploadH5ad(file);
    return;
  }
  const text = await file.text();
  loadDataset(text, file.name);
});

controls.loadExample.addEventListener("click", () => {
  loadDataset(makeExampleDataset(), "example_pbmc_like.csv");
});

controls.loadServerDataset.addEventListener("click", async () => {
  const datasetId = controls.serverDataset.value;
  if (!datasetId) return;
  await loadServerDataset(datasetId);
});

controls.runAnalysis.addEventListener("click", runAnalysis);
controls.markerCluster.addEventListener("change", renderMarkerTable);
controls.groupA.addEventListener("change", renderDifferentialTable);
controls.groupB.addEventListener("change", renderDifferentialTable);
controls.colorBy.addEventListener("change", schedulePlotRender);
controls.geneOverlay.addEventListener("change", schedulePlotRender);
controls.pcaTab.addEventListener("click", () => setActiveView("pca"));
controls.globalTab.addEventListener("click", () => setActiveView("global"));
controls.globalMetric.addEventListener("change", runGlobalClustering);
controls.globalLinkage.addEventListener("change", runGlobalClustering);
controls.runGlobalClustering.addEventListener("click", runGlobalClustering);
if (controls.multimodalForm) {
  controls.multimodalForm.addEventListener("submit", async (event) => {
    event.preventDefault();
    await runMultimodalQuery(controls.multimodalInput.value.trim());
  });
}

controls.chatForm.addEventListener("submit", async (event) => {
  event.preventDefault();
  const message = controls.chatInput.value.trim();
  if (!message) return;
  addChat("user", message);
  controls.chatInput.value = "";
  addChat("assistant", await answerQuestion(message));
});

function loadDataset(text, fileName) {
  const parsed = parseTable(text);
  if (!parsed.rows.length) {
    addChat("assistant", "I could not find expression values in that file. Try a CSV or TSV with genes as rows and cell types as columns.");
    return;
  }

  loadParsedDataset(parsed, fileName);
}

function loadParsedDataset(parsed, fileName) {
  state.fileName = fileName;
  state.rows = parsed.rows;
  state.columns = parsed.columns;
  state.serverProjection = parsed.projection || [];
  state.serverClusters = parsed.clusters || [];
  state.geneMetadata = parsed.geneMetadata || {};
  state.globalClustering = null;
  state.globalHighlightIndices = new Set();
  if (parsed.geneColumns && parsed.metadataColumns) {
    state.genes = parsed.geneColumns;
    state.metadata = parsed.metadataColumns;
  } else {
    classifyColumns();
  }
  setDefaultColumns();
  populateSelectors();
  runAnalysis();
}

function setActiveView(view) {
  state.activeView = view;
  $("pcaView").classList.toggle("active", view === "pca");
  $("globalView").classList.toggle("active", view === "global");
  controls.pcaTab.classList.toggle("active", view === "pca");
  controls.globalTab.classList.toggle("active", view === "global");
  controls.pcaTab.setAttribute("aria-selected", String(view === "pca"));
  controls.globalTab.setAttribute("aria-selected", String(view === "global"));
  if (view === "pca") schedulePlotRender();
  if (view === "global") scheduleGlobalClusteringRender();
}

async function uploadH5ad(file) {
  $("analysisStatus").textContent = "Reading H5AD";
  addChat("assistant", `Uploading ${file.name}. I will load a fast browser-ready H5AD view with readable gene names.`);
  const formData = new FormData();
  formData.append("file", file);

  try {
    const response = await fetch(apiUrl("/api/upload-h5ad"), {
      method: "POST",
      body: formData,
    });

    if (!response.ok) {
      const detail = await response.text();
      throw new Error(detail || `Upload failed with status ${response.status}`);
    }

    const parsed = await response.json();
    loadParsedDataset(parsed, file.name);
  } catch (error) {
    const message = explainH5adUploadError(error);
    $("analysisStatus").textContent = message.status;
    addChat("assistant", message.text);
  }
}

function explainH5adUploadError(error) {
  const detail = String(error?.message || "");
  if (/no tunnel here/i.test(detail)) {
    return {
      status: "Tunnel expired",
      text:
        "That H5AD request went to an expired temporary tunnel, so the upload never reached Cell Portal. " +
        "Use the permanent Render URL or the current localhost page, not an old lhr.life/loca.lt tunnel link.",
    };
  }
  if (location.hostname.endsWith("onrender.com")) {
    return {
      status: "Hosted H5AD upload failed",
      text:
        "The Render free deployment is live, but it has no persistent dataset disk and large H5AD uploads may exceed free-tier limits. " +
        `For full H5AD datasets, use the local app or a paid Render/VPS deployment with datasets mounted on the server. Details: ${detail}`,
    };
  }
  if (location.protocol === "file:") {
    return {
      status: "Use public link",
      text:
        "This page was opened as a local file, so I routed the upload to the public Cell Portal backend. " +
        "If the upload still fails, open https://celegans-single-cell-rnaseq-explorer.vercel.app directly and try again. " +
        `Details: ${detail}`,
    };
  }
  if (location.hostname.endsWith("vercel.app")) {
    return {
      status: "Hosted H5AD upload failed",
      text:
        "The public Vercel backend is available for H5AD parsing, but very large H5AD files can exceed serverless upload or execution limits. " +
        `Try a smaller processed H5AD view or run the local backend for full-size files. Details: ${detail}`,
    };
  }
  return {
    status: "H5AD backend unavailable",
    text:
      "I can accept .h5ad files when the Python backend is running. Start it with the README instructions, then try again. " +
      `Details: ${detail}`,
  };
}

async function refreshServerDatasets() {
  try {
    const response = await fetch(apiUrl("/api/datasets"));
    if (!response.ok) throw new Error("Dataset registry is unavailable");
    const datasets = await response.json();
    controls.serverDataset.innerHTML = "";
    if (!datasets.length) {
      fillSelect(controls.serverDataset, ["No downloaded H5AD files"], "No downloaded H5AD files");
      controls.loadServerDataset.disabled = true;
      return;
    }
    controls.loadServerDataset.disabled = false;
    for (const dataset of datasets) {
      const option = document.createElement("option");
      option.value = dataset.id;
      option.textContent = dataset.label;
      controls.serverDataset.appendChild(option);
    }
  } catch (error) {
    const label = location.protocol === "file:" ? "Open public website for H5AD" : "Backend not ready";
    fillSelect(controls.serverDataset, [label], label);
    controls.loadServerDataset.disabled = true;
    console.warn("Dataset backend unavailable", error);
  }
}

async function loadServerDataset(datasetId) {
  $("analysisStatus").textContent = "Loading downloaded H5AD";
  addChat("assistant", "Loading the server-side H5AD view directly. This skips the browser re-upload step.");
  try {
    const response = await fetch(apiUrl(`/api/datasets/${encodeURIComponent(datasetId)}`));
    if (!response.ok) {
      const detail = await response.text();
      throw new Error(detail || `Load failed with status ${response.status}`);
    }
    const parsed = await response.json();
    loadParsedDataset(parsed, parsed.label || datasetId);
  } catch (error) {
    $("analysisStatus").textContent = "Dataset load failed";
    addChat("assistant", `I could not load that downloaded dataset. Details: ${error.message}`);
  }
}

function parseTable(text) {
  const firstLine = text.split(/\r?\n/).find((line) => line.trim().length);
  const delimiter = firstLine && firstLine.includes("\t") ? "\t" : ",";
  const lines = text.trim().split(/\r?\n/).filter(Boolean);
  const columns = splitLine(lines[0], delimiter).map((value) => value.trim());
  const rawRows = lines.slice(1).map((line) => splitLine(line, delimiter));
  if (isGeneRowMatrix(columns, rawRows)) {
    return transposeGeneRowMatrix(columns, rawRows);
  }
  const rows = rawRows.map((values) => Object.fromEntries(columns.map((column, index) => [column, coerce(values[index])])));
  return { columns, rows };
}

function isGeneRowMatrix(columns, rawRows) {
  if (columns.length < 2 || !rawRows.length) return false;
  const firstHeader = columns[0].toLowerCase();
  const firstColumnLooksLikeGene = /gene|symbol|feature|transcript|id/.test(firstHeader);
  const numericValueCount = rawRows.flatMap((row) => row.slice(1)).filter((value) => Number.isFinite(Number(String(value).trim()))).length;
  const totalValueCount = rawRows.reduce((sum, row) => sum + Math.max(0, row.length - 1), 0);
  const numericShare = totalValueCount ? numericValueCount / totalValueCount : 0;
  const firstValuesAreNames = rawRows.filter((row) => row[0] && !Number.isFinite(Number(String(row[0]).trim()))).length / rawRows.length;
  return firstHeader === "" || firstColumnLooksLikeGene || (numericShare >= numericThreshold && firstValuesAreNames >= numericThreshold);
}

function transposeGeneRowMatrix(columns, rawRows) {
  const cellColumns = columns.slice(1).filter(Boolean);
  const genes = makeUniqueNames(rawRows.map((row, index) => String(row[0] || `gene_${index + 1}`).trim()));
  const rows = cellColumns.map((cellName, cellIndex) => {
    const row = {
      cell_id: cellName,
      cell_type: cellName,
    };
    genes.forEach((gene, geneIndex) => {
      row[gene] = coerce(rawRows[geneIndex][cellIndex + 1]);
    });
    return row;
  });
  return {
    columns: ["cell_id", "cell_type", ...genes],
    metadataColumns: ["cell_id", "cell_type"],
    geneColumns: genes,
    rows,
  };
}

function makeUniqueNames(names) {
  const counts = new Map();
  return names.map((name) => {
    const cleanName = name || "unnamed_gene";
    const count = (counts.get(cleanName) || 0) + 1;
    counts.set(cleanName, count);
    return count === 1 ? cleanName : `${cleanName}_${count}`;
  });
}

function splitLine(line, delimiter) {
  const values = [];
  let current = "";
  let quoted = false;
  for (let i = 0; i < line.length; i += 1) {
    const char = line[i];
    if (char === '"') quoted = !quoted;
    else if (char === delimiter && !quoted) {
      values.push(current);
      current = "";
    } else current += char;
  }
  values.push(current);
  return values;
}

function coerce(value = "") {
  const trimmed = String(value).trim();
  if (trimmed === "") return "";
  const number = Number(trimmed);
  return Number.isFinite(number) ? number : trimmed;
}

async function responseErrorMessage(response, fallback) {
  const text = await response.text();
  if (!text) return `${fallback} with status ${response.status}`;
  try {
    const parsed = JSON.parse(text);
    return parsed.detail || parsed.message || text;
  } catch {
    return text;
  }
}

function classifyColumns() {
  state.genes = [];
  state.metadata = [];
  for (const column of state.columns) {
    const numericShare = state.rows.filter((row) => typeof row[column] === "number").length / state.rows.length;
    if (numericShare >= numericThreshold) state.genes.push(column);
    else state.metadata.push(column);
  }
}

function setDefaultColumns() {
  state.cellIdColumn = state.columns.find((column) => /cell|barcode|id/i.test(column)) || state.columns[0];
  state.groupColumn =
    state.metadata.find((column) => /condition|sample|group|type|cluster/i.test(column) && column !== state.cellIdColumn) ||
    state.metadata.find((column) => column !== state.cellIdColumn) ||
    "cluster";
}

function populateSelectors() {
  fillSelect(controls.cellIdColumn, state.columns, state.cellIdColumn);
  fillSelect(controls.groupColumn, [...state.metadata, "cluster"].filter(unique), state.groupColumn);
  fillSelect(controls.colorBy, [...state.metadata, "cluster"].filter(unique), "cluster");
  fillSelect(controls.geneOverlay, ["none", ...state.genes], state.genes[0] || "none");
}

function fillSelect(select, options, selected) {
  select.innerHTML = "";
  for (const option of options) {
    const element = document.createElement("option");
    element.value = option;
    element.textContent = option;
    element.selected = option === selected;
    select.appendChild(element);
  }
}

function unique(value, index, values) {
  return value && values.indexOf(value) === index;
}

function runAnalysis() {
  if (!state.rows.length || !state.genes.length) {
    $("analysisStatus").textContent = "Waiting for expression matrix";
    return;
  }

  state.cellIdColumn = controls.cellIdColumn.value;
  state.groupColumn = controls.groupColumn.value;
  state.projection = state.serverProjection.length === state.rows.length ? state.serverProjection : projectCells();
  if (state.serverClusters.length === state.rows.length) {
    state.rows.forEach((row, index) => {
      row.cluster = state.serverClusters[index];
    });
  } else {
    state.clusters = kmeans(state.projection, Math.min(6, Math.max(2, Math.round(Math.sqrt(state.rows.length / 6)))));
    state.rows.forEach((row, index) => {
      row.cluster = `Cluster ${state.clusters[index] + 1}`;
    });
  }
  fillSelect(controls.colorBy, [...state.metadata, "cluster"].filter(unique), controls.colorBy.value || "cluster");
  fillSelect(controls.markerCluster, clusterLabels(), clusterLabels()[0]);
  fillSelect(controls.groupA, groupValues(), groupValues()[0]);
  fillSelect(controls.groupB, groupValues(), groupValues()[1] || groupValues()[0]);
  computeMarkers();
  renderAll();
  if (state.activeView === "global") scheduleGlobalClusteringRender();
  addChat("assistant", `Loaded ${state.fileName}. Ask for markers, clusters, expression of a gene, or a group comparison.`);
}

function projectCells() {
  const matrix = normalizedMatrix();
  const centered = centerColumns(matrix);
  const pc1 = powerIteration(centered, 18);
  const scores1 = centered.map((row) => dot(row, pc1));
  const deflated = centered.map((row, i) => row.map((value, j) => value - scores1[i] * pc1[j]));
  const pc2 = powerIteration(deflated, 18);
  const scores2 = deflated.map((row) => dot(row, pc2));
  return scores1.map((x, index) => ({ x, y: scores2[index] || 0 }));
}

function normalizedMatrix() {
  const totals = state.rows.map((row) => state.genes.reduce((sum, gene) => sum + Math.max(0, Number(row[gene]) || 0), 0) || 1);
  return state.rows.map((row, rowIndex) =>
    state.genes.map((gene) => Math.log1p(((Number(row[gene]) || 0) / totals[rowIndex]) * 10000))
  );
}

function centerColumns(matrix) {
  const means = state.genes.map((_, index) => mean(matrix.map((row) => row[index])));
  return matrix.map((row) => row.map((value, index) => value - means[index]));
}

function powerIteration(matrix, iterations) {
  let vector = Array.from({ length: state.genes.length }, (_, index) => (index % 2 === 0 ? 1 : -1));
  vector = normalize(vector);
  for (let i = 0; i < iterations; i += 1) {
    const mv = matrix.map((row) => dot(row, vector));
    const mtmv = vector.map((_, columnIndex) => dot(mv, matrix.map((row) => row[columnIndex])));
    vector = normalize(mtmv);
  }
  return vector;
}

function kmeans(points, k) {
  const centroids = Array.from({ length: k }, (_, index) => points[Math.floor((index * points.length) / k)] || { x: 0, y: 0 });
  let labels = Array(points.length).fill(0);
  for (let iteration = 0; iteration < 20; iteration += 1) {
    labels = points.map((point) => closestCentroid(point, centroids));
    for (let cluster = 0; cluster < k; cluster += 1) {
      const members = points.filter((_, index) => labels[index] === cluster);
      if (members.length) {
        centroids[cluster] = {
          x: mean(members.map((point) => point.x)),
          y: mean(members.map((point) => point.y)),
        };
      }
    }
  }
  return labels;
}

function closestCentroid(point, centroids) {
  let best = 0;
  let bestDistance = Infinity;
  centroids.forEach((centroid, index) => {
    const distance = (point.x - centroid.x) ** 2 + (point.y - centroid.y) ** 2;
    if (distance < bestDistance) {
      best = index;
      bestDistance = distance;
    }
  });
  return best;
}

function computeMarkers() {
  state.markers = new Map();
  for (const cluster of clusterLabels()) {
    const inCluster = state.rows.filter((row) => row.cluster === cluster);
    const outCluster = state.rows.filter((row) => row.cluster !== cluster);
    const markers = state.genes
      .map((gene) => {
        const a = mean(inCluster.map((row) => Number(row[gene]) || 0));
        const b = mean(outCluster.map((row) => Number(row[gene]) || 0));
        return { gene, a, b, fold: Math.log2((a + 0.1) / (b + 0.1)) };
      })
      .sort((left, right) => right.fold - left.fold)
      .slice(0, 20);
    state.markers.set(cluster, markers);
  }
}

function renderAll() {
  $("datasetTitle").textContent = state.fileName;
  $("analysisStatus").textContent = "Analysis ready";
  $("metricCells").textContent = formatNumber(state.rows.length);
  $("metricGenes").textContent = formatNumber(state.genes.length);
  $("metricUmis").textContent = formatNumber(Math.round(median(state.rows.map((row) => state.genes.reduce((sum, gene) => sum + (Number(row[gene]) || 0), 0)))));
  $("metricClusters").textContent = clusterLabels().length;
  schedulePlotRender();
  scheduleGlobalClusteringRender();
  renderMarkerTable();
  renderDifferentialTable();
  renderSelectedCell(state.selectedCellIndex ?? 0);
}

function schedulePlotRender() {
  requestAnimationFrame(() => {
    void renderPlot().then(() => {
      const plot = $("pcaPlot");
      if (window.Plotly && plot?.data) Plotly.Plots.resize(plot);
    });
  });
}

async function renderPlot() {
  const plot = $("pcaPlot");
  if (!state.projection.length) {
    plot.innerHTML = '<div class="plot-empty">Upload a dataset to render cells.</div>';
    return;
  }

  if (!window.Plotly) {
    plot.innerHTML = '<div class="plot-empty">Plotly is still loading. Refresh once if the plot does not appear.</div>';
    return;
  }

  const colorBy = controls.colorBy.value;
  const overlayGene = controls.geneOverlay.value;
  const colorByGene = state.genes.includes(colorBy) ? colorBy : "";
  const gene = colorByGene || (overlayGene && overlayGene !== "none" ? overlayGene : "");
  const useGene = Boolean(gene) && colorBy !== "cell_type" && colorBy !== "assigned_cell_type" && colorBy !== "cell_type_group";
  const values = useGene ? state.rows.map((row) => Number(row[gene]) || 0) : state.rows.map((row) => row[colorBy] || "NA");
  const categories = useGene ? [] : [...new Set(values.map(String))];
  const markerColors = useGene ? values : values.map((value) => palette[categories.indexOf(String(value)) % palette.length]);
  const hasGlobalHighlight = state.globalHighlightIndices.size > 0;
  plot.querySelectorAll(".plot-empty").forEach((element) => element.remove());
  const hoverText = state.rows.map((row, index) => {
    const point = state.projection[index];
    const topGenes = topExpressedGenes(row, 3).map((item) => `${item.gene}: ${item.value.toFixed(2)}`).join("<br>");
    return [
      `<b>${cellLabel(row)}</b>`,
      `Cluster: ${row.cluster || "NA"}`,
      `Cell type: ${cellTypeLabel(row)}`,
      `PC1: ${point.x.toFixed(3)}`,
      `PC2: ${point.y.toFixed(3)}`,
      topGenes ? `Top genes:<br>${topGenes}` : "",
    ]
      .filter(Boolean)
      .join("<br>");
  });

  const trace = {
    type: "scattergl",
    mode: "markers",
    x: state.projection.map((point) => point.x),
    y: state.projection.map((point) => point.y),
    customdata: state.rows.map((_, index) => index),
    text: hoverText,
    hovertemplate: "%{text}<extra></extra>",
    marker: {
      color: markerColors,
      colorscale: useGene ? [[0, "#dbe7eb"], [1, "#ba4a68"]] : undefined,
      showscale: useGene,
      colorbar: useGene ? { title: gene, thickness: 12 } : undefined,
      line: {
        color: hasGlobalHighlight || state.multimodalHighlightIds.size
          ? state.rows.map((row, index) => (state.globalHighlightIndices.has(index) || state.multimodalHighlightIds.has(normalizeNeuronId(cellLabel(row))) ? "#111827" : "rgba(29,39,51,0.18)"))
          : "rgba(29,39,51,0.25)",
        width:
          hasGlobalHighlight || state.multimodalHighlightIds.size
            ? state.rows.map((row, index) => (state.globalHighlightIndices.has(index) || state.multimodalHighlightIds.has(normalizeNeuronId(cellLabel(row))) ? 2 : 0.4))
            : 0.5,
      },
      opacity:
        hasGlobalHighlight || state.multimodalHighlightIds.size
          ? state.rows.map((row, index) => (state.globalHighlightIndices.has(index) || state.multimodalHighlightIds.has(normalizeNeuronId(cellLabel(row))) ? 0.98 : 0.2))
          : 0.86,
      size: state.rows.length > 1000 ? 5 : 9,
    },
  };

  const layout = {
    margin: { l: 52, r: useGene ? 58 : 18, t: 8, b: 44 },
    paper_bgcolor: "#fbfcfc",
    plot_bgcolor: "#fbfcfc",
    dragmode: "pan",
    hovermode: "closest",
    xaxis: { title: "PC1", zeroline: false, gridcolor: "#e5ecef" },
    yaxis: { title: "PC2", zeroline: false, gridcolor: "#e5ecef" },
  };

  await Plotly.react(plot, [trace], layout, {
    responsive: true,
    displaylogo: false,
    modeBarButtonsToRemove: ["lasso2d", "select2d"],
  });
  Plotly.Plots.resize(plot);
  if (plot.removeAllListeners) plot.removeAllListeners("plotly_click");
  plot.on("plotly_click", (event) => {
    const index = event.points?.[0]?.customdata;
    if (Number.isInteger(index)) {
      state.selectedCellIndex = index;
      renderSelectedCell(index);
    }
  });
  renderLegend(useGene, categories, gene);
}

function scheduleGlobalClusteringRender() {
  if (state.activeView !== "global") return;
  requestAnimationFrame(() => {
    void runGlobalClustering(false);
  });
}

async function runGlobalClustering(force = true) {
  const status = $("globalStatus");
  if (!state.rows.length || !state.genes.length) {
    status.textContent = "Load a dataset first";
    $("dendrogramPlot").innerHTML = '<div class="plot-empty">Load neurons to compute a dendrogram.</div>';
    $("globalHeatmap").innerHTML = '<div class="plot-empty">Load neurons to compute pairwise similarity.</div>';
    return;
  }
  if (!window.Plotly) {
    status.textContent = "Plotly is still loading";
    return;
  }

  const payload = globalClusteringPayload();
  if (payload.matrix.length < 2) {
    status.textContent = "Need at least two cells";
    return;
  }
  const cacheKey = payload.cacheKey;
  if (!force && globalClusterCache.has(cacheKey)) {
    state.globalClustering = globalClusterCache.get(cacheKey);
    renderGlobalClustering();
    return;
  }

  status.textContent = `Computing ${payload.metric} similarity for ${payload.matrix.length} profiles and ${payload.genes.length} informative genes`;
  controls.runGlobalClustering.disabled = true;
  try {
    const response = await fetch(apiUrl("/api/global-clustering"), {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    if (!response.ok) {
      throw new Error(await responseErrorMessage(response, "Global clustering failed"));
    }
    state.globalClustering = await response.json();
    globalClusterCache.set(cacheKey, state.globalClustering);
    renderGlobalClustering();
  } catch (error) {
    status.textContent = `Global clustering failed: ${error.message}`;
    $("globalHeatmap").innerHTML = `<div class="plot-empty">${error.message}</div>`;
  } finally {
    controls.runGlobalClustering.disabled = false;
  }
}

function globalClusteringPayload() {
  const limit = Math.min(state.rows.length, maxGlobalClusterCells);
  const rows = state.rows.slice(0, limit);
  const metric = controls.globalMetric.value || "pearson";
  const method = controls.globalLinkage.value || "average";
  const genes = selectGlobalClusteringGenes(rows);
  const matrix = rows.map((row) => genes.map((gene) => Number(row[gene]) || 0));
  return {
    metric,
    method,
    matrix,
    genes,
    labels: rows.map(cellLabel),
    classes: rows.map(cellTypeLabel),
    cacheKey: stableGlobalCacheKey(rows, genes, metric, method),
  };
}

function selectGlobalClusteringGenes(rows) {
  if (state.genes.length <= maxGlobalClusterGenes) return state.genes;
  return state.genes
    .map((gene) => {
      const values = rows.map((row) => Number(row[gene]) || 0);
      const avg = mean(values);
      const variance = mean(values.map((value) => (value - avg) ** 2));
      const detected = values.filter((value) => value > 0).length / Math.max(1, values.length);
      return { gene, score: variance * (0.25 + detected) };
    })
    .sort((left, right) => right.score - left.score)
    .slice(0, maxGlobalClusterGenes)
    .map((item) => item.gene);
}

function stableGlobalCacheKey(rows, genes, metric, method) {
  const total = rows.reduce((sum, row) => sum + genes.reduce((geneSum, gene) => geneSum + (Number(row[gene]) || 0), 0), 0);
  return [state.fileName, rows.length, genes.length, metric, method, total.toFixed(3), genes.slice(0, 8).join(",")].join("|");
}

function renderGlobalClustering() {
  const result = state.globalClustering;
  if (!result) return;
  $("globalStatus").textContent = `${result.summary.cells} profiles, ${result.summary.genes} genes, ${result.summary.clusterCount} groups`;
  renderDendrogram(result);
  renderSimilarityHeatmap(result);
}

function renderDendrogram(result) {
  const traces = result.dendrogram.icoord.map((xValues, index) => ({
    type: "scatter",
    mode: "lines",
    x: xValues.map((value) => (value - 5) / 10),
    y: result.dendrogram.dcoord[index],
    line: { color: "#8a96a3", width: 1.4 },
    hoverinfo: "skip",
    showlegend: false,
  }));
  Plotly.react(
    $("dendrogramPlot"),
    traces,
    {
      margin: { l: 54, r: 18, t: 8, b: 6 },
      paper_bgcolor: "#fbfcfc",
      plot_bgcolor: "#fbfcfc",
      xaxis: { showticklabels: false, range: [-0.5, result.labels.length - 0.5], fixedrange: false },
      yaxis: { title: "Distance", zeroline: false, gridcolor: "#e5ecef" },
    },
    { responsive: true, displaylogo: false }
  );
}

function renderSimilarityHeatmap(result) {
  const labels = result.labels;
  const customdata = result.similarity.map((row, rowIndex) =>
    row.map((_, colIndex) => ({
      rowIndex: result.order[rowIndex],
      colIndex: result.order[colIndex],
      cluster: result.clusters[rowIndex],
      rowLabel: labels[rowIndex],
      colLabel: labels[colIndex],
    }))
  );
  const heatmap = {
    type: "heatmap",
    z: result.similarity,
    x: labels,
    y: labels,
    customdata,
    colorscale: [
      [0, "#1946d2"],
      [0.5, "#ffffff"],
      [1, "#f01818"],
    ],
    zmin: -1,
    zmax: 1,
    colorbar: { title: result.metric === "cosine" ? "cosine" : "corr", thickness: 12 },
    hovertemplate: "%{customdata.rowLabel} vs %{customdata.colLabel}<br>Similarity: %{z:.3f}<extra></extra>",
  };
  const layout = {
    margin: { l: 74, r: 22, t: 8, b: 86 },
    paper_bgcolor: "#fbfcfc",
    plot_bgcolor: "#fbfcfc",
    dragmode: "pan",
    xaxis: { tickangle: 60, tickfont: { size: labels.length > 80 ? 7 : 9 }, automargin: true },
    yaxis: { autorange: "reversed", tickfont: { size: labels.length > 80 ? 7 : 9 }, automargin: true },
    shapes: globalBlockShapes(result),
  };
  const heatmapElement = $("globalHeatmap");
  Plotly.react(heatmapElement, [heatmap], layout, {
    responsive: true,
    displaylogo: false,
    scrollZoom: true,
    modeBarButtonsToRemove: ["lasso2d", "select2d"],
  });
  if (heatmapElement.removeAllListeners) heatmapElement.removeAllListeners("plotly_click");
  heatmapElement.on("plotly_click", (event) => {
    const data = event.points?.[0]?.customdata;
    if (!data) return;
    highlightGlobalCluster(data.cluster);
    state.selectedCellIndex = data.rowIndex;
    renderSelectedCell(data.rowIndex);
  });
}

function globalBlockShapes(result) {
  const shapes = [];
  result.blocks.forEach((block, index) => {
    const color = palette[index % palette.length];
    shapes.push({
      type: "rect",
      xref: "x",
      yref: "y",
      x0: block.start - 0.5,
      x1: block.end + 0.5,
      y0: block.start - 0.5,
      y1: block.end + 0.5,
      line: { color, width: 2 },
      fillcolor: "rgba(0,0,0,0)",
    });
  });
  return shapes;
}

function renderClusterBlocks(result) {
  const container = $("clusterBlocks");
  if (!container) return;
  container.innerHTML = "";
  result.blocks.forEach((block, index) => {
    const button = document.createElement("button");
    button.className = "cluster-chip";
    button.type = "button";
    button.innerHTML = `<i style="background:${palette[index % palette.length]}"></i><div><strong>${block.label}</strong><br><span>Cluster ${block.cluster}, n=${block.count}</span></div>`;
    button.addEventListener("click", () => highlightGlobalCluster(block.cluster));
    container.appendChild(button);
  });
}

function highlightGlobalCluster(clusterId) {
  const result = state.globalClustering;
  if (!result) return;
  const selected = result.order.filter((originalIndex, orderedIndex) => result.clusters[orderedIndex] === clusterId);
  state.globalHighlightIndices = new Set(selected);
  $("globalStatus").textContent = `Highlighted cluster ${clusterId} (${selected.length} profiles)`;
  schedulePlotRender();
}

async function loadMultimodalIndex() {
  if (!$("multimodalSummary")) return;
  try {
    const response = await fetch(apiUrl("/api/multimodal/neurons"));
    if (!response.ok) throw new Error("Multimodal API unavailable");
    state.multimodalNeurons = await response.json();
    renderMultimodalViews(state.multimodalNeurons.slice(0, 6).map((item) => item.neuron_id));
    $("multimodalSummary").textContent = `${state.multimodalNeurons.length} multimodal neuron records indexed`;
  } catch (error) {
    $("multimodalSummary").textContent = `Multimodal index unavailable: ${error.message}`;
  }
}

async function runMultimodalQuery(query) {
  if (!$("multimodalSummary")) return;
  if (!query) return;
  $("multimodalSummary").textContent = "Running multimodal query";
  try {
    const result = await executeMultimodalQuery(query);
    $("multimodalSummary").textContent = `${result.kind.replaceAll("_", " ")}: ${result.highlighted_neurons.length} highlighted neurons`;
    renderMultimodalResults(result);
    renderMultimodalViews(result.highlighted_neurons);
    schedulePlotRender();
  } catch (error) {
    $("multimodalSummary").textContent = `Query failed: ${error.message}`;
  }
}

async function executeMultimodalQuery(query) {
  const response = await fetch(apiUrl("/api/multimodal/query"), {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      query,
      dataset: {
        fileName: state.fileName,
        cells: state.rows.length,
        genes: state.genes.length,
      },
    }),
  });
  if (!response.ok) throw new Error(await response.text());
  const result = await response.json();
  state.multimodalHighlightIds = new Set(result.highlighted_neurons.map(normalizeNeuronId));
  return result;
}

function renderMultimodalResults(result) {
  const container = $("multimodalResults");
  const rows = Array.isArray(result.result) ? result.result : [result.result];
  container.innerHTML = rows
    .slice(0, 12)
    .map((item) => {
      const label = item.neuron_id || item.source || item.target || result.kind;
      const detail = Object.entries(item)
        .filter(([key]) => key !== "neighbors")
        .map(([key, value]) => `${key}: ${Array.isArray(value) ? value.join(" > ") : value}`)
        .join(" | ");
      return `<div class="result-row"><div><strong>${label}</strong>${detail}</div></div>`;
    })
    .join("");
}

function renderMultimodalViews(highlightedIds = []) {
  const highlighted = new Set(highlightedIds.map(normalizeNeuronId));
  const neurons = state.multimodalNeurons;
  if (!neurons.length || !window.Plotly) return;
  renderConnectomeGraph(neurons, highlighted);
  renderSpatialMap(neurons, highlighted);
  renderLineageTree(neurons, highlighted);
}

function renderConnectomeGraph(neurons, highlighted) {
  const ids = neurons.map((item) => item.neuron_id);
  const positions = new Map(ids.map((id, index) => [id, circularPoint(index, ids.length, 1)]));
  const edgeX = [];
  const edgeY = [];
  neurons.forEach((source) => {
    const start = positions.get(source.neuron_id);
    source.connectome_edges.forEach((edge) => {
      const end = positions.get(edge.target);
      if (!start || !end) return;
      edgeX.push(start.x, end.x, null);
      edgeY.push(start.y, end.y, null);
    });
  });
  const nodeTrace = {
    type: "scatter",
    mode: "markers+text",
    x: ids.map((id) => positions.get(id).x),
    y: ids.map((id) => positions.get(id).y),
    text: ids,
    textposition: "top center",
    marker: {
      color: ids.map((id) => (highlighted.has(id) ? "#ba4a68" : "#1f7a8c")),
      size: ids.map((id) => (highlighted.has(id) ? 16 : 10)),
      line: { color: "#ffffff", width: 1 },
    },
    hovertemplate: "%{text}<extra>connectome</extra>",
  };
  Plotly.react(
    $("connectomeGraph"),
    [{ type: "scatter", mode: "lines", x: edgeX, y: edgeY, line: { color: "#c6d0d7", width: 1 }, hoverinfo: "skip" }, nodeTrace],
    multimodalLayout("Connectome"),
    { responsive: true, displaylogo: false }
  );
}

function renderSpatialMap(neurons, highlighted) {
  Plotly.react(
    $("spatialMap"),
    [
      {
        type: "scatter",
        mode: "markers+text",
        x: neurons.map((item) => item.spatial_coordinates.x),
        y: neurons.map((item) => item.spatial_coordinates.y),
        text: neurons.map((item) => item.neuron_id),
        textposition: "top center",
        marker: {
          color: neurons.map((item) => (highlighted.has(item.neuron_id) ? "#ba4a68" : "#4d8d62")),
          size: neurons.map((item) => (highlighted.has(item.neuron_id) ? 15 : 9)),
        },
        hovertemplate: "%{text}<br>x=%{x:.1f}<br>y=%{y:.1f}<extra>spatial</extra>",
      },
    ],
    multimodalLayout("Spatial map"),
    { responsive: true, displaylogo: false }
  );
}

function renderLineageTree(neurons, highlighted) {
  const x = [];
  const y = [];
  const text = [];
  neurons.forEach((item, index) => {
    x.push(item.lineage_path.length);
    y.push(index);
    text.push(item.neuron_id);
  });
  Plotly.react(
    $("lineageTree"),
    [
      {
        type: "scatter",
        mode: "markers+text",
        x,
        y,
        text,
        textposition: "middle right",
        marker: {
          color: text.map((id) => (highlighted.has(id) ? "#ba4a68" : "#d88742")),
          size: text.map((id) => (highlighted.has(id) ? 15 : 9)),
        },
        hovertemplate: "%{text}<br>lineage depth=%{x}<extra>lineage</extra>",
      },
    ],
    multimodalLayout("Lineage depth"),
    { responsive: true, displaylogo: false }
  );
}

function multimodalLayout(title) {
  return {
    title: { text: title, font: { size: 13 } },
    margin: { l: 28, r: 16, t: 34, b: 28 },
    paper_bgcolor: "#fbfcfc",
    plot_bgcolor: "#fbfcfc",
    xaxis: { zeroline: false, showgrid: false },
    yaxis: { zeroline: false, showgrid: false, showticklabels: false },
  };
}

function circularPoint(index, total, radius) {
  const angle = (index / Math.max(1, total)) * Math.PI * 2;
  return { x: Math.cos(angle) * radius, y: Math.sin(angle) * radius };
}

function normalizeNeuronId(value) {
  return String(value || "").trim().toUpperCase();
}

function renderSelectedCell(index) {
  if (!state.rows.length || index === null || index === undefined || !state.rows[index]) return;
  const row = state.rows[index];
  const point = state.projection[index] || { x: 0, y: 0 };
  const topGenes = topExpressedGenes(row, 6);
  $("selectedCellTitle").textContent = cellLabel(row);
  $("selectedCellCluster").textContent = row.cluster || "NA";
  $("selectedCellType").textContent = cellTypeLabel(row);
  $("selectedCellCoords").textContent = `${point.x.toFixed(3)} / ${point.y.toFixed(3)}`;
  $("selectedCellGenes").textContent = topGenes.map((item) => `${item.gene} ${item.value.toFixed(2)}`).join(", ") || "No expression values";
  $("selectedCellBio").textContent = explainSelectedCell(row, point, topGenes);
}

function topExpressedGenes(row, limit = 5) {
  return state.genes
    .map((gene) => ({ gene, value: Number(row[gene]) || 0 }))
    .sort((left, right) => right.value - left.value)
    .slice(0, limit)
    .filter((item) => item.value > 0);
}

function explainSelectedCell(row, point, topGenes) {
  const cellType = cellTypeLabel(row);
  const cluster = row.cluster || "unassigned cluster";
  const markerGenes = (state.markers.get(cluster) || []).slice(0, 5).map((item) => item.gene);
  const direction = `${point.x >= 0 ? "positive" : "negative"} PC1 and ${point.y >= 0 ? "positive" : "negative"} PC2`;
  const topGeneText = topGenes.map((item) => item.gene).join(", ") || "no high-expression genes in the selected view";
  const markerText = markerGenes.length ? ` Cluster marker genes include ${markerGenes.join(", ")}.` : "";
  return `This profile represents ${cellType} in ${cluster}. It sits in the ${direction} region of the PCA plot because its expression profile is dominated by ${topGeneText}; these genes pull the cell toward nearby profiles with similar expression patterns.${markerText} In a Scanpy or Seurat workflow, these PC1/PC2 coordinates correspond to the object's PCA embedding, and the explanation is based on the displayed expression values plus cluster-level marker genes.`;
}

function cellTypeLabel(row) {
  return String(row.cell_type || row.assigned_cell_type || row.cell_type_group || row.annotation || row.type || row.sample || "unspecified");
}

function renderLegend(useGene, categories, gene) {
  const legend = $("plotLegend");
  legend.innerHTML = "";
  if (useGene) {
    legend.innerHTML = `<span><i style="background:#dbe7eb"></i>low ${gene}</span><span><i style="background:#ba4a68"></i>high ${gene}</span>`;
    return;
  }
  categories.slice(0, 8).forEach((category, index) => {
    const item = document.createElement("span");
    item.innerHTML = `<i style="background:${palette[index % palette.length]}"></i>${category}`;
    legend.appendChild(item);
  });
}

function renderMarkerTable() {
  const cluster = controls.markerCluster.value || clusterLabels()[0];
  const markers = state.markers.get(cluster) || [];
  $("markerTable").innerHTML = markers
    .map((row) => `<tr><td>${row.gene}</td><td>${row.a.toFixed(2)}</td><td>${row.b.toFixed(2)}</td><td>${row.fold.toFixed(2)}</td></tr>`)
    .join("");
}

function renderDifferentialTable() {
  const groupColumn = controls.groupColumn.value || state.groupColumn;
  const aName = controls.groupA.value;
  const bName = controls.groupB.value;
  const groupA = state.rows.filter((row) => String(row[groupColumn]) === aName || row.cluster === aName);
  const groupB = state.rows.filter((row) => String(row[groupColumn]) === bName || row.cluster === bName);
  const rows = state.genes
    .map((gene) => {
      const a = mean(groupA.map((row) => Number(row[gene]) || 0));
      const b = mean(groupB.map((row) => Number(row[gene]) || 0));
      return { gene, a, b, fold: Math.log2((a + 0.1) / (b + 0.1)) };
    })
    .sort((left, right) => Math.abs(right.fold) - Math.abs(left.fold))
    .slice(0, 20);
  $("deTable").innerHTML = rows
    .map((row) => `<tr><td>${row.gene}</td><td>${row.a.toFixed(2)}</td><td>${row.b.toFixed(2)}</td><td>${row.fold.toFixed(2)}</td></tr>`)
    .join("");
}

async function answerQuestion(message) {
  const intent = analyzeQuery(message);
  const ambiguity = ambiguityNote(intent);
  const response = await executeQueryIntent(intent, message);
  return ambiguity ? `${ambiguity} ${response}` : response;
}

function isMultimodalQuestion(query) {
  const intent = analyzeQuery(query);
  return intent.modalities.some((modality) => ["connectome", "spatial", "lineage", "multimodal"].includes(modality.name));
}

function analyzeQuery(message) {
  const normalized = normalizeName(message);
  const tokens = new Set(normalized.split(/\s+/).filter(Boolean));
  const genes = extractGenes(message);
  const neurons = extractNeuronEntities(message);
  const cluster = clusterLabels().find((label) => normalized.includes(normalizeName(label))) || matchClusterNumber(normalized);
  const modalities = scoreSemanticDimensions(tokens, normalized, {
    transcriptomics: ["transcript", "expression", "express", "gene", "marker", "cluster", "pca", "umap", "similar"],
    connectome: ["connectome", "synapse", "synaptic", "partner", "edge", "weight", "gap", "junction", "connected", "neighbor"],
    spatial: ["spatial", "anatomy", "anatomical", "near", "nearest", "close", "proximal", "coordinate", "position", "location", "distance", "region"],
    lineage: ["lineage", "ancestor", "descendant", "parent", "development", "developmental", "tree"],
    annotation: ["describe", "annotation", "homology", "ortholog", "disease", "phenotype", "pathway", "wormbase", "interact"],
  });
  const operations = scoreSemanticDimensions(tokens, normalized, {
    rank: ["highest", "top", "max", "most", "strongest", "enriched", "rank"],
    compare: ["compare", "versus", "vs", "different", "differential", "between"],
    visualize: ["show", "plot", "color", "highlight", "display", "map"],
    summarize: ["summary", "summarize", "summarise", "overview", "many", "count", "size", "describe"],
    retrieve: ["find", "which", "what", "where", "who", "list", "near"],
    annotate: ["describe", "annotation", "homology", "ortholog", "disease", "phenotype", "pathway", "wormbase"],
  });
  const relationCount = modalities.filter((item) => item.score > 0).length;
  if (relationCount > 1) modalities.unshift({ name: "multimodal", score: relationCount });
  return {
    raw: message,
    normalized,
    genes,
    neurons,
    cluster,
    limit: parseResultLimit(normalized),
    modalities: modalities.filter((item) => item.score > 0).sort((a, b) => b.score - a.score),
    operations: operations.filter((item) => item.score > 0).sort((a, b) => b.score - a.score),
  };
}

function scoreSemanticDimensions(tokens, normalized, dimensions) {
  return Object.entries(dimensions).map(([name, words]) => {
    const score = words.reduce((sum, word) => {
      const normalizedWord = normalizeName(word);
      const tokenHit = tokens.has(normalizedWord);
      const stemHit = [...tokens].some((token) => token.startsWith(normalizedWord) || normalizedWord.startsWith(token));
      const phraseHit = normalized.includes(normalizedWord);
      return sum + (tokenHit ? 2 : 0) + (!tokenHit && stemHit ? 1 : 0) + (!tokenHit && phraseHit ? 1 : 0);
    }, 0);
    return { name, score };
  });
}

async function executeQueryIntent(intent, message) {
  if (shouldUseMultimodalEngine(intent)) return await answerMultimodalQuestion(message);
  if (!state.rows.length) return "Upload or load a dataset first, then I can query transcriptomic values from it. I can still answer connectome, spatial, or lineage questions from the multimodal index.";

  const topOperation = intent.operations[0]?.name || "retrieve";
  const topModality = intent.modalities[0]?.name || "transcriptomics";
  const gene = intent.genes[0];

  if (topModality === "annotation" || (topOperation === "annotate" && (gene || extractGeneLikeToken(message)))) {
    const targetGene = gene || extractGeneLikeToken(message);
    if (!targetGene) return "I interpret this as a gene annotation question, but I need a gene symbol or WormBase ID to answer it.";
    return await answerWormBaseQuestion(targetGene, intent.normalized);
  }
  if (topOperation === "rank" && intent.genes.length) return answerHighestExpression(intent.genes, intent.limit);
  if ((topOperation === "visualize" || topModality === "transcriptomics") && gene && /express|gene|plot|show|color|display/.test(intent.normalized)) {
    controls.geneOverlay.value = gene;
    schedulePlotRender();
    const positive = state.rows.filter((row) => Number(row[gene]) > 0).length;
    return `${gene} is detected in ${positive} of ${state.rows.length} cells. I colored the projection by ${gene}.`;
  }
  if (/marker|signature/.test(intent.normalized) || (topModality === "transcriptomics" && intent.cluster && topOperation !== "compare")) {
    const label = intent.cluster || controls.markerCluster.value || clusterLabels()[0];
    controls.markerCluster.value = label;
    renderMarkerTable();
    const markers = (state.markers.get(label) || []).slice(0, 5).map((row) => row.gene).join(", ");
    return `For ${label}, the strongest displayed marker genes are ${markers || "not available"}. I updated the marker table.`;
  }
  if (topOperation === "compare") {
    renderDifferentialTable();
    const first = $("deTable tr td");
    return first ? `I interpreted this as a group comparison. The strongest displayed changing gene is ${first.textContent}.` : "I interpreted this as a comparison, but I need at least two groups to compare.";
  }
  if (topOperation === "summarize") {
    return `This loaded dataset has ${state.rows.length} cells/profiles, ${state.genes.length} numeric gene columns, and ${clusterLabels().length} computed clusters.`;
  }
  if (intent.normalized.includes("cluster")) {
    const counts = clusterLabels().map((label) => `${label}: ${state.rows.filter((row) => row.cluster === label).length}`).join("; ");
    return `Cluster sizes are ${counts}.`;
  }
  return fallbackIntentAnswer(intent);
}

function shouldUseMultimodalEngine(intent) {
  if (!$("multimodalSummary")) return false;
  const names = intent.modalities.map((item) => item.name);
  return names.includes("multimodal") || names.includes("connectome") || names.includes("spatial") || names.includes("lineage");
}

function ambiguityNote(intent) {
  const closeOperations = closeScoredNames(intent.operations);
  const closeModalities = closeScoredNames(intent.modalities.filter((item) => item.name !== "multimodal"));
  const notes = [];
  if (closeOperations.length > 1) notes.push(`I see more than one possible operation (${closeOperations.join(", ")}), so I chose the one most supported by the query.`);
  if (closeModalities.length > 1) notes.push(`The question touches multiple data modalities (${closeModalities.join(", ")}), so I treated it as an integrated query.`);
  return notes.join(" ");
}

function closeScoredNames(items) {
  if (!items.length) return [];
  const best = items[0].score;
  return items.filter((item) => item.score > 0 && best - item.score <= 1).map((item) => item.name);
}

function fallbackIntentAnswer(intent) {
  const modalities = intent.modalities.map((item) => item.name).join(", ") || "transcriptomics";
  const operations = intent.operations.map((item) => item.name).join(", ") || "retrieve";
  return `I parsed this as operation(s): ${operations}; modality/modalities: ${modalities}. I need a more specific target such as a gene, neuron ID, cluster, comparison groups, or relationship to execute a precise query.`;
}

function extractNeuronEntities(message) {
  const upper = message.toUpperCase();
  return state.multimodalNeurons
    .map((item) => item.neuron_id)
    .filter((id) => new RegExp(`(^|[^A-Z0-9])${escapeRegex(id)}([^A-Z0-9]|$)`).test(upper));
}

async function answerMultimodalQuestion(message) {
  if (!$("multimodalSummary")) {
    return "Multimodal querying is temporarily hidden in this version. I can still answer questions about the loaded transcriptomic dataset.";
  }
  try {
    const result = await executeMultimodalQuery(message);
    $("multimodalSummary").textContent = `${result.kind.replaceAll("_", " ")}: ${result.highlighted_neurons.length} highlighted neurons`;
    renderMultimodalResults(result);
    renderMultimodalViews(result.highlighted_neurons);
    schedulePlotRender();
    return formatMultimodalChatAnswer(result);
  } catch (error) {
    return `I could not run the multimodal query right now. Details: ${error.message}`;
  }
}

function formatMultimodalChatAnswer(result) {
  const highlighted = result.highlighted_neurons || [];
  const resultRows = Array.isArray(result.result) ? result.result : [result.result];
  const topRows = resultRows.slice(0, 5).map(formatMultimodalResultRow).filter(Boolean);
  const datasetContext = transcriptomicContextForNeurons(highlighted);
  const intro = `I interpreted this as a ${result.kind.replaceAll("_", " ")} query and highlighted ${highlighted.length} neuron(s) across the PCA, connectome, spatial, and lineage views.`;
  const findings = topRows.length ? ` Key results: ${topRows.join("; ")}.` : "";
  return `${intro}${findings}${datasetContext}`;
}

function formatMultimodalResultRow(item) {
  if (!item || typeof item !== "object") return "";
  if (item.source && item.target) {
    return `${item.source}-${item.target} (${item.relationship || "relationship"}, weight ${Number(item.total_weight || 0).toFixed(1)})`;
  }
  if (item.neuron_id && item.distance !== undefined) {
    const distance = item.distance === null ? "not spatially indexed" : `${Number(item.distance).toFixed(2)} units`;
    const weight = item.total_weight !== undefined ? `, synaptic weight ${Number(item.total_weight).toFixed(1)}` : "";
    return `${item.neuron_id} (${distance}${weight})`;
  }
  if (item.neuron_id && item.ancestors) {
    return `${item.neuron_id} lineage ${item.ancestors.join(" > ")}; similar lineage: ${(item.similar_lineage || []).join(", ") || "none"}`;
  }
  if (item.neuron_id && item.neighbors) {
    return `${item.neuron_id} has ${item.degree || item.neighbors.length} connected neighbor(s)`;
  }
  return Object.entries(item)
    .slice(0, 4)
    .map(([key, value]) => `${key}: ${Array.isArray(value) ? value.join(", ") : value}`)
    .join(", ");
}

function transcriptomicContextForNeurons(neuronIds) {
  if (!state.rows.length) {
    return " No expression dataset is currently loaded, so this answer uses the multimodal neuron index only.";
  }
  const matches = neuronIds
    .map((neuronId) => ({ neuronId, rows: matchingDatasetRows(neuronId).slice(0, 3) }))
    .filter((item) => item.rows.length);
  if (!matches.length) {
    return ` I did not find exact matching neuron IDs in the loaded expression view (${state.fileName}), so transcriptomic highlighting uses the multimodal IDs only.`;
  }
  const summaries = matches.slice(0, 6).map(({ neuronId, rows }) => {
    const clusters = [...new Set(rows.map((row) => row.cluster || cellTypeLabel(row)).filter(Boolean))].slice(0, 3).join(", ");
    const topGenes = topExpressedGenes(rows[0], 3).map((item) => item.gene).join(", ");
    return `${neuronId}: ${clusters || "unclustered"}${topGenes ? `; top displayed genes ${topGenes}` : ""}`;
  });
  const similarity = transcriptomicSimilaritySummary(matches.map((item) => item.neuronId));
  return ` Transcriptomic context from ${state.fileName}: ${summaries.join("; ")}.${similarity}`;
}

function matchingDatasetRows(neuronId) {
  const normalized = normalizeNeuronId(neuronId);
  return state.rows.filter((row) => {
    const candidates = [cellLabel(row), cellTypeLabel(row), row.neuron_id, row.cell_id, row.annotation, row.assigned_cell_type].map(normalizeNeuronId);
    return candidates.some((value) => value === normalized || new RegExp(`(^|[^A-Z0-9])${escapeRegex(normalized)}([^A-Z0-9]|$)`).test(value));
  });
}

function transcriptomicSimilaritySummary(neuronIds) {
  const matchedRows = neuronIds.flatMap((id) => matchingDatasetRows(id).slice(0, 1));
  if (matchedRows.length < 2) return "";
  const pairs = [];
  for (let i = 0; i < matchedRows.length; i += 1) {
    for (let j = i + 1; j < matchedRows.length; j += 1) {
      pairs.push({
        a: cellLabel(matchedRows[i]),
        b: cellLabel(matchedRows[j]),
        value: pearsonForRows(matchedRows[i], matchedRows[j]),
      });
    }
  }
  const best = pairs.sort((left, right) => right.value - left.value)[0];
  return best ? ` Strongest matched transcriptomic similarity: ${best.a}-${best.b} r=${best.value.toFixed(2)}.` : "";
}

function pearsonForRows(left, right) {
  const a = state.genes.map((gene) => Number(left[gene]) || 0);
  const b = state.genes.map((gene) => Number(right[gene]) || 0);
  const meanA = mean(a);
  const meanB = mean(b);
  const numerator = a.reduce((sum, value, index) => sum + (value - meanA) * (b[index] - meanB), 0);
  const denomA = Math.sqrt(a.reduce((sum, value) => sum + (value - meanA) ** 2, 0));
  const denomB = Math.sqrt(b.reduce((sum, value) => sum + (value - meanB) ** 2, 0));
  return denomA && denomB ? numerator / (denomA * denomB) : 0;
}

function isGeneInfoQuestion(query) {
  return /wormbase|describe|description|homolog|homology|ortholog|disease|phenotype|pathway|location|interact|interaction|where is .*located|what is .* gene|tell me about/.test(query);
}

async function answerWormBaseQuestion(gene, query) {
  const metadata = state.geneMetadata[gene] || {};
  if (!metadata.id || !metadata.id.startsWith("WBGene")) {
    const localDescription = metadata.description ? ` Local dataset description: ${metadata.description}` : "";
    return `I found ${gene} in the loaded expression matrix, but I do not have its WormBase WBGene identifier.${localDescription} Load one of the H5AD datasets or use a gene that came from a WormBase-backed H5AD file.`;
  }

  try {
    const response = await fetch(apiUrl(`/api/wormbase/gene/${encodeURIComponent(gene)}?gene_id=${encodeURIComponent(metadata.id)}`));
    if (!response.ok) {
      const detail = await response.text();
      throw new Error(detail || `status ${response.status}`);
    }
    const info = await response.json();
    return formatWormBaseAnswer(info, query);
  } catch (error) {
    return `I could not retrieve WormBase annotations for ${gene} right now. ${metadata.description || ""} Details: ${error.message}`;
  }
}

function formatWormBaseAnswer(info, query) {
  const sections = [];
  const header = `${info.gene} (${info.id})`;

  if (/describe|description|what is|tell me about|wormbase/.test(query) && info.description) {
    sections.push(`Description: ${info.description}`);
  }
  if (/homolog|homology|ortholog|wormbase/.test(query) && info.homology?.length) {
    sections.push(`Homology: ${info.homology.join("; ")}`);
  }
  if (/disease|wormbase/.test(query)) {
    const diseaseText = [info.humanDiseaseRelevance, ...(info.humanDiseases || [])].filter(Boolean).slice(0, 8).join("; ");
    if (diseaseText) sections.push(`Human disease: ${diseaseText}`);
  }
  if (/interact|interaction|wormbase/.test(query) && info.interactingGenes?.length) {
    sections.push(`Interacting genes: ${info.interactingGenes.slice(0, 8).join(", ")}`);
  }
  if (/location|located|where|wormbase/.test(query)) {
    const cellular = info.location?.cellular?.slice(0, 6).join(", ");
    const genomic = info.location?.genomic;
    const locationText = [genomic && `genomic ${genomic}`, cellular && `cellular ${cellular}`].filter(Boolean).join("; ");
    if (locationText) sections.push(`Location: ${locationText}`);
  }
  if (/pathway|process|wormbase/.test(query) && info.pathways?.length) {
    sections.push(`Pathways/processes: ${info.pathways.slice(0, 8).join(", ")}`);
  }
  if (/phenotype|wormbase/.test(query) && info.phenotypes?.length) {
    sections.push(`Phenotypes: ${info.phenotypes.slice(0, 8).join(", ")}`);
  }
  if (/expression|expressed|where is .*expressed|wormbase/.test(query)) {
    const expressionText = [...(info.expression?.tissues || []), ...(info.expression?.notes || [])].slice(0, 6).join("; ");
    if (expressionText) sections.push(`WormBase expression: ${expressionText}`);
  }

  if (!sections.length && info.description) sections.push(`Description: ${info.description}`);
  if (!sections.length) sections.push("No matching WormBase annotation fields were returned for that query.");
  return `${header}. ${sections.join(" ")}`;
}

function answerHighestExpression(genes, limit) {
  const ranked = rankCellsByExpression(genes, limit);
  if (!ranked.length) return `I could not find expression values for ${genes.join(", ")} in this dataset.`;

  controls.geneOverlay.value = genes[0];
  schedulePlotRender();

  const geneText = genes.length === 1 ? genes[0] : genes.join(", ");
  const label = state.metadata.includes("cell_type") ? "cell types/cells" : "cells";
  const results = ranked
    .map((row, index) => `${index + 1}. ${row.label} (${row.scoreLabel})`)
    .join("; ");
  return `Highest ${geneText} expression is in these ${label}: ${results}. I also colored the projection by ${genes[0]}.`;
}

function rankCellsByExpression(genes, limit = 5) {
  return state.rows
    .map((row) => {
      const values = genes.map((geneName) => Number(row[geneName]) || 0);
      const score = values.reduce((sum, value) => sum + value, 0) / values.length;
      const scoreLabel =
        genes.length === 1
          ? `${genes[0]}=${score.toFixed(2)}`
          : `mean=${score.toFixed(2)}; ${genes.map((geneName, index) => `${geneName}=${values[index].toFixed(2)}`).join(", ")}`;
      return {
        label: cellLabel(row),
        score,
        scoreLabel,
      };
    })
    .sort((left, right) => right.score - left.score)
    .slice(0, limit)
    .filter((row) => row.score > 0);
}

function extractGenes(message) {
  const normalizedMessage = normalizeName(message);
  return state.genes.filter((geneName) => {
    const normalizedGene = normalizeName(geneName);
    return new RegExp(`(^|[^a-z0-9])${escapeRegex(normalizedGene)}([^a-z0-9]|$)`).test(normalizedMessage);
  });
}

function extractGeneLikeToken(message) {
  const token = message.match(/\b[a-z]{2,5}-\d+(?:\.\d+)?\b/i) || message.match(/\bWBGene\d+\b/i);
  return token ? token[0] : "";
}

function normalizeName(value) {
  return String(value).toLowerCase().replace(/[^a-z0-9]+/g, " ").trim();
}

function escapeRegex(value) {
  return value.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function parseResultLimit(query) {
  const match = query.match(/top\s+(\d+)|(\d+)\s+(cells|cell types|types)/);
  const value = Number(match?.[1] || match?.[2] || 5);
  return Math.max(1, Math.min(20, value || 5));
}

function cellLabel(row) {
  const preferredColumns = [state.cellIdColumn, "cell_type", "assigned_cell_type", "cell", "cell_id", "barcode", "sample"];
  const labels = preferredColumns
    .filter(unique)
    .map((column) => row[column])
    .filter((value) => value !== undefined && value !== "");
  return String(labels[0] || "unknown cell");
}

function matchClusterNumber(query) {
  const match = query.match(/cluster\s*(\d+)/);
  if (!match) return null;
  return `Cluster ${match[1]}`;
}

function addChat(role, text) {
  const bubble = document.createElement("div");
  bubble.className = `bubble ${role}`;
  bubble.textContent = text;
  $("chatLog").appendChild(bubble);
  $("chatLog").scrollTop = $("chatLog").scrollHeight;
}

function groupValues() {
  const column = controls.groupColumn.value || state.groupColumn;
  return [...new Set(state.rows.map((row) => String(row[column] || row.cluster || "NA")))].filter(Boolean);
}

function clusterLabels() {
  return [...new Set(state.rows.map((row) => row.cluster).filter(Boolean))].sort();
}

function mean(values) {
  return values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : 0;
}

function median(values) {
  if (!values.length) return 0;
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
}

function dot(a, b) {
  return a.reduce((sum, value, index) => sum + value * b[index], 0);
}

function normalize(vector) {
  const length = Math.sqrt(dot(vector, vector)) || 1;
  return vector.map((value) => value / length);
}

function extent(values) {
  const min = Math.min(...values);
  const max = Math.max(...values);
  return min === max ? [min - 1, max + 1] : [min, max];
}

function scale(value, domain, range) {
  return range[0] + ((value - domain[0]) / (domain[1] - domain[0])) * (range[1] - range[0]);
}

function interpolateColor(t) {
  const clamped = Math.max(0, Math.min(1, t));
  const start = [219, 231, 235];
  const end = [186, 74, 104];
  const rgb = start.map((value, index) => Math.round(value + (end[index] - value) * clamped));
  return `rgb(${rgb.join(",")})`;
}

function formatNumber(value) {
  return new Intl.NumberFormat().format(value);
}

function makeExampleDataset() {
  const genes = ["CD3D", "IL7R", "MS4A1", "CD79A", "LYZ", "S100A8", "NKG7", "GNLY", "PPBP", "FCGR3A"];
  const cellTypes = ["T_cells", "B_cells", "Monocytes", "NK_cells", "Platelets"];
  const header = ["gene", ...cellTypes];
  const rows = [header.join(",")];
  const programs = {
    T_cells: ["CD3D", "IL7R"],
    B_cells: ["MS4A1", "CD79A"],
    Monocytes: ["LYZ", "S100A8", "FCGR3A"],
    NK_cells: ["NKG7", "GNLY"],
    Platelets: ["PPBP"],
  };
  genes.forEach((gene, geneIndex) => {
    const values = cellTypes.map((cellType, cellIndex) => {
      const base = ((geneIndex + 2) * (cellIndex + 3)) % 4;
      return programs[cellType].includes(gene) ? base + 10 : base + 1;
    });
    rows.push([gene, ...values].join(","));
  });
  return rows.join("\n");
}

schedulePlotRender();
refreshServerDatasets();
if ($("multimodalSummary")) loadMultimodalIndex();
