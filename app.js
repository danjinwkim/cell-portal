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
  markers: new Map(),
};

const palette = ["#1f7a8c", "#d88742", "#ba4a68", "#4d8d62", "#7353ba", "#247ba0", "#c44536", "#5b8e7d"];
const numericThreshold = 0.9;

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

async function uploadH5ad(file) {
  $("analysisStatus").textContent = "Reading H5AD";
  addChat("assistant", `Uploading ${file.name}. I will load a fast browser-ready H5AD view with readable gene names.`);
  const formData = new FormData();
  formData.append("file", file);

  try {
    const response = await fetch("/api/upload-h5ad", {
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
    $("analysisStatus").textContent = "H5AD backend unavailable";
    addChat(
      "assistant",
      `I can accept .h5ad files when the Python backend is running. Start it with the README instructions, then try again. Details: ${error.message}`
    );
  }
}

async function refreshServerDatasets() {
  try {
    const response = await fetch("/api/datasets");
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
  } catch {
    fillSelect(controls.serverDataset, ["Backend not ready"], "Backend not ready");
    controls.loadServerDataset.disabled = true;
  }
}

async function loadServerDataset(datasetId) {
  $("analysisStatus").textContent = "Loading downloaded H5AD";
  addChat("assistant", "Loading the server-side H5AD view directly. This skips the browser re-upload step.");
  try {
    const response = await fetch(`/api/datasets/${encodeURIComponent(datasetId)}`);
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
      line: { color: "rgba(29,39,51,0.25)", width: 0.5 },
      opacity: 0.86,
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
  if (!state.rows.length) return "Upload or load a dataset first, then I can query it.";
  const query = message.toLowerCase();
  const cluster = clusterLabels().find((label) => query.includes(label.toLowerCase())) || matchClusterNumber(query);
  const genes = extractGenes(message);
  const gene = genes[0];

  if (isGeneInfoQuestion(query)) {
    const targetGene = gene || extractGeneLikeToken(message);
    if (!targetGene) return "Tell me which C. elegans gene you want WormBase information for, for example: describe daf-2.";
    return await answerWormBaseQuestion(targetGene, query);
  }

  if (/(highest|high|most|top|max|maximum|enriched|strongest).*(express|expression)|which .*express|express.*highest/.test(query) && genes.length) {
    return answerHighestExpression(genes, parseResultLimit(query));
  }

  if (/marker|top gene|signature/.test(query)) {
    const label = cluster || controls.markerCluster.value || clusterLabels()[0];
    controls.markerCluster.value = label;
    renderMarkerTable();
    const markers = (state.markers.get(label) || []).slice(0, 5).map((row) => row.gene).join(", ");
    return `Top marker genes for ${label}: ${markers}. I also updated the marker table.`;
  }

  if (/express|expression|show|plot|color/.test(query) && gene) {
    controls.geneOverlay.value = gene;
    schedulePlotRender();
    const positive = state.rows.filter((row) => Number(row[gene]) > 0).length;
    return `${gene} is detected in ${positive} of ${state.rows.length} cells. I colored the projection by ${gene}.`;
  }

  if (/compare|differential|versus| vs /.test(query)) {
    renderDifferentialTable();
    const first = $("deTable tr td");
    return first ? `I ran the current group comparison. The strongest changing gene is ${first.textContent}.` : "I need at least two groups to compare.";
  }

  if (/summary|overview|how many/.test(query)) {
    return `This dataset has ${state.rows.length} cells, ${state.genes.length} numeric gene columns, and ${clusterLabels().length} computed clusters.`;
  }

  if (/cluster/.test(query)) {
    const counts = clusterLabels().map((label) => `${label}: ${state.rows.filter((row) => row.cluster === label).length}`).join("; ");
    return `Cluster sizes are ${counts}. Ask for marker genes for a specific cluster to drill in.`;
  }

  return "I can answer questions about summaries, clusters, marker genes, gene expression overlays, highest-expressing cells, WormBase gene annotations, and group comparisons.";
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
    const response = await fetch(`/api/wormbase/gene/${encodeURIComponent(gene)}?gene_id=${encodeURIComponent(metadata.id)}`);
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
