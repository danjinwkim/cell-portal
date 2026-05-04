from __future__ import annotations

import csv
import io
import json
import math
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.error import URLError
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from scipy.io import loadmat


USER_AGENT = "CellPortal/0.2 multimodal-ingest"
SOURCE_TIMEOUT_SECONDS = 25


@dataclass(frozen=True)
class SourceSpec:
    name: str
    url: str
    modality: str
    parser: str
    required: bool = False


DEFAULT_SOURCES = [
    SourceSpec(
        "OpenWorm CElegansNeuroML edge list",
        "https://raw.githubusercontent.com/openworm/CElegansNeuroML/master/herm_full_edgelist.csv",
        "connectome",
        "edge_table",
    ),
    SourceSpec(
        "OpenWorm Varshney NeuronConnectFormatted",
        "https://raw.githubusercontent.com/openworm/CElegansNeuroML/master/NeuronConnectFormatted.xlsx",
        "connectome",
        "table_or_matrix",
    ),
    SourceSpec(
        "OpenWorm CElegansNeuronTables",
        "https://raw.githubusercontent.com/openworm/CElegansNeuroML/master/CElegansNeuronTables.xls",
        "connectome_spatial",
        "table_or_matrix",
    ),
    SourceSpec(
        "WormAtlas neuronal wiring page",
        "https://www.wormatlas.org/neuronalwiring.html",
        "connectome_spatial_lineage",
        "wormatlas_page",
    ),
    SourceSpec(
        "Figshare C. elegans lineage tree",
        "https://figshare.com/articles/dataset/C_elegans_lineage_tree/22491361",
        "lineage",
        "figshare_landing",
    ),
]


def build_unified_dataset(data_dir: Path, sources: list[SourceSpec] | None = None) -> dict[str, Any]:
    """Fetch and normalize public multimodal neuron datasets into a JSON-ready graph index."""
    data_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = data_dir / "raw"
    raw_dir.mkdir(exist_ok=True)
    neurons: dict[str, dict[str, Any]] = {}
    provenance: list[dict[str, Any]] = []

    for source in sources or DEFAULT_SOURCES:
        try:
            fetched = fetch_source(source, raw_dir)
            parsed = parse_source(source, fetched)
            merge_records(neurons, parsed)
            provenance.append({"name": source.name, "url": source.url, "status": "ok", "records": len(parsed)})
        except Exception as exc:  # noqa: BLE001 - ingestion should degrade source-by-source
            provenance.append({"name": source.name, "url": source.url, "status": "failed", "error": str(exc)})
            if source.required:
                raise

    if not neurons:
        neurons = fallback_demo_records()
        provenance.append({"name": "built-in demo fallback", "url": "", "status": "ok", "records": len(neurons)})

    normalized = [finalize_record(record) for record in neurons.values()]
    normalized.sort(key=lambda item: item["id"])
    payload = {
        "schema": "cell-portal-multimodal-v1",
        "neuron_count": len(normalized),
        "edge_count": sum(len(item["connections"]) for item in normalized),
        "neurons": normalized,
        "provenance": provenance,
    }
    (data_dir / "unified_neurons.json").write_text(json.dumps(payload, indent=2))
    return payload


def fetch_source(source: SourceSpec, raw_dir: Path) -> Path:
    target = raw_dir / safe_filename(source.name, source.url)
    if target.exists() and target.stat().st_size:
        return target
    request = Request(source.url, headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(request, timeout=SOURCE_TIMEOUT_SECONDS) as response:
            target.write_bytes(response.read())
    except URLError as exc:
        raise RuntimeError(f"Could not fetch {source.url}: {exc}") from exc
    return target


def parse_source(source: SourceSpec, path: Path) -> list[dict[str, Any]]:
    if source.parser == "edge_table":
        return parse_edge_table(read_table(path))
    if source.parser == "table_or_matrix":
        return parse_workbook_or_table(path)
    if source.parser == "wormatlas_page":
        return parse_wormatlas_page(path)
    if source.parser == "figshare_landing":
        return parse_figshare_page(path)
    if path.suffix.lower() == ".mat":
        return parse_mat_positions(path)
    return []


def parse_workbook_or_table(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for table in read_all_tables(path):
        records.extend(parse_edge_table(table))
        records.extend(parse_spatial_table(table))
        records.extend(parse_lineage_table(table))
        if not records:
            records.extend(parse_adjacency_matrix(table))
    return records


def parse_wormatlas_page(path: Path) -> list[dict[str, Any]]:
    html = path.read_text(errors="ignore")
    records: list[dict[str, Any]] = []
    for href in re.findall(r'href=["\']([^"\']+)["\']', html, flags=re.I):
        if not re.search(r"(Neuron.*Fixed|Lineage|wiring|connect).*\.(csv|xls|xlsx)", href, flags=re.I):
            continue
        url = href if href.startswith("http") else f"https://www.wormatlas.org/{href.lstrip('/')}"
        try:
            local = fetch_source(SourceSpec(f"WormAtlas discovered {Path(href).name}", url, "discovered", "table_or_matrix"), path.parent)
            records.extend(parse_workbook_or_table(local))
        except Exception:
            continue
    return records


def parse_figshare_page(path: Path) -> list[dict[str, Any]]:
    html = path.read_text(errors="ignore")
    records: list[dict[str, Any]] = []
    for url in sorted(set(re.findall(r'https://ndownloader\.figshare\.com/files/\d+', html))):
        try:
            local = fetch_source(SourceSpec(f"Figshare file {url.rsplit('/', 1)[-1]}", url, "lineage", "table_or_matrix"), path.parent)
            if local.suffix.lower() == ".gml":
                records.extend(parse_gml_lineage(local))
            elif zipfile.is_zipfile(local):
                records.extend(parse_zip_archive(local))
            else:
                records.extend(parse_workbook_or_table(local))
        except Exception:
            continue
    return records


def parse_zip_archive(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if not re.search(r"\.(csv|tsv|txt|gml|mat)$", name, flags=re.I):
                continue
            extracted = path.parent / safe_filename(name, name)
            extracted.write_bytes(archive.read(name))
            if name.lower().endswith(".gml"):
                records.extend(parse_gml_lineage(extracted))
            elif name.lower().endswith(".mat"):
                records.extend(parse_mat_positions(extracted))
            else:
                records.extend(parse_workbook_or_table(extracted))
    return records


def read_all_tables(path: Path) -> list[pd.DataFrame]:
    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        return list(pd.read_excel(path, sheet_name=None).values())
    return [read_table(path)]


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    if suffix in (".xlsx", ".xls"):
        return pd.read_excel(path)
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.read_csv(path, sep=None, engine="python")


def parse_edge_table(table: pd.DataFrame) -> list[dict[str, Any]]:
    columns = normalized_columns(table)
    source_col = find_column(columns, "source", "pre", "from", "neuron1", "n1", "sender")
    target_col = find_column(columns, "target", "post", "to", "neuron2", "n2", "receiver")
    if not source_col or not target_col:
        return []
    type_col = find_column(columns, "type", "synapse_type", "connection_type")
    weight_col = find_column(columns, "weight", "synapses", "count", "number", "synaptic_weight")
    records = []
    for _, row in table.iterrows():
        source = normalize_neuron_name(row.get(source_col))
        target = normalize_neuron_name(row.get(target_col))
        if not source or not target:
            continue
        kind = normalize_connection_type(row.get(type_col, "chemical"))
        records.append(
            {
                "id": source,
                "connections": [
                    {
                        "target": target,
                        "type": kind,
                        "weight": float_or_zero(row.get(weight_col, 1) if weight_col else 1),
                    }
                ],
            }
        )
        records.append({"id": target})
    return records


def parse_adjacency_matrix(table: pd.DataFrame) -> list[dict[str, Any]]:
    if table.shape[0] < 2 or table.shape[1] < 2:
        return []
    first = table.columns[0]
    targets = [normalize_neuron_name(col) for col in table.columns[1:]]
    records = []
    for _, row in table.iterrows():
        source = normalize_neuron_name(row.get(first))
        if not source:
            continue
        edges = []
        for col, target in zip(table.columns[1:], targets):
            weight = float_or_zero(row.get(col))
            if target and weight > 0:
                edges.append({"target": target, "type": "chemical", "weight": weight})
        if edges:
            records.append({"id": source, "connections": edges})
    return records


def parse_spatial_table(table: pd.DataFrame) -> list[dict[str, Any]]:
    columns = normalized_columns(table)
    neuron_col = find_column(columns, "neuron", "neuron_id", "cell", "cell_name", "name")
    x_col = find_column(columns, "x", "ap", "position", "cellbody_position")
    y_col = find_column(columns, "y", "dv")
    z_col = find_column(columns, "z", "lr")
    if not neuron_col or not x_col:
        return []
    records = []
    for _, row in table.iterrows():
        neuron_id = normalize_neuron_name(row.get(neuron_col))
        if not neuron_id:
            continue
        records.append(
            {
                "id": neuron_id,
                "coordinates": {
                    "x": float_or_none(row.get(x_col)),
                    "y": float_or_none(row.get(y_col)) if y_col else None,
                    "z": float_or_none(row.get(z_col)) if z_col else None,
                },
            }
        )
    return records


def parse_lineage_table(table: pd.DataFrame) -> list[dict[str, Any]]:
    columns = normalized_columns(table)
    neuron_col = find_column(columns, "neuron", "neuron_id", "cell", "cell_name", "name")
    parent_col = find_column(columns, "parent", "parent_cell", "mother")
    path_col = find_column(columns, "lineage", "lineage_path", "ancestry", "path")
    if not neuron_col or not (parent_col or path_col):
        return []
    records = []
    for _, row in table.iterrows():
        neuron_id = normalize_neuron_name(row.get(neuron_col))
        if not neuron_id:
            continue
        lineage_path = split_lineage(row.get(path_col)) if path_col else []
        parent = str(row.get(parent_col) or "").strip() if parent_col else (lineage_path[-1] if lineage_path else "")
        records.append({"id": neuron_id, "lineage": {"parent": parent, "lineage_path": lineage_path}})
    return records


def parse_mat_positions(path: Path) -> list[dict[str, Any]]:
    mat = loadmat(path, squeeze_me=True, struct_as_record=False)
    records = []
    names = first_array(mat, "names", "neurons", "labels", "cell_names")
    coords = first_array(mat, "coords", "coordinates", "positions", "xyz")
    if names is None or coords is None:
        return records
    coords = np.asarray(coords)
    for name, coord in zip(np.ravel(names), coords):
        neuron_id = normalize_neuron_name(name)
        values = np.ravel(coord).astype(float)
        if neuron_id and len(values):
            records.append({"id": neuron_id, "coordinates": {"x": float(values[0]), "y": float(values[1]) if len(values) > 1 else None, "z": float(values[2]) if len(values) > 2 else None}})
    return records


def parse_gml_lineage(path: Path) -> list[dict[str, Any]]:
    try:
        import networkx as nx
    except ImportError:
        return []
    graph = nx.read_gml(path)
    records = []
    for node in graph.nodes:
        neuron_id = normalize_neuron_name(node)
        if not neuron_id:
            continue
        predecessors = [str(item) for item in graph.predecessors(node)]
        ancestors = [str(item) for item in nx.ancestors(graph, node)]
        records.append({"id": neuron_id, "lineage": {"parent": predecessors[0] if predecessors else "", "lineage_path": ancestors}})
    return records


def merge_records(neurons: dict[str, dict[str, Any]], records: list[dict[str, Any]]) -> None:
    for record in records:
        neuron_id = normalize_neuron_name(record.get("id"))
        if not neuron_id:
            continue
        target = neurons.setdefault(neuron_id, empty_record(neuron_id))
        if record.get("type"):
            target["type"] = record["type"]
        if record.get("coordinates"):
            target["coordinates"] = merge_coordinates(target.get("coordinates", {}), record["coordinates"])
        if record.get("lineage"):
            target["lineage"].update({key: value for key, value in record["lineage"].items() if value not in ("", None, [])})
        for edge in record.get("connections", []):
            edge_target = normalize_neuron_name(edge.get("target"))
            if not edge_target:
                continue
            neurons.setdefault(edge_target, empty_record(edge_target))
            upsert_edge(target["connections"], edge_target, edge.get("type", "chemical"), float_or_zero(edge.get("weight")))


def empty_record(neuron_id: str) -> dict[str, Any]:
    return {"id": neuron_id, "type": "", "coordinates": {}, "lineage": {}, "connections": []}


def finalize_record(record: dict[str, Any]) -> dict[str, Any]:
    lineage = record.get("lineage") or {}
    path = lineage.get("lineage_path") or []
    root = lineage.get("lineage_root") or (path[0] if path else "")
    return {
        "id": record["id"],
        "type": record.get("type") or "unannotated",
        "coordinates": {
            "x": record.get("coordinates", {}).get("x"),
            "y": record.get("coordinates", {}).get("y"),
            "z": record.get("coordinates", {}).get("z"),
        },
        "lineage": {
            "parent": lineage.get("parent", ""),
            "lineage_depth": int(lineage.get("lineage_depth") or len(path)),
            "lineage_root": root,
            "lineage_path": path,
        },
        "connections": sorted(record.get("connections", []), key=lambda item: item["weight"], reverse=True),
    }


def upsert_edge(edges: list[dict[str, Any]], target: str, kind: str, weight: float) -> None:
    kind = normalize_connection_type(kind)
    for edge in edges:
        if edge["target"] == target and edge["type"] == kind:
            edge["weight"] += weight
            return
    edges.append({"target": target, "type": kind, "weight": weight})


def query_unified_dataset(dataset: dict[str, Any], text: str) -> dict[str, Any]:
    index = {item["id"]: item for item in dataset.get("neurons", [])}
    query = text.lower()
    anchor = resolve_anchor(index, text)
    needs_spatial = any(term in query for term in ("near", "nearest", "spatial", "distance", "coordinate", "position"))
    needs_lineage = any(term in query for term in ("lineage", "ancestor", "parent", "related"))
    needs_connectome = any(term in query for term in ("connect", "synap", "partner", "gap", "edge", "weight"))
    radius = parse_radius(query, default=15.0)
    min_weight = parse_min_weight(query, default=5.0 if "strong" in query else 1.0)
    steps = []
    candidates = set(index)

    if needs_spatial and anchor:
        steps.append(f"Find neurons within radius {radius:g} of {anchor} using Euclidean coordinate distance.")
        candidates &= {row["id"] for row in nearby(index, anchor, radius)}
    if needs_lineage and anchor:
        steps.append(f"Keep neurons sharing a lineage ancestor with {anchor}.")
        candidates &= set(shared_lineage(index, anchor))
    if needs_connectome and anchor:
        steps.append(f"Rank remaining neurons by chemical + gap-junction weight connected to {anchor}; minimum weight {min_weight:g}.")
        connected = {edge["target"]: edge for edge in index.get(anchor, {}).get("connections", []) if edge.get("weight", 0) >= min_weight}
        if needs_spatial or needs_lineage:
            candidates &= set(connected)
        else:
            candidates = set(connected)

    rows = []
    for neuron_id in candidates:
        if neuron_id == anchor:
            continue
        row = index[neuron_id]
        distance = distance_between(index, anchor, neuron_id) if anchor else None
        connected_weight = connection_weight(index, anchor, neuron_id) if anchor else 0.0
        rows.append(
            {
                "neuron_id": neuron_id,
                "type": row.get("type", ""),
                "distance": None if distance is None or math.isinf(distance) else round(distance, 3),
                "shared_lineage_depth": shared_lineage_depth(index.get(anchor), row) if anchor else 0,
                "connection_weight": round(connected_weight, 3),
            }
        )
    rows.sort(key=lambda item: (-item["connection_weight"], item["distance"] is None, item["distance"] or 999999, item["neuron_id"]))
    cypher = cypher_plan(anchor, needs_spatial, needs_lineage, needs_connectome, radius, min_weight)
    return {
        "modalities": [name for name, enabled in (("spatial", needs_spatial), ("lineage", needs_lineage), ("connectome", needs_connectome)) if enabled],
        "anchor": anchor,
        "query_plan": steps or ["No specific modality keywords were found; return indexed neuron metadata."],
        "structured_query": cypher,
        "results": rows[:20],
        "result_format": {"neuron_id": "string", "type": "string", "distance": "float|null", "shared_lineage_depth": "int", "connection_weight": "float"},
    }


def nearby(index: dict[str, dict[str, Any]], anchor: str, radius: float) -> list[dict[str, Any]]:
    return [
        {"id": neuron_id, "distance": distance_between(index, anchor, neuron_id)}
        for neuron_id in index
        if neuron_id != anchor and distance_between(index, anchor, neuron_id) <= radius
    ]


def shared_lineage(index: dict[str, dict[str, Any]], anchor: str) -> list[str]:
    source = index.get(anchor)
    return [neuron_id for neuron_id, row in index.items() if neuron_id != anchor and shared_lineage_depth(source, row) > 0]


def distance_between(index: dict[str, dict[str, Any]], left_id: str, right_id: str) -> float:
    left = (index.get(left_id) or {}).get("coordinates") or {}
    right = (index.get(right_id) or {}).get("coordinates") or {}
    if left.get("x") is None or right.get("x") is None:
        return math.inf
    return math.sqrt(sum(((left.get(axis) or 0) - (right.get(axis) or 0)) ** 2 for axis in ("x", "y", "z")))


def connection_weight(index: dict[str, dict[str, Any]], source: str, target: str) -> float:
    return sum(edge.get("weight", 0) for edge in (index.get(source) or {}).get("connections", []) if edge.get("target") == target)


def shared_lineage_depth(left: dict[str, Any] | None, right: dict[str, Any] | None) -> int:
    if not left or not right:
        return 0
    a = left.get("lineage", {}).get("lineage_path", [])
    b = right.get("lineage", {}).get("lineage_path", [])
    count = 0
    for one, two in zip(a, b):
        if one != two:
            break
        count += 1
    return count


def cypher_plan(anchor: str | None, spatial: bool, lineage: bool, connectome: bool, radius: float, min_weight: float) -> str:
    where = []
    if spatial:
        where.append(f"distance(n.coordinates, anchor.coordinates) <= {radius:g}")
    if lineage:
        where.append("any(a IN n.lineage.lineage_path WHERE a IN anchor.lineage.lineage_path)")
    if connectome:
        where.append(f"edge.weight >= {min_weight:g}")
    match = "MATCH (anchor:Neuron {id: $anchor})"
    if connectome:
        match += "-[edge:SYNAPSE|GAP_JUNCTION]-(n:Neuron)"
    else:
        match += ", (n:Neuron)"
    return f"{match}\nWHERE {' AND '.join(where) if where else 'true'}\nRETURN n.id, n.type, n.coordinates, n.lineage, coalesce(edge.weight, 0) AS connection_weight\nORDER BY connection_weight DESC"


def resolve_anchor(index: dict[str, dict[str, Any]], text: str) -> str | None:
    for token in re.findall(r"\b[A-Za-z]{2,5}[LR]?\d*\b", text):
        normalized = normalize_neuron_name(token)
        if normalized in index:
            return normalized
        if normalized.endswith(("L", "R")) and normalized[:-1] in index:
            return normalized[:-1]
    return next(iter(index), None)


def normalized_columns(table: pd.DataFrame) -> dict[str, str]:
    return {normalize_column(column): column for column in table.columns}


def find_column(columns: dict[str, str], *names: str) -> str | None:
    wanted = {normalize_column(name) for name in names}
    for normalized, original in columns.items():
        if normalized in wanted or any(name in normalized for name in wanted):
            return original
    return None


def normalize_column(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def normalize_neuron_name(value: Any) -> str:
    text = str(value or "").strip().upper()
    text = re.sub(r"[^A-Z0-9]", "", text)
    return text if re.match(r"^[A-Z]{1,5}[LR]?\d*$", text) else ""


def normalize_connection_type(value: Any) -> str:
    text = str(value or "").lower()
    if "gap" in text or "elect" in text or "junction" in text:
        return "gap_junction"
    return "chemical"


def float_or_zero(value: Any) -> float:
    parsed = float_or_none(value)
    return parsed if parsed is not None else 0.0


def float_or_none(value: Any) -> float | None:
    try:
        if value is None or pd.isna(value):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def split_lineage(value: Any) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [part for part in re.split(r"\s*(?:>|/|,|;|\|)\s*", str(value).strip()) if part]


def parse_radius(query: str, default: float) -> float:
    match = re.search(r"(?:within|radius|near)\s+(\d+(?:\.\d+)?)", query)
    return float(match.group(1)) if match else default


def parse_min_weight(query: str, default: float) -> float:
    match = re.search(r"(?:weight|strong(?:ly)?|synapses?)\s*(?:>=|over|above)?\s*(\d+(?:\.\d+)?)", query)
    return float(match.group(1)) if match else default


def first_array(mapping: dict[str, Any], *names: str) -> Any:
    for name in names:
        if name in mapping:
            return mapping[name]
    return None


def merge_coordinates(left: dict[str, Any], right: dict[str, Any]) -> dict[str, Any]:
    return {axis: right.get(axis) if right.get(axis) is not None else left.get(axis) for axis in ("x", "y", "z")}


def safe_filename(name: str, url: str) -> str:
    suffix = Path(url.split("?", 1)[0]).suffix or ".html"
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_")[:90]
    return f"{stem}{suffix}" if not stem.endswith(suffix) else stem


def fallback_demo_records() -> dict[str, dict[str, Any]]:
    from multimodal_store import demo_neurons

    records: dict[str, dict[str, Any]] = {}
    for neuron in demo_neurons().values():
        records[neuron.neuron_id] = {
            "id": neuron.neuron_id,
            "type": neuron.neuron_class,
            "coordinates": neuron.spatial_coordinates,
            "lineage": {"parent": neuron.parent_cell, "lineage_path": neuron.lineage_path},
            "connections": [
                {"target": edge.target, "type": "chemical", "weight": edge.chemical_weight}
                for edge in neuron.connectome_edges
                if edge.chemical_weight
            ]
            + [
                {"target": edge.target, "type": "gap_junction", "weight": edge.gap_weight}
                for edge in neuron.connectome_edges
                if edge.gap_weight
            ],
        }
    return records
