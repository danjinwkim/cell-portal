from __future__ import annotations

import csv
import os
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


DEFAULT_DATA_DIR = Path("/tmp/cell_portal_multimodal") if os.environ.get("VERCEL") else Path(__file__).parent / "multimodal_data"
DATA_DIR = Path(os.environ.get("MULTIMODAL_DATA_DIR", DEFAULT_DATA_DIR))


@dataclass
class ConnectomeEdge:
    target: str
    chemical_weight: float = 0.0
    gap_weight: float = 0.0


@dataclass
class NeuronRecord:
    neuron_id: str
    neuron_class: str = "unannotated"
    cluster: str = ""
    transcriptomic_neighbors: list[str] = field(default_factory=list)
    connectome_edges: list[ConnectomeEdge] = field(default_factory=list)
    spatial_coordinates: dict[str, float] = field(default_factory=dict)
    lineage_path: list[str] = field(default_factory=list)
    parent_cell: str = ""


class MultimodalStore:
    def __init__(self, data_dir: Path = DATA_DIR) -> None:
        self.data_dir = data_dir
        self.neurons: dict[str, NeuronRecord] = {}
        self.load()

    def load(self) -> None:
        self.neurons = demo_neurons()
        self.load_unified_json()
        self.load_neuron_metadata()
        self.load_connectome_edges()
        self.load_spatial_coordinates()
        self.load_lineage()

    def load_unified_json(self) -> None:
        path = self.data_dir / "unified_neurons.json"
        if not path.exists():
            return
        try:
            payload = json.loads(path.read_text())
        except json.JSONDecodeError:
            return
        for row in payload.get("neurons", []):
            record = self.ensure(row.get("id", ""))
            record.neuron_class = row.get("type") or record.neuron_class
            coords = row.get("coordinates") or {}
            if coords.get("x") is not None:
                record.spatial_coordinates = {axis: float(coords.get(axis) or 0) for axis in ("x", "y", "z")}
            lineage = row.get("lineage") or {}
            record.parent_cell = lineage.get("parent") or record.parent_cell
            record.lineage_path = lineage.get("lineage_path") or record.lineage_path
            for edge in row.get("connections", []):
                target = normalize_id(edge.get("target"))
                if not target:
                    continue
                self.ensure(target)
                if edge.get("type") == "gap_junction":
                    record.connectome_edges.append(ConnectomeEdge(target, gap_weight=float(edge.get("weight") or 0)))
                else:
                    record.connectome_edges.append(ConnectomeEdge(target, chemical_weight=float(edge.get("weight") or 0)))

    def all_neurons(self) -> list[dict[str, Any]]:
        return [serialize_neuron(record) for record in sorted(self.neurons.values(), key=lambda item: item.neuron_id)]

    def unified_dataset(self) -> dict[str, Any]:
        neurons = []
        for record in sorted(self.neurons.values(), key=lambda item: item.neuron_id):
            lineage_path = record.lineage_path or []
            connections = []
            for edge in record.connectome_edges:
                if edge.chemical_weight:
                    connections.append({"target": edge.target, "type": "chemical", "weight": edge.chemical_weight})
                if edge.gap_weight:
                    connections.append({"target": edge.target, "type": "gap_junction", "weight": edge.gap_weight})
            neurons.append(
                {
                    "id": record.neuron_id,
                    "type": record.neuron_class,
                    "coordinates": {
                        "x": record.spatial_coordinates.get("x"),
                        "y": record.spatial_coordinates.get("y"),
                        "z": record.spatial_coordinates.get("z"),
                    },
                    "lineage": {
                        "parent": record.parent_cell,
                        "lineage_depth": len(lineage_path),
                        "lineage_root": lineage_path[0] if lineage_path else "",
                        "lineage_path": lineage_path,
                    },
                    "connections": connections,
                }
            )
        return {
            "schema": "cell-portal-multimodal-v1",
            "neuron_count": len(neurons),
            "edge_count": sum(len(item["connections"]) for item in neurons),
            "neurons": neurons,
        }

    def get_neuron(self, neuron_id: str) -> dict[str, Any] | None:
        record = self.neurons.get(self.resolve_id(neuron_id))
        return serialize_neuron(record) if record else None

    def connectivity(self, neuron_id: str, min_weight: float = 0.0) -> dict[str, Any]:
        record = self.neurons.get(self.resolve_id(neuron_id))
        if not record:
            return {"neuron_id": neuron_id, "neighbors": [], "degree": 0, "weighted_degree": 0}
        neighbors = [
            {
                "target": edge.target,
                "chemical_weight": edge.chemical_weight,
                "gap_weight": edge.gap_weight,
                "total_weight": edge.chemical_weight + edge.gap_weight,
            }
            for edge in record.connectome_edges
            if edge.chemical_weight + edge.gap_weight >= min_weight
        ]
        return {
            "neuron_id": record.neuron_id,
            "neighbors": sorted(neighbors, key=lambda item: item["total_weight"], reverse=True),
            "degree": len(neighbors),
            "weighted_degree": round(sum(item["total_weight"] for item in neighbors), 3),
        }

    def nearest(self, neuron_id: str, limit: int = 8) -> dict[str, Any]:
        source = self.neurons.get(self.resolve_id(neuron_id))
        if not source or not source.spatial_coordinates:
            return {"neuron_id": neuron_id, "neighbors": []}
        distances = []
        for record in self.neurons.values():
            if record.neuron_id == source.neuron_id or not record.spatial_coordinates:
                continue
            distances.append(
                {
                    "neuron_id": record.neuron_id,
                    "distance": round(euclidean(source.spatial_coordinates, record.spatial_coordinates), 3),
                    "neuron_class": record.neuron_class,
                }
            )
        return {"neuron_id": source.neuron_id, "neighbors": sorted(distances, key=lambda item: item["distance"])[:limit]}

    def lineage(self, neuron_id: str) -> dict[str, Any]:
        record = self.neurons.get(self.resolve_id(neuron_id))
        if not record:
            return {"neuron_id": neuron_id, "ancestors": [], "descendants": [], "similar_lineage": []}
        prefix = record.lineage_path[: max(1, len(record.lineage_path) - 1)]
        similar = [
            item.neuron_id
            for item in self.neurons.values()
            if item.neuron_id != record.neuron_id and item.lineage_path[: len(prefix)] == prefix
        ]
        descendants = [item.neuron_id for item in self.neurons.values() if record.neuron_id in item.lineage_path]
        return {
            "neuron_id": record.neuron_id,
            "parent_cell": record.parent_cell,
            "ancestors": record.lineage_path,
            "descendants": descendants,
            "similar_lineage": similar,
        }

    def natural_language_query(self, text: str) -> dict[str, Any]:
        query = text.lower()
        entities = [name for name in self.neurons if re.search(rf"\b{re.escape(name.lower())}\b", query)]
        anchor = self.resolve_id(entities[0] if entities else infer_anchor(query))
        if "transcription" in query and ("connected" in query or "synap" in query):
            rows = self.transcriptionally_similar_and_connected()
            highlighted = sorted({item["source"] for item in rows} | {item["target"] for item in rows})
            return response("transcriptomic_connectome_overlap", rows, highlighted)
        spatial_terms = ("near", "close", "proximal", "spatial", "position", "coordinate", "distance")
        if any(term in query for term in spatial_terms) and ("synap" in query or "connect" in query):
            result = self.nearby_connected(anchor or "AVA")
            highlighted = [anchor or "AVA", *[item["neuron_id"] for item in result]]
            return response("spatial_connectivity", result, highlighted)
        if any(term in query for term in spatial_terms):
            result = self.nearest(anchor or "AVA")
            return response("spatial_nearest", result["neighbors"], [anchor] + [item["neuron_id"] for item in result["neighbors"]])
        if "lineage" in query or "ancestor" in query or "descendant" in query:
            result = self.lineage(anchor or "AVA")
            return response("lineage", result, [result["neuron_id"], *result.get("similar_lineage", [])])
        if "connect" in query or "synap" in query or "gap junction" in query or "partner" in query:
            result = self.connectivity(anchor or "AVA", min_weight=1)
            highlighted = [result["neuron_id"], *[item["target"] for item in result["neighbors"][:8]]]
            return response("connectivity", result, highlighted)
        result = self.connectivity(anchor or "AVA", min_weight=1)
        return response("summary", result, [result["neuron_id"], *[item["target"] for item in result["neighbors"][:5]]])

    def transcriptionally_similar_and_connected(self) -> list[dict[str, Any]]:
        rows = []
        for record in self.neurons.values():
            transcriptomic = set(record.transcriptomic_neighbors)
            for edge in record.connectome_edges:
                if edge.target in transcriptomic:
                    rows.append(
                        {
                            "source": record.neuron_id,
                            "target": edge.target,
                            "total_weight": edge.chemical_weight + edge.gap_weight,
                            "relationship": "transcriptionally similar and synaptically connected",
                        }
                    )
        return sorted(rows, key=lambda item: item["total_weight"], reverse=True)

    def nearby_connected(self, neuron_id: str, limit: int = 8) -> list[dict[str, Any]]:
        anchor = normalize_id(neuron_id)
        nearest = {item["neuron_id"]: item["distance"] for item in self.nearest(anchor, limit=limit * 2)["neighbors"]}
        connected = {item["target"]: item for item in self.connectivity(anchor, min_weight=1)["neighbors"]}
        rows = []
        for target, edge in connected.items():
            rows.append(
                {
                    "neuron_id": target,
                    "distance": nearest.get(target, None),
                    "total_weight": edge["total_weight"],
                    "chemical_weight": edge["chemical_weight"],
                    "gap_weight": edge["gap_weight"],
                    "relationship": f"near {anchor} and synaptically connected",
                }
            )
        return sorted(rows, key=lambda item: (item["distance"] is None, item["distance"] or 999, -item["total_weight"]))[:limit]

    def load_neuron_metadata(self) -> None:
        path = self.data_dir / "neurons.csv"
        if not path.exists():
            return
        for row in csv.DictReader(path.open()):
            record = self.ensure(row.get("neuron_id", ""))
            record.neuron_class = row.get("neuron_class") or record.neuron_class
            record.cluster = row.get("cluster") or record.cluster

    def load_connectome_edges(self) -> None:
        path = self.data_dir / "connectome_edges.csv"
        if not path.exists():
            return
        for row in csv.DictReader(path.open()):
            source = self.ensure(row.get("source", ""))
            target = normalize_id(row.get("target", ""))
            if target:
                self.ensure(target)
                source.connectome_edges.append(
                    ConnectomeEdge(target, float(row.get("chemical_weight") or 0), float(row.get("gap_weight") or 0))
                )

    def load_spatial_coordinates(self) -> None:
        path = self.data_dir / "spatial_coordinates.csv"
        if not path.exists():
            return
        for row in csv.DictReader(path.open()):
            record = self.ensure(row.get("neuron_id", ""))
            record.spatial_coordinates = {axis: float(row.get(axis) or 0) for axis in ("x", "y", "z")}

    def load_lineage(self) -> None:
        path = self.data_dir / "lineage.json"
        if not path.exists():
            return
        for row in json.loads(path.read_text()):
            record = self.ensure(row.get("neuron_id", ""))
            record.parent_cell = row.get("parent_cell") or ""
            record.lineage_path = row.get("lineage_path") or []

    def ensure(self, neuron_id: str) -> NeuronRecord:
        normalized = normalize_id(neuron_id)
        if normalized not in self.neurons:
            self.neurons[normalized] = NeuronRecord(normalized)
        return self.neurons[normalized]

    def resolve_id(self, neuron_id: str | None) -> str:
        normalized = normalize_id(neuron_id)
        if normalized in self.neurons:
            return normalized
        if normalized.endswith(("L", "R")) and normalized[:-1] in self.neurons:
            return normalized[:-1]
        return normalized


def response(kind: str, result: Any, highlighted: list[str]) -> dict[str, Any]:
    return {
        "kind": kind,
        "result": result,
        "highlighted_neurons": [normalize_id(item) for item in highlighted if item],
    }


def serialize_neuron(record: NeuronRecord) -> dict[str, Any]:
    return {
        "neuron_id": record.neuron_id,
        "neuron_class": record.neuron_class,
        "cluster": record.cluster,
        "transcriptomic_neighbors": record.transcriptomic_neighbors,
        "connectome_edges": [edge.__dict__ for edge in record.connectome_edges],
        "spatial_coordinates": record.spatial_coordinates,
        "lineage_path": record.lineage_path,
        "parent_cell": record.parent_cell,
    }


def normalize_id(value: str | None) -> str:
    return str(value or "").strip().upper()


def infer_anchor(query: str) -> str:
    match = re.search(r"\b[A-Z]{2,4}[LR]?\b", query.upper())
    return match.group(0) if match else "AVA"


def euclidean(left: dict[str, float], right: dict[str, float]) -> float:
    return math.sqrt(sum((left.get(axis, 0) - right.get(axis, 0)) ** 2 for axis in ("x", "y", "z")))


def demo_neurons() -> dict[str, NeuronRecord]:
    rows = [
        ("AVA", "interneuron", "Cluster 1", (-12, 4, 2), ["AB", "ABpl", "ABplp"], "ABplp", ["AVB", "AVE"]),
        ("AVB", "interneuron", "Cluster 1", (-10, 3, 2), ["AB", "ABpr", "ABprp"], "ABprp", ["AVA", "PVC"]),
        ("AVE", "interneuron", "Cluster 1", (-11, 6, 1), ["AB", "ABpl", "ABplp"], "ABplp", ["AVA", "AIB"]),
        ("AIB", "interneuron", "Cluster 2", (-16, 8, 1), ["AB", "ABal", "ABalp"], "ABalp", ["AVE", "AIY"]),
        ("AIY", "interneuron", "Cluster 2", (-18, 9, 1), ["AB", "ABal", "ABalp"], "ABalp", ["AIB", "ASE"]),
        ("ASE", "sensory neuron", "Cluster 3", (-24, 12, 0), ["AB", "ABal", "ABala"], "ABala", ["AIY", "AWA"]),
        ("AWA", "sensory neuron", "Cluster 3", (-26, 15, 1), ["AB", "ABar", "ABara"], "ABara", ["ASE", "AWC"]),
        ("AWC", "sensory neuron", "Cluster 3", (-28, 14, 1), ["AB", "ABar", "ABara"], "ABara", ["AWA", "AIY"]),
        ("PVC", "interneuron", "Cluster 1", (4, -1, 0), ["AB", "ABpr", "ABprp"], "ABprp", ["AVB", "DB"]),
        ("VB", "motor neuron", "Cluster 4", (10, -4, 0), ["P1", "EMS", "MS"], "MS", ["DB", "AVA"]),
        ("DB", "motor neuron", "Cluster 4", (12, -6, 0), ["P1", "EMS", "MS"], "MS", ["VB", "AVB"]),
    ]
    neurons = {}
    for neuron_id, neuron_class, cluster, coords, lineage, parent, transcriptomic in rows:
        neurons[neuron_id] = NeuronRecord(
            neuron_id=neuron_id,
            neuron_class=neuron_class,
            cluster=cluster,
            transcriptomic_neighbors=transcriptomic,
            spatial_coordinates={"x": coords[0], "y": coords[1], "z": coords[2]},
            lineage_path=lineage,
            parent_cell=parent,
        )
    edges = {
        "AVA": [("AVB", 18, 3), ("AVE", 12, 2), ("VB", 9, 0)],
        "AVB": [("AVA", 14, 3), ("PVC", 6, 0), ("DB", 4, 1)],
        "AVE": [("AVA", 10, 2), ("AIB", 8, 1)],
        "AIB": [("AIY", 11, 2), ("AVE", 7, 1)],
        "AIY": [("AIB", 8, 2), ("ASE", 13, 0), ("AWC", 5, 0)],
        "ASE": [("AIY", 10, 0), ("AWA", 3, 1)],
        "AWA": [("AWC", 7, 2), ("ASE", 4, 1)],
        "AWC": [("AWA", 7, 2), ("AIY", 5, 0)],
        "PVC": [("AVB", 6, 0), ("DB", 5, 1)],
        "VB": [("DB", 11, 1), ("AVA", 7, 0)],
        "DB": [("VB", 11, 1), ("AVB", 4, 1)],
    }
    for source, targets in edges.items():
        neurons[source].connectome_edges = [ConnectomeEdge(target, chemical, gap) for target, chemical, gap in targets]
    return neurons


STORE = MultimodalStore()
