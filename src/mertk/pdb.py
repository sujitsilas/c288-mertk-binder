"""RCSB assembly discovery and downloads with recorded provenance."""

import hashlib
import json
import re
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
TAM_ACCESSIONS = {"MERTK": "Q12866", "AXL": "P30530", "TYRO3": "Q06418"}


def session():
    client = requests.Session()
    retry = Retry(
        total=4,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
    )
    client.mount("https://", HTTPAdapter(max_retries=retry))
    client.headers["User-Agent"] = "c288-mertk-binder/0.1 (academic research)"
    return client


def terminal(attribute, operator, value):
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {"attribute": attribute, "operator": operator, "value": value},
    }


def assembly_query(resolution=3.0):
    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                terminal(
                    "rcsb_assembly_info.polymer_entity_instance_count_protein",
                    "greater_or_equal",
                    2,
                ),
                terminal(
                    "rcsb_entry_info.selected_polymer_entity_types", "exact_match", "Protein (only)"
                ),
                terminal("rcsb_entry_info.resolution_combined", "less_or_equal", resolution),
                terminal("exptl.method", "in", ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]),
            ],
        },
        "return_type": "assembly",
        "request_options": {
            "paginate": {"start": 0, "rows": 100},
            "sort": [{"sort_by": "rcsb_id", "direction": "asc"}],
        },
    }


def discover(client, query, limit):
    if limit < 1:
        raise ValueError("limit must be positive")
    identifiers = []
    while len(identifiers) < limit:
        query["request_options"]["paginate"] = {
            "start": len(identifiers),
            "rows": min(100, limit - len(identifiers)),
        }
        response = client.post(SEARCH_URL, json=query, timeout=60)
        response.raise_for_status()
        if response.status_code == 204:
            break
        hits = response.json().get("result_set", [])
        if not hits:
            break
        identifiers.extend(hit["identifier"] for hit in hits)
        if len(identifiers) >= response.json().get("total_count", 0):
            break
    return identifiers


def download_assembly(client, identifier, output):
    if not re.fullmatch(r"[A-Za-z0-9]{4}-[0-9]+", identifier):
        raise ValueError(f"Expected a legacy PDB assembly identifier such as 1BRS-1: {identifier}")
    pdb_id, assembly_id = identifier.upper().split("-")
    url = f"https://files.rcsb.org/download/{pdb_id}-assembly{assembly_id}.cif.gz"
    path = Path(output) / f"{pdb_id}-{assembly_id}.cif.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        temporary = path.with_suffix(".part")
        with client.get(url, stream=True, timeout=120) as response:
            response.raise_for_status()
            with temporary.open("wb") as handle:
                for chunk in response.iter_content(1024 * 1024):
                    handle.write(chunk)
        temporary.replace(path)
    with path.open("rb") as handle:
        sha256 = hashlib.file_digest(handle, "sha256").hexdigest()
    return {"assembly": identifier.upper(), "path": str(path), "url": url, "sha256": sha256}


def download_tam_seeds(client, output):
    """Include full sequences and annotated domains so fragments can match the holdout."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    fasta = []
    for gene, accession in TAM_ACCESSIONS.items():
        response = client.get(f"https://rest.uniprot.org/uniprotkb/{accession}.json", timeout=60)
        response.raise_for_status()
        record = response.json()
        (output / f"{accession}.json").write_text(json.dumps(record, indent=2) + "\n")
        sequence = record["sequence"]["value"]
        fasta.append(f">TAM_{gene}_full\n{sequence}\n")
        for i, feature in enumerate(record.get("features", [])):
            is_domain = feature["type"] == "Domain"
            extracellular = (
                feature["type"] == "Topological domain"
                and "extracellular" in feature.get("description", "").lower()
            )
            if not (is_domain or extracellular):
                continue
            start = feature["location"]["start"].get("value")
            end = feature["location"]["end"].get("value")
            if start and end:
                fasta.append(f">TAM_{gene}_region{i}\n{sequence[start - 1 : end]}\n")
    path = output / "tam_seeds.fasta"
    path.write_text("".join(fasta))
    return path
