"""Group both sides of each complex, all assemblies of an entry, and TAM homologs."""

import hashlib


def assign_splits(records, clusters, validation=0.1, test=0.1, seed=42):
    if validation < 0 or test < 0 or validation + test >= 1:
        raise ValueError("Invalid split fractions")
    parents = {}

    def root(key):
        parents.setdefault(key, key)
        while parents[key] != key:
            parents[key] = parents[parents[key]]
            key = parents[key]
        return key

    def union(left, right):
        left, right = root(left), root(right)
        parents[max(left, right)] = min(left, right)

    if not any(key.startswith("TAM_") for key in clusters):
        raise ValueError("Cluster input must include the TAM holdout seed sequences")
    for row in records:
        for key in (row["binder_key"], row["target_key"]):
            if key not in clusters:
                raise ValueError(f"Missing cluster assignment: {key}")
            union("cluster:" + clusters[key], "entry:" + row["pdb_id"])
    tam_roots = {
        root("cluster:" + cluster) for key, cluster in clusters.items() if key.startswith("TAM_")
    }
    results = []
    for row in records:
        component = root("entry:" + row["pdb_id"])
        score = int(hashlib.sha256(f"{seed}:{component}".encode()).hexdigest()[:8], 16) / 2**32
        split = "test" if score < test else "val" if score < test + validation else "train"
        if component in tam_roots:
            split = "application_holdout"
        if clusters[row["binder_key"]] == clusters[row["target_key"]]:
            split = "excluded_homologous_pair"
        results.append({**row, "component": component, "split": split})
    return results
