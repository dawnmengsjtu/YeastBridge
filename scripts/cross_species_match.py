#!/usr/bin/env python3
"""Leakage-resistant cross-species perturbation-to-yeast matching benchmark.

The benchmark deliberately makes no target-family assumption.  It projects signed
human/mouse and yeast signatures into an explicitly supplied conserved-module
space, compares three locked arms, and emits a GO/NO_GO routing decision.

Only NumPy and the Python standard library are required.  No network access or
large pretrained model is used.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


SCHEMA_VERSION = "1.0"
# `human_only` is the competition-facing arm name.  Its configured scope may
# include mouse rows as an additional higher-organism reference panel.
ARM_NAMES = ("human_only", "yeast_only", "joint")
REQUIRED_META_COLUMNS = {
    "signature_id",
    "species",
    "function_label",
    "perturbation_genes",
    "split_group",
    "split",
}
REQUIRED_PROJECTION_COLUMNS = {
    "species",
    "gene",
    "module_id",
    "evidence",
    "weight",
}


class BenchmarkError(ValueError):
    """Raised when inputs violate the predeclared benchmark protocol."""


@dataclass(frozen=True)
class SignatureMeta:
    signature_id: str
    species: str
    function_label: str
    perturbation_genes: tuple[str, ...]
    split_group: str
    split: str
    assay_type: str
    source: str

    @property
    def perturbation_key(self) -> str:
        return f"{self.species}|{'+'.join(self.perturbation_genes)}"


@dataclass(frozen=True)
class ProjectionEdge:
    module_id: str
    effective_weight: float
    evidence: str


def _normalise_name(value: str) -> str:
    return value.strip()


def _split_genes(value: str) -> tuple[str, ...]:
    genes = sorted(
        {gene.strip() for gene in re.split(r"[;,|+]", value) if gene.strip()}
    )
    if not genes:
        raise BenchmarkError("perturbation_genes must contain at least one gene")
    return tuple(genes)


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise BenchmarkError(f"missing input file: {path}")
    handle = (
        gzip.open(path, "rt", encoding="utf-8-sig", newline="")
        if path.suffix.lower() == ".gz"
        else path.open("r", encoding="utf-8-sig", newline="")
    )
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise BenchmarkError(f"TSV has no header: {path}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise BenchmarkError(f"TSV has no data rows: {path}")
    return list(reader.fieldnames), rows


def _parse_float(value: str, context: str) -> float:
    if value is None or not value.strip() or value.strip().upper() in {"NA", "NAN"}:
        return 0.0
    try:
        parsed = float(value)
    except ValueError as exc:
        raise BenchmarkError(f"non-numeric value in {context}: {value!r}") from exc
    if not math.isfinite(parsed):
        raise BenchmarkError(f"non-finite value in {context}: {value!r}")
    return parsed


def _parse_matrix_float(value: str | None, context: str) -> float:
    """Parse an observed matrix value while preserving unmeasured entries.

    A blank/NA cell means "not observed" and is represented as NaN internally.
    It is deliberately different from a measured zero effect.
    """

    if value is None or not value.strip() or value.strip().upper() in {
        "NA",
        "NAN",
        "NULL",
        ".",
    }:
        return float("nan")
    try:
        parsed = float(value)
    except ValueError as exc:
        raise BenchmarkError(f"non-numeric value in {context}: {value!r}") from exc
    if not math.isfinite(parsed):
        raise BenchmarkError(f"non-finite observed value in {context}: {value!r}")
    return parsed


def load_matrix(path: Path) -> tuple[list[str], dict[str, np.ndarray]]:
    # Wide transcriptomic matrices (GSE42528 is 1,484 x 6,170) must be read one
    # row at a time.  Routing them through ``_read_tsv`` would first retain
    # millions of Python strings in a list of dictionaries and can consume far
    # more memory than the numeric matrix itself.
    if not path.is_file():
        raise BenchmarkError(f"missing input file: {path}")
    handle = (
        gzip.open(path, "rt", encoding="utf-8-sig", newline="")
        if path.suffix.lower() == ".gz"
        else path.open("r", encoding="utf-8-sig", newline="")
    )
    with handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            columns = next(reader)
        except StopIteration as exc:
            raise BenchmarkError(f"matrix has no header: {path}") from exc
        if not columns:
            raise BenchmarkError(f"matrix has no header: {path}")
        if columns[0] != "signature_id":
            raise BenchmarkError(f"first matrix column must be signature_id: {path}")
        genes = columns[1:]
        if not genes or len(genes) != len(set(genes)):
            raise BenchmarkError(f"matrix genes must be non-empty and unique: {path}")
        matrix: dict[str, np.ndarray] = {}
        for row_number, values in enumerate(reader, start=2):
            if len(values) != len(columns):
                raise BenchmarkError(
                    f"matrix row width differs from header at {path}:{row_number}"
                )
            signature_id = _normalise_name(values[0])
            if not signature_id:
                raise BenchmarkError(f"blank signature_id at {path}:{row_number}")
            if signature_id in matrix:
                raise BenchmarkError(f"duplicate signature_id {signature_id!r} in {path}")
            matrix[signature_id] = np.asarray(
                [
                    _parse_matrix_float(
                        value, f"{path}:{row_number}:{gene}"
                    )
                    for gene, value in zip(genes, values[1:], strict=True)
                ],
                dtype=np.float64,
            )
    if not matrix:
        raise BenchmarkError(f"matrix has no data rows: {path}")
    return genes, matrix


def load_metadata(path: Path, source: str) -> list[SignatureMeta]:
    columns, rows = _read_tsv(path)
    missing = REQUIRED_META_COLUMNS.difference(columns)
    if missing:
        raise BenchmarkError(f"metadata {path} is missing columns: {sorted(missing)}")
    records: list[SignatureMeta] = []
    seen: set[str] = set()
    for row_number, row in enumerate(rows, start=2):
        signature_id = _normalise_name(row["signature_id"])
        if not signature_id or signature_id in seen:
            raise BenchmarkError(
                f"blank or duplicate signature_id {signature_id!r} at {path}:{row_number}"
            )
        seen.add(signature_id)
        species = _normalise_name(row["species"]).lower()
        label = _normalise_name(row["function_label"])
        split_group = _normalise_name(row["split_group"])
        split = _normalise_name(row["split"]).lower()
        if not species or not label or not split_group:
            raise BenchmarkError(
                f"species, function_label and split_group are required at {path}:{row_number}"
            )
        if split not in {"", "train", "test"}:
            raise BenchmarkError(
                f"split must be train, test or blank at {path}:{row_number}; got {split!r}"
            )
        records.append(
            SignatureMeta(
                signature_id=signature_id,
                species=species,
                function_label=label,
                perturbation_genes=_split_genes(row["perturbation_genes"]),
                split_group=split_group,
                split=split,
                assay_type=_normalise_name(row.get("assay_type", "unspecified"))
                or "unspecified",
                source=source,
            )
        )
    return records


def validate_matrix_metadata(
    matrix: Mapping[str, np.ndarray], metadata: Sequence[SignatureMeta], path_label: str
) -> None:
    matrix_ids = set(matrix)
    metadata_ids = {record.signature_id for record in metadata}
    if matrix_ids != metadata_ids:
        missing_meta = sorted(matrix_ids - metadata_ids)[:10]
        missing_matrix = sorted(metadata_ids - matrix_ids)[:10]
        raise BenchmarkError(
            f"{path_label} matrix/metadata IDs differ; missing metadata={missing_meta}, "
            f"missing matrix={missing_matrix}"
        )


def assign_locked_splits(
    records: Sequence[SignatureMeta], seed: int, test_fraction: float
) -> list[SignatureMeta]:
    """Use explicit splits or assign deterministic label-stratified group splits."""

    any_explicit = any(record.split for record in records)
    all_explicit = all(record.split for record in records)
    if any_explicit and not all_explicit:
        raise BenchmarkError(
            "split is partially specified; either lock every row or leave every row blank"
        )

    group_labels: dict[str, set[str]] = defaultdict(set)
    group_splits: dict[str, set[str]] = defaultdict(set)
    for record in records:
        group_labels[record.split_group].add(record.function_label)
        if record.split:
            group_splits[record.split_group].add(record.split)
    conflicting_splits = {
        group: splits for group, splits in group_splits.items() if len(splits) != 1
    }
    if conflicting_splits:
        raise BenchmarkError(
            f"split_group appears in both train and test: {conflicting_splits}"
        )

    if all_explicit:
        return list(records)
    if not (0.0 < test_fraction < 1.0):
        raise BenchmarkError("test_fraction must be strictly between 0 and 1")

    label_groups: dict[str, list[str]] = defaultdict(list)
    for group, labels in group_labels.items():
        if len(labels) != 1:
            raise BenchmarkError(
                "automatic label-stratified splitting cannot assign a split_group "
                "containing multiple function_label values. Supply an explicit "
                "component-grouped split instead; conflicts="
                f"{group}: {sorted(labels)}"
            )
        label_groups[next(iter(labels))].append(group)
    assignments: dict[str, str] = {}
    for label, groups in sorted(label_groups.items()):
        if len(groups) < 2:
            raise BenchmarkError(
                f"function_label {label!r} needs at least two split_groups for automatic split"
            )
        ranked = sorted(
            groups,
            key=lambda group: hashlib.sha256(
                f"{seed}|{label}|{group}".encode("utf-8")
            ).hexdigest(),
        )
        n_test = min(len(groups) - 1, max(1, round(len(groups) * test_fraction)))
        test_groups = set(ranked[:n_test])
        for group in groups:
            assignments[group] = "test" if group in test_groups else "train"

    return [
        SignatureMeta(**{**asdict(record), "split": assignments[record.split_group]})
        for record in records
    ]


def _vector_hash(vector: np.ndarray) -> str:
    values = vector.astype(np.float64)
    observed = np.isfinite(values)
    rounded = np.round(np.where(observed, values, 0.0), decimals=12)
    digest = hashlib.sha256()
    digest.update(observed.tobytes())
    digest.update(rounded.tobytes())
    return digest.hexdigest()


def validate_no_leakage(
    records: Sequence[SignatureMeta],
    raw_vectors: Mapping[str, np.ndarray],
) -> dict[str, Any]:
    """Reject group, exact perturbation-combination, and exact-vector leakage."""

    checks: dict[str, dict[str, set[str]]] = {
        "split_group": defaultdict(set),
        "perturbation_combination": defaultdict(set),
        "exact_signature": defaultdict(set),
    }
    for record in records:
        checks["split_group"][record.split_group].add(record.split)
        checks["perturbation_combination"][record.perturbation_key].add(record.split)
        vector_key = (
            f"{record.species}|{_vector_hash(raw_vectors[record.signature_id])}"
        )
        checks["exact_signature"][vector_key].add(record.split)

    leaks: dict[str, list[str]] = {}
    for name, values in checks.items():
        crossed = sorted(key for key, splits in values.items() if len(splits) > 1)
        if crossed:
            leaks[name] = crossed[:20]
    if leaks:
        raise BenchmarkError(
            "train/test leakage detected. Replicates, doses, cell lines, orthologous KOs and "
            f"gene combinations must share one split_group. Details: {leaks}"
        )

    return {
        "status": "PASS",
        "n_split_groups": len(checks["split_group"]),
        "n_species_specific_perturbation_combinations": len(
            checks["perturbation_combination"]
        ),
        "n_species_specific_exact_signatures": len(checks["exact_signature"]),
        "checks": [
            "split_group_not_crossed",
            "exact_species_gene_combination_not_crossed",
            "exact_species_signature_not_crossed",
        ],
        "user_obligation": (
            "split_group must unite every replicate/dose/cell-line and cross-species "
            "orthologous version of one perturbation or gene combination"
        ),
    }


def load_projection(
    path: Path,
    allowed_evidence: set[str],
    evidence_weights: Mapping[str, float],
) -> tuple[dict[tuple[str, str], list[ProjectionEdge]], list[str], dict[str, int]]:
    columns, rows = _read_tsv(path)
    missing = REQUIRED_PROJECTION_COLUMNS.difference(columns)
    if missing:
        raise BenchmarkError(f"projection {path} is missing columns: {sorted(missing)}")

    collapsed: dict[tuple[str, str, str], list[Any]] = {}
    evidence_counts: Counter[str] = Counter()
    for row_number, row in enumerate(rows, start=2):
        species = _normalise_name(row["species"]).lower()
        gene = _normalise_name(row["gene"])
        module_id = _normalise_name(row["module_id"])
        evidence = _normalise_name(row["evidence"]).lower()
        if not species or not gene or not module_id:
            raise BenchmarkError(f"blank projection identifier at {path}:{row_number}")
        if evidence not in allowed_evidence:
            raise BenchmarkError(
                f"unsupported evidence {evidence!r} at {path}:{row_number}; "
                f"allowed={sorted(allowed_evidence)}"
            )
        weight = _parse_float(row["weight"], f"{path}:{row_number}:weight")
        if weight <= 0:
            raise BenchmarkError(
                f"projection weight must be positive at {path}:{row_number}"
            )
        evidence_scale = float(evidence_weights.get(evidence, 1.0))
        if not math.isfinite(evidence_scale) or evidence_scale <= 0:
            raise BenchmarkError(f"invalid evidence weight for {evidence!r}")
        key = (species, gene, module_id)
        if key not in collapsed:
            collapsed[key] = [0.0, set()]
        collapsed[key][0] += weight * evidence_scale
        collapsed[key][1].add(evidence)
        evidence_counts[evidence] += 1

    projection: dict[tuple[str, str], list[ProjectionEdge]] = defaultdict(list)
    modules: set[str] = set()
    for (species, gene, module_id), (weight, evidences) in collapsed.items():
        modules.add(module_id)
        projection[(species, gene)].append(
            ProjectionEdge(
                module_id=module_id,
                effective_weight=float(weight),
                evidence="+".join(sorted(evidences)),
            )
        )
    for key in projection:
        projection[key] = sorted(projection[key], key=lambda edge: edge.module_id)
    return dict(projection), sorted(modules), dict(sorted(evidence_counts.items()))


def project_signatures(
    genes: Sequence[str],
    raw_vectors: Mapping[str, np.ndarray],
    records: Sequence[SignatureMeta],
    projection: Mapping[tuple[str, str], Sequence[ProjectionEdge]],
    modules: Sequence[str],
    min_nonzero_coverage: float,
    require_bidirectional: bool,
) -> tuple[dict[str, np.ndarray], list[dict[str, Any]]]:
    if not (0.0 <= min_nonzero_coverage <= 1.0):
        raise BenchmarkError("min_nonzero_coverage must be between 0 and 1")
    module_index = {module: index for index, module in enumerate(modules)}
    projected: dict[str, np.ndarray] = {}
    coverage_rows: list[dict[str, Any]] = []
    for record in records:
        raw = raw_vectors[record.signature_id]
        observed = np.isfinite(raw)
        if not np.any(observed):
            raise BenchmarkError(
                f"signature {record.signature_id!r} has no observed features"
            )
        observed_raw = raw[observed]
        if require_bidirectional and not (
            np.any(observed_raw > 0) and np.any(observed_raw < 0)
        ):
            raise BenchmarkError(
                f"signature {record.signature_id!r} is not signed bidirectionally; "
                "provide positive and negative effects or disable require_bidirectional"
            )
        numerator = np.zeros(len(modules), dtype=np.float64)
        denominator = np.zeros(len(modules), dtype=np.float64)
        total_abs = float(np.abs(observed_raw).sum())
        mapped_abs = 0.0
        mapped_observed_features = 0
        mapped_nonzero_features = 0
        nonzero_features = int(np.count_nonzero(observed_raw))
        for gene, value in zip(genes, raw, strict=True):
            if not math.isfinite(float(value)):
                continue
            edges = projection.get((record.species, gene), ())
            if not edges:
                continue
            mapped_observed_features += 1
            mapped_nonzero_features += int(value != 0)
            mapped_abs += abs(float(value))
            edge_total = sum(abs(edge.effective_weight) for edge in edges)
            for edge in edges:
                weight = edge.effective_weight / edge_total
                index = module_index[edge.module_id]
                numerator[index] += float(value) * weight
                denominator[index] += abs(weight)
        coverage = mapped_abs / total_abs if total_abs else 0.0
        if coverage < min_nonzero_coverage:
            raise BenchmarkError(
                f"signature {record.signature_id!r} projection coverage {coverage:.3f} "
                f"is below min_nonzero_coverage={min_nonzero_coverage:.3f}"
            )
        vector = np.divide(
            numerator,
            denominator,
            out=np.zeros_like(numerator),
            where=denominator > 0,
        )
        norm = float(np.linalg.norm(vector))
        if norm == 0:
            raise BenchmarkError(
                f"signature {record.signature_id!r} projects to a zero module vector"
            )
        vector /= norm
        projected[record.signature_id] = vector
        coverage_rows.append(
            {
                "signature_id": record.signature_id,
                "source": record.source,
                "species": record.species,
                "split": record.split,
                "observed_features": int(np.count_nonzero(observed)),
                "missing_features": int(len(raw) - np.count_nonzero(observed)),
                "observed_fraction": float(np.mean(observed)),
                "nonzero_features": nonzero_features,
                "mapped_observed_features": mapped_observed_features,
                "mapped_nonzero_features": mapped_nonzero_features,
                "absolute_effect_coverage": coverage,
                "nonzero_projected_modules": int(np.count_nonzero(vector)),
            }
        )
    return projected, coverage_rows


def fit_centroids(
    records: Sequence[SignatureMeta],
    projected: Mapping[str, np.ndarray],
    classes: Sequence[str],
) -> np.ndarray:
    centroids: list[np.ndarray] = []
    for label in classes:
        vectors = [
            projected[record.signature_id]
            for record in records
            if record.split == "train" and record.function_label == label
        ]
        if not vectors:
            raise BenchmarkError(f"no training signature for function_label={label!r}")
        centroid = np.mean(np.stack(vectors), axis=0)
        norm = float(np.linalg.norm(centroid))
        if norm == 0:
            raise BenchmarkError(f"zero centroid for function_label={label!r}")
        centroids.append(centroid / norm)
    return np.stack(centroids)


def _rank_order(scores: np.ndarray, classes: Sequence[str]) -> list[int]:
    return sorted(
        range(len(classes)), key=lambda index: (-scores[index], classes[index])
    )


def evaluate_scores(
    scores: np.ndarray,
    true_labels: Sequence[str],
    classes: Sequence[str],
    top_k: int,
) -> tuple[dict[str, float], list[str], np.ndarray, list[int]]:
    class_index = {label: index for index, label in enumerate(classes)}
    predictions: list[str] = []
    correct: list[float] = []
    reciprocal_ranks: list[float] = []
    top_k_correct: list[float] = []
    ranks: list[int] = []
    for row, truth in zip(scores, true_labels, strict=True):
        order = _rank_order(row, classes)
        prediction = classes[order[0]]
        predictions.append(prediction)
        if truth not in class_index:
            # In the preregistered open-set route, a test label absent from one
            # or both training arms remains in the all-test primary endpoint.
            # It is not added to the prediction universe and every arm receives
            # the same conservative zero credit.
            correct.append(0.0)
            reciprocal_ranks.append(0.0)
            top_k_correct.append(0.0)
            ranks.append(0)
        else:
            truth_rank = order.index(class_index[truth]) + 1
            correct.append(float(prediction == truth))
            reciprocal_ranks.append(1.0 / truth_rank)
            top_k_correct.append(float(truth_rank <= min(top_k, len(classes))))
            ranks.append(truth_rank)

    per_class_f1: list[float] = []
    # Unknown truth labels get F1=0 because the train-defined classifier cannot
    # emit them; retaining them prevents macro-F1 from hiding the open-set cost.
    evaluation_labels = sorted(set(classes).union(true_labels))
    for label in evaluation_labels:
        tp = sum(p == label and t == label for p, t in zip(predictions, true_labels))
        fp = sum(p == label and t != label for p, t in zip(predictions, true_labels))
        fn = sum(p != label and t == label for p, t in zip(predictions, true_labels))
        denominator = 2 * tp + fp + fn
        per_class_f1.append((2 * tp / denominator) if denominator else 0.0)
    metrics = {
        "accuracy": float(np.mean(correct)),
        "macro_f1": float(np.mean(per_class_f1)),
        "mrr": float(np.mean(reciprocal_ranks)),
        f"top_{min(top_k, len(classes))}_accuracy": float(np.mean(top_k_correct)),
    }
    return metrics, predictions, np.asarray(correct), ranks


def _aggregate_by_group(
    groups: Sequence[str], values: np.ndarray
) -> tuple[list[str], np.ndarray]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for group, value in zip(groups, values, strict=True):
        grouped[group].append(float(value))
    ordered = sorted(grouped)
    return ordered, np.asarray([np.mean(grouped[group]) for group in ordered])


def paired_inference(
    joint_correct: np.ndarray,
    comparator_correct: np.ndarray,
    groups: Sequence[str],
    seed: int,
    bootstrap_replicates: int,
    permutation_replicates: int,
) -> dict[str, Any]:
    joint_groups, joint_values = _aggregate_by_group(groups, joint_correct)
    comparator_groups, comparator_values = _aggregate_by_group(
        groups, comparator_correct
    )
    if joint_groups != comparator_groups:
        raise BenchmarkError("paired inference groups are not aligned")
    differences = joint_values - comparator_values
    observed = float(np.mean(differences))
    rng = np.random.default_rng(seed)

    if bootstrap_replicates < 100:
        raise BenchmarkError("bootstrap_replicates must be at least 100")
    n_groups = len(differences)
    indices = rng.integers(0, n_groups, size=(bootstrap_replicates, n_groups))
    bootstrap = differences[indices].mean(axis=1)
    ci_low, ci_high = np.quantile(bootstrap, [0.025, 0.975])

    if n_groups <= 16:
        signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=n_groups)))
        permuted = (signs * differences).mean(axis=1)
        permutation_p = float(np.mean(permuted >= observed - 1e-12))
        permutation_method = f"exact_sign_flip_{2**n_groups}_assignments"
    else:
        if permutation_replicates < 100:
            raise BenchmarkError("permutation_replicates must be at least 100")
        signs = rng.choice((-1.0, 1.0), size=(permutation_replicates, n_groups))
        permuted = (signs * differences).mean(axis=1)
        permutation_p = float(
            (1 + np.count_nonzero(permuted >= observed - 1e-12))
            / (permutation_replicates + 1)
        )
        permutation_method = f"monte_carlo_sign_flip_{permutation_replicates}"

    return {
        "n_paired_split_groups": n_groups,
        "delta_accuracy_joint_minus_comparator": observed,
        "bootstrap_ci95": [float(ci_low), float(ci_high)],
        "bootstrap_replicates": bootstrap_replicates,
        "permutation_p_one_sided": permutation_p,
        "permutation_method": permutation_method,
        "unit": "split_group",
    }


def make_decision(
    arm_metrics: Mapping[str, Mapping[str, float]],
    paired_tests: Mapping[str, Mapping[str, Any]],
    alpha: float,
    min_accuracy_gain: float,
    min_test_groups: int,
    tie_break_route: str,
) -> dict[str, Any]:
    failures: list[str] = []
    for comparator in ("human_only", "yeast_only"):
        test = paired_tests[f"joint_vs_{comparator}"]
        delta = float(test["delta_accuracy_joint_minus_comparator"])
        ci_low = float(test["bootstrap_ci95"][0])
        p_value = float(test["permutation_p_one_sided"])
        n_groups = int(test["n_paired_split_groups"])
        if n_groups < min_test_groups:
            failures.append(
                f"{comparator}: n_groups={n_groups} < min_test_groups={min_test_groups}"
            )
        if delta <= min_accuracy_gain:
            failures.append(
                f"{comparator}: delta={delta:.4f} is not > {min_accuracy_gain:.4f}"
            )
        if ci_low <= min_accuracy_gain:
            failures.append(
                f"{comparator}: CI_low={ci_low:.4f} is not > {min_accuracy_gain:.4f}"
            )
        if p_value > alpha:
            failures.append(f"{comparator}: paired p={p_value:.4g} > alpha={alpha:.4g}")

    if not failures:
        return {
            "status": "GO_JOINT",
            "selected_route": "joint",
            "joint_claim_allowed": True,
            "reason": (
                "joint beat both locked single-source arms with the predeclared paired gate"
            ),
            "failed_gate_checks": [],
        }

    higher_accuracy = float(arm_metrics["human_only"]["accuracy"])
    yeast_accuracy = float(arm_metrics["yeast_only"]["accuracy"])
    if higher_accuracy > yeast_accuracy:
        selected = "human_only"
    elif yeast_accuracy > higher_accuracy:
        selected = "yeast_only"
    else:
        selected = tie_break_route
    return {
        "status": "NO_GO_JOINT_ROUTE_SINGLE",
        "selected_route": selected,
        "joint_claim_allowed": False,
        "reason": (
            "joint did not clear every locked comparison; route to the strongest single "
            "arm and do not claim cross-species gain"
        ),
        "failed_gate_checks": failures,
    }


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"cannot serialise {type(value).__name__}")


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
        os.replace(temp_name, path)
    except BaseException:
        try:
            os.unlink(temp_name)
        except FileNotFoundError:
            pass
        raise


def _tsv_text(fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> str:
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fieldnames})
    return buffer.getvalue()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_path(config_dir: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = config_dir / path
    return path.resolve()


def _top_module_contributions(
    query: np.ndarray, candidate: np.ndarray, modules: Sequence[str], limit: int = 5
) -> str:
    contributions = query * candidate
    order = sorted(
        range(len(modules)),
        key=lambda index: (-abs(contributions[index]), modules[index]),
    )
    return ";".join(
        f"{modules[index]}:{contributions[index]:+.4f}"
        for index in order[:limit]
        if contributions[index] != 0
    )


def run_benchmark(
    config_path: Path,
    output_override: Path | None = None,
    bootstrap_override: int | None = None,
    permutation_override: int | None = None,
    seed_override: int | None = None,
) -> dict[str, Any]:
    with config_path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    config_dir = config_path.parent.resolve()
    seed = int(
        seed_override if seed_override is not None else config.get("seed", 20260901)
    )

    inputs = config["inputs"]
    input_paths = {
        key: _resolve_path(config_dir, inputs[key])
        for key in (
            "higher_matrix",
            "higher_metadata",
            "yeast_matrix",
            "yeast_metadata",
            "projection",
        )
    }
    configured_output = _resolve_path(config_dir, config["output_dir"])
    output_dir = (output_override or configured_output).resolve()

    higher_genes, higher_matrix = load_matrix(input_paths["higher_matrix"])
    yeast_genes, yeast_matrix = load_matrix(input_paths["yeast_matrix"])
    higher_metadata = load_metadata(input_paths["higher_metadata"], "higher")
    yeast_metadata = load_metadata(input_paths["yeast_metadata"], "yeast")
    validate_matrix_metadata(higher_matrix, higher_metadata, "higher")
    validate_matrix_metadata(yeast_matrix, yeast_metadata, "yeast")
    all_ids = [record.signature_id for record in higher_metadata + yeast_metadata]
    if len(all_ids) != len(set(all_ids)):
        raise BenchmarkError(
            "signature_id values must be globally unique across both sources"
        )

    protocol = config.get("protocol", {})
    higher_species = {
        species.lower()
        for species in protocol.get("higher_species", ["human", "mouse"])
    }
    yeast_species = {
        species.lower() for species in protocol.get("yeast_species", ["yeast"])
    }
    invalid_higher = sorted(
        {record.species for record in higher_metadata}.difference(higher_species)
    )
    invalid_yeast = sorted(
        {record.species for record in yeast_metadata}.difference(yeast_species)
    )
    if invalid_higher or invalid_yeast:
        raise BenchmarkError(
            f"species do not match protocol; higher={invalid_higher}, yeast={invalid_yeast}"
        )

    all_metadata = assign_locked_splits(
        higher_metadata + yeast_metadata,
        seed=seed,
        test_fraction=float(config.get("split", {}).get("test_fraction", 0.3)),
    )
    metadata_by_id = {record.signature_id: record for record in all_metadata}
    higher_metadata = [
        metadata_by_id[record.signature_id] for record in higher_metadata
    ]
    yeast_metadata = [metadata_by_id[record.signature_id] for record in yeast_metadata]

    raw_all = {**higher_matrix, **yeast_matrix}
    leakage_audit = validate_no_leakage(all_metadata, raw_all)

    projection_config = config.get("projection", {})
    allowed_evidence = {
        evidence.lower()
        for evidence in projection_config.get(
            "allowed_evidence", ["ortholog", "go", "conserved_module"]
        )
    }
    evidence_weights = {
        str(key).lower(): float(value)
        for key, value in projection_config.get("evidence_weights", {}).items()
    }
    projection, modules, evidence_counts = load_projection(
        input_paths["projection"], allowed_evidence, evidence_weights
    )
    if not modules:
        raise BenchmarkError("projection defines no conserved modules")

    min_coverage = float(projection_config.get("min_nonzero_coverage", 0.5))
    require_bidirectional = bool(
        projection_config.get("require_bidirectional_signature", True)
    )
    higher_projected, higher_coverage = project_signatures(
        higher_genes,
        higher_matrix,
        higher_metadata,
        projection,
        modules,
        min_coverage,
        require_bidirectional,
    )
    yeast_projected, yeast_coverage = project_signatures(
        yeast_genes,
        yeast_matrix,
        yeast_metadata,
        projection,
        modules,
        min_coverage,
        require_bidirectional,
    )
    projected = {**higher_projected, **yeast_projected}

    test_records = [record for record in higher_metadata if record.split == "test"]
    if not test_records:
        raise BenchmarkError("locked test split has no higher-organism signatures")
    higher_train_labels = {
        record.function_label for record in higher_metadata if record.split == "train"
    }
    yeast_train_labels = {
        record.function_label for record in yeast_metadata if record.split == "train"
    }
    shared_train_labels = higher_train_labels.intersection(yeast_train_labels)
    configured_classes = protocol.get("benchmark_function_labels")
    if configured_classes is None:
        # The class universe is derived only from locked training data.  In
        # particular, held-out labels cannot remove difficult distractor classes.
        classes = sorted(shared_train_labels)
    else:
        classes = sorted({_normalise_name(str(label)) for label in configured_classes})
        absent_from_training = sorted(set(classes).difference(shared_train_labels))
        if absent_from_training:
            raise BenchmarkError(
                "predeclared benchmark_function_labels lack training examples in one or "
                f"both sources: {absent_from_training}"
            )
    if len(classes) < 2:
        raise BenchmarkError(
            "training split needs at least two shared function_label classes"
        )
    unseen_test_labels = sorted(
        {record.function_label for record in test_records}.difference(classes)
    )
    unseen_policy = str(protocol.get("unseen_test_label_policy", "error"))
    if unseen_policy not in {"error", "score_incorrect_open_set"}:
        raise BenchmarkError(
            "protocol.unseen_test_label_policy must be error or "
            "score_incorrect_open_set"
        )
    if unseen_test_labels and unseen_policy == "error":
        raise BenchmarkError(
            "test contains labels outside the training-defined class universe: "
            f"{unseen_test_labels}"
        )

    higher_centroids = fit_centroids(higher_metadata, projected, classes)
    yeast_centroids = fit_centroids(yeast_metadata, projected, classes)
    query_matrix = np.stack([projected[record.signature_id] for record in test_records])
    higher_scores = query_matrix @ higher_centroids.T
    yeast_scores = query_matrix @ yeast_centroids.T
    model_config = config.get("model", {})
    joint_weight_higher = float(model_config.get("joint_weight_higher", 0.5))
    if not (0.0 <= joint_weight_higher <= 1.0):
        raise BenchmarkError("joint_weight_higher must be between 0 and 1")
    joint_scores = (
        joint_weight_higher * higher_scores + (1.0 - joint_weight_higher) * yeast_scores
    )
    score_by_arm = {
        "human_only": higher_scores,
        "yeast_only": yeast_scores,
        "joint": joint_scores,
    }

    top_k = int(model_config.get("top_k", 3))
    true_labels = [record.function_label for record in test_records]
    arm_metrics: dict[str, dict[str, float]] = {}
    predictions_by_arm: dict[str, list[str]] = {}
    correct_by_arm: dict[str, np.ndarray] = {}
    ranks_by_arm: dict[str, list[int]] = {}
    for arm in ARM_NAMES:
        metrics, predictions, correct, ranks = evaluate_scores(
            score_by_arm[arm], true_labels, classes, top_k
        )
        arm_metrics[arm] = metrics
        predictions_by_arm[arm] = predictions
        correct_by_arm[arm] = correct
        ranks_by_arm[arm] = ranks
    known_test_indices = [
        index for index, label in enumerate(true_labels) if label in classes
    ]
    known_only_metrics: dict[str, dict[str, float]] = {}
    if known_test_indices:
        known_truth = [true_labels[index] for index in known_test_indices]
        for arm in ARM_NAMES:
            metrics, _, _, _ = evaluate_scores(
                score_by_arm[arm][known_test_indices], known_truth, classes, top_k
            )
            known_only_metrics[arm] = metrics

    inference = config.get("inference", {})
    bootstrap_replicates = int(
        bootstrap_override
        if bootstrap_override is not None
        else inference.get("bootstrap_replicates", 5000)
    )
    permutation_replicates = int(
        permutation_override
        if permutation_override is not None
        else inference.get("permutation_replicates", 20000)
    )
    groups = [record.split_group for record in test_records]
    paired_tests = {
        f"joint_vs_{comparator}": paired_inference(
            correct_by_arm["joint"],
            correct_by_arm[comparator],
            groups,
            seed + offset,
            bootstrap_replicates,
            permutation_replicates,
        )
        for offset, comparator in enumerate(("human_only", "yeast_only"), start=1)
    }
    tie_break_route = str(inference.get("no_go_tie_break", "yeast_only"))
    if tie_break_route not in {"human_only", "yeast_only"}:
        raise BenchmarkError("no_go_tie_break must be human_only or yeast_only")
    alpha = float(inference.get("alpha", 0.05))
    min_accuracy_gain = float(inference.get("min_accuracy_gain", 0.0))
    min_test_groups = int(inference.get("min_test_groups", 8))
    decision = make_decision(
        arm_metrics,
        paired_tests,
        alpha=alpha,
        min_accuracy_gain=min_accuracy_gain,
        min_test_groups=min_test_groups,
        tie_break_route=tie_break_route,
    )

    lock_rows = [
        {
            "signature_id": record.signature_id,
            "source": record.source,
            "species": record.species,
            "function_label": record.function_label,
            "split_group": record.split_group,
            "perturbation_key": record.perturbation_key,
            "split": record.split,
        }
        for record in sorted(all_metadata, key=lambda item: item.signature_id)
    ]
    lock_fields = [
        "signature_id",
        "source",
        "species",
        "function_label",
        "split_group",
        "perturbation_key",
        "split",
    ]
    lock_text = _tsv_text(lock_fields, lock_rows)
    lock_path = output_dir / "split_lock.tsv"
    enforce_existing_lock = bool(
        config.get("split", {}).get("enforce_existing_lock", True)
    )
    fingerprints = {
        key: {"path": str(path), "sha256": _sha256_file(path)}
        for key, path in input_paths.items()
    }
    protocol_lock_payload = {
        "schema_version": SCHEMA_VERSION,
        "seed": seed,
        "script_sha256": _sha256_file(Path(__file__).resolve()),
        "config_sha256": _sha256_file(config_path),
        "input_sha256": {
            key: value["sha256"] for key, value in sorted(fingerprints.items())
        },
        "split_lock_sha256": hashlib.sha256(lock_text.encode("utf-8")).hexdigest(),
        "class_universe": classes,
        "protocol": protocol,
        "projection": projection_config,
        "model": model_config,
        "inference": {
            "bootstrap_replicates": bootstrap_replicates,
            "permutation_replicates": permutation_replicates,
            "alpha": alpha,
            "min_accuracy_gain": min_accuracy_gain,
            "min_test_groups": min_test_groups,
            "no_go_tie_break": tie_break_route,
        },
    }
    protocol_lock_text = (
        json.dumps(protocol_lock_payload, ensure_ascii=False, indent=2, sort_keys=True)
        + "\n"
    )
    protocol_lock_path = output_dir / "protocol_lock.json"
    if lock_path.exists() and enforce_existing_lock:
        existing = lock_path.read_text(encoding="utf-8")
        if existing != lock_text:
            raise BenchmarkError(
                f"existing split lock changed at {lock_path}; use a new versioned output_dir"
            )
    if protocol_lock_path.exists() and enforce_existing_lock:
        existing = protocol_lock_path.read_text(encoding="utf-8")
        if existing != protocol_lock_text:
            raise BenchmarkError(
                "existing analysis protocol or input fingerprint changed at "
                f"{protocol_lock_path}; use a new versioned output_dir"
            )
    _atomic_write_text(lock_path, lock_text)
    _atomic_write_text(protocol_lock_path, protocol_lock_text)

    prediction_rows: list[dict[str, Any]] = []
    for row_index, record in enumerate(test_records):
        row: dict[str, Any] = {
            "signature_id": record.signature_id,
            "species": record.species,
            "split_group": record.split_group,
            "perturbation_genes": "+".join(record.perturbation_genes),
            "true_function_label": record.function_label,
            "truth_in_training_class_universe": int(
                record.function_label in classes
            ),
        }
        for arm in ARM_NAMES:
            row[f"{arm}_prediction"] = predictions_by_arm[arm][row_index]
            row[f"{arm}_correct"] = int(correct_by_arm[arm][row_index])
            row[f"{arm}_truth_rank"] = ranks_by_arm[arm][row_index]
        prediction_rows.append(row)
    prediction_fields = list(prediction_rows[0])
    _atomic_write_text(
        output_dir / "predictions.tsv", _tsv_text(prediction_fields, prediction_rows)
    )

    match_top_k = int(model_config.get("match_top_k", 5))
    yeast_train = [record for record in yeast_metadata if record.split == "train"]
    yeast_train_matrix = np.stack(
        [projected[record.signature_id] for record in yeast_train]
    )
    match_rows: list[dict[str, Any]] = []
    for query_record, query_vector in zip(test_records, query_matrix, strict=True):
        similarities = yeast_train_matrix @ query_vector
        order = sorted(
            range(len(yeast_train)),
            key=lambda index: (-similarities[index], yeast_train[index].signature_id),
        )
        for rank, index in enumerate(order[:match_top_k], start=1):
            candidate = yeast_train[index]
            match_rows.append(
                {
                    "query_signature_id": query_record.signature_id,
                    "query_species": query_record.species,
                    "query_perturbation_genes": "+".join(
                        query_record.perturbation_genes
                    ),
                    "query_function_label_eval_only": query_record.function_label,
                    "rank": rank,
                    "yeast_signature_id": candidate.signature_id,
                    "yeast_assay_type": candidate.assay_type,
                    "yeast_perturbation_genes": "+".join(candidate.perturbation_genes),
                    "yeast_function_label_eval_only": candidate.function_label,
                    "cosine_similarity": float(similarities[index]),
                    "same_function_label_eval_only": int(
                        query_record.function_label == candidate.function_label
                    ),
                    "top_module_contributions": _top_module_contributions(
                        query_vector,
                        projected[candidate.signature_id],
                        modules,
                    ),
                }
            )
    match_fields = list(match_rows[0])
    _atomic_write_text(output_dir / "matches.tsv", _tsv_text(match_fields, match_rows))

    coverage_rows = higher_coverage + yeast_coverage
    coverage_fields = list(coverage_rows[0])
    _atomic_write_text(
        output_dir / "projection_coverage.tsv",
        _tsv_text(coverage_fields, coverage_rows),
    )

    metrics_payload: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "decision": decision,
        "arms": arm_metrics,
        "paired_tests": paired_tests,
        "secondary_known_only": {
            "role": "secondary_descriptive_only_not_used_for_GO_NO_GO",
            "n_test_signatures": len(known_test_indices),
            "arms": known_only_metrics,
        },
        "evaluation": {
            "query_source": "higher_test",
            "n_test_signatures": len(test_records),
            "n_test_split_groups": len(set(groups)),
            "function_labels": classes,
            "training_class_universe": classes,
            "unseen_test_labels": unseen_test_labels,
            "unseen_test_signatures": sum(
                record.function_label in unseen_test_labels for record in test_records
            ),
            "primary_endpoint_scope": "all_higher_test_signatures",
            "primary_metric": "accuracy",
            "paired_unit": "split_group",
        },
        "protocol": {
            "target_preselection": "none",
            "human_only_arm_scope": sorted(higher_species),
            "projection_space": "explicit_signed_conserved_modules",
            "allowed_projection_evidence": sorted(allowed_evidence),
            "evidence_weights": evidence_weights,
            "joint_weight_higher_fixed_before_test": joint_weight_higher,
            "test_labels_used_for_scoring_features": False,
            "class_universe_source": (
                "predeclared_config"
                if configured_classes is not None
                else "intersection_of_locked_higher_and_yeast_training_labels"
            ),
            "unseen_test_label_policy": unseen_policy,
            "unknown_truth_credit": {
                "accuracy": 0.0,
                "top_k": 0.0,
                "mrr": 0.0,
                "included_in_primary_all_test_endpoint": True,
            },
            "label_role": "benchmark ground truth and training prototypes only",
            "no_go_tie_break": tie_break_route,
        },
        "projection_audit": {
            "n_modules": len(modules),
            "modules": modules,
            "evidence_row_counts": evidence_counts,
            "minimum_absolute_effect_coverage": min_coverage,
            "observed_minimum_absolute_effect_coverage": min(
                row["absolute_effect_coverage"] for row in coverage_rows
            ),
        },
        "leakage_audit": leakage_audit,
        "split_lock": {
            "path": str(lock_path),
            "sha256": hashlib.sha256(lock_text.encode("utf-8")).hexdigest(),
            "enforced_on_rerun": enforce_existing_lock,
        },
        "protocol_lock": {
            "path": str(protocol_lock_path),
            "sha256": hashlib.sha256(protocol_lock_text.encode("utf-8")).hexdigest(),
            "enforced_on_rerun": enforce_existing_lock,
        },
        "random_seed": seed,
        "input_fingerprints": fingerprints,
    }
    _atomic_write_text(
        output_dir / "metrics.json",
        json.dumps(
            metrics_payload,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
            default=_json_default,
        )
        + "\n",
    )

    summary_lines = [
        "# Cross-species match benchmark",
        "",
        f"- Decision: **{decision['status']}**",
        f"- Selected route: **{decision['selected_route']}**",
        f"- Joint gain claim allowed: **{decision['joint_claim_allowed']}**",
        f"- Locked higher-organism test signatures: {len(test_records)}",
        f"- Independent split groups: {len(set(groups))}",
        "",
        "## Three-arm accuracy",
        "",
    ]
    summary_lines.extend(
        f"- {arm}: {arm_metrics[arm]['accuracy']:.4f}" for arm in ARM_NAMES
    )
    summary_lines.extend(
        [
            "",
            "## Interpretation",
            "",
            decision["reason"] + ".",
            "This is a computational matching benchmark, not evidence that a candidate "
            "drug works in yeast or that yeast predicts human efficacy. Wet-lab throughput "
            "and biological transfer must be measured separately.",
            "",
        ]
    )
    if decision["failed_gate_checks"]:
        summary_lines.append("## Failed joint gate checks")
        summary_lines.append("")
        summary_lines.extend(f"- {item}" for item in decision["failed_gate_checks"])
        summary_lines.append("")
    _atomic_write_text(output_dir / "RUN_SUMMARY.md", "\n".join(summary_lines))
    return metrics_payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", required=True, type=Path, help="JSON protocol config"
    )
    parser.add_argument(
        "--output-dir", type=Path, help="versioned output directory override"
    )
    parser.add_argument("--bootstrap-replicates", type=int)
    parser.add_argument("--permutation-replicates", type=int)
    parser.add_argument("--seed", type=int)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = run_benchmark(
            args.config.resolve(),
            output_override=args.output_dir,
            bootstrap_override=args.bootstrap_replicates,
            permutation_override=args.permutation_replicates,
            seed_override=args.seed,
        )
    except (BenchmarkError, KeyError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "decision": result["decision"],
                "arms": result["arms"],
            },
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
