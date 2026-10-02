"""Strain-level HIP/HOP responses from processed barcode expression.

The E-MTAB-2391 design pools thousands of deletion strains into one culture,
so every array measures every barcode probe set at once:

- vehicle control arrays are the SDRF rows without any compound identity;
- the raw response of a probe set is its normalised expression minus the
  median of vehicle arrays from the same processing batch (the leading date
  token of the array filename); arrays whose batch has no vehicle fall back
  to the global vehicle median;
- probe sets are aggregated to strain ORFs by averaging their sets;
- outputs carry the exact method note so no RMA-grade claim can be smuggled
  in: this is log2 + quantile normalisation + median polish + vehicle
  contrast.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from typing import Any

import numpy as np

from .workspace import exclusive_json_write, project_root, resolve_in_root, sha256_file


class HiphopResponseError(ValueError):
    """Raised when HIP/HOP response computation cannot stay auditable."""


_METHOD = (
    "barcode response = median-polished log2 quantile-normalised expression "
    "(no background model) minus batch-matched vehicle median"
)


def _batch_of(array_name: str) -> str:
    return str(array_name).split("_")[0] if "_" in str(array_name) else str(array_name)


def build_hiphop_response(
    root: str,
    *,
    expr_path: str,
    sdrf_path: str,
    output_npz: str,
    output_json: str,
) -> dict[str, Any]:
    base = project_root(root)
    expr_file = resolve_in_root(base, expr_path, must_exist=True)
    sdrf_file = resolve_in_root(base, sdrf_path, must_exist=True)

    with gzip.open(expr_file, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        set_ids: list[str] = []
        values: list[list[float]] = []
        for row in reader:
            if not row:
                continue
            set_ids.append(row[0])
            values.append([float(item) for item in row[1:]])
    # Some writers emit a leading empty header cell for the row-name column,
    # some do not; align on the data-row width instead of trusting either.
    data_width = len(values[0]) if values else 0
    if len(header) == data_width + 1:
        header_aligned = [""] + list(header)
    elif len(header) == data_width:
        header_aligned = [""]
        header_aligned.extend(str(item) for item in header)
    else:
        raise HiphopResponseError(
            f"header width {len(header)} incompatible with data width {data_width}"
        )
    arrays = [name.strip() for name in header_aligned[1:]]
    expr = np.asarray(values, dtype=np.float64)
    if expr.shape[0] != len(set_ids) or expr.shape[1] != len(arrays):
        raise HiphopResponseError("expression matrix shape does not match its header")

    with sdrf_file.open(encoding="utf-8-sig", newline="") as handle:
        sdrf_rows = list(csv.DictReader(handle, delimiter="\t"))
    metadata = {str(row["Array Data File"]).strip(): row for row in sdrf_rows}
    missing_meta = [name for name in arrays if name not in metadata]
    if missing_meta:
        raise HiphopResponseError(
            f"{len(missing_meta)} processed arrays lack SDRF rows, e.g. "
            f"{missing_meta[:3]}"
        )

    vehicle_names = [
        name
        for name in arrays
        if not str(metadata[name].get("Comment[inchikey]", "")).strip()
    ]
    treated_names = [name for name in arrays if name not in set(vehicle_names)]
    if not vehicle_names or not treated_names:
        raise HiphopResponseError("SDRF lacks vehicle or treated arrays")

    vehicle_index = np.array([arrays.index(name) for name in vehicle_names])
    global_vehicle = np.median(expr[:, vehicle_index], axis=1)

    vehicle_set = set(vehicle_names)
    batch_vehicle = {
        batch: np.array(
            [
                arrays.index(name)
                for name in vehicle_names
                if _batch_of(name) == batch
            ]
        )
        for batch in {_batch_of(n) for n in vehicle_names}
    }
    fallback_batches: list[str] = []

    response = np.empty_like(expr)
    for column, name in enumerate(arrays):
        batch = _batch_of(name)
        index = batch_vehicle.get(batch)
        if index is None or index.size == 0 or name in vehicle_set:
            reference = global_vehicle
            if batch not in batch_vehicle and name not in vehicle_set:
                fallback_batches.append(batch)
        else:
            reference = np.median(expr[:, index], axis=1)
        response[:, column] = expr[:, column] - reference

    # aggregate probe sets to strain ORFs (mean over the sets of one strain)
    strains: dict[str, list[int]] = {}
    for row, set_id in enumerate(set_ids):
        orf = str(set_id).split(".")[0].strip()
        strains.setdefault(orf, []).append(row)
    ordered_strains = sorted(strains)
    strain_matrix = np.vstack(
        [response[strains[orf], :].mean(axis=0) for orf in ordered_strains]
    )

    inchikey_of = []
    dose_of = []
    unit_of = []
    for name in arrays:
        row = metadata[name]
        key = str(row.get("Comment[inchikey]", "")).strip().upper()
        inchikey_of.append(key)
        dose_of.append(str(row.get("Factor Value[dose]", "")).strip())
        unit_of.append(str(row.get("Unit[concentration unit]", "")).strip())

    npz_path = resolve_in_root(base, output_npz)
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    with open(npz_path, "wb") as handle:
        np.savez_compressed(
            handle,
            strain_orfs=np.array(ordered_strains),
            array_names=np.array(arrays),
            compound_inchikeys=np.array(inchikey_of),
            doses=np.array(dose_of),
            dose_units=np.array(unit_of),
            is_vehicle=np.array([n in set(vehicle_names) for n in arrays]),
            z_score=strain_matrix.astype(np.float32),
            method=np.array([_METHOD]),
        )

    treated_keys = sorted({k for k in inchikey_of if k})
    document = {
        "schema_version": "yeastbridge.discovery.hiphop-response.v1",
        "method": _METHOD,
        "arrays_processed": len(arrays),
        "vehicle_arrays": len(vehicle_names),
        "treated_arrays": len(treated_names),
        "probe_sets": len(set_ids),
        "strains": len(ordered_strains),
        "compounds_with_identity": len(treated_keys),
        "batches_without_matched_vehicle": sorted(set(fallback_batches)),
        "vehicle_definition": "SDRF rows with an empty Comment[inchikey]",
        "sources": {
            "expression_path": expr_path,
            "expression_sha256": sha256_file(expr_file),
            "sdrf_path": sdrf_path,
            "sdrf_sha256": sha256_file(sdrf_file),
            "output_npz_sha256": sha256_file(npz_path),
        },
        "claim_boundary": (
            "Barcode fitness contrasts rank strain sensitivity within this "
            "public dataset; they are not target identities and do not by "
            "themselves validate any compound."
        ),
    }
    exclusive_json_write(resolve_in_root(base, output_json), document)
    return document


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m yeastbridge_vs.discovery_release.hiphop_response"
    )
    parser.add_argument("--root", default=".")
    parser.add_argument("--expr", required=True)
    parser.add_argument("--sdrf", required=True)
    parser.add_argument("--output-npz", required=True)
    parser.add_argument("--output-json", required=True)
    arguments = parser.parse_args(argv)
    build_hiphop_response(
        arguments.root,
        expr_path=arguments.expr,
        sdrf_path=arguments.sdrf,
        output_npz=arguments.output_npz,
        output_json=arguments.output_json,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
