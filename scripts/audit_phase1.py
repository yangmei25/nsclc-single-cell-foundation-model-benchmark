#!/usr/bin/env python3
"""Validate GSE205335 Phase I metadata and write auditable summary tables."""

from __future__ import annotations

import sys
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.data import load_cell_metadata, load_clinical_metadata, normalized_sample_key  # noqa: E402


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def child_text(element: ET.Element, name: str) -> str:
    child = next((item for item in element if local_name(item.tag) == name), None)
    return "" if child is None or child.text is None else child.text.strip()


def parse_geo_miniml() -> pd.DataFrame:
    archive = ROOT / "data/metadata/source/GSE205335_family.xml.tgz"
    with tarfile.open(archive, "r:gz") as tar:
        member = next(item for item in tar.getmembers() if item.name.endswith(".xml"))
        stream = tar.extractfile(member)
        if stream is None:
            raise ValueError("Could not extract GEO MINiML XML")
        root = ET.parse(stream).getroot()

    records: list[dict[str, str]] = []
    for sample in (item for item in root.iter() if local_name(item.tag) == "Sample"):
        record = {
            "geo_accession": child_text(sample, "Accession"),
            "geo_title": child_text(sample, "Title"),
            "geo_description": child_text(sample, "Description").splitlines()[0],
        }
        channel = next(
            (item for item in sample if local_name(item.tag) == "Channel"), None
        )
        if channel is not None:
            record["geo_source"] = child_text(channel, "Source")
            for item in channel:
                if local_name(item.tag) == "Characteristics":
                    record[f"geo_{item.attrib.get('tag', 'characteristic')}"] = (
                        item.text or ""
                    ).strip()
        records.append(record)
    frame = pd.DataFrame(records)
    if len(frame) != 33:
        raise ValueError(f"Expected 33 GEO samples, found {len(frame)}")
    return frame


def add_row(
    rows: list[dict[str, object]],
    category: str,
    metric: str,
    value: object,
    source: str,
    notes: str = "",
) -> None:
    rows.append(
        {
            "category": category,
            "metric": metric,
            "value": value,
            "source": source,
            "notes": notes,
        }
    )


def main() -> int:
    tables = ROOT / "results/tables"
    metadata_dir = ROOT / "data/metadata"
    tables.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    cells = load_cell_metadata(ROOT).reset_index(drop=True)
    clinical = load_clinical_metadata(ROOT)
    geo = parse_geo_miniml()

    clinical["sample_base_key"] = clinical["Sample"].map(normalized_sample_key)
    geo["sample_base_key"] = geo["geo_description"].map(normalized_sample_key)
    samples = clinical.merge(geo, on="sample_base_key", how="left", validate="one_to_one")
    if samples["geo_accession"].eq("").any() or samples["geo_accession"].isna().any():
        raise ValueError("Clinical samples did not map completely to GEO samples")

    cell_counts = cells.groupby("sample_key").size().rename("n_cells")
    samples = samples.merge(
        cell_counts, left_on="sample_key", right_index=True, how="left", validate="one_to_one"
    )
    if samples["n_cells"].isna().any() or int(samples["n_cells"].sum()) != len(cells):
        raise ValueError("Cell annotations did not map completely to clinical samples")
    samples["n_cells"] = samples["n_cells"].astype(int)

    samples["sampling_date_parsed"] = pd.to_datetime(samples["Sampling date"])
    samples["ici_start_date_parsed"] = pd.to_datetime(samples["Immunotherapy start date"])
    samples["timepoint_derived"] = "same_day_or_unknown"
    samples.loc[
        samples["sampling_date_parsed"] < samples["ici_start_date_parsed"],
        "timepoint_derived",
    ] = "pre"
    samples.loc[
        samples["sampling_date_parsed"] > samples["ici_start_date_parsed"],
        "timepoint_derived",
    ] = "post"

    cell_clinical = cells.merge(
        samples[
            [
                "sample_key",
                "Patient",
                "Class",
                "RECIST",
                "Response",
                "Tissue origin",
                "timepoint_derived",
            ]
        ],
        on="sample_key",
        how="left",
        validate="many_to_one",
    )
    if cell_clinical["Patient"].eq("").any() or cell_clinical["Patient"].isna().any():
        raise ValueError("Some cells lack patient metadata after joining")

    core_samples = samples[samples["Class"].eq("Core")].copy()
    if core_samples.groupby("Patient")["Response"].nunique().max() > 1:
        raise ValueError("A core patient has inconsistent response labels")
    if core_samples.groupby("Patient")["RECIST"].nunique().max() > 1:
        raise ValueError("A core patient has inconsistent RECIST labels")
    core_patients = core_samples.drop_duplicates("Patient")

    matrix_shape = pd.read_csv(
        ROOT / "data/processed/sparse_export/matrix_shape.csv"
    ).iloc[0]
    n_genes = int(matrix_shape["n_genes"])
    n_matrix_cells = int(matrix_shape["n_cells"])
    if n_matrix_cells != len(cells):
        raise ValueError(
            f"Matrix/annotation cell-count mismatch: {n_matrix_cells} versus {len(cells)}"
        )

    rows: list[dict[str, object]] = []
    add_row(rows, "dataset", "geo_accession", "GSE205335", "GEO series record")
    add_row(rows, "dataset", "publication_doi", "10.7554/eLife.98366", "eLife/PubMed")
    add_row(rows, "dataset", "n_cells", len(cells), "GEO cell-identity file")
    add_row(rows, "dataset", "n_genes", n_genes, "Processed UMI matrix audit")
    add_row(rows, "dataset", "n_samples", samples["Sample"].nunique(), "GEO + eLife supplement")
    add_row(rows, "dataset", "n_patients", samples["Patient"].nunique(), "GEO + eLife supplement")
    add_row(rows, "clinical_core", "n_samples", len(core_samples), "eLife Supplementary file 1")
    add_row(rows, "clinical_core", "n_patients", core_patients["Patient"].nunique(), "eLife Supplementary file 1")
    for label, count in core_patients["Response"].value_counts().items():
        add_row(rows, "clinical_core", f"patients_{label.lower().replace('-', '_')}", count, "eLife Supplementary file 1")
    for label, count in core_patients["RECIST"].value_counts().items():
        add_row(rows, "clinical_core", f"patients_recist_{label.lower()}", count, "eLife Supplementary file 1")
    add_row(
        rows,
        "clinical_core",
        "response_definition",
        "Responder=PR; Non-responder=SD or PD",
        "eLife Methods: Clinical outcomes",
        "No complete responses were observed.",
    )
    add_row(
        rows,
        "metadata",
        "explicit_malignant_status",
        "not available",
        "GEO cell-identity file",
        "Epithelial/tumor annotations exist, but no authoritative per-cell malignant flag is supplied.",
    )
    add_row(
        rows,
        "metadata",
        "treatment_timepoint",
        "derived from sampling and ICI start dates",
        "eLife Supplementary file 1",
        "All core samples are pretreatment; derived status must remain documented rather than presented as a source column.",
    )
    pd.DataFrame(rows).to_csv(tables / "dataset_summary.csv", index=False)

    samples.to_csv(metadata_dir / "sample_metadata.csv", index=False)
    (
        cell_clinical.groupby(["Patient", "Class", "Response"], dropna=False)
        .size()
        .rename("n_cells")
        .reset_index()
        .to_csv(tables / "cells_per_patient.csv", index=False)
    )

    core_cells = cell_clinical[cell_clinical["Class"].eq("Core")]
    lineage_rows = []
    for lineage, group in cell_clinical.groupby("lineage.total", dropna=False):
        core_group = core_cells[core_cells["lineage.total"].eq(lineage)]
        lineage_rows.append(
            {
                "cell_type": lineage,
                "n_cells_all": len(group),
                "n_patients_all": group["Patient"].nunique(),
                "n_cells_core": len(core_group),
                "n_patients_core": core_group["Patient"].nunique(),
                "n_responder_patients": core_group.loc[
                    core_group["Response"].eq("Responder"), "Patient"
                ].nunique(),
                "n_nonresponder_patients": core_group.loc[
                    core_group["Response"].eq("Non-responder"), "Patient"
                ].nunique(),
            }
        )
    pd.DataFrame(lineage_rows).sort_values("n_cells_all", ascending=False).to_csv(
        tables / "cell_type_patient_coverage.csv", index=False
    )

    missingness = pd.DataFrame(
        {
            "field": samples.columns,
            "n_missing_or_blank": [
                int(samples[column].isna().sum() + samples[column].astype(str).eq("").sum())
                for column in samples.columns
            ],
        }
    )
    missingness.to_csv(tables / "metadata_missingness.csv", index=False)

    report = f"""# GSE205335 Phase I Feasibility Report

## Verified scope

- Full atlas: {len(cells):,} cells from {samples['Sample'].nunique()} samples and {samples['Patient'].nunique()} patients.
- Response-analysis core: {len(core_samples)} pretreatment samples from {core_patients['Patient'].nunique()} patients.
- Patient-level response groups: {int((core_patients['Response'] == 'Responder').sum())} responders (PR) and {int((core_patients['Response'] == 'Non-responder').sum())} nonresponders (SD/PD).
- Processed matrix: {n_genes:,} genes x {n_matrix_cells:,} cells; raw UMI counts in a sparse R matrix.

## Feasibility conclusion

The dataset is feasible for conventional single-cell analysis and exploratory,
patient-aware response association. It is not large enough for a credible
high-capacity clinical prediction model: the effective response sample size is
11 patients. Foundation-model representations may later be compared with simple
baselines, but only with grouped validation and strong uncertainty reporting.

## Material ambiguities and limitations

- GEO has no linked citation; the associated publication was verified through
  Europe PMC/PubMed and the paper's data-availability statement.
- Response analysis uses only the paper-defined core cohort, not all 26 patients.
- Response is PR versus SD/PD under RECIST 1.1; this is not the same endpoint as
  pathological response used by some NSCLC studies.
- Multiple samples come from some patients, so sample rows are not independent.
- Tissue sites are heterogeneous, especially metastatic lymph node, lung,
  pleural effusion, and liver.
- Treatments include different PD-1/PD-L1 regimens and lines of therapy.
- `malignant_status` is not supplied as an explicit authoritative per-cell field.
- Cell annotations are published labels and should be reused, with marker-based
  validation rather than unnecessary re-annotation.
- The paper's reported response classifier was based on in-sample performance;
  it must not be treated as externally validated performance.
"""
    (ROOT / "results/tables/phase1_feasibility_report.md").write_text(report)

    print(f"Cells: {len(cells):,}")
    print(f"Samples: {samples['Sample'].nunique()}")
    print(f"Patients: {samples['Patient'].nunique()}")
    print(f"Core samples/patients: {len(core_samples)}/{core_patients['Patient'].nunique()}")
    print(core_patients["Response"].value_counts().to_string())
    print("Phase I metadata audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
