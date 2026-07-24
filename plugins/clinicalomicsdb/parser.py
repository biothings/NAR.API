"""
parser.py — ClinicalOmicsDB (trials.linkedomics.org) plugin for pending.api

Data source: 64 per-treatment-arm bulk JSON downloads from the ClinicalOmicsDB
REST API (table/study/gene/{study} endpoint), one file per clinical-trial
treatment arm. Each downloaded file is a JSON array of per-gene significance
statistics (p-value, AUROC, FDR) computed for that arm's responder vs.
non-responder transcriptomic comparison.

Study-level clinical metadata (disease, treatment, sample sizes, NCT ID,
PubMed ID, GEO series, raw-data download link) is NOT re-fetched at parse
time. It is bundled as a small static reference file, ``study_metadata.json``,
shipped alongside this parser. That file was built once, ahead of time, from
the same API's info/{study} endpoint for the identical 64 study identifiers
used in manifest.json's data_url list (see design_rationale.md for the exact
provenance and why this avoids a network dependency + a filename collision
at parse time).

Each output document represents one clinical-trial treatment arm, with the
full per-gene significance-statistics table nested as a list
(``gene_stats``), mirroring the disgenet-style "associatedWith list" pattern
for merged multi-row records (see references/production-plugin-examples.md
Pattern 2).
"""

import glob
import json
import logging
import os

from biothings.utils.dataload import dict_sweep, unlist

logger = logging.getLogger(__name__)


def _load_study_metadata():
    """Load the bundled static per-study metadata reference file."""
    meta_path = os.path.join(os.path.dirname(__file__), "study_metadata.json")
    with open(meta_path, "r", encoding="utf-8") as f:
        return json.load(f)


def _to_bool(value):
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in ("true", "1", "yes"):
        return True
    if value in ("false", "0", "no"):
        return False
    return None


def _to_int(value):
    try:
        if value is None or value == "":
            return None
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _to_float(value):
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _parse_gene_stats(rows):
    """Convert the raw table/study/gene/{study} rows to the output schema."""
    stats = []
    for row in rows:
        gene = row.get("analyte")
        if not gene:
            continue
        stats.append(
            {
                "gene": gene,
                "p_value": _to_float(row.get("p")),
                "auroc": _to_float(row.get("auc")),
                "fdr": _to_float(row.get("fdr")),
                "sorted_p": _to_float(row.get("sorted_p")),
                "sorted_fdr": _to_float(row.get("sorted_fdr")),
            }
        )
    return stats


def load_data(data_folder):
    """Parse ClinicalOmicsDB per-study gene-significance JSON files and yield
    one BioThings-compatible document per clinical-trial treatment arm.
    """
    study_metadata = _load_study_metadata()

    # The Hub's HTTPDumper names downloaded files after the URL basename, e.g.
    # ".../table/study/gene/GSE14764.csv" -> "GSE14764.csv" (JSON content
    # despite the .csv extension inherited from the source study identifier).
    infiles = sorted(glob.glob(os.path.join(data_folder, "*.csv")))
    assert infiles, f"No downloaded study files found in {data_folder}"

    seen_ids = set()
    for infile in infiles:
        study_key = os.path.basename(infile)
        meta = study_metadata.get(study_key)
        if meta is None:
            logger.warning("No bundled metadata for study file %s — skipping", study_key)
            continue

        study_id = os.path.splitext(study_key)[0]
        if study_id in seen_ids:
            logger.warning("Duplicate study_id %s — skipping", study_id)
            continue
        seen_ids.add(study_id)

        with open(infile, "r", encoding="utf-8") as f:
            try:
                rows = json.load(f)
            except json.JSONDecodeError:
                logger.error("Could not parse JSON from %s — skipping", infile)
                continue

        gene_stats = _parse_gene_stats(rows)
        if not gene_stats:
            logger.warning("No gene stats parsed for %s — skipping", study_id)
            continue

        treatment_raw = meta.get("treatment") or ""
        treatment = [t.strip() for t in treatment_raw.split(",") if t.strip()]

        doc = {
            "_id": study_id,
            "clinicalomicsdb": {
                "study_id": study_id,
                "geo_series": meta.get("series"),
                "disease": meta.get("disease"),
                "subtype": meta.get("subtype"),
                "adjuvant": _to_bool(meta.get("adjuvant")),
                "treatment": treatment,
                "response_eval": meta.get("response_eval"),
                "sample_size": _to_int(meta.get("sample_size")),
                "responder_size": _to_int(meta.get("responder_size")),
                "non_responder_size": _to_int(meta.get("non_responder_size")),
                "download_url": meta.get("download_url"),
                "xrefs": {
                    "geo": meta.get("series"),
                    "clinicaltrials_gov": meta.get("clinical_trial_id") or None,
                    "pubmed": meta.get("pub_med_id") or None,
                },
                "gene_stats": gene_stats,
            },
        }
        yield dict_sweep(unlist(doc), [None, "", []])
