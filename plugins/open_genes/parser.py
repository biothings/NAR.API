"""
Parser for the Open Genes datasource (https://open-genes.com/).

Open Genes has no directly-fetchable bulk file: its /download page is a
JS-rendered SPA (confirmed via `curl -sIL` returning `content-type: text/html`
for both the download page and a guessed `open_genes_sql_dump.zip` path — no
static file exists at the canonical domain). The confirmed-live REST API
(`/api/gene/search`) does, however, support a single request with a
generous `pageSize` that returns the entire gene collection (2,405 genes)
in one JSON response (~3.7 MB) — see manifest.json `dumper.data_url`. This
lets the plugin follow the standard manifest-first bulk-download strategy
even though the underlying source is API-backed.

Because the dumped filename can vary depending on how the Hub's HTTP
dumper names a query-string URL, this parser globs `data_folder` for any
file and picks the first one whose contents parse as JSON with an
`items` list, rather than assuming a specific filename.
"""

import glob
import json
import os

from biothings.utils.dataload import dict_sweep, unlist


def _find_source_file(data_folder):
    """Locate the dumped Open Genes API response inside data_folder."""
    candidates = sorted(
        p for p in glob.glob(os.path.join(data_folder, "*")) if os.path.isfile(p)
    )
    assert candidates, f"No files found in {data_folder}"

    for path in candidates:
        try:
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        if isinstance(data, dict) and isinstance(data.get("items"), list):
            return path, data
    raise AssertionError(
        f"No file in {data_folder} parsed as an Open Genes API JSON response "
        f"(expected a dict with an 'items' list). Files seen: {candidates}"
    )


def _clean_list(records, keep_keys=None):
    """Return a list of dicts, optionally restricted to keep_keys, dropping empties."""
    out = []
    for rec in records or []:
        if not isinstance(rec, dict):
            continue
        if keep_keys:
            rec = {k: rec[k] for k in keep_keys if k in rec}
        rec = dict_sweep(unlist(rec), [None, "", [], {}])
        if rec:
            out.append(rec)
    return out


def load_data(data_folder):
    """Parse Open Genes API data and yield one BioThings document per gene.

    _id is the NCBI (Entrez) Gene ID (`ncbiId`), matching MyGene.info-style
    gene identifiers and the primary_key verified in the site inspection.
    """
    _, data = _find_source_file(data_folder)
    items = data.get("items", [])

    n_yielded = 0
    n_skipped_no_id = 0

    for gene in items:
        ncbi_id = gene.get("ncbiId")
        if not ncbi_id:
            n_skipped_no_id += 1
            continue

        confidence = gene.get("confidenceLevel") or {}
        origin = gene.get("origin") or {}
        family_origin = gene.get("familyOrigin") or {}

        doc = {
            "_id": str(ncbi_id),
            "open_genes": {
                "symbol": gene.get("symbol"),
                "xrefs": {
                    "ncbigene": ncbi_id,
                    "ensembl": gene.get("ensembl"),
                    "uniprot": gene.get("uniprot"),
                },
                "confidence_level": confidence.get("name"),
                "origin": {
                    "phylum": origin.get("phylum"),
                    "age_million_years": origin.get("age"),
                },
                "family_origin": {
                    "phylum": family_origin.get("phylum"),
                    "age_million_years": family_origin.get("age"),
                },
                "disease_categories": _clean_list(
                    gene.get("diseaseCategories"),
                    keep_keys=["icdCode", "icdCategoryName"],
                ),
                "diseases": _clean_list(
                    gene.get("diseases"),
                    keep_keys=["icdCode", "name", "icdName"],
                ),
                "aging_mechanisms": [
                    m.get("name") for m in (gene.get("agingMechanisms") or []) if m.get("name")
                ],
                "functional_clusters": [
                    c.get("name") for c in (gene.get("functionalClusters") or []) if c.get("name")
                ],
                "protein_classes": [
                    c.get("name") for c in (gene.get("proteinClasses") or []) if c.get("name")
                ],
                "longevity_associations": [
                    c.get("name") for c in (gene.get("commentCause") or []) if c.get("name")
                ],
                "expression_change": gene.get("expressionChange"),
                "methylation_correlation": gene.get("methylationCorrelation"),
            },
        }

        doc = dict_sweep(unlist(doc), [None, "", [], {}])
        n_yielded += 1
        yield doc

    import logging

    logging.info(
        "open_genes parser: %d documents yielded, %d rows skipped (missing ncbiId)",
        n_yielded,
        n_skipped_no_id,
    )
