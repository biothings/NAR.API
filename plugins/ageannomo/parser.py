"""
Parser for AgeAnnoMO — Multi-omics of Animal Ageing (NAR 2024 Database Issue).

Ingests four independent per-hallmark bulk files from the AgeAnnoMO GitHub
companion repo, each covering a distinct entity type:

  1. Differential expression.xlsx  (Genomic instability)      -> gene_expression
  2. Differential protein.xlsx     (Loss of proteostasis)     -> protein_expression
  3. Differential metabolite.xlsx  (Dysregulated metabolism)  -> metabolite
  4. Lifespan regulators.xlsx      (Lifespan regulators)      -> lifespan_regulator

There is no single global primary key across AgeAnnoMO — each hallmark
category uses a different composite key (dataset/species/tissue + gene
symbol, or dataset/tissue + UniProt/PubChem identifier, or gene symbol
alone). The parser therefore builds a composite, entity-type-prefixed `_id`
per row so all four entity types can coexist in one pending.api collection
without collision.
"""
import os
import re

import pandas as pd
from biothings.utils.dataload import dict_sweep, unlist

FILES = {
    # The biothings HTTPDumper saves files under the literal (percent-encoded)
    # basename of the source URL, so these must match the %20-encoded names
    # actually written to disk rather than the decoded "Differential expression.xlsx".
    "gene_expression": "Differential%20expression.xlsx",
    "protein_expression": "Differential%20protein.xlsx",
    "metabolite": "Differential%20metabolite.xlsx",
    "lifespan_regulator": "Lifespan%20regulators.xlsx",
}


def _slug(value):
    """Normalize a free-text field into a compact, _id-safe token."""
    if value is None:
        return "na"
    s = str(value).strip().lower()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^a-z0-9_.\-]", "", s)
    return s or "na"


def _to_float(value):
    try:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_int(value):
    try:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return None
        return int(value)
    except (TypeError, ValueError):
        return None


def _to_bool(value):
    if isinstance(value, bool):
        return value
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    s = str(value).strip().lower()
    if s in ("true", "1", "yes", "up"):
        return True
    if s in ("false", "0", "no", "down"):
        return False
    return None


def _read_excel(data_folder, filename, sheet_name=0):
    path = os.path.join(data_folder, filename)
    assert os.path.exists(path), f"Expected file not found: {path}"
    df = pd.read_excel(path, sheet_name=sheet_name, engine="openpyxl")
    # pandas' default "str" dtype (since pandas 3.0) silently keeps float NaN
    # in string columns when using `.where(pd.notnull(df), None)` directly;
    # casting to object first makes the None replacement actually stick.
    df = df.astype(object).where(pd.notnull(df), None)
    return df


def _parse_gene_expression(data_folder, seen_ids):
    """Genomic instability hallmark: age-related differential gene expression."""
    df = _read_excel(data_folder, FILES["gene_expression"], sheet_name="Sheet1")
    for row in df.to_dict(orient="records"):
        symbol = row.get("symbol")
        if not symbol:
            continue
        dataset_id = row.get("number")
        species = row.get("species")
        tissue = row.get("tissue")
        group = row.get("group")

        _id = "amoexpr:{}:{}:{}:{}:{}".format(
            _slug(dataset_id), _slug(species), _slug(symbol), _slug(tissue), _slug(group)
        )
        if _id in seen_ids:
            continue
        seen_ids.add(_id)

        doc = {
            "_id": _id,
            "ageannomo": {
                "entity_type": "gene_expression",
                "dataset_id": dataset_id,
                "species": species,
                "tissue": tissue,
                "comparison_group": group,
                "gene": {
                    "symbol": symbol,
                    "species_specific_gene": _to_bool(row.get("species_specific_gene")),
                },
                "statistics": {
                    "pvalue": _to_float(row.get("P.Value")),
                    "fdr": _to_float(row.get("FDR")),
                    "logfc": _to_float(row.get("logFC")),
                    "direction": row.get("Up/Down"),
                },
            },
        }
        yield doc


def _parse_protein_expression(data_folder, seen_ids):
    """Loss-of-proteostasis hallmark: age-related differential protein abundance."""
    df = _read_excel(data_folder, FILES["protein_expression"])
    for row in df.to_dict(orient="records"):
        uniprot = row.get("Uniprot entry")
        if not uniprot:
            continue
        dataset_id = row.get("Dataset ID")
        animal = row.get("Animal")
        tissue = row.get("Tissue")
        category = row.get("Category")

        _id = "amoprot:{}:{}:{}:{}".format(
            _slug(dataset_id), _slug(uniprot), _slug(tissue), _slug(category)
        )
        if _id in seen_ids:
            continue
        seen_ids.add(_id)

        pmid = _to_int(row.get("Pubmed"))
        doc = {
            "_id": _id,
            "ageannomo": {
                "entity_type": "protein_expression",
                "dataset_id": dataset_id,
                "species": animal,
                "tissue": tissue,
                "comparison_group": category,
                "protein": {
                    "name": row.get("Name"),
                    "xrefs": {"uniprot": uniprot},
                },
                "statistics": {
                    "pvalue": _to_float(row.get("P.Value")),
                    "fdr": _to_float(row.get("adj.P.Val")),
                    "logfc": _to_float(row.get("logFC")),
                    "direction": row.get("Up_or_down"),
                },
                "project_id": row.get("ProjectID"),
                "pubmed": pmid,
            },
        }
        yield doc


def _parse_metabolite(data_folder, seen_ids):
    """Dysregulated-metabolism hallmark: age-related differential metabolites."""
    df = _read_excel(data_folder, FILES["metabolite"])
    for row in df.to_dict(orient="records"):
        pubchem_cid = _to_int(row.get("Id"))
        if pubchem_cid is None:
            continue
        dataset_id = row.get("Dataset ID")
        animal = row.get("Animal")
        tissue = row.get("Tissue")

        _id = "amomet:{}:{}:{}".format(_slug(dataset_id), pubchem_cid, _slug(tissue))
        if _id in seen_ids:
            continue
        seen_ids.add(_id)

        pmid = _to_int(row.get("Pubmed"))
        doc = {
            "_id": _id,
            "ageannomo": {
                "entity_type": "metabolite",
                "dataset_id": dataset_id,
                "species": animal,
                "tissue": tissue,
                "age_group": row.get("Age group"),
                "metabolite": {
                    "name": row.get("Name"),
                    "molecular_formula": row.get("Molecular Formula"),
                    "molecular_weight": _to_float(row.get("Molecular weight")),
                    "xrefs": {"pubchem_cid": pubchem_cid},
                },
                "statistics": {
                    "plsda_vip": _to_float(row.get("plsda_vip")),
                },
                "is_age_related_metabolite": _to_bool(row.get("is_AgeRelatedMetabolite")),
                "project_id": row.get("Project ID"),
                "pubmed": pmid,
            },
        }
        yield doc


def _parse_lifespan_regulators(data_folder, seen_ids):
    """Lifespan-regulators hallmark: gene-lifespan correlation (mouse GenAge/MLS)."""
    df = _read_excel(data_folder, FILES["lifespan_regulator"], sheet_name="lifespans")
    for row in df.to_dict(orient="records"):
        gene = row.get("Gene")
        if not gene:
            continue

        _id = "amolife:{}".format(_slug(gene))
        if _id in seen_ids:
            continue
        seen_ids.add(_id)

        doc = {
            "_id": _id,
            "ageannomo": {
                "entity_type": "lifespan_regulator",
                "gene": {"symbol": gene},
                "statistics": {
                    "r_value": _to_float(row.get("R_value")),
                    "pvalue": _to_float(row.get("p_value")),
                },
                "category": row.get("Category"),
            },
        }
        yield doc


def load_data(data_folder):
    """Parse all four AgeAnnoMO hallmark files and yield BioThings-compatible documents.

    Each source row becomes one document. `_id` is a composite string prefixed by
    entity type (amoexpr/amoprot/amomet/amolife) to keep the four otherwise
    independent per-hallmark keyspaces collision-free within a single collection.
    Exact-duplicate composite keys within a file (a small number of true
    duplicate rows in the source spreadsheets) are dropped, keeping the first
    occurrence.
    """
    seen_ids = set()

    for doc in _parse_gene_expression(data_folder, seen_ids):
        yield dict_sweep(unlist(doc), [None])

    for doc in _parse_protein_expression(data_folder, seen_ids):
        yield dict_sweep(unlist(doc), [None])

    for doc in _parse_metabolite(data_folder, seen_ids):
        yield dict_sweep(unlist(doc), [None])

    for doc in _parse_lifespan_regulators(data_folder, seen_ids):
        yield dict_sweep(unlist(doc), [None])
