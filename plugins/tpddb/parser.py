"""
TPDdb parser — Targeted Protein Degraders Database

Parses 11 TSV files from tpddb.idrblab.net:
  - 6 compound main tables (PROTAC, MG, LYTAC, ATTEC, AUTAC, AUTOTAC)
  - 3 activity tables (PROTAC_activity, MG_activity, Lysosome-based_TPD_activity)
  - 1 disease table (TPD_Related_Diseases)
  - 1 PDB table (TPD_PDB)

Primary key: TPD ID (e.g. TPD-P58QBR) — unique across all modalities.
Target API: pending.api
"""

import os
import csv
from collections import defaultdict
from biothings.utils.dataload import dict_sweep, unlist


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Maps filename stem → TPD modality label
_MODALITY_MAP = {
    "PROTAC_main_table": "PROTAC",
    "MG_main_table": "Molecular Glue",
    "LYTAC_main_table": "LYTAC",
    "ATTEC_main_table": "ATTEC",
    "AUTAC_main_table": "AUTAC",
    "AUTOTAC_main_table": "AUTOTAC",
}

# Activity file stems
_ACTIVITY_FILES = {
    "PROTAC_activity",
    "MG_activity",
    "Lysosome-based_TPD_activity",
}


def _split_multi(value, sep=";"):
    """Split a delimited string, strip whitespace, drop empty/dot values."""
    if not value or value.strip() in (".", "N/A", ""):
        return []
    return [v.strip() for v in value.split(sep) if v.strip() and v.strip() != "."]


def _split_target_ids(value):
    """Split target IDs (UniProt), which use '/' as delimiter."""
    return _split_multi(value, sep="/")


def _clean_val(value):
    """Return None for placeholder values, else the stripped string."""
    if value is None:
        return None
    v = value.strip()
    return None if v in (".", "N/A", "") else v


def _parse_pubchem_synonyms(raw):
    """Parse semicolon-delimited synonym string, extract ChEMBL/CAS if present."""
    if not raw or raw.strip() in (".", ""):
        return None, []
    parts = [p.strip() for p in raw.split(";") if p.strip() and p.strip() != "."]
    chembl_ids = [p for p in parts if p.upper().startswith("CHEMBL")]
    cas_ids = [p for p in parts if p[:4].replace("-", "").isdigit() and "-" in p]
    return parts if parts else None, {"chembl": chembl_ids or None, "cas": cas_ids or None}


# ---------------------------------------------------------------------------
# Phase 1: Load compound main tables
# ---------------------------------------------------------------------------

def _load_main_tables(data_folder):
    """
    Load all 6 modality main tables. Returns dict: tpd_id -> compound doc.
    Handles the common schema (all files) and LYTAC-specific extra columns.
    """
    compounds = {}

    for stem, modality in _MODALITY_MAP.items():
        fpath = os.path.join(data_folder, stem + ".txt")
        if not os.path.exists(fpath):
            continue

        with open(fpath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                tpd_id = _clean_val(row.get("TPD ID"))
                if not tpd_id:
                    continue
                # Skip duplicate header rows embedded mid-file (data quality issue in PROTAC_main_table)
                if tpd_id == "TPD ID":
                    continue
                if tpd_id in compounds:
                    # Should not happen — TPD IDs are unique across modalities
                    continue

                smiles = _clean_val(row.get("SMILES"))
                raw_synonyms = row.get("PubChem synonyms", "")
                synonyms, xref_parts = _parse_pubchem_synonyms(raw_synonyms)

                target_symbols = _split_multi(row.get("Target Symbol", ""))
                target_ids = _split_target_ids(row.get("Target ID", ""))
                ligase = _clean_val(row.get("Ligase"))
                source = _clean_val(row.get("Source"))
                formula = _clean_val(row.get("Fomula"))  # typo in source data
                subtype = _clean_val(row.get("Subtype"))  # MG only

                # LYTAC-specific fields
                lytac_target = _clean_val(row.get("Lytac_target"))
                ltr = _clean_val(row.get("Lysosome-targeting receptors"))
                linker_type = _clean_val(row.get("Linker type"))
                lytac_linker = _clean_val(row.get("Lytac_Linker"))

                doc = {
                    "tpd_id": tpd_id,
                    "name": _clean_val(row.get("TPD NAME")),
                    "modality": modality,
                    "smiles": smiles,
                    "formula": formula,
                    "subtype": subtype,
                    "synonyms": synonyms,
                    "xrefs": {
                        "chembl": xref_parts.get("chembl") if xref_parts else None,
                        "cas": xref_parts.get("cas") if xref_parts else None,
                    },
                    "target": {
                        "symbols": target_symbols if target_symbols else None,
                        "uniprot_ids": target_ids if target_ids else None,
                    },
                    "ligase": ligase,
                    "source": source,
                }

                # Add LYTAC-specific fields only if present
                if lytac_target or ltr or linker_type or lytac_linker:
                    doc["lytac"] = {
                        "target": lytac_target,
                        "lysosome_targeting_receptor": ltr,
                        "linker_type": linker_type,
                        "linker": lytac_linker,
                    }

                compounds[tpd_id] = doc

    return compounds


# ---------------------------------------------------------------------------
# Phase 2: Load activity tables
# ---------------------------------------------------------------------------

def _load_activities(data_folder):
    """
    Load all 3 activity tables. Returns dict: tpd_id -> list of activity records.
    """
    activities = defaultdict(list)

    for stem in _ACTIVITY_FILES:
        fpath = os.path.join(data_folder, stem + ".txt")
        if not os.path.exists(fpath):
            continue

        with open(fpath, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                tpd_id = _clean_val(row.get("TPD ID"))
                if not tpd_id:
                    continue

                activity_type = _clean_val(row.get("Activity Type"))
                activity_val = _clean_val(row.get("Activity"))
                cell_line = _clean_val(row.get("Cell Line"))
                target_symbols = _split_multi(row.get("Target Symbols", ""))
                target_ids = _split_target_ids(row.get("Target IDs", ""))

                if not activity_type and not activity_val:
                    continue

                record = {
                    "activity_type": activity_type,
                    "activity": activity_val,
                    "cell_line": cell_line,
                }
                if target_symbols:
                    record["target_symbols"] = target_symbols
                if target_ids:
                    record["target_uniprot_ids"] = target_ids

                activities[tpd_id].append(record)

    return activities


# ---------------------------------------------------------------------------
# Phase 3: Load disease table
# ---------------------------------------------------------------------------

def _load_diseases(data_folder):
    """
    Load TPD_Related_Diseases.txt. Returns dict: tpd_id -> list of disease records.
    One row per disease per compound.
    """
    diseases = defaultdict(list)

    fpath = os.path.join(data_folder, "TPD_Related_Diseases.txt")
    if not os.path.exists(fpath):
        return diseases

    with open(fpath, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tpd_id = _clean_val(row.get("TPD_ID"))
            if not tpd_id:
                continue

            disease_name = _clean_val(row.get("Disease_Names"))
            icd11 = _clean_val(row.get("ICD-11"))
            source = _clean_val(row.get("Source"))
            orpha_ncit = _clean_val(row.get("OrphaID/Cellosaurus verified NCIt code"))

            if not disease_name:
                continue

            record = {
                "name": disease_name,
                "icd11": icd11,
                "source": source,
            }
            if orpha_ncit:
                record["orpha_or_ncit_id"] = orpha_ncit

            diseases[tpd_id].append(record)

    return diseases


# ---------------------------------------------------------------------------
# Phase 4: Load PDB table
# ---------------------------------------------------------------------------

def _load_pdb(data_folder):
    """
    Load TPD_PDB.txt. Returns dict: tpd_id -> list of PDB IDs.
    Note: Some entries have multiple PDB IDs separated by '/'.
    """
    pdb_map = defaultdict(list)

    fpath = os.path.join(data_folder, "TPD_PDB.txt")
    if not os.path.exists(fpath):
        return pdb_map

    with open(fpath, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tpd_id = _clean_val(row.get("TPD_ID"))
            pdb_raw = _clean_val(row.get("Complex_PDB_ID"))
            if not tpd_id or not pdb_raw:
                continue
            pdbs = _split_multi(pdb_raw, sep="/")
            pdb_map[tpd_id].extend(pdbs)

    return pdb_map


# ---------------------------------------------------------------------------
# Main parser
# ---------------------------------------------------------------------------

def load_data(data_folder):
    """
    Parse TPDdb data files and yield BioThings-compatible documents.

    _id = TPD ID (e.g. TPD-P58QBR) — unique proprietary identifier.
    All TPD-specific data nested under top-level 'tpddb' key.
    """
    # Load all sub-tables
    compounds = _load_main_tables(data_folder)
    activities = _load_activities(data_folder)
    diseases = _load_diseases(data_folder)
    pdb_map = _load_pdb(data_folder)

    seen_ids = set()

    for tpd_id, compound in compounds.items():
        if tpd_id in seen_ids:
            continue
        seen_ids.add(tpd_id)

        # Attach activity records
        activity_list = activities.get(tpd_id)
        # Attach disease records
        disease_list = diseases.get(tpd_id)
        # Attach PDB IDs
        pdb_ids = pdb_map.get(tpd_id)

        if activity_list:
            compound["activities"] = activity_list
        if disease_list:
            compound["diseases"] = disease_list
        if pdb_ids:
            compound["pdb_ids"] = pdb_ids

        doc = {
            "_id": tpd_id,
            "tpddb": compound,
        }

        doc = dict_sweep(unlist(doc), [None])
        yield doc
