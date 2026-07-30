"""
PersADE parser — personalized adverse drug event association database.

This plugin ingests the available mtADENet dataset (downloaded 2023-06-27),
which is the predecessor to the full PersADE database described in the 2026
NAR paper. The full PersADE dataset (4M personalized associations with
demographic stratification) is not yet available for bulk download.

Primary entity: drug compound.
_id strategy: InChIKey (MyChem.info standard), from Global_Drug file.
Structure: one document per unique InChIKey, with nested drug metadata
and a list of ADE associations grouped from Global_associations file.

Files consumed:
  Global_Drug_20230627.xlsx         — drug metadata: InChIKey, MeSH ID, CASRN,
                                       SMILES, xrefs (PubChem, ChEMBL, etc.),
                                       ATC code, drug type (8,802 drugs)
  Global_ADE_20230627.xlsx          — ADE metadata: MeSH ID, Tree Number, Name
                                       (4,380 unique ADEs)
  Global_associations_20230627.xlsx — drug-ADE pairs: InChIKey + MeSH ID
                                       (461,848 associations)

Data gap warning: The 2026 NAR PersADE paper describes 4,061,772 personalized
associations with demographic stratification (age, sex, route, dose). These
are NOT in the available download files (which are from the 2023 mtADENet
version). This plugin ingests the available data only.
"""

import os
import glob
from collections import defaultdict
from biothings.utils.dataload import dict_sweep, unlist

try:
    import openpyxl
except ImportError:
    raise ImportError("openpyxl is required: pip install openpyxl")


# ─── helpers ────────────────────────────────────────────────────────────────

def _safe_float(val):
    if val is None:
        return None
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _or_none(val):
    if val is None:
        return None
    s = str(val).strip()
    return s if s and s not in (".", "None", "N/A", "-", "nan", "") else None


def _find_file(data_folder, pattern):
    """Find a file by glob pattern in data_folder."""
    matches = glob.glob(os.path.join(data_folder, pattern))
    if matches:
        return matches[0]
    return None


def _read_xlsx(path):
    """Yield rows as dicts from an xlsx file (first sheet)."""
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows = iter(ws.iter_rows(values_only=True))
    # First row is header
    try:
        header = [str(h).strip() if h is not None else f"col_{i}"
                  for i, h in enumerate(next(rows))]
    except StopIteration:
        wb.close()
        return
    for row in rows:
        yield dict(zip(header, row))
    wb.close()


def _parse_xref_links(link_str):
    """
    Parse pipe-delimited xref string like:
    'PubChem$46943432|BindingDB$50365463|ChEBI$95082|ChEMBL$CHEMBL1232461|...'
    Returns a dict: {source_lower: id}.
    """
    if not link_str or str(link_str).strip() in ("", "None", "nan"):
        return {}
    xrefs = {}
    for part in str(link_str).split("|"):
        part = part.strip()
        if "$" in part:
            src, val = part.split("$", 1)
            src_key = src.lower().replace("-", "_")
            val = val.strip()
            if val:
                xrefs[src_key] = val
    return xrefs


# ─── index builders ──────────────────────────────────────────────────────────

def _build_ade_index(data_folder):
    """Build MeSH ID → ADE metadata dict from Global_ADE file."""
    ade_index = {}
    fpath = _find_file(data_folder, "Global_ADE_*.xlsx")
    if not fpath:
        return ade_index
    for row in _read_xlsx(fpath):
        mesh_id = _or_none(row.get("MeSH ID") or row.get("MeSH_ID") or row.get("MeSH"))
        if not mesh_id:
            continue
        ade_index[mesh_id] = dict_sweep({
            "mesh_id": mesh_id,
            "tree_number": _or_none(row.get("Tree Number") or row.get("Tree_Number")),
            "name": _or_none(row.get("Name") or row.get("name") or row.get("ADE Name")),
        }, [None])
    return ade_index


def _build_association_index(data_folder, ade_index):
    """
    Build InChIKey → list of ADE records from Global_associations file.
    Each ADE record has mesh_id, tree_number, name.
    """
    assoc_index = defaultdict(list)
    fpath = _find_file(data_folder, "Global_associations_*.xlsx")
    if not fpath:
        return assoc_index

    for row in _read_xlsx(fpath):
        inchikey = _or_none(
            row.get("InChI Key") or row.get("InChIKey") or row.get("inchikey")
        )
        mesh_id = _or_none(row.get("MeSH ID") or row.get("MeSH_ID") or row.get("MeSH"))
        if not inchikey or not mesh_id:
            continue

        ade_meta = ade_index.get(mesh_id, {"mesh_id": mesh_id})
        assoc_index[inchikey].append(ade_meta)

    return assoc_index


# ─── main loader ─────────────────────────────────────────────────────────────

def load_data(data_folder):
    """
    Parse PersADE/mtADENet data and yield BioThings-compatible documents.

    One document per unique InChIKey from Global_Drug file.
    Each document has drug metadata + list of ADE associations.

    Note: InChIKey column is named 'InChI Key' (with space) in Global_Drug
    and 'InChI Key' (with space) in Global_associations.
    """
    # Build supporting indices
    ade_index = _build_ade_index(data_folder)
    assoc_index = _build_association_index(data_folder, ade_index)

    # Main loop: iterate Global_Drug file
    drug_path = _find_file(data_folder, "Global_Drug_*.xlsx")
    assert drug_path, (
        "Global_Drug_*.xlsx not found in data_folder: " + data_folder
    )

    seen_ids = set()

    for row in _read_xlsx(drug_path):
        inchikey = _or_none(
            row.get("InChI Key") or row.get("InChIKey") or row.get("inchikey")
        )
        if not inchikey:
            continue
        if inchikey in seen_ids:
            continue
        seen_ids.add(inchikey)

        mesh_drug_id = _or_none(
            row.get("MeSH ID") or row.get("MeSH_ID") or row.get("Drug MeSH ID")
        )
        casrn = _or_none(row.get("CASRN") or row.get("CAS") or row.get("CASrn"))
        # CASRN sometimes has leading '$' from source formatting
        if casrn and casrn.startswith("$"):
            casrn = casrn[1:]

        drug_name = _or_none(row.get("Drug name") or row.get("drug_name") or row.get("Name"))
        smiles = _or_none(row.get("SMILES") or row.get("smiles"))
        formula = _or_none(row.get("Formula") or row.get("formula"))
        mw = _safe_float(row.get("MW") or row.get("mw") or row.get("Molecular Weight"))
        logp = _safe_float(row.get("logP") or row.get("logp"))
        atc_code = _or_none(row.get("ATC Code") or row.get("atc_code") or row.get("ATC"))
        drug_type = _or_none(row.get("Type") or row.get("drug_type") or row.get("Drug Type"))

        # Parse cross-references from pipe-delimited Link field
        link_str = row.get("Link") or row.get("link") or row.get("Links") or ""
        xref_dict = _parse_xref_links(link_str)

        # Also extract individual xref columns if present
        if not xref_dict.get("pubchem"):
            pc = _or_none(row.get("PubChem ID") or row.get("pubchem_id"))
            if pc:
                xref_dict["pubchem"] = pc
        if not xref_dict.get("chembl"):
            ch = _or_none(row.get("ChEMBL ID") or row.get("chembl_id"))
            if ch:
                xref_dict["chembl"] = ch

        xrefs = dict_sweep(xref_dict, [None]) if xref_dict else None

        # Resolve CASRN as separate xref
        if casrn and xrefs:
            xrefs["cas"] = casrn
        elif casrn:
            xrefs = {"cas": casrn}

        # Get ADE associations
        ade_associations = assoc_index.get(inchikey, []) or None

        doc = {
            "_id": inchikey,
            "persade": dict_sweep({
                "inchikey": inchikey,
                "mesh_drug_id": mesh_drug_id,
                "name": drug_name,
                "smiles": smiles,
                "properties": dict_sweep({
                    "formula": formula,
                    "mw": mw,
                    "logp": logp,
                }, [None]) or None,
                "atc_code": atc_code,
                "drug_type": drug_type,
                "xrefs": xrefs,
                "ade_count": len(ade_associations) if ade_associations else None,
                "ade_associations": ade_associations,
            }, [None])
        }
        doc = dict_sweep(unlist(doc), [None])
        yield doc
