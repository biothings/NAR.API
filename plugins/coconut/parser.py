import os
import csv
import glob
import logging

from biothings.utils.dataload import dict_sweep, unlist

logger = logging.getLogger(__name__)

# Pipe-delimited multi-value fields
MULTI_VALUE_FIELDS = {"organisms", "collections", "dois", "synonyms", "cas"}

FLOAT_FIELDS = {
    "molecular_weight", "exact_molecular_weight", "alogp",
    "topological_polar_surface_area", "np_likeness", "qed_drug_likeliness",
    "fractioncsp3", "van_der_walls_volume",
}

INT_FIELDS = {
    "total_atom_count", "heavy_atom_count", "rotatable_bond_count",
    "hydrogen_bond_acceptors", "hydrogen_bond_donors",
    "hydrogen_bond_acceptors_lipinski", "hydrogen_bond_donors_lipinski",
    "lipinski_rule_of_five_violations", "aromatic_rings_count",
    "formal_charge", "number_of_minimal_rings", "annotation_level",
}

BOOL_FIELDS = {
    "contains_sugar", "contains_ring_sugars", "contains_linear_sugars",
    "np_classifier_is_glycoside",
}


def _to_float(val):
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def _to_int(val):
    try:
        return int(float(val))
    except (ValueError, TypeError):
        return None


def _to_bool(val):
    if not val or not str(val).strip():
        return None
    return str(val).strip().lower() == "true"


def _split_pipe(val):
    if not val:
        return None
    parts = [v.strip() for v in val.split("|") if v.strip()]
    return parts if parts else None


def load_data(data_folder):
    """Parse COCONUT full CSV and yield BioThings-compatible documents.

    Each document is keyed by InChIKey and contains identifiers, molecular
    properties, chemical classification, NP classifier annotations, organism
    sources, collection provenance, and cross-references.
    """
    seen_ids = set()
    csv_files = sorted(glob.glob(os.path.join(data_folder, "coconut_csv*.csv")))
    if not csv_files:
        csv_files = sorted(glob.glob(os.path.join(data_folder, "*.csv")))
    assert csv_files, f"No CSV files found in {data_folder}"

    infile = csv_files[0]
    logger.info("Parsing %s", os.path.basename(infile))

    with open(infile, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            inchi_key = (row.get("standard_inchi_key") or "").strip()
            if not inchi_key or len(inchi_key) < 10:
                continue
            if inchi_key in seen_ids:
                continue
            seen_ids.add(inchi_key)

            identifiers = {
                "coconut_id": (row.get("identifier") or "").strip() or None,
                "inchi_key": inchi_key,
                "inchi": (row.get("standard_inchi") or "").strip() or None,
                "smiles": (row.get("canonical_smiles") or "").strip() or None,
                "name": (row.get("name") or "").strip() or None,
                "iupac_name": (row.get("iupac_name") or "").strip() or None,
                "molecular_formula": (row.get("molecular_formula") or "").strip() or None,
            }

            properties = {}
            for field in FLOAT_FIELDS:
                val = _to_float(row.get(field))
                if val is not None:
                    properties[field] = val
            for field in INT_FIELDS:
                val = _to_int(row.get(field))
                if val is not None:
                    properties[field] = val
            for field in BOOL_FIELDS:
                raw = (row.get(field) or "").strip()
                if raw:
                    properties[field] = _to_bool(raw)

            classification = {
                "chemical_class": (row.get("chemical_class") or "").strip() or None,
                "chemical_sub_class": (row.get("chemical_sub_class") or "").strip() or None,
                "chemical_super_class": (row.get("chemical_super_class") or "").strip() or None,
                "direct_parent": (row.get("direct_parent_classification") or "").strip() or None,
            }

            np_classifier = {
                "pathway": (row.get("np_classifier_pathway") or "").strip() or None,
                "superclass": (row.get("np_classifier_superclass") or "").strip() or None,
                "class": (row.get("np_classifier_class") or "").strip() or None,
                "is_glycoside": _to_bool(row.get("np_classifier_is_glycoside")),
            }

            xrefs = {}
            cas = _split_pipe(row.get("cas"))
            dois = _split_pipe(row.get("dois"))
            if cas:
                xrefs["cas"] = cas
            if dois:
                xrefs["doi"] = dois

            doc = {
                "_id": inchi_key,
                "coconut": {
                    **identifiers,
                    "properties": properties or None,
                    "murcko_framework": (row.get("murcko_framework") or "").strip() or None,
                    "classification": classification,
                    "np_classifier": np_classifier,
                    "organisms": _split_pipe(row.get("organisms")),
                    "collections": _split_pipe(row.get("collections")),
                    "synonyms": _split_pipe(row.get("synonyms")),
                    "xrefs": xrefs or None,
                },
            }

            doc = dict_sweep(unlist(doc), [None, ""])
            yield doc
