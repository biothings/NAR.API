import csv
import os

from biothings.utils.dataload import dict_sweep, unlist

_FLOAT_FIELDS = {"mw", "tpsa", "fp3", "logp"}
_INT_FIELDS = {"hba", "hbd", "rb", "violates_ro5"}

_SUBLIBRARY_MAP = {
    "bioactives": "bioactives",
    "fragments": "fragments",
    "nuisance_set": "nuisance",
    "academic": "academic",
    "diverse_library": "diverse",
}

_FILE_ORDER = [
    "bioactives.csv",
    "fragments.csv",
    "nuisance_set.csv",
    "academic.csv",
    "diverse_library.csv",
]


def load_data(data_folder):
    seen_ids = set()

    for fname in _FILE_ORDER:
        fpath = os.path.join(data_folder, fname)
        if not os.path.exists(fpath):
            continue

        stem = fname.replace(".csv", "")
        sub_library = _SUBLIBRARY_MAP.get(stem, stem)

        with open(fpath, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                inchi_key = row.get("inchikey", "").strip()
                if not inchi_key:
                    continue
                if inchi_key in seen_ids:
                    continue
                seen_ids.add(inchi_key)

                props = {}
                for field in _FLOAT_FIELDS:
                    val = row.get(field, "").strip()
                    if val:
                        try:
                            props[field] = float(val)
                        except ValueError:
                            pass
                for field in _INT_FIELDS:
                    val = row.get(field, "").strip()
                    if val:
                        try:
                            props[field] = int(float(val))
                        except ValueError:
                            pass

                xrefs = {}
                if v := row.get("pubchem", "").strip():
                    xrefs["pubchem"] = v
                if v := row.get("chembl", "").strip():
                    xrefs["chembl"] = v
                if v := row.get("zinc", "").strip():
                    xrefs["zinc"] = v

                doc = {
                    "_id": inchi_key,
                    "ecbd": {
                        "eos_id": row.get("eos", "").strip() or None,
                        "inchikey": inchi_key,
                        "inchi": row.get("inchi", "").strip() or None,
                        "smiles": row.get("smiles", "").strip() or None,
                        "formula": row.get("formula", "").strip() or None,
                        "sub_library": sub_library,
                        "properties": props if props else None,
                        "xrefs": xrefs if xrefs else None,
                    },
                }
                doc = dict_sweep(unlist(doc), [None])
                yield doc
