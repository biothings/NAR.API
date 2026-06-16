import csv
import os

from biothings.utils.dataload import dict_sweep, unlist

# Column indices for getData.php?organism=9606 (no header row, 28 columns)
# Derived from official release file header (Apr2026_release.txt)
_COL = {
    "ENTITYA":           0,
    "TYPEA":             1,
    "IDA":               2,
    "DATABASEA":         3,
    "ENTITYB":           4,
    "TYPEB":             5,
    "IDB":               6,
    "DATABASEB":         7,
    "EFFECT":            8,
    "MECHANISM":         9,
    "RESIDUE":          10,
    "SEQUENCE":         11,
    "TAX_ID":           12,
    "CELL_DATA":        13,
    "TISSUE_DATA":      14,
    "MODULATOR_COMPLEX":15,
    "TARGET_COMPLEX":   16,
    "MODIFICATIONA":    17,
    "MODASEQ":          18,
    "MODIFICATIONB":    19,
    "MODBSEQ":          20,
    "PMID":             21,
    "DIRECT":           22,
    "NOTES":            23,
    "ANNOTATOR":        24,
    "SENTENCE":         25,
    "SIGNOR_ID":        26,
    "SIGNOR_SCORE":     27,
}


def _get(row, name):
    idx = _COL[name]
    val = row[idx].strip() if idx < len(row) else ""
    return val if val else None


def _split_bto(val):
    if not val:
        return None
    parts = [v.strip() for v in val.split(";") if v.strip()]
    return parts if len(parts) > 1 else parts[0]


def load_data(data_folder):
    infile = os.path.join(data_folder, "signor_human_interactions.tsv")
    if not os.path.exists(infile):
        # biothings dumper saves with URL-derived filename
        candidates = [f for f in os.listdir(data_folder) if f.startswith("getData")]
        if not candidates:
            raise FileNotFoundError(
                f"SIGNOR data file not found in {data_folder}. "
                "Expected 'signor_human_interactions.tsv' or 'getData.php?organism=9606'."
            )
        infile = os.path.join(data_folder, candidates[0])

    seen_ids = set()

    with open(infile, "r", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 27:
                continue

            signor_id = _get(row, "SIGNOR_ID")
            if not signor_id or not signor_id.startswith("SIGNOR-"):
                continue
            if signor_id in seen_ids:
                continue
            seen_ids.add(signor_id)

            # Entity A
            entity_a = {}
            if v := _get(row, "ENTITYA"):
                entity_a["name"] = v
            if v := _get(row, "TYPEA"):
                entity_a["type"] = v
            if v := _get(row, "IDA"):
                entity_a["id"] = v
            if v := _get(row, "DATABASEA"):
                entity_a["database"] = v
            if v := _get(row, "MODIFICATIONA"):
                entity_a["modification"] = v
            if v := _get(row, "MODASEQ"):
                entity_a["modification_sequence"] = v
            if v := _get(row, "MODULATOR_COMPLEX"):
                entity_a["complex_id"] = v

            # Entity B
            entity_b = {}
            if v := _get(row, "ENTITYB"):
                entity_b["name"] = v
            if v := _get(row, "TYPEB"):
                entity_b["type"] = v
            if v := _get(row, "IDB"):
                entity_b["id"] = v
            if v := _get(row, "DATABASEB"):
                entity_b["database"] = v
            if v := _get(row, "MODIFICATIONB"):
                entity_b["modification"] = v
            if v := _get(row, "MODBSEQ"):
                entity_b["modification_sequence"] = v
            if v := _get(row, "TARGET_COMPLEX"):
                entity_b["complex_id"] = v

            interaction = {
                "signor_id": signor_id,
                "entity_a": entity_a if entity_a else None,
                "entity_b": entity_b if entity_b else None,
                "effect": _get(row, "EFFECT"),
                "mechanism": _get(row, "MECHANISM"),
            }

            if v := _get(row, "RESIDUE"):
                interaction["residue"] = v
            if v := _get(row, "SEQUENCE"):
                interaction["sequence"] = v

            if v := _get(row, "TAX_ID"):
                interaction["tax_id"] = v

            if v := _get(row, "CELL_DATA"):
                interaction["cell_data"] = _split_bto(v)
            if v := _get(row, "TISSUE_DATA"):
                interaction["tissue_data"] = _split_bto(v)

            if v := _get(row, "PMID"):
                interaction["pmid"] = v

            if v := _get(row, "DIRECT"):
                interaction["direct"] = v == "t"

            if v := _get(row, "NOTES"):
                interaction["notes"] = v
            if v := _get(row, "SENTENCE"):
                interaction["sentence"] = v
            if v := _get(row, "ANNOTATOR"):
                interaction["annotator"] = v

            if v := _get(row, "SIGNOR_SCORE"):
                try:
                    interaction["score"] = float(v)
                except ValueError:
                    pass

            doc = {
                "_id": signor_id,
                "signor": interaction,
            }
            doc = dict_sweep(unlist(doc), [None])
            yield doc
