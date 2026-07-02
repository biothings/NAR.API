import os
import csv
import glob
import collections

from biothings.utils.dataload import dict_sweep, unlist


_GENE_TYPE_MAP = {
    "protein coding": "protein_coding",
}

_INTERACTION_TYPE_MAP = {
    "BSJ supported": "BSJ_supported",
    "nonBSJ supported": "nonBSJ_supported",
}


def _extract_rar(data_folder):
    """Extract all.rar → all.txt in data_folder; return path to all.txt."""
    rar_path = os.path.join(data_folder, "all.rar")
    txt_path = os.path.join(data_folder, "all.txt")
    if not os.path.exists(txt_path):
        import rarfile
        with rarfile.RarFile(rar_path) as rf:
            rf.extract("all.txt", path=data_folder)
    return txt_path


def _parse_pvalue(raw):
    try:
        return float(raw)
    except (ValueError, TypeError):
        return None


def load_data(data_folder):
    """Parse CircTarget circRNA-target RNA interactions; yield one doc per circRNA."""
    txt_path = _extract_rar(data_folder)
    assert os.path.exists(txt_path), f"Expected file not found: {txt_path}"

    groups = collections.defaultdict(list)

    with open(txt_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            circ_id = row.get("circRNA ID", "").strip()
            if not circ_id:
                continue
            gene_type = row.get("Gene type", "").strip()
            gene_type = _GENE_TYPE_MAP.get(gene_type, gene_type)

            interaction_type = row.get("Interaction type", "").strip()
            interaction_type = _INTERACTION_TYPE_MAP.get(interaction_type, interaction_type)

            p_val = _parse_pvalue(row.get("P-value", ""))
            try:
                chimeric = int(row.get("Chimeric read count", ""))
            except (ValueError, TypeError):
                chimeric = None

            interaction = {
                "target": {
                    "ensembl_id": row.get("Ensembl ID", "").strip() or None,
                    "gene_name": row.get("Gene name", "").strip() or None,
                    "gene_type": gene_type or None,
                },
                "chimeric_read_count": chimeric,
                "p_value": p_val,
                "cell_line": row.get("Cell line/Tissue", "").strip() or None,
                "species": row.get("Species", "").strip() or None,
                "interaction_type": interaction_type or None,
                "detected_method": row.get("Detected method", "").strip() or None,
            }
            groups[circ_id].append(interaction)

    for circ_id, interactions in groups.items():
        species_set = {i["species"] for i in interactions if i.get("species")}
        species = list(species_set)[0] if len(species_set) == 1 else list(species_set)

        doc = {
            "_id": circ_id,
            "circtarget": {
                "circrna_id": circ_id,
                "species": species,
                "interaction_count": len(interactions),
                "interactions": interactions,
            },
        }
        doc = dict_sweep(unlist(doc), [None])
        yield doc
