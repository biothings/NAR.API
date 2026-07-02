"""
Chem(Pro)2 parser — atlas of chemoproteomic probes labelling human proteins.

Primary entity: chemoproteomic probe (CPP) or competitor compound.
_id strategy: InChIKey (MyChem.info standard).
Structure: one document per unique InChIKey, with nested probe metadata and
a list of experiment records linking the probe to protein targets.

Files consumed:
  general_probe.txt        — probe structures + chemical properties
  general_competitor.txt   — competitor compound structures + properties
  chemoproteomics_experiment.txt — experiment-level probe/competitor-target links
  general_cell.txt         — cell system metadata (joined by cell_id)
  general_target_*.txt     — target metadata (joined by target_id)
"""

import os
import csv
from collections import defaultdict
from biothings.utils.dataload import dict_sweep, unlist

# ─── helpers ────────────────────────────────────────────────────────────────

def _safe_float(val):
    try:
        return float(val) if val and val.strip() not in (".", "", "None") else None
    except (ValueError, AttributeError):
        return None


def _safe_int(val):
    try:
        return int(val) if val and val.strip() not in (".", "", "None") else None
    except (ValueError, AttributeError):
        return None


def _or_none(val):
    if val is None:
        return None
    s = str(val).strip()
    return s if s and s not in (".", "None", "N/A") else None


def _read_tsv(path, delimiter="\t"):
    """Yield rows as dicts from a TSV file, skipping empty lines."""
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for row in reader:
            yield row


# ─── target index builder ────────────────────────────────────────────────────

def _build_target_index(data_folder):
    """Build a dict mapping target_id → minimal target info from all target files."""
    target_info = {}
    target_files = [
        "general_target_enzyme.txt",
        "general_target_gpcr.txt",
        "general_target_immunoglobulin.txt",
        "general_target_transporter.txt",
        "general_target_other.txt",
    ]
    for fname in target_files:
        fpath = os.path.join(data_folder, fname)
        if not os.path.exists(fpath):
            continue
        for row in _read_tsv(fpath):
            tid = row.get("targetid", "").strip()
            if not tid:
                continue
            target_info[tid] = {
                "target_id": tid,
                "name": _or_none(row.get("targetname")),
                "gene_symbol": _or_none(row.get("genename")),
                "uniprot_id": _or_none(row.get("uniprotid")),
                "entrez_gene_id": _safe_int(row.get("geneid")),
                "bioclass": _or_none(row.get("bioclass")),
                "hgnc_id": _or_none(row.get("hgncid")),
                "chembl_id": _or_none(row.get("chemblid")),
                "ensembl_id": _or_none(row.get("ensemblid")),
            }
    return target_info


def _build_cell_index(data_folder):
    """Build a dict mapping cell_id → cell info."""
    cell_info = {}
    fpath = os.path.join(data_folder, "general_cell.txt")
    if not os.path.exists(fpath):
        return cell_info
    for row in _read_tsv(fpath):
        cid = row.get("cell_id", "").strip()
        if not cid:
            continue
        cell_info[cid] = {
            "cell_id": cid,
            "cell_name": _or_none(row.get("cellname_full")),
            "model_type": _or_none(row.get("model_type")),
            "disease": _or_none(row.get("cell_diseases_name")),
            "tissue": _or_none(row.get("cell_tissue")),
            "species": _or_none(row.get("Species")),
            "cellosaurus_accession": _or_none(row.get("cellosaurus_accession")),
        }
    return cell_info


def _build_experiment_index(data_folder, target_index, cell_index):
    """
    Build a dict mapping probeid → list of experiment records.
    Also returns competitor_experiments: competitorid → list of experiment records.
    """
    probe_experiments = defaultdict(list)
    competitor_experiments = defaultdict(list)

    fpath = os.path.join(data_folder, "chemoproteomics_experiment.txt")
    if not os.path.exists(fpath):
        return probe_experiments, competitor_experiments

    for row in _read_tsv(fpath):
        probeid = _or_none(row.get("probeid"))
        cpid = _or_none(row.get("cpid"))  # competitor id
        method_id = _or_none(row.get("method_id"))
        ref_id = _or_none(row.get("referenceid"))
        criteria = _or_none(row.get("criteria"))
        probe_conc = _or_none(row.get("probe_concentration"))
        cp_conc = _or_none(row.get("cp_concentration"))
        exp_method = _or_none(row.get("experiment_method"))
        quant_method = _or_none(row.get("Quantitative Method"))
        cell_id = _or_none(row.get("cell_id"))

        cell_data = dict_sweep(cell_index.get(cell_id, {"cell_id": cell_id}), [None]) if cell_id else None

        exp_record = dict_sweep({
            "method_id": method_id,
            "reference_id": ref_id,
            "criteria": criteria,
            "probe_concentration": probe_conc,
            "cp_concentration": cp_conc,
            "experiment_method": exp_method,
            "quantitative_method": quant_method,
            "cell": cell_data,
        }, [None])

        if probeid:
            probe_experiments[probeid].append(exp_record)
        if cpid:
            competitor_experiments[cpid].append(exp_record)

    return probe_experiments, competitor_experiments


# ─── chemical property extractor ─────────────────────────────────────────────

def _extract_chem_props(row):
    return dict_sweep({
        "mw": _safe_float(row.get("mw")),
        "mf": _or_none(row.get("mf")),
        "polar_area": _safe_float(row.get("polararea")),
        "complexity": _safe_float(row.get("complexity")),
        "xlogp": _safe_float(row.get("xlogp")),
        "heavy_atom_count": _safe_int(row.get("heavycnt")),
        "hbond_donor": _safe_int(row.get("hbonddonor")),
        "hbond_acceptor": _safe_int(row.get("hbondacc")),
        "rotatable_bonds": _safe_int(row.get("rotbonds")),
    }, [None])


# ─── main loader ─────────────────────────────────────────────────────────────

def load_data(data_folder):
    """Parse Chem(Pro)2 data and yield BioThings-compatible documents."""

    # Build supporting indices
    target_index = _build_target_index(data_folder)
    cell_index = _build_cell_index(data_folder)
    probe_experiments, competitor_experiments = _build_experiment_index(
        data_folder, target_index, cell_index
    )

    seen_ids = set()

    # ── Part 1: Probes ────────────────────────────────────────────────────────
    probe_path = os.path.join(data_folder, "general_probe.txt")
    assert os.path.exists(probe_path), f"Expected file not found: {probe_path}"

    for row in _read_tsv(probe_path):
        inchikey = _or_none(row.get("inchikey"))
        if not inchikey:
            continue
        if inchikey in seen_ids:
            continue
        seen_ids.add(inchikey)

        probeid = _or_none(row.get("probeid"))

        doc = {
            "_id": inchikey,
            "chemprob": {
                "probe_id": probeid,
                "name": _or_none(row.get("name_show")),
                "probe_type": _or_none(row.get("type")),
                "entity_type": "probe",
                "inchi": _or_none(row.get("inchi")),
                "smiles": _or_none(row.get("isosmiles")),
                "iupac_name": _or_none(row.get("iupacname")),
                "properties": _extract_chem_props(row),
                "xrefs": dict_sweep({
                    "pubchem": _or_none(row.get("pubchemid")),
                    "synonyms": _or_none(row.get("cmpdsynonym")),
                }, [None]),
                "experiments": probe_experiments.get(probeid, []) or None,
            }
        }
        doc = dict_sweep(unlist(doc), [None])
        yield doc

    # ── Part 2: Competitors ───────────────────────────────────────────────────
    comp_path = os.path.join(data_folder, "general_competitor.txt")
    assert os.path.exists(comp_path), f"Expected file not found: {comp_path}"

    for row in _read_tsv(comp_path):
        inchikey = _or_none(row.get("inchikey"))
        if not inchikey:
            continue
        if inchikey in seen_ids:
            # Probe already created doc — competitor data is secondary for same compound
            continue
        seen_ids.add(inchikey)

        cpid = _or_none(row.get("cpid"))

        doc = {
            "_id": inchikey,
            "chemprob": {
                "competitor_id": cpid,
                "name": _or_none(row.get("name_competitor")),
                "entity_type": "competitor",
                "inchi": _or_none(row.get("inchi")),
                "smiles": _or_none(row.get("isosmiles")),
                "iupac_name": _or_none(row.get("iupacname")),
                "properties": _extract_chem_props(row),
                "xrefs": dict_sweep({
                    "pubchem": _or_none(row.get("pubchemid")),
                    "synonyms": _or_none(row.get("cmpdsynonym")),
                }, [None]),
                "experiments": competitor_experiments.get(cpid, []) or None,
            }
        }
        doc = dict_sweep(unlist(doc), [None])
        yield doc
