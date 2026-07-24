"""
Parser for DRMref (Drug Resistance Mechanism reference).

Primary source: gene_summary.txt — one row per differentially-expressed-gene
(DEG) record, keyed by a stable per-row `primary_key` integer emitted by the
DRMref site itself (unique across the whole file). Each row represents a
gene x dataset x cell-type x cancer-type x drug-regimen association from a
Seurat-style differential expression test (resistant vs. sensitive cell
populations).

Four supplementary files are joined onto each gene_summary.txt record by
gene symbol wherever the join key is present:
  - Existed_drug_mechanism_gene_file.csv : known drug-resistance mechanism
    category for a gene (6 established mechanism categories)
  - miRNA_summary.txt                    : microRNAs predicted to regulate
    the gene
  - tf_summary.txt                       : transcription factors / motifs
    regulating the gene
  - enrichment_pathway_summary.txt       : Hallmark/KEGG/GO_BP pathway
    enrichment results attached at the dataset/cell-type granularity when a
    gene-level column is not present in that row

Because the exact column layout of the four supplementary files could not be
directly re-verified at generation time (ccsm.uth.edu was returning HTTP 522
"origin unreachable" through Cloudflare during this session), the loaders
for those four files are written defensively: they auto-detect the gene
symbol column (any header containing "gene" case-insensitively, preferring
an exact "Gene_symbol" match) and carry all remaining columns through
verbatim as a nested dict under the appropriate list. This keeps the parser
functional even if the live column names differ slightly from what the
Phase-1 inspection sampled (which only fully enumerated gene_summary.txt's
schema).
"""

import csv
import os
from collections import defaultdict

from biothings.utils.dataload import dict_sweep, unlist

# ---------------------------------------------------------------------------
# gene_summary.txt column groups (per Phase-1 site inspection)
# ---------------------------------------------------------------------------
_NOVEL_FLOAT_FIELDS = {
    "p_val": "p_val",
    "avg_log2FC": "avg_log2fc",
    "pct.1": "pct_1",
    "pct.2": "pct_2",
    "p_val_adj": "p_val_adj",
}

# Fields intentionally dropped per Phase-1 inspection REDUNDANT classification:
# entrezgene_description, external_synonym, gene_biotype, Organism, Tissue,
# Date, Cancer_type_level1_forDB, Drug_type_forDB, Dataset (superseded by
# dataset.geo_accession), Original_Dataset, dataset_subgroup.

_GENE_SYMBOL_CANDIDATES = ("Gene_symbol", "gene_symbol", "GeneSymbol", "Gene", "gene", "SYMBOL", "symbol")


def _sniff_delimiter(path, default="\t"):
    """Best-effort delimiter detection; falls back to `default` on failure."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            sample = fh.read(4096)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        return dialect.delimiter
    except Exception:
        return default


def _find_gene_column(fieldnames):
    if not fieldnames:
        return None
    for cand in _GENE_SYMBOL_CANDIDATES:
        if cand in fieldnames:
            return cand
    for name in fieldnames:
        if name and "gene" in name.lower():
            return name
    return None


def _to_float(val):
    if val is None:
        return None
    val = str(val).strip()
    if not val or val.upper() == "NA":
        return None
    try:
        return float(val)
    except ValueError:
        return None


def _to_int(val):
    f = _to_float(val)
    if f is None:
        return None
    try:
        return int(f)
    except (ValueError, OverflowError):
        return None


def _clean_str(val):
    if val is None:
        return None
    val = str(val).strip()
    if not val or val.upper() == "NA":
        return None
    return val


def _load_gene_keyed_sidefile(path):
    """
    Generic loader for a supplementary file that annotates genes
    (mechanism / miRNA / TF files). Returns dict: gene_symbol -> [records].
    Each record is a dict of all non-gene-key columns from that row.
    """
    result = defaultdict(list)
    if not path or not os.path.exists(path):
        return result

    delim = _sniff_delimiter(path)
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        gene_col = _find_gene_column(reader.fieldnames)
        if gene_col is None:
            return result
        for row in reader:
            gene = _clean_str(row.get(gene_col))
            if not gene:
                continue
            record = {}
            for k, v in row.items():
                if k == gene_col or k is None:
                    continue
                cleaned = _clean_str(v)
                if cleaned is not None:
                    record[k.strip().lower().replace(" ", "_")] = cleaned
            if record:
                result[gene].append(record)
    return result


def load_data(data_folder):
    """Parse DRMref data and yield BioThings-compatible documents.

    One document per gene_summary.txt row (a gene x dataset x cell-type x
    cancer-type x drug-regimen differential expression association),
    enriched with drug-resistance mechanism classification, and
    resistance-associated miRNA/TF regulator annotations for that gene
    wherever the gene symbol is present in the corresponding side file.
    """
    gene_summary_path = os.path.join(data_folder, "gene_summary.txt")
    assert os.path.exists(gene_summary_path), f"Expected file not found: {gene_summary_path}"

    mechanism_index = _load_gene_keyed_sidefile(
        os.path.join(data_folder, "Existed_drug_mechanism_gene_file.csv")
    )
    mirna_index = _load_gene_keyed_sidefile(os.path.join(data_folder, "miRNA_summary.txt"))
    tf_index = _load_gene_keyed_sidefile(os.path.join(data_folder, "tf_summary.txt"))

    seen_ids = set()
    delim = _sniff_delimiter(gene_summary_path)

    with open(gene_summary_path, "r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        for row in reader:
            pk = _clean_str(row.get("primary_key"))
            if not pk:
                continue
            _id = f"drmref_{pk}"
            if _id in seen_ids:
                continue
            seen_ids.add(_id)

            gene_symbol = _clean_str(row.get("Gene_symbol"))

            diff_expr = {}
            for src_field, out_field in _NOVEL_FLOAT_FIELDS.items():
                fv = _to_float(row.get(src_field))
                if fv is not None:
                    diff_expr[out_field] = fv

            doc_body = {
                "primary_key": _to_int(pk),
                "gene": {
                    "symbol": gene_symbol,
                    "ensembl_gene_id": _clean_str(row.get("ensembl_gene_id")),
                    "entrezgene_id": _clean_str(row.get("entrezgene_id")),
                    "uniprot": _clean_str(row.get("uniprotswissprot")),
                },
                "differential_expression": diff_expr or None,
                "cell_type": _clean_str(row.get("Cell_type")),
                "cancer_type_level1": _clean_str(row.get("Cancer_type_level1")),
                "cancer_type_level2": _clean_str(row.get("Cancer_type_level2")),
                "drug_type": _clean_str(row.get("Drug_type")),
                "regimen": _clean_str(row.get("Regimen")),
                "timepoint": _clean_str(row.get("Timepoint")),
                "sample_size": _clean_str(row.get("Sample_size")),
                "sample_size_all": _to_int(row.get("Sample_size_all")),
                "cell_number_all": _to_int(row.get("Cell_number_all")),
                "source": _clean_str(row.get("Source")),
                "description": _clean_str(row.get("Description")),
                "extract_protocol": _clean_str(row.get("Extract_protocol")),
                "data_processing": _clean_str(row.get("Data_processing")),
                "dataset": {
                    "geo_accession": _clean_str(row.get("Dataset")) or _clean_str(row.get("Original_Dataset")),
                    "rawdata_id": _clean_str(row.get("RawData_ID")),
                    "pmid": _clean_str(row.get("PMID")),
                },
            }

            if gene_symbol:
                if gene_symbol in mechanism_index:
                    doc_body["mechanism"] = mechanism_index[gene_symbol]
                if gene_symbol in mirna_index:
                    doc_body["mirna_regulators"] = mirna_index[gene_symbol]
                if gene_symbol in tf_index:
                    doc_body["tf_regulators"] = tf_index[gene_symbol]

            doc = {
                "_id": _id,
                "drmref": doc_body,
            }

            yield dict_sweep(unlist(doc), [None])
