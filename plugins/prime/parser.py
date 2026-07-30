import os
import csv
from biothings.utils.dataload import dict_sweep, unlist

# Columns that map to structured sub-objects in the document.
# Any CSV column NOT in this set is treated as a sparse host/clinical field.
_KNOWN_COLS = {
    "Run", "Experiment", "BioProject", "BioSample", "SRA_Study",
    "Sample_Name", "Participant_Id", "Project_name",
    "AvgSpotLen", "Bases", "Instrument", "Library_Layout",
    "Sequencing_Type", "Variable_Region",
    "Primer_Cut", "Primer_Par", "Denoise_Par", "Collapse",
    "Sequencing_Quality", "Filter_Pass",
    "Collection_Date", "Country", "Continent",
    "Body_Site", "Systems", "Phenotype", "Study_Group",
    "Time_series", "Comparison", "Matched",
    "Release_Date", "Create_Date",
}


def _safe_int(val):
    if not val or not str(val).strip():
        return None
    try:
        return int(str(val).strip())
    except (ValueError, TypeError):
        return None


def _parse_bool(val):
    if not val:
        return None
    return str(val).strip().lower() in ("yes", "true", "1")


def load_data(data_folder):
    """Parse PRIME samples_metadata.csv and yield one document per SRA run."""
    infile = os.path.join(data_folder, "samples_metadata.csv")
    assert os.path.exists(infile), f"Expected file not found: {infile}"

    # utf-8-sig strips a leading BOM (﻿) that Zenodo CSVs sometimes carry,
    # which would otherwise corrupt the first column name and silently skip all rows.
    with open(infile, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        # Normalize fieldnames once so _KNOWN_COLS lookup is BOM/whitespace safe.
        if reader.fieldnames:
            reader.fieldnames = [fn.strip() for fn in reader.fieldnames]
        for row in reader:
            run_id = row.get("Run", "").strip()
            if not run_id:
                continue

            # Sparse clinical/host fields: any column not in _KNOWN_COLS
            host = {}
            for col, val in row.items():
                if col not in _KNOWN_COLS and val and str(val).strip():
                    host[col.strip().lower()] = str(val).strip()

            doc = {
                "_id": run_id,
                "prime": {
                    "run": run_id,
                    "experiment": row.get("Experiment") or None,
                    "project": {
                        "bioproject": row.get("BioProject") or None,
                        "sra_study": row.get("SRA_Study") or None,
                        "name": row.get("Project_name") or None,
                    },
                    "sample": {
                        "biosample": row.get("BioSample") or None,
                        "name": row.get("Sample_Name") or None,
                        "participant_id": row.get("Participant_Id") or None,
                    },
                    "sequencing": {
                        "instrument": row.get("Instrument") or None,
                        "library_layout": row.get("Library_Layout") or None,
                        "type": row.get("Sequencing_Type") or None,
                        "variable_region": row.get("Variable_Region") or None,
                        "avg_spot_len": _safe_int(row.get("AvgSpotLen")),
                        "bases": _safe_int(row.get("Bases")),
                    },
                    "processing": {
                        "primer_cut": row.get("Primer_Cut") or None,
                        "primer_params": row.get("Primer_Par") or None,
                        "denoise_params": row.get("Denoise_Par") or None,
                        "taxonomy_db": row.get("Collapse") or None,
                        "sequencing_quality": row.get("Sequencing_Quality") or None,
                        "filter_pass": row.get("Filter_Pass") or None,
                    },
                    "geography": {
                        "country": row.get("Country") or None,
                        "continent": row.get("Continent") or None,
                    },
                    "phenotype": {
                        "body_site": row.get("Body_Site") or None,
                        "system": row.get("Systems") or None,
                        "phenotype_category": row.get("Phenotype") or None,
                        "study_group": row.get("Study_Group") or None,
                    },
                    "flags": {
                        "time_series": _parse_bool(row.get("Time_series")),
                        "comparison": _parse_bool(row.get("Comparison")),
                        "matched": _parse_bool(row.get("Matched")),
                    },
                    "collection_date": row.get("Collection_Date") or None,
                    "release_date": row.get("Release_Date") or None,
                    "create_date": row.get("Create_Date") or None,
                },
            }

            if host:
                doc["prime"]["host"] = host

            doc = dict_sweep(unlist(doc), [None])
            yield doc
