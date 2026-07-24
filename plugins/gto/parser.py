import os
import glob
import pandas as pd
from biothings.utils.dataload import dict_sweep, unlist


# Column name normalization map for clinical data file
_CLINICAL_RENAME = {
    "GTOID": "gtoid",
    "Year": "year",
    "Trial_ID": "trial_id",
    "link": "trial_link",
    "Country": "country",
    "Phase": "phase",
    "Status": "status",
    "Title": "title",
    "Major_Therapy_Category": "major_therapy_category",
    "Therapy_Type": "therapy_type",
    "Treatment": "treatment",
    "Location_Approved": "location_approved",
    "Co_Treatment": "co_treatment",
    "Altered_Gene": "altered_gene",
    "Target/Therapeutic_Gene": "target_gene",
    "Generation": "generation",
    "Vector": "vector",
    "Construct": "construct",
    "Vector_Type": "vector_type",
    "Transgene/Inserted Gene": "transgene",
    "Regulatory_Element": "regulatory_element",
    "Viral_Genome_Modification": "viral_genome_modification",
    "Vector_Production_Method": "vector_production_method",
    "Additional_Feature": "additional_feature",
    "Administration": "administration",
    "Dose": "dose",
    "Remark": "remark",
    "Disease_group": "disease_group",
    "Disease": "disease",
    "HLA": "hla",
    "Ex/In_Vivo": "ex_in_vivo",
    "Donor_Type": "donor_type",
    "Pts": "pts",
    "Age": "age",
    "Activation": "activation",
    "Lymph_depletion": "lymph_depletion",
    "Adverse_Reactions": "adverse_reactions",
    "CR": "cr",
    "PR": "pr",
    "SD": "sd",
    "Outcome": "outcome",
    "Company-Sponsor/Collab": "sponsor",
    "Other_IDs": "other_ids",
    "References": "references",
    "Ref_Link": "ref_link",
    "GTDID": "gtdid",
}

# Column name normalization map for indication file
_INDICATION_RENAME = {
    "study_ID": "gtoid",
    "trial_ID": "trial_id",
    "disease": "disease",
    "doid": "doid",
    "disease_ID": "disease_umls",
    "xref": "xref",
    "synonyms": "synonyms",
    "definition": "definition",
    "disgenet_disease_name": "disgenet_name",
    "disease_type": "disease_type",
    "diseaseclassmsh": "mesh_class_id",
    "diseaseclassnamemsh": "mesh_class_name",
    "hpoclassid": "hpo_class_id",
    "hpoclassname": "hpo_class_name",
    "doclassid": "do_class_id",
    "doclassname": "do_class_name",
    "umlssemantictypeid": "umls_semantic_type_id",
    "umlssemantictypename": "umls_semantic_type_name",
}


def _extract_mondo(xref_str):
    """Extract MONDO ID from the pipe-delimited xref string."""
    if not xref_str:
        return None
    for part in str(xref_str).split("|"):
        part = part.strip()
        if part.startswith("MONDO:"):
            return part
    return None


def _extract_mesh(xref_str):
    """Extract MeSH ID from the pipe-delimited xref string."""
    if not xref_str:
        return None
    for part in str(xref_str).split("|"):
        part = part.strip()
        if part.startswith("MESH:") or part.startswith("MSH:"):
            return part.replace("MSH:", "MESH:")
    return None


def _extract_omim(xref_str):
    """Extract OMIM ID from the pipe-delimited xref string."""
    if not xref_str:
        return None
    for part in str(xref_str).split("|"):
        part = part.strip()
        if part.startswith("OMIM:") or part.startswith("MIM:"):
            return part.replace("MIM:", "OMIM:")
    return None


def _split_pipe(value, sep="|"):
    """Split a pipe-delimited string into a list, stripping empties."""
    if pd.isna(value) or value is None:
        return None
    parts = [p.strip() for p in str(value).split(sep) if p.strip()]
    return parts if parts else None


def load_data(data_folder):
    """Parse GTO clinical trial and indication data, yield BioThings-compatible documents.

    Joins GTO_clinic_data.xlsx with indication data on GTOID to enrich each
    clinical trial record with DOID/MONDO/MeSH disease cross-references.
    Primary key: GTOID (e.g., GTC0001) — one document per clinical trial record.
    """
    # Locate files — biothings-cli downloads them with their Content-Disposition filenames
    # The indication file may be named 'indication' (no extension) or 'GTO_indication.xlsx'
    clinic_file = None
    indication_file = None

    for f in glob.glob(os.path.join(data_folder, "*")):
        basename = os.path.basename(f).lower()
        if "clinic_data" in basename or "clinical_data" in basename:
            clinic_file = f
        elif "indication" in basename:
            indication_file = f

    # Fallback: check for the exact Content-Disposition filename pattern
    if clinic_file is None:
        candidate = os.path.join(data_folder, "GTO_clinic_data.xlsx")
        if os.path.exists(candidate):
            clinic_file = candidate
    if indication_file is None:
        candidate = os.path.join(data_folder, "GTO_indication.xlsx")
        if os.path.exists(candidate):
            indication_file = candidate

    assert clinic_file is not None and os.path.exists(clinic_file), \
        f"GTO clinical data file not found in {data_folder}. Expected GTO_clinic_data.xlsx"

    # Load clinical data
    df_clinic = pd.read_excel(clinic_file, dtype=str)
    # pandas' default "str" dtype (since pandas 3.0) silently keeps float NaN
    # in string columns when using `.where(pd.notnull(df), None)` directly;
    # casting to object first makes the None replacement actually stick.
    df_clinic = df_clinic.astype(object).where(pd.notnull(df_clinic), None)
    df_clinic.rename(columns=_CLINICAL_RENAME, inplace=True)

    # Load indication data (optional — enrich if available)
    indication_map = {}  # gtoid -> list of indication rows
    if indication_file and os.path.exists(indication_file):
        df_indication = pd.read_excel(indication_file, dtype=str)
        df_indication = df_indication.astype(object).where(pd.notnull(df_indication), None)
        df_indication.rename(columns=_INDICATION_RENAME, inplace=True)
        for row in df_indication.to_dict(orient="records"):
            gid = row.get("gtoid")
            if not gid:
                continue
            if gid not in indication_map:
                indication_map[gid] = []
            indication_map[gid].append(row)

    seen_ids = set()

    for row in df_clinic.to_dict(orient="records"):
        gtoid = row.get("gtoid")
        if not gtoid:
            continue
        gtoid = str(gtoid).strip()
        if gtoid in seen_ids:
            continue
        seen_ids.add(gtoid)

        # Build base document
        gto_doc = {
            "gtoid": gtoid,
            "year": _to_int(row.get("year")),
            "trial_id": row.get("trial_id"),
            "trial_link": row.get("trial_link"),
            "country": row.get("country"),
            "phase": row.get("phase"),
            "status": row.get("status"),
            "title": row.get("title"),
            "major_therapy_category": row.get("major_therapy_category"),
            "therapy_type": row.get("therapy_type"),
            "treatment": row.get("treatment"),
            "location_approved": row.get("location_approved"),
            "co_treatment": row.get("co_treatment"),
            "altered_gene": row.get("altered_gene"),
            "target_gene": row.get("target_gene"),
            "generation": row.get("generation"),
            "vector": row.get("vector"),
            "construct": row.get("construct"),
            "vector_type": row.get("vector_type"),
            "transgene": row.get("transgene"),
            "regulatory_element": row.get("regulatory_element"),
            "viral_genome_modification": row.get("viral_genome_modification"),
            "vector_production_method": row.get("vector_production_method"),
            "additional_feature": row.get("additional_feature"),
            "administration": row.get("administration"),
            "dose": row.get("dose"),
            "disease_group": row.get("disease_group"),
            "disease": row.get("disease"),
            "hla": row.get("hla"),
            "ex_in_vivo": row.get("ex_in_vivo"),
            "donor_type": row.get("donor_type"),
            "pts": _to_int(row.get("pts")),
            "age": row.get("age"),
            "activation": row.get("activation"),
            "lymph_depletion": row.get("lymph_depletion"),
            "adverse_reactions": row.get("adverse_reactions"),
            "outcome": row.get("outcome"),
            "sponsor": row.get("sponsor"),
            "other_ids": row.get("other_ids"),
            "references": row.get("references"),
            "ref_link": row.get("ref_link"),
            "gtdid": row.get("gtdid"),
        }

        # Add indication cross-references if available
        indications = indication_map.get(gtoid, [])
        if indications:
            ind = indications[0]  # primary indication for this trial
            doid = ind.get("doid")
            xref_str = ind.get("xref")
            gto_doc["disease_xrefs"] = {
                "doid": doid,
                "mondo": _extract_mondo(xref_str),
                "mesh": _extract_mesh(xref_str),
                "omim": _extract_omim(xref_str),
                "umls": ind.get("disease_umls"),
                "disease_type": ind.get("disease_type"),
                "mesh_class": ind.get("mesh_class_name"),
                "hpo_class": ind.get("hpo_class_name"),
            }

        doc = {
            "_id": gtoid,
            "gto": gto_doc,
        }
        doc = dict_sweep(unlist(doc), [None, "None", "nan"])
        yield doc


def _to_int(value):
    """Convert string to int, returning None if not convertible."""
    if value is None:
        return None
    try:
        v = str(value).strip()
        if not v or v.lower() in ("none", "nan", ""):
            return None
        return int(float(v))
    except (ValueError, TypeError):
        return None
