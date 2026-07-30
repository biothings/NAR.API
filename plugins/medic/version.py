"""
version.py — MeDIC release string

Fetches the latest release tags from the GitHub APIs for the three MeDIC repos
and returns the combined version string.

Returns: "{drug_list_tag}_{indication_list_tag}" e.g. "v2.3.1_v1.4.1"
Fallback: current date as YYYYMMDD.
"""


def get_release(self):
    import requests
    from datetime import datetime

    drug_tag = None
    indication_tag = None

    # Strategy 1: Query GitHub API for latest release tags
    try:
        resp1 = requests.get(
            "https://api.github.com/repos/everycure-org/matrix-drug-list/releases/latest",
            timeout=30,
            headers={"Accept": "application/vnd.github+json"},
        )
        if resp1.status_code == 200:
            drug_tag = resp1.json().get("tag_name")
    except Exception:
        pass

    try:
        resp2 = requests.get(
            "https://api.github.com/repos/everycure-org/matrix-indication-list/releases/latest",
            timeout=30,
            headers={"Accept": "application/vnd.github+json"},
        )
        if resp2.status_code == 200:
            indication_tag = resp2.json().get("tag_name")
    except Exception:
        pass

    if drug_tag and indication_tag:
        return f"{drug_tag}_{indication_tag}"
    elif drug_tag:
        return drug_tag
    elif indication_tag:
        return indication_tag

    # Fallback: today's date
    return datetime.utcnow().strftime("%Y%m%d")
