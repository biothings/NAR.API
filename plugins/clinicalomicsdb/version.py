def get_release(self):
    """Return the ClinicalOmicsDB (trials.linkedomics.org) API release string.

    Primary source: the API's own published OpenAPI 3.1 spec
    (https://trials.linkedomics.org/trials_api.json), which carries an
    `info.version` field. Falls back to a Last-Modified header check on one
    of the bulk-downloaded per-study endpoints, then to today's date.
    """
    import requests
    from datetime import datetime, timezone

    try:
        resp = requests.get(
            "https://trials.linkedomics.org/trials_api.json",
            timeout=30,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        resp.raise_for_status()
        spec = resp.json()
        api_version = spec.get("info", {}).get("version")
        if api_version:
            return f"api-v{api_version}"
    except Exception:
        pass

    try:
        resp = requests.head(
            "https://trials.linkedomics.org/api/table/study/gene/GSE14764.csv",
            timeout=30,
            allow_redirects=True,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        last_modified = resp.headers.get("Last-Modified")
        if last_modified:
            dt = datetime.strptime(last_modified, "%a, %d %b %Y %H:%M:%S %Z")
            return dt.strftime("%Y%m%d")
    except Exception:
        pass

    return datetime.now(timezone.utc).strftime("%Y%m%d")
