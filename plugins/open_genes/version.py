def get_release(self):
    """Return a release string for the Open Genes datasource.

    Open Genes has no dedicated version/release endpoint. Its REST API
    (https://open-genes.com/api/gene/search) does, however, report a
    per-gene `timestamp.changed` (unix epoch of last curation edit) and
    an `options.objTotal` gene count. We fetch the full result set once
    and derive a release string from the most recent per-gene edit
    timestamp plus the current total gene count, e.g. "2405-genes-2026-06-30".
    This changes whenever genes are added/removed or any record is edited,
    which is the closest available proxy for a dataset release version.
    """
    import requests
    from datetime import datetime, timezone

    url = "https://open-genes.com/api/gene/search"
    resp = requests.get(url, params={"pageSize": 5000, "page": 1}, timeout=30)
    resp.raise_for_status()
    data = resp.json()

    total = data.get("options", {}).get("objTotal")
    items = data.get("items", [])

    latest_epoch = 0
    for item in items:
        ts = (item.get("timestamp") or {}).get("changed")
        if isinstance(ts, (int, float)) and ts > latest_epoch:
            latest_epoch = ts

    if latest_epoch:
        latest_date = datetime.fromtimestamp(latest_epoch, tz=timezone.utc).strftime("%Y-%m-%d")
    else:
        latest_date = "unknown-date"

    if total is not None:
        return f"{total}-genes-{latest_date}"
    return latest_date if latest_date != "unknown-date" else None
