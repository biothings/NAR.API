def get_release(self):
    """Determine the current release string for the AgeAnnoMO GitHub companion repo.

    Strategy:
    1. Query the GitHub Commits API for the most recent commit touching any of the
       four data files this plugin ingests. Return the newest commit date as YYYYMMDD.
    2. Fall back to the repo's latest tagged GitHub release (v1.0, published 2023-09-30).
    3. Fall back to a hardcoded date matching the known v1.0 snapshot.
    """
    import requests

    paths = [
        "Genomic instability/Differential expression.xlsx",
        "Loss of proteostasis/Differential protein.xlsx",
        "Dysregulated metabolism in aging/Differential metabolite.xlsx",
        "Lifespan regulators/Lifespan regulators.xlsx",
    ]
    latest_date = None
    try:
        for path in paths:
            resp = requests.get(
                "https://api.github.com/repos/vikkihuangkexin/AgeAnnoMO/commits",
                params={"path": path, "per_page": 1},
                timeout=30,
            )
            resp.raise_for_status()
            data = resp.json()
            if data:
                commit_date = data[0]["commit"]["committer"]["date"]  # e.g. 2023-08-16T05:35:36Z
                date_str = commit_date.split("T")[0].replace("-", "")
                if latest_date is None or date_str > latest_date:
                    latest_date = date_str
    except Exception:
        latest_date = None

    if latest_date:
        return latest_date

    try:
        resp = requests.get(
            "https://api.github.com/repos/vikkihuangkexin/AgeAnnoMO/releases",
            timeout=30,
        )
        resp.raise_for_status()
        releases = resp.json()
        if releases:
            tag = releases[0].get("tag_name", "v1.0")
            published = releases[0].get("published_at", "")
            date_str = published.split("T")[0].replace("-", "") if published else ""
            return f"{tag}_{date_str}" if date_str else tag
    except Exception:
        pass

    return "v1.0_20230930"
