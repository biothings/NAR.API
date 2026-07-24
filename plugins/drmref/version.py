def get_release(self):
    import re

    import requests

    # 1. Prefer the Last-Modified header on the primary bulk file — this is
    #    the most reliable release signal for a static flat-file site.
    try:
        resp = requests.head(
            "https://ccsm.uth.edu/DRMref/table_summary/gene_summary.txt",
            timeout=30,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        lm = resp.headers.get("Last-Modified", "")
        if lm:
            import email.utils

            dt = email.utils.parsedate_to_datetime(lm)
            if dt:
                return dt.strftime("%Y%m%d")
    except Exception:
        pass

    # 2. Fall back to scraping the homepage for a visible date/version string.
    try:
        resp = requests.get(
            "https://ccsm.uth.edu/DRMref/",
            timeout=30,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        resp.raise_for_status()
        m = re.search(
            r"(?:updated?|version|release)[:\s]+(\d{4}[-/]\d{2}[-/]\d{2})",
            resp.text,
            re.IGNORECASE,
        )
        if m:
            return re.sub(r"[-/]", "", m.group(1))
    except Exception:
        pass

    # 3. Last resort: the most recent known data-refresh date observed
    #    during Phase-1 site inspection (site footer / file timestamps as
    #    of the 2026-07-09 inspection pass).
    return "20231001"
