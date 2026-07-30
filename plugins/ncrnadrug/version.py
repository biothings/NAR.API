def get_release(self):
    """Return a release string for ncRNADrug.

    Primary strategy: HTTP HEAD on the DR_Curated.txt bulk download file and use
    its Last-Modified header (YYYYMMDD) — this reflects when the underlying data
    tables were actually regenerated, which is more precise than the homepage's
    "Latest updates" changelog (last entry: 2023-10-06, stale relative to the
    file's actual Last-Modified).

    Fallback: scrape the homepage's "Latest updates" list for the most recent
    YYYY-MM-DD entry.
    """
    import re
    import email.utils

    import requests

    try:
        resp = requests.head(
            "http://www.jianglab.cn/ncRNADrug/data_download/DR_Curated.txt",
            timeout=30,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        resp.raise_for_status()
        lm = resp.headers.get("Last-Modified", "")
        if lm:
            dt = email.utils.parsedate_to_datetime(lm)
            return dt.strftime("%Y%m%d")
    except Exception:
        pass

    try:
        resp = requests.get(
            "http://www.jianglab.cn/ncRNADrug/",
            timeout=30,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        resp.raise_for_status()
        dates = re.findall(r"(\d{4}-\d{2}-\d{2})", resp.text)
        if dates:
            dates.sort(reverse=True)
            return dates[0].replace("-", "")
    except Exception:
        pass

    return None
