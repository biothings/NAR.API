def get_release(self):
    import re
    import requests

    try:
        resp = requests.get(
            "https://ecbd.eu/download/",
            timeout=30,
            verify=False,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        resp.raise_for_status()
        # Look for a date pattern like "Updated: 2026-05-01" or a version string
        m = re.search(r'(?:updated?|version|release)[:\s]+(\d{4}[-/]\d{2}[-/]\d{2})', resp.text, re.IGNORECASE)
        if m:
            return re.sub(r'[-/]', '', m.group(1))
    except Exception:
        pass

    # Fallback: Last-Modified header on bioactives.csv
    try:
        r2 = requests.head(
            "https://ecbd.eu/static/core/compounds/bioactives.csv",
            timeout=30,
            verify=False,
        )
        lm = r2.headers.get("Last-Modified", "")
        if lm:
            import email.utils
            t = email.utils.parsedate_to_datetime(lm)
            return t.strftime("%Y%m%d")
    except Exception:
        pass

    return None
