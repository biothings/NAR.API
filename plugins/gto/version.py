def get_release(self):
    """Return the GTO database release string.

    Fetches the GTO download page and extracts any visible date/version marker.
    Falls back to a Last-Modified header check on the clinical data file.
    """
    import requests
    import re
    from datetime import datetime

    # Try to get Last-Modified from the clinical data download
    try:
        resp = requests.head(
            "http://www.inbirg.com/gto/download/clinical_data_download",
            timeout=30,
            allow_redirects=True,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        last_modified = resp.headers.get("Last-Modified")
        if last_modified:
            # Parse HTTP date: "Mon, 01 Jan 2024 00:00:00 GMT"
            dt = datetime.strptime(last_modified, "%a, %d %b %Y %H:%M:%S %Z")
            return dt.strftime("%Y%m%d")
    except Exception:
        pass

    # Fallback: try to find a version string on the GTO homepage
    try:
        resp = requests.get(
            "http://www.inbirg.com/gto/home",
            timeout=30,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        text = resp.text
        # Look for year patterns like "2024" or "2025" near version/update keywords
        match = re.search(r"(?:version|update|release|last\s+updated?)[^<]*(\d{4})", text, re.IGNORECASE)
        if match:
            return match.group(1)
    except Exception:
        pass

    # Last resort: return the paper publication year as a static fallback
    return "20250106"
