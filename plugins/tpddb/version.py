"""
version.py — TPDdb release string

Fetches the download page at tpddb.idrblab.net/download and extracts
the "Last updated: YYYY-MM-DD" date displayed on the page.

Fallback: HTTP HEAD Last-Modified on the PROTAC main table file → YYYYMMDD.
"""


def get_release(self):
    import re
    import requests
    from datetime import datetime

    # Strategy 1: parse "Last updated: YYYY-MM-DD" from the download page HTML
    try:
        resp = requests.get(
            "https://tpddb.idrblab.net/download",
            timeout=30,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        if resp.status_code == 200:
            # Look for patterns like "Last updated: 2025-08-31" or "Updated: Aug 31, 2025"
            match = re.search(
                r"(?:last\s+updated|updated)[:\s]+(\d{4}-\d{2}-\d{2})",
                resp.text,
                re.IGNORECASE,
            )
            if match:
                date_str = match.group(1).replace("-", "")
                return date_str  # e.g. "20250831"

            # Alternative format: "31 August 2025" or "August 31, 2025"
            match2 = re.search(
                r"(\d{1,2})\s+(January|February|March|April|May|June|July|August|September|October|November|December)\s+(\d{4})",
                resp.text,
                re.IGNORECASE,
            )
            if match2:
                try:
                    d = datetime.strptime(
                        f"{match2.group(1)} {match2.group(2)} {match2.group(3)}",
                        "%d %B %Y",
                    )
                    return d.strftime("%Y%m%d")
                except ValueError:
                    pass
    except Exception:
        pass

    # Strategy 2: HTTP HEAD Last-Modified on PROTAC main table
    try:
        head = requests.head(
            "https://tpddb.idrblab.net/sites/files/tpd_download/PROTAC_main_table.txt",
            timeout=30,
            headers={"User-Agent": "Mozilla/5.0"},
            allow_redirects=True,
        )
        lm = head.headers.get("Last-Modified")
        if lm:
            from email.utils import parsedate
            t = parsedate(lm)
            if t:
                return datetime(*t[:6]).strftime("%Y%m%d")
    except Exception:
        pass

    return None
