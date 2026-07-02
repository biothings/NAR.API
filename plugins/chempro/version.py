def get_release(self):
    """Return Chem(Pro)2 release string from the download page last-updated date."""
    import requests
    import re

    try:
        url = "https://chemprosquare.idrblab.net/download"
        resp = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=30, verify=False)
        resp.raise_for_status()
        # Look for last updated date in the page
        match = re.search(r"Last\s+[Uu]pdated[:\s]+([A-Za-z]+\s+\d{1,2},?\s+\d{4}|\d{4}-\d{2}-\d{2})", resp.text)
        if match:
            date_str = match.group(1)
            return date_str.replace(" ", "").replace(",", "")
    except Exception:
        pass

    # Fallback: Last-Modified header on general_probe.txt
    try:
        r = requests.head(
            "https://chemprosquare.idrblab.net/sites/default/files/download/general_probe.txt",
            headers={"User-Agent": "Mozilla/5.0"},
            timeout=30,
            verify=False,
        )
        lm = r.headers.get("Last-Modified", "")
        if lm:
            from email.utils import parsedate
            from datetime import datetime
            parsed = parsedate(lm)
            if parsed:
                dt = datetime(*parsed[:6])
                return dt.strftime("%Y%m%d")
    except Exception:
        pass

    return "20240923"  # Last known update date
