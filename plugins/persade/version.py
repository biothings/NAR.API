def get_release(self):
    """Return PersADE release string from Last-Modified header on the drug file."""
    import requests
    import re

    try:
        r = requests.head(
            "https://persade.idrblab.net/download/mtADENet/Global_Drug_20230627.xlsx",
            headers={"User-Agent": "Mozilla/5.0"},
            timeout=30,
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

    # Fallback: extract date from filename pattern in download URL
    try:
        resp = requests.get(
            "https://persade.idrblab.net/resource.php",
            headers={"User-Agent": "Mozilla/5.0"},
            timeout=30,
        )
        match = re.search(r"Global_Drug_(\d{8})\.xlsx", resp.text)
        if match:
            return match.group(1)
    except Exception:
        pass

    return "20230627"  # Date in available download filenames
