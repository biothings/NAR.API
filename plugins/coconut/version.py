def get_release(self):
    import re
    import requests
    # Parse version from the canonical COCONUT download page
    resp = requests.get(
        "https://coconut.naturalproducts.net/download",
        timeout=30,
        headers={"User-Agent": "Mozilla/5.0"}
    )
    resp.raise_for_status()
    # Page contains "Version: May 2026" — extract and convert to YYYYMM
    m = re.search(r"Version:\s*(\w+)\s+(\d{4})", resp.text)
    if m:
        month_str, year = m.group(1), m.group(2)
        months = {
            "january": "01", "february": "02", "march": "03", "april": "04",
            "may": "05", "june": "06", "july": "07", "august": "08",
            "september": "09", "october": "10", "november": "11", "december": "12"
        }
        mm = months.get(month_str.lower(), "00")
        return f"{year}{mm}"  # e.g. "202605"
    # Fallback: extract date from CSV download filenames on the page
    dates = re.findall(r'coconut_csv[^"]*-(\d{2})-(\d{4})\.zip', resp.text)
    if dates:
        mm, yyyy = dates[-1]
        return f"{yyyy}{mm}"
    return None
