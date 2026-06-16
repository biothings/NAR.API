def get_release(self):
    import re
    import requests

    # getLatestRelease.php returns Content-Disposition: filename="Apr2026_release.txt"
    resp = requests.head(
        "https://signor.uniroma2.it/releases/getLatestRelease.php",
        timeout=30,
        allow_redirects=True,
    )
    resp.raise_for_status()
    cd = resp.headers.get("Content-Disposition", "")
    m = re.search(r'filename="(\w{3})(\d{4})_release\.txt"', cd)
    if m:
        months = {
            "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
            "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
            "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12",
        }
        return f"{m.group(2)}{months.get(m.group(1), '00')}"

    # Fallback: scrape downloads page for most recent release link
    resp2 = requests.get("https://signor.uniroma2.it/downloads.php", timeout=30)
    m2 = re.search(r'/releases/(\w{3})(\d{4})_release\.txt', resp2.text)
    if m2:
        months = {
            "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
            "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
            "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12",
        }
        return f"{m2.group(2)}{months.get(m2.group(1), '00')}"

    return None
