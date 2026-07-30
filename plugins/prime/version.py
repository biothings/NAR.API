def get_release(self):
    import requests
    r = requests.get("https://zenodo.org/api/records/15711237", timeout=30)
    r.raise_for_status()
    data = r.json()
    pub_date = data.get("metadata", {}).get("publication_date", "")
    if pub_date:
        return pub_date.replace("-", "")
    # Fallback: use last-modified header on the metadata file
    r2 = requests.head(
        "https://zenodo.org/records/15711237/files/samples_metadata.csv?download=1",
        timeout=30,
        allow_redirects=True,
    )
    lm = r2.headers.get("Last-Modified", "")
    if lm:
        from email.utils import parsedate
        from datetime import datetime
        parsed = parsedate(lm)
        if parsed:
            return datetime(*parsed[:6]).strftime("%Y%m%d")
    return None
