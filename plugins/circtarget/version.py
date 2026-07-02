def get_release(self):
    import requests
    from email.utils import parsedate_to_datetime

    url = "https://circtarget.cn/static/download/all.rar"
    try:
        resp = requests.head(url, timeout=30, allow_redirects=True)
        resp.raise_for_status()
        last_modified = resp.headers.get("Last-Modified")
        if last_modified:
            dt = parsedate_to_datetime(last_modified)
            return dt.strftime("%Y%m%d")
    except Exception:
        pass
    return None
