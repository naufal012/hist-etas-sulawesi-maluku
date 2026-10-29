#!/usr/bin/env python
"""
Fetch an FDSN event extract for the region in region.conf.

The USGS FDSN endpoint refuses any query matching more than 20000 events, so
requests are issued in yearly chunks and concatenated. Chunking by time (not
by magnitude) keeps each slice a complete, unbiased sample of its interval,
which matters because this catalogue feeds a completeness analysis.

Usage:
    python fetch_catalog.py <region.conf> <out.csv> [--minmag M] [--chunk-years N]
"""
import urllib.request, urllib.parse, argparse, datetime, sys, time


def load(path):
    d = {}
    for line in open(path, encoding="utf-8"):
        line = line.split("#")[0].strip()
        if line and "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def fetch(params, retries=3):
    url = "https://earthquake.usgs.gov/fdsnws/event/1/query?" + urllib.parse.urlencode(params)
    for a in range(retries):
        try:
            return urllib.request.urlopen(url, timeout=300).read().decode()
        except Exception as e:
            if a == retries - 1:
                raise
            print(f"    retry {a+1} after {e}", file=sys.stderr)
            time.sleep(5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("region"); ap.add_argument("outfile")
    ap.add_argument("--minmag", type=float, default=None)
    ap.add_argument("--chunk-years", type=int, default=1)
    a = ap.parse_args()

    c = load(a.region)
    xmin, ymin = float(c["XMIN"]), float(c["YMIN"])
    xmax, ymax = xmin + float(c["TX"]), ymin + float(c["TY"])
    minmag = a.minmag if a.minmag is not None else float(c.get("FETCH_MINMAG", 4.5))
    t0 = datetime.date.fromisoformat(c["CATALOG_START"])
    t1 = datetime.date.fromisoformat(c["CATALOG_END"])

    header, rows = None, []
    s = t0
    while s < t1:
        e = min(datetime.date(s.year + a.chunk_years, 1, 1), t1)
        p = dict(format="csv", starttime=s.isoformat(), endtime=e.isoformat(),
                 minlatitude=ymin, maxlatitude=ymax,
                 minlongitude=xmin, maxlongitude=xmax,
                 minmagnitude=minmag, orderby="time-asc")
        txt = fetch(p).strip().splitlines()
        if txt:
            if header is None:
                header = txt[0]
            body = txt[1:] if txt[0].startswith("time,") else txt
            rows.extend(body)
            print(f"  {s} .. {e}: {len(body)} events")
        s = e

    with open(a.outfile, "w", encoding="utf-8", newline="\n") as f:
        f.write(header + "\n")
        f.write("\n".join(rows) + "\n")
    print(f"wrote {a.outfile}: {len(rows)} events, M>={minmag}")


if __name__ == "__main__":
    main()
