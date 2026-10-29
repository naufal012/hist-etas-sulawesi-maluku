#!/usr/bin/env python
"""
Convert an FDSN event CSV (USGS / ISC / BMKG export) into the 10-column
whitespace-separated format that HIST-PPM's tsei2etas reads:

    year month day hour min sec lon lat depth mag

tsei2etas uses a list-directed read, so only field ORDER matters, not columns.

Usage:
    python make_hypo.py <in.csv> <out.ts> [--region-conf region.conf]

Filters applied:
  * type == "earthquake"   (drops explosions, quarry blasts, landslides)
  * rows with missing time / lat / lon / mag are dropped and counted
  * events outside the region box in region.conf are dropped (if conf given)
Events are sorted by origin time ascending -- tsei2etas takes the FIRST record
as its time origin, so an unsorted file silently produces negative times.
"""
import csv, sys, os, datetime, argparse


def read_region(path):
    if not path or not os.path.exists(path):
        return None
    r = {}
    for line in open(path, encoding="utf-8"):
        line = line.split("#")[0].strip()
        if not line or "=" not in line:
            continue
        k, v = line.split("=", 1)
        r[k.strip()] = v.strip()
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--region-conf", default=None)
    a = ap.parse_args()

    reg = read_region(a.region_conf)
    if reg:
        xmin = float(reg["XMIN"]); ymin = float(reg["YMIN"])
        xmax = xmin + float(reg["TX"]); ymax = ymin + float(reg["TY"])
    else:
        xmin = ymin = -1e9; xmax = ymax = 1e9

    rows, drop_type, drop_bad, drop_box = [], 0, 0, 0
    with open(a.infile, newline="", encoding="utf-8") as fh:
        for rec in csv.DictReader(fh):
            if (rec.get("type") or "earthquake").strip() != "earthquake":
                drop_type += 1
                continue
            try:
                t = rec["time"].strip().replace("Z", "")
                # FDSN times are ISO8601, with or without fractional seconds
                dt = datetime.datetime.fromisoformat(t)
                lon = float(rec["longitude"]); lat = float(rec["latitude"])
                mag = float(rec["mag"])
                dep = float(rec["depth"]) if rec.get("depth") not in (None, "") else 0.0
            except (ValueError, KeyError, TypeError):
                drop_bad += 1
                continue
            if not (xmin <= lon <= xmax and ymin <= lat <= ymax):
                drop_box += 1
                continue
            sec = dt.second + dt.microsecond / 1e6
            rows.append((dt, dt.year, dt.month, dt.day, dt.hour, dt.minute,
                         sec, lon, lat, dep, mag))

    rows.sort(key=lambda r: r[0])

    with open(a.outfile, "w", encoding="utf-8", newline="\n") as out:
        for (_, y, mo, d, h, mi, s, lon, lat, dep, mag) in rows:
            out.write(f"{y:4d} {mo:02d} {d:02d} {h:02d} {mi:02d} {s:6.2f} "
                      f"{lon:9.4f} {lat:9.4f} {dep:7.2f} {mag:4.1f}\n")

    print(f"wrote {a.outfile}: {len(rows)} events")
    if rows:
        print(f"  time  : {rows[0][0]}  ..  {rows[-1][0]}")
        print(f"  lon   : {min(r[7] for r in rows):.3f} .. {max(r[7] for r in rows):.3f}")
        print(f"  lat   : {min(r[8] for r in rows):.3f} .. {max(r[8] for r in rows):.3f}")
        print(f"  mag   : {min(r[10] for r in rows):.1f} .. {max(r[10] for r in rows):.1f}")
        span = (rows[-1][0] - rows[0][0]).total_seconds() / 86400.0
        print(f"  span  : {span:.2f} days")
    print(f"  dropped: {drop_type} non-earthquake, {drop_bad} malformed, "
          f"{drop_box} outside region box")


if __name__ == "__main__":
    main()
