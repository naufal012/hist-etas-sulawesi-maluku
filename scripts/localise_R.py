#!/usr/bin/env python
"""
Rewrite the hardcoded Japan map origin in the package's R plotting scripts.

delone2.out stores coordinates RELATIVE to (xmin, ymin), and every plotting
script adds the origin back before drawing. All five scripts hardcode Japan's
origin (128, 30), and delo2d-bvalues.R additionally hardcodes Japan's 21x17
degree extent. Run unpatched on another region they produce a plot with the
right values at entirely the wrong coordinates -- which looks completely
plausible unless you read the axis labels.

Usage:  python localise_R.py <script.R> <region.conf>   (edits in place)
"""
import re, sys, os


def load(path):
    d = {}
    for line in open(path, encoding="utf-8"):
        line = line.split("#")[0].strip()
        if line and "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def main():
    rfile, region = sys.argv[1], sys.argv[2]
    c = load(region)
    xmin, ymin = float(c["XMIN"]), float(c["YMIN"])
    tx, ty = float(c["TX"]), float(c["TY"])

    src = open(rfile, encoding="utf-8").read()
    orig = src
    n = 0

    def sub(pat, rep, s):
        nonlocal n
        s2, k = re.subn(pat, rep, s)
        n += k
        return s2

    # lon0 <- 128.   /  lon0<-128  /  lon0<- 128.
    src = sub(r"(lon0\s*<-\s*)[-\d.]+", rf"\g<1>{xmin}", src)
    src = sub(r"(lat0\s*<-\s*)[-\d.]+", rf"\g<1>{ymin}", src)

    # lon2 = lon0+21. ; lat2 = lat0+17.   (Japan's extent, in delo2d-bvalues.R)
    src = sub(r"(lon2\s*=\s*lon0\s*\+\s*)[-\d.]+", rf"\g<1>{tx}", src)
    src = sub(r"(lat2\s*=\s*lat0\s*\+\s*)[-\d.]+", rf"\g<1>{ty}", src)

    # dataPoint[,2]=dataPoint[,2]+128.   (delone-plot.R)
    src = sub(r"(dataPoint\[\s*,\s*2\s*\]\s*=\s*dataPoint\[\s*,\s*2\s*\]\s*\+\s*)[-\d.]+",
              rf"\g<1>{xmin}", src)
    src = sub(r"(dataPoint\[\s*,\s*3\s*\]\s*=\s*dataPoint\[\s*,\s*3\s*\]\s*\+\s*)[-\d.]+",
              rf"\g<1>{ymin}", src)

    if src != orig:
        open(rfile, "w", encoding="utf-8", newline="\n").write(src)
    print(f"   localised {os.path.basename(rfile)}: {n} coordinate(s) "
          f"-> origin ({xmin}, {ymin}), extent {tx}x{ty}")


if __name__ == "__main__":
    main()
