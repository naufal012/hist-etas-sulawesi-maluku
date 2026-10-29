"""Set up the smoothing-sensitivity scan (stage 15) for hist-etas5pa-fixw.

Each scan point s gets <run>/s15/s<tag>/ with delone2.out, the fixed-weight
program, a conf in the package's init=1 form, and start.upda: the ABIC
solution's coefficients with its five weights multiplied by s.  The program
then performs one inner solve at those weights (src/hist-etas5pa-fixw.f).

    python scripts/setup_smoothing_scan.py --run run-mc50-ss
    python scripts/setup_smoothing_scan.py --run run-mc50-ss --tags 1_32 --from-tag 1_16

--from-tag warm-starts from another scan point's solution instead of the ABIC
one (used for s = 1/32, whose solve from the ABIC solution diverged).  Then, in
each directory:  ./hist-etas5pa-fixw.exe > hist-etas5pa.prt
"""
from __future__ import annotations

import argparse
import os
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCALES = {"4": 4, "2": 2, "1": 1, "1_2": 1 / 2, "1_4": 1 / 4, "1_8": 1 / 8,
          "1_16": 1 / 16, "1_32": 1 / 32, "1_64": 1 / 64}


def read_upda(path):
    lines = open(path).read().split("\n")
    w = [float(v) for v in lines[0].split()] + [float(lines[1].split()[0])]
    t = lines[1].split()
    return w, float(t[1]), int(t[2]), float(t[3]), lines[2:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="run-mc50-ss")
    ap.add_argument("--tags", nargs="*", default=list(SCALES))
    ap.add_argument("--from-tag", help="warm-start from s15/s<tag> instead of s10")
    a = ap.parse_args()
    run = os.path.join(ROOT, a.run)
    w_hat, abic, m, w01, _ = read_upda(os.path.join(run, "s10", "hist-etas5pa.upda"))
    src = (os.path.join(run, "s15", "s" + a.from_tag, "hist-etas5pa.upda") if a.from_tag
           else os.path.join(run, "s10", "hist-etas5pa.upda"))
    _, _, _, _, body = read_upda(src)
    # catalogue, region and time lines of the fitted run's own conf
    with open(os.path.join(run, "s10", "hist-etas5pa.conf")) as fh:
        head = "".join(fh.readlines()[:3])
    for tag in a.tags:
        s = SCALES[tag]
        d = os.path.join(run, "s15", "s" + tag)
        os.makedirs(d, exist_ok=True)
        shutil.copy(os.path.join(run, "s10", "delone2.out"), d)
        shutil.copy(os.path.join(ROOT, "bin", "hist-etas5pa-fixw.exe"), d)
        ws = [s * v for v in w_hat]
        with open(os.path.join(d, "start.upda"), "w") as fh:
            fh.write("".join(f"{v:20.12e}" for v in ws[:4]) + "\n")
            fh.write(f"{ws[4]:20.12e}{abic:20.12e}{m:10d}{'':10s}{w01:20.12e}\n")
            fh.write("\n".join(body))
        with open(os.path.join(d, "hist-etas5pa.conf"), "w") as fh:
            fh.write(head + "1\n0\n1\nstart.upda\n"
                     f"0. {w01:.1f} " + " ".join(f"{v:.6e}" for v in ws) + "\n"
                     "0               !ipr\n0.1d-3  0.1d-3  !tau1,tau2(davidn)\n"
                     "0.1d-3  0.1d-3  !eps1,eps2(davidn)\n0 1.d0 0.5d-0    !nhesapp,dist,eps(simplex)\n")
        print("set up", d, "s =", s)


if __name__ == "__main__":
    main()
