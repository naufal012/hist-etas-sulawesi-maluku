"""Conditional-intensity snapshots from a fit that ends before the instants.

The fitted run (--fit) supplies the tessellation and the frozen fields.  The
full catalogue (--full) supplies the events observed after the fit ended;
they enter the history as they occur, with parameters interpolated on the
fitted tessellation exactly as in oos_eval.py.  Nothing is refitted.  This is
the package's forecast mode (intensityday2magnod5pa.f with a work.para
catalogue) computed in Python, and writes img1-<tag>.data in that program's
format so intensity_plot.R draws it unchanged.

Like the Fortran program, lambda sums over every earlier event (no search
radius) and mu is linear in the value on the triangles.

Validation: with --fit and --full the same run, the grids must reproduce the
Fortran's img1-*.data in <run>/s13.

    python scripts/prospective_fig4.py --fit run-mc50 --full run-mc50 --validate
    python scripts/prospective_fig4.py --fit run-mc50-pre78 --full run-mc50 \
        --cut "2026-06-07 12:00"
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import etas_lik as E  # noqa: E402
import oos_eval as O  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402

RSCRIPT = r"C:/Program Files/R/R-4.5.1/bin/Rscript.exe"


def read_snapshots(path):
    out = []
    with open(path) as fh:
        for line in fh:
            tag, day, label = line.rstrip("\n").split("\t")
            out.append((tag, float(day), label))
    return out


def intensity_grid(cat, fe, mu_at, gx, gy, t, chunk=200):
    """log10 lambda(t, gx, gy) given every event of `cat` with t_j < t."""
    lam = mu_at(gx, gy)
    idx = np.nonzero(cat.t < t)[0]
    for s in range(0, len(idx), chunk):
        j = idx[s:s + chunk]
        sc = np.exp(-fe.a[j] * cat.mag[j])[:, None]
        u = (gx[None, :] - cat.x[j, None]) * sc
        v = (gy[None, :] - cat.y[j, None]) * sc
        Q = cat.qa[j, None] * u ** 2 - 2 * cat.qb[j, None] * u * v + cat.qc[j, None] * v ** 2
        lam = lam + np.sum(fe.K[j, None] / (t - cat.t[j, None] + fe.c) ** fe.p[j, None]
                           * (Q + fe.d) ** (-fe.q[j, None]), axis=0)
    return np.log10(lam)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--full", default="run-mc50")
    ap.add_argument("--cut", default=None, help="UTC datetime the fitted catalogue ends")
    ap.add_argument("--hypo", default="catalog/hypo.ts")
    ap.add_argument("--snap", action="append", default=[],
                    help="DAY:TAG:LABEL; default the full run's s13/snapshots.txt")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--figs", default="figs-mc50")
    a = ap.parse_args()
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    fit_run, full_run = os.path.join(root, a.fit), os.path.join(root, a.full)

    cfit = E.read_conf(os.path.join(fit_run, "s10", "hist-etas5pa.conf"))
    cfull = E.read_conf(os.path.join(full_run, "s10", "hist-etas5pa.conf"))
    fit = E.read_delone2(os.path.join(fit_run, "s10", "delone2.out"), cfit)
    full = E.read_delone2(os.path.join(full_run, "s10", "delone2.out"), cfull)
    f5 = E.read_5pa(os.path.join(fit_run, "s10", "hist-etas5pa.upda"))

    if a.snap:
        snaps = [(s.split(":", 2)[1], float(s.split(":")[0]), s.split(":", 2)[2]) for s in a.snap]
    else:
        snaps = read_snapshots(os.path.join(full_run, "s13", "snapshots.txt"))
    T2 = max(d for _, d, _ in snaps) + 1e-6

    if a.cut:
        origin = O.catalogue_origin(os.path.join(root, a.hypo))
        T1 = (dt.datetime.fromisoformat(a.cut) - origin).total_seconds() / 86400.0
        cat = O.combined_catalog(fit, full, T1, max(T2, cfull["tz"]))
        fe = O.event_fields(fit, f5, cat)
        print(f"cut {a.cut} = day {T1:.4f}; fitted events {fit.n}, "
              f"later events available {cat.n - fit.n}")
    else:
        cat, fe = fit, f5.at_events(fit.n)

    tri = mtri.Triangulation(fit.node_x, fit.node_y, fit.tri)
    mu_interp = mtri.LinearTriInterpolator(tri, f5.mu)

    def mu_at(x, y):
        return np.asarray(mu_interp(x, y).filled(np.nan))

    node = np.loadtxt(os.path.join(full_run, "s13", "node.dat"))
    lon, lat = node[:, 0], node[:, 1]
    gx = (lon - fit.xmin) * fit.cosfac
    gy = lat - fit.ymin
    pix = (14.0 / 210) * fit.cosfac * (12.0 / 180)   # node.conf: 210 x 180 over the box

    out = os.path.join(fit_run, "s13p" if not a.validate else "s13v")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "snapshots.txt"), "w") as fh:
        for tag, day, label in snaps:
            z = intensity_grid(cat, fe, mu_at, gx, gy, day)
            np.savetxt(os.path.join(out, f"img1-{tag}.data"),
                       np.column_stack([np.full(len(z), day - cfit["tz"]), lon, lat, z]),
                       fmt=["%11.5f", "%9.4f", "%9.4f", "%10.5f"], delimiter="")
            fh.write(f"{tag}\t{day}\t{label}\n")
            total = np.sum(10 ** z) * pix
            k = np.argmax(z)
            msg = (f"{tag:4s} day {day:9.3f}  region total {total:7.3f} /d  "
                   f"peak {10 ** z[k]:8.3f} at ({lon[k]:.2f}E, {lat[k]:.2f}N)")
            if a.validate:
                ref = np.loadtxt(os.path.join(full_run, "s13", f"img1-{tag}.data"))
                msg += f"   max |d log10| vs Fortran {np.max(np.abs(ref[:, 3] - z)):.2e}"
            print(msg)
    if a.validate:
        return

    # draw with the house script: full-catalogue epicentres, fitted run's origin
    s13 = os.path.join(fit_run, "s13")
    shutil.copy(os.path.join(full_run, "s10", "delone2.out"), out)
    shutil.copy(os.path.join(s13, "hist-etas5pa-forecast.conf"), out)
    for f in ("drawmap.R", "mapkit.R", "intensity_plot.R"):
        shutil.copy(os.path.join(root, "scripts", f), out)
    os.makedirs(os.path.join(out, "maps"), exist_ok=True)
    for f in ("w1.pp.data", "w1.br.data"):
        shutil.copy(os.path.join(root, "maps", f), os.path.join(out, "maps"))
    r = subprocess.run([RSCRIPT, "intensity_plot.R"], cwd=out, capture_output=True, text=True)
    print(r.stdout.strip() or r.stderr.strip())
    if r.returncode == 0:
        dst = os.path.join(root, a.figs, "fig6-intensity-prospective.pdf")
        shutil.copy(os.path.join(out, "intensity.pdf"), dst)
        print("wrote", dst)


if __name__ == "__main__":
    main()
