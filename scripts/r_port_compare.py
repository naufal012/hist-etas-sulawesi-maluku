"""Numbers of the paper that no pipeline script writes, computed the same way
for the Fortran run (run-mc50) and the R-port run (run-mc50-R), so the
comparison draft can quote both. Run after r_port_post.R, residuals.py and
tectonic_context.py.

    python scripts/r_port_compare.py   ->  run-mc50-R/compare.json
"""
import json
import os

import numpy as np

import etas_lik as E

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIX = (14 / 210) * (12 / 180)


def one(run):
    r = os.path.join(ROOT, run)
    out = {}
    f5 = E.read_5pa(os.path.join(r, "s10", "hist-etas5pa.upda"))
    n = 2480
    fields = {"mu": f5.mu[:n], "K0": f5.K[:n], "alpha": 2 * f5.a[:n],
              "p": f5.p[:n], "q": f5.q[:n]}
    out["vertices"] = {k: dict(min=float(v.min()), median=float(np.median(v)),
                               max=float(v.max()), ratio=float(v.max() / v.min()))
                       for k, v in fields.items()}
    # conditional intensity snapshots
    snaps = {}
    for tag in ("pre", "d1", "end"):
        d = np.loadtxt(os.path.join(r, "s13", f"img1-{tag}.data"))
        lam = 10 ** d[:, 3]
        i = int(lam.argmax())
        snaps[tag] = dict(rate=float(lam.sum() * PIX), peak=float(lam[i]),
                          at=[float(d[i, 1]), float(d[i, 2])])
    out["intensity"] = snaps
    # fields on the grid against tectonic classes
    px = np.loadtxt(os.path.join(r, "s13", "interp5pa.pixel"))
    cv = np.loadtxt(os.path.join(r, "tect", "covariates.pixel"), skiprows=1)
    inh = (px[:, 2] > -998) & (cv[:, 2] > 0)
    p = px[:, 5]
    cut = np.quantile(p[inh], 0.05)
    low = inh & (p < cut)
    ev = np.genfromtxt(os.path.join(r, "tect", "events-class.txt"), names=True, dtype=None,
                       encoding=None)
    # event -> pixel
    ix = np.clip(((ev["lon"] - 118) / (14 / 210)).astype(int), 0, 209)
    iy = np.clip(((ev["lat"] + 6) / (12 / 180)).astype(int), 0, 179)
    pid = ix * 180 + iy
    in_low = low[pid]
    intra = ev["class"] == "intraslab"
    mak = inh & (px[:, 0] < 119.5) & (px[:, 1] > -4) & (px[:, 1] < 0)
    out["p_low5"] = dict(cut=float(cut), n_events=int(in_low.sum()),
                         intraslab_share=float(intra[in_low].mean()),
                         intraslab_share_all=float(intra.mean()),
                         makassar_median_p=float(np.median(p[mak])))
    # K0 by event class (field at the event's pixel)
    K = px[pid, 3]
    out["K0_median_by_class"] = {c: float(np.median(K[ev["class"] == c]))
                                 for c in np.unique(ev["class"])}
    tect = json.load(open(os.path.join(r, "tect", "tectonic-summary.json")))
    out["spearman"] = tect["spearman_in_hull"]
    out["K0_top10_none"] = tect["top10_pixel_slab_context"]["K0"].get("none")
    res = os.path.join(r, "s14", "residuals.json")
    if os.path.exists(res):
        out["residuals"] = json.load(open(res))
    return out


if __name__ == "__main__":
    both = {"fortran": one("run-mc50"), "r": one("run-mc50-R")}
    with open(os.path.join(ROOT, "run-mc50-R", "compare.json"), "w") as fh:
        json.dump(both, fh, indent=1)
    for k in ("fortran", "r"):
        o = both[k]
        print("==", k)
        for f, v in o["vertices"].items():
            print(f"  {f:6s} {v['min']:.3g} {v['median']:.4g} {v['max']:.3g} ratio {v['ratio']:.4g}")
        for t, v in o["intensity"].items():
            print(f"  {t:4s} rate {v['rate']:.3f} peak {v['peak']:.3g} at {v['at']}")
        print("  p low5", {a: round(b, 3) for a, b in o["p_low5"].items()})
        print("  K0 med by class", {a: f"{b:.2g}" for a, b in o["K0_median_by_class"].items()})
        print("  K0 top10 none", o["K0_top10_none"])
        print("  rho mu", o["spearman"]["log10_mu"], "\n  rho K0", o["spearman"]["log10_K0"])
