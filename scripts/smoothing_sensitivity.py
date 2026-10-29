"""Smoothing-sensitivity diagnostic for the hist-etas5pa fit (stage 15).

ABIC chooses the five smoothness weights w.  This scales all five together by
s and refits the coefficients at the fixed weights s*w_hat (one inner solve,
src/hist-etas5pa-fixw.f, warm-started from the ABIC optimum), then asks the
residuals of residuals.py what the choice costs:

  * ABIC and the exact log-likelihood at each s;
  * temporal residuals -- KS of tau_i = Lambda(t_i) against uniform, and the
    largest |N(t) - Lambda(t)|;
  * Voronoi residuals -- KS against Gamma(3.569, 3.569), and the Spearman
    correlation between the normal score z_i and the cell area, the
    cluster-scale bias seen at the ABIC choice (negative: the model expects
    too many events in large sparse cells and too few in small clustered ones).

    python scripts/smoothing_sensitivity.py --run run-mc50-ss
reads <run>/s15/s*/hist-etas5pa.upda and .prt, writes <run>/s15/summary.json
and figs-mc50-ss/fig13-smoothing.pdf.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import etas_lik as E  # noqa: E402
from residuals import (COLOURS, INK, MUTED, VORONOI_SHAPE, WIDTH, _style,  # noqa: E402
                       cell_compensator, voronoi_assign)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def scale_of(tag):
    """s4 -> 4, s1_8 -> 1/8."""
    a = tag[1:].split("_")
    return float(a[0]) / (float(a[1]) if len(a) > 1 else 1.0)


def load_catalog(run):
    conf = E.read_conf(os.path.join(run, "s10", "hist-etas5pa.conf"))
    cat = E.read_delone2(os.path.join(run, "s10", "delone2.out"), conf)
    return cat, E.neighbour_pairs(cat, conf["bi2"])


def cell_areas(cat, assign, nv, h=0.01):
    xs = np.arange(h / 2, cat.tx, h)
    ys = np.arange(h / 2, cat.ty, h)
    X, Y = np.meshgrid(xs, ys)
    return np.bincount(assign(X.ravel(), Y.ravel()), minlength=nv) * h * h


def diagnose(args):
    run, path, n_mc = args
    cat, pr = load_catalog(run)
    if path.endswith(".json"):        # the mk and stationary references
        f = E.constant_fields(list(json.load(open(path))["theta"].values()), cat.n)
    elif "hist-etas-mk" in path:
        f = E.read_mk(path)
    else:
        f = E.read_5pa(path)
    mass = E.region_mass(cat, f.at_events(cat.n), nseg=128)
    ll = E.loglik(cat, f, pr, "exact", mass=mass)
    tau = E.compensator_at_events(cat, f, mass)
    lam_T = ll["trig_int"] + ll["bg_int"]
    n = len(tau)
    assign, nv = voronoi_assign(cat)
    r = cell_compensator(cat, f, assign, nv, n_mc=n_mc)
    z = stats.norm.ppf(1 - stats.gamma.cdf(r, VORONOI_SHAPE, scale=1 / VORONOI_SHAPE))
    area = cell_areas(cat, assign, nv)
    q = np.quantile(area, [0.25, 0.75])
    return dict(
        n_target=int(n),
        loglik_exact=float(ll["loglik"]),
        background_share=float(ll["bg_int"] / lam_T),
        ks_uniform_p=float(stats.kstest(tau / lam_T, "uniform").pvalue),
        ks_gaps_p=float(stats.kstest(np.diff(np.concatenate([[0.0], tau])), "expon").pvalue),
        max_abs_N_minus_Lambda=float(np.max(np.abs(np.arange(1, n + 1) - tau))),
        voronoi_sum=float(r.sum()),
        voronoi_ks_p=float(stats.kstest(r, "gamma", args=(VORONOI_SHAPE, 0, 1 / VORONOI_SHAPE)).pvalue),
        spearman_z_area=float(stats.spearmanr(z, area)[0]),
        z_mean_small_cells=float(z[area <= q[0]].mean()),
        z_mean_large_cells=float(z[area >= q[1]].mean()),
    )


def read_abic(prt):
    vals = [float(m.group(1)) for m in re.finditer(r"abic=\s*([-0-9.E+]+)", open(prt).read())]
    return vals[-1] if vals else None


def figure(rows, refs, out):
    s = np.array([r["scale"] for r in rows])
    fig, axes = plt.subplots(1, 3, figsize=(WIDTH, 0.36 * WIDTH))
    c = COLOURS["5pa"]
    base = next(r for r in rows if r["scale"] == 1)
    n = base["n_target"]
    panels = (
        ("abic", "ABIC $-$ ABIC at $\\hat{\\mathbf{w}}$", "(a) ABIC", base["abic"]),
        ("spearman_z_area", "Spearman $\\rho$(score, cell area)", "(b) spatial residuals", 0.0),
        ("max_abs_N_minus_Lambda", "largest $|N(t)-\\Lambda(t)|$", "(c) temporal residuals", 0.0),
    )
    for ax, (key, ylab, title, off) in zip(axes, panels):
        y = np.array([r[key] for r in rows]) - off
        ax.plot(s, y, color=c, lw=1.2, marker="o", ms=3)
        ax.axvline(1, color=MUTED, lw=0.6, ls=":")
        if key == "max_abs_N_minus_Lambda":
            ax.axhline(1.36 * np.sqrt(n), color=MUTED, lw=0.9, ls="--")
            ax.annotate("95 % K–S band", (s.min(), 1.36 * np.sqrt(n)), fontsize=6, color=MUTED,
                        ha="left", va="bottom", xytext=(0, 1.5), textcoords="offset points")
        if key == "spearman_z_area":
            for k, ls in (("mk", "--"), ("st", ":")):
                ax.axhline(refs[k]["spearman_z_area"], color=COLOURS[k], lw=1.0, ls=ls)
                ax.annotate({"mk": "hist-etas-mk", "st": "stationary"}[k],
                            (s.max(), refs[k]["spearman_z_area"]), fontsize=6, color=COLOURS[k],
                            ha="right", va="bottom", xytext=(0, 1.5), textcoords="offset points")
        ax.set_xscale("log", base=2)
        ax.set_xticks(s)
        # label every other power of two so the fractions do not collide
        ax.set_xticklabels([(f"{v:g}" if v >= 1 else f"1/{1 / v:g}")
                            if round(np.log2(v)) % 2 == 0 else "" for v in s], fontsize=6)
        ax.set_xlabel("weight scale $s$", fontsize=7, color=INK)
        ax.set_ylabel(ylab, fontsize=6.5, color=INK)
        ax.set_title(title, fontsize=7.5, color=INK, loc="left")
        _style(ax)
        ax.tick_params(axis="x", labelsize=6)
    fig.tight_layout(pad=0.3, w_pad=0.9)
    fig.savefig(out)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="run-mc50-ss")
    ap.add_argument("--figs", default="figs-mc50-ss")
    ap.add_argument("--n-mc", type=int, default=4000)
    ap.add_argument("--plot-only", action="store_true",
                    help="redraw the figure from s15/summary.json without rescoring")
    a = ap.parse_args()
    run = os.path.join(ROOT, a.run)
    s15 = os.path.join(run, "s15")
    if a.plot_only:
        with open(os.path.join(s15, "summary.json")) as fh:
            sm = json.load(fh)
        n = int(load_catalog(run)[0].target().sum())
        for r in sm["rows"]:
            r.setdefault("n_target", n)
        figure(sm["rows"], sm["references"], os.path.join(ROOT, a.figs, "fig13-smoothing.pdf"))
        return
    tags = sorted((d for d in os.listdir(s15) if d.startswith("s") and
                   os.path.exists(os.path.join(s15, d, "hist-etas5pa.upda"))), key=scale_of)
    jobs = [(run, os.path.join(s15, t, "hist-etas5pa.upda"), a.n_mc) for t in tags]
    jobs += [(run, os.path.join(run, "s09", "hist-etas-mk.upda"), a.n_mc),
             (run, os.path.join(run, "s14", "stationary.json"), a.n_mc)]
    with ProcessPoolExecutor(max_workers=min(len(jobs), 12)) as ex:
        res = list(ex.map(diagnose, jobs))
    rows = []
    for t, r in zip(tags, res):
        r.update(scale=scale_of(t), tag=t, abic=read_abic(os.path.join(s15, t, "hist-etas5pa.prt")))
        rows.append(r)
        print(f"s={r['scale']:<8g} ABIC {r['abic']:.2f}  lnL {r['loglik_exact']:.1f}  "
              f"rho {r['spearman_z_area']:+.3f}  KS_t p {r['ks_uniform_p']:.3f}  "
              f"maxN-L {r['max_abs_N_minus_Lambda']:.1f}  KS_v p {r['voronoi_ks_p']:.2g}  "
              f"zS {r['z_mean_small_cells']:+.2f} zL {r['z_mean_large_cells']:+.2f}  "
              f"bg {r['background_share']:.3f}")
    refs = {"mk": res[-2], "st": res[-1]}
    for k, r in refs.items():
        print(k, {x: round(v, 4) for x, v in r.items()})
    with open(os.path.join(s15, "summary.json"), "w") as fh:
        json.dump(dict(rows=rows, references=refs), fh, indent=2)
    figure(rows, refs, os.path.join(ROOT, a.figs, "fig13-smoothing.pdf"))
    print("wrote", os.path.join(s15, "summary.json"))


if __name__ == "__main__":
    main()
