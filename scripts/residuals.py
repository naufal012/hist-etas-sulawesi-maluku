"""Residual analysis of the fitted ETAS models (stage 14).

For the stationary, hist-etas-mk and hist-etas5pa fits:

  * temporal residuals -- random time change tau_i = Lambda(t_i), Ogata (1988);
    under the model the tau_i are a unit-rate Poisson process;
  * Voronoi residuals -- Bray et al. (2014): Lambda over [S,T] x (Voronoi cell
    of each target event), approximately Gamma(3.569, 3.569) under the model;
  * the expected number of survivors of Schoenberg (2003) thinning, which
    decides whether a thinned-point K-function is usable at all.

The stationary model is refitted here rather than taken from st-etas, whose
compensator places the events at absolute lon/lat against a box [0,TX]x[0,TY]
(see etas_lik.read_out3).  The refit maximises the same likelihood the
hierarchical programs use, with every field held constant -- the phi == 0
member of the hist-etas family, so all three log-likelihoods are comparable.

    python scripts/residuals.py --region region-mc50.conf
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from scipy import stats
from scipy.optimize import minimize
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import etas_lik as E  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# categorical slots 1-3 of the reference palette, validated all-pairs
COLOURS = {"st": "#2a78d6", "mk": "#eb6834", "5pa": "#1baf7a"}
STYLES = {"st": ":", "mk": "--", "5pa": "-"}
LABELS = {"st": "stationary", "mk": "hist-etas-mk", "5pa": "hist-etas5pa"}
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#d9d8d4"
DIVERGING = LinearSegmentedColormap.from_list(
    "bnr", ["#1d5fae", "#6f9fd8", "#f0efec", "#e38a86", "#b3292b"])
VORONOI_SHAPE = 3.569  # Bray et al. (2014), after Hinde & Miles (1980)


def region_conf(path: str) -> dict:
    out = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#")[0].strip()
            if "=" in line:
                k, v = line.split("=", 1)
                out[k.strip()] = v.strip()
    return out


# ---------------------------------------------------------------- stationary refit
def fit_stationary(cat, pr, theta0):
    """MLE of the seven constants under the hist likelihood (fortran mode)."""
    names = ["mu", "K", "c", "alpha", "p", "d", "q"]

    def nll(z):
        th = np.exp(z)
        if abs(th[4] - 1) < 1e-9 or th[6] <= 1.0:
            return 1e30
        f = E.constant_fields(th, cat.n)
        with np.errstate(all="ignore"):
            v = -E.loglik(cat, f, pr, "fortran")["loglik"]
        return v if np.isfinite(v) else 1e30

    z0 = np.log(theta0)
    r = minimize(nll, z0, method="Nelder-Mead",
                 options=dict(maxiter=6000, maxfev=12000, xatol=1e-7, fatol=1e-7, adaptive=True))
    r = minimize(nll, r.x, method="L-BFGS-B", options=dict(maxiter=500))
    r = minimize(nll, r.x, method="Nelder-Mead",
                 options=dict(maxiter=6000, xatol=1e-9, fatol=1e-9, adaptive=True))
    th = np.exp(r.x)
    return dict(zip(names, th.tolist())), -r.fun


# ---------------------------------------------------------------- Voronoi residuals
def sample_kernel(cat, fe, j, n, rng):
    """n draws from event j's spatial kernel over the whole plane (scaled deg)."""
    q, d, a = fe.q[j], fe.d, fe.a[j]
    u = rng.random(n)
    s = d * ((1 - u) ** (1 / (1 - q)) - 1)          # |w|^2, radial CDF inverted
    th = 2 * np.pi * rng.random(n)
    w1, w2 = np.sqrt(s) * np.cos(th), np.sqrt(s) * np.sin(th)
    l11 = np.sqrt(cat.qa[j])
    l21 = -cat.qb[j] / l11
    l22 = np.sqrt(cat.qc[j] - l21 ** 2)
    v2 = w2 / l22
    v1 = (w1 - l21 * v2) / l11
    sc = np.exp(a * cat.mag[j])
    return cat.x[j] + sc * v1, cat.y[j] + sc * v2


def cell_compensator(cat, f, assign, ncell, n_mc=4000, grid=0.02, seed=1):
    """Lambda([S,T] x C_k) for a partition of the box.

    assign(x, y) maps scaled coordinates to a cell index in [0, ncell).
    Background: mu linear on the Delaunay triangles, midpoint rule on a grid.
    Triggering: n_mc draws from each event's kernel, weighted by its mass and
    its Omori integral over [S, T]; draws outside the box are dropped.
    """
    rng = np.random.default_rng(seed)
    fe = f.at_events(cat.n)
    lam = np.zeros(ncell)

    gx = np.arange(grid / 2, cat.tx, grid)
    gy = np.arange(grid / 2, cat.ty, grid)
    GX, GY = np.meshgrid(gx, gy)
    if f.const_mu:
        mu = np.full(GX.shape, fe.mu[0])
    else:
        tri = mtri.Triangulation(cat.node_x, cat.node_y, cat.tri)
        mu = mtri.LinearTriInterpolator(tri, f.mu)(GX, GY).filled(np.nan)
        if np.isnan(mu).any():
            raise RuntimeError("background grid falls outside the triangulation")
    lam += np.bincount(assign(GX.ravel(), GY.ravel()), weights=mu.ravel() * grid * grid,
                       minlength=ncell) * (cat.tz - cat.tstart)

    T = (E.omori_int(cat.tz - cat.t, fe.c, fe.p)
         - E.omori_int(cat.tstart - cat.t, fe.c, fe.p))
    total = np.pi * np.exp(2 * fe.a * cat.mag) * fe.d ** (1 - fe.q) / (fe.q - 1)
    for j in range(cat.n):
        x, y = sample_kernel(cat, fe, j, n_mc, rng)
        inside = (x >= 0) & (x <= cat.tx) & (y >= 0) & (y <= cat.ty)
        if inside.any():
            lam += np.bincount(assign(x[inside], y[inside]), minlength=ncell) * (
                fe.K[j] * T[j] * total[j] / n_mc)
    return lam


def voronoi_assign(cat):
    tgt = np.nonzero(cat.target())[0]
    tree = cKDTree(np.column_stack([cat.x[tgt], cat.y[tgt]]))
    return (lambda x, y: tree.query(np.column_stack([x, y]))[1]), len(tgt)


def grid_assign(cat, h):
    nx = int(round(cat.tx / h))
    ny = int(round(cat.ty / h))

    def assign(x, y):
        i = np.clip((x / cat.tx * nx).astype(int), 0, nx - 1)
        k = np.clip((y / cat.ty * ny).astype(int), 0, ny - 1)
        return k * nx + i
    return assign, nx * ny, (nx, ny)


def deviance_by_cell(cat, pr, f1, f2, lam1, lam2, assign, ncell):
    """Clements et al. (2011): per cell, log-likelihood of model 1 minus model 2."""
    tgt = cat.target()
    l1 = sum(E.intensity_at_events(cat, f1.at_events(cat.n), pr))[tgt]
    l2 = sum(E.intensity_at_events(cat, f2.at_events(cat.n), pr))[tgt]
    k = assign(cat.x[tgt], cat.y[tgt])
    return np.bincount(k, weights=np.log(l1) - np.log(l2), minlength=ncell) - (lam1 - lam2)


# ---------------------------------------------------------------- figures
# Authored at the width they are included at -- the Springer Nature text block,
# 372 pt = 5.17 in -- so the type sizes below are the printed sizes.  Notes on
# shading and colour belong in the captions, not in the figures.
WIDTH = 5.15


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=6.5, length=2.5)
    ax.grid(color=GRID, lw=0.4)
    ax.set_axisbelow(True)


def _map_axes(ax, ext, left=True):
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])
    # stop short of the right edge, where the label would meet the next panel's
    ax.set_xticks(np.arange(np.ceil(ext[0] / 4) * 4, ext[1] - 0.5, 4))
    ax.set_yticks(np.arange(np.ceil(ext[2] / 3) * 3, ext[3] + 0.1, 3))
    ax.tick_params(labelsize=6.5, colors=MUTED, labelcolor=INK, length=2)
    if not left:
        ax.set_yticklabels([])


def fig_temporal(res, cat, year0, out, width=WIDTH):
    fig, axes = plt.subplots(1, 2, figsize=(width, 0.46 * width))
    ax = axes[0]
    n = len(next(iter(res.values()))["tau"])
    lim = max(max(r["tau"][-1] for r in res.values()), n) * 1.02
    xs = np.linspace(0, lim, 200)
    for z, a_ in ((1.63, 0.10), (1.36, 0.16)):
        ax.fill_between(xs, xs - z * np.sqrt(n), xs + z * np.sqrt(n),
                        color=MUTED, alpha=a_, lw=0)
    ax.plot(xs, xs, color=INK, lw=0.6)
    for key, r in res.items():
        ax.step(r["tau"], np.arange(1, n + 1), where="post", color=COLOURS[key],
                ls=STYLES[key], lw=1.2, label=LABELS[key])
    ax.legend(fontsize=6.5, frameon=False, loc="upper left", handlelength=2.2)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("transformed time $\\tilde t = \\Lambda(t)$", fontsize=7.5, color=INK)
    ax.set_ylabel("cumulative number of events", fontsize=7.5, color=INK)
    ax.set_title("(a) residual process", fontsize=8, color=INK, loc="left")
    _style(ax)

    ax = axes[1]
    tgt = cat.target()
    yrs = year0 + cat.t[tgt] / 365.25
    band = 1.36 * np.sqrt(n)
    ax.axhspan(-1.63 * np.sqrt(n), 1.63 * np.sqrt(n), color=MUTED, alpha=0.10, lw=0)
    ax.axhspan(-band, band, color=MUTED, alpha=0.16, lw=0)
    ax.axhline(0, color=INK, lw=0.6)
    for key, r in res.items():
        ax.plot(yrs, np.arange(1, n + 1) - r["tau"], color=COLOURS[key],
                ls=STYLES[key], lw=1.2, label=LABELS[key])
    ax.set_xlabel("year", fontsize=7.5, color=INK)
    ax.set_ylabel("$N(t) - \\Lambda(t)$", fontsize=7.5, color=INK)
    ax.set_title("(b) observed minus expected", fontsize=8, color=INK, loc="left")
    _style(ax)
    fig.tight_layout(pad=0.3, w_pad=1.2)
    fig.savefig(out)
    plt.close(fig)


def load_coast(maps_dir):
    pb = np.loadtxt(os.path.join(maps_dir, "w1.pp.data")).reshape(-1, 2)
    br = np.loadtxt(os.path.join(maps_dir, "w1.br.data")).astype(int)
    ends = np.cumsum(br)
    return [pb[e - n:e] for n, e in zip(br, ends)]


def fig_voronoi(vres, cat, coast, out, width=WIDTH):
    keys = list(vres)
    nk = len(keys)
    fig = plt.figure(figsize=(width, 0.36 * width))
    gs = fig.add_gridspec(1, nk + 1, width_ratios=[1] * nk + [0.06],
                          left=0.06, right=0.90, bottom=0.12, top=0.90, wspace=0.08)
    axes = [fig.add_subplot(gs[0, i]) for i in range(nk)]
    tgt = np.nonzero(cat.target())[0]
    lon = cat.xmin + cat.x[tgt] / cat.cosfac
    lat = cat.ymin + cat.y[tgt]
    ext = (cat.xmin, cat.xmin + cat.tx / cat.cosfac, cat.ymin, cat.ymin + cat.ty)
    # pixel map of the cell each pixel belongs to
    px = np.linspace(ext[0], ext[1], 700)
    py = np.linspace(ext[2], ext[3], 600)
    PX, PY = np.meshgrid(px, py)
    tree = cKDTree(np.column_stack([cat.x[tgt], cat.y[tgt]]))
    _, cell = tree.query(np.column_stack([((PX - cat.xmin) * cat.cosfac).ravel(),
                                          (PY - cat.ymin).ravel()]))
    cell = cell.reshape(PX.shape)
    for i, (ax, key) in enumerate(zip(axes, keys)):
        u = stats.gamma.cdf(vres[key], VORONOI_SHAPE, scale=1 / VORONOI_SHAPE)
        # sign chosen so red = more events than the model expects
        z = stats.norm.ppf(np.clip(1 - u, 1e-6, 1 - 1e-6))
        im = ax.imshow(z[cell], origin="lower", extent=ext, cmap=DIVERGING,
                       vmin=-3, vmax=3, interpolation="nearest", aspect="equal")
        for seg in coast:
            ax.plot(seg[:, 1], seg[:, 0], color=INK, lw=0.3)
        ax.plot(lon, lat, ".", ms=0.5, color=INK, alpha=0.4, mew=0)
        _map_axes(ax, ext, left=(i == 0))
        ax.set_title(f"({'abc'[i]}) {LABELS[key]}", fontsize=7.5, color=INK, loc="left")
    cax = fig.add_subplot(gs[0, nk])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("normal score", fontsize=6.5, color=INK)
    cb.ax.tick_params(labelsize=6, length=2)
    fig.savefig(out)
    plt.close(fig)


def fig_deviance(dev, shape, cat, coast, out, width=WIDTH):
    """Per-cell log-likelihood differences, one panel per model pair."""
    nx, ny = shape
    ext = (cat.xmin, cat.xmin + cat.tx / cat.cosfac, cat.ymin, cat.ymin + cat.ty)
    vmax = max(np.abs(v).max() for v in dev.values())
    norm = matplotlib.colors.SymLogNorm(linthresh=2.0, vmin=-vmax, vmax=vmax, base=10)
    nd = len(dev)
    fig = plt.figure(figsize=(width, 0.44 * width))
    gs = fig.add_gridspec(1, nd + 1, width_ratios=[1] * nd + [0.05],
                          left=0.06, right=0.89, bottom=0.09, top=0.84, wspace=0.08)
    axes = [fig.add_subplot(gs[0, i]) for i in range(nd)]
    tgt = np.nonzero(cat.target())[0]
    for i, (ax, ((k1, k2), v)) in enumerate(zip(axes, dev.items())):
        im = ax.imshow(v.reshape(ny, nx), origin="lower", extent=ext,
                       cmap=DIVERGING.reversed(), norm=norm,
                       interpolation="nearest", aspect="equal")
        for seg in coast:
            ax.plot(seg[:, 1], seg[:, 0], color=INK, lw=0.3)
        ax.plot(cat.xmin + cat.x[tgt] / cat.cosfac, cat.ymin + cat.y[tgt], ".",
                ms=0.5, color=INK, alpha=0.35, mew=0)
        _map_axes(ax, ext, left=(i == 0))
        ax.set_title(f"({'ab'[i]}) {LABELS[k1]} vs {LABELS[k2]}\ntotal {v.sum():+.1f} nats",
                     fontsize=7.5, color=INK, loc="left")
    cax = fig.add_subplot(gs[0, nd])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("nats per 1° cell", fontsize=6.5, color=INK)
    cb.ax.tick_params(labelsize=6, length=2)
    fig.savefig(out)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", default="region-mc50.conf")
    ap.add_argument("--n-mc", type=int, default=4000)
    ap.add_argument("--refit", action="store_true",
                    help="refit the stationary model even if s14/stationary.json exists")
    a = ap.parse_args()
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    rc = region_conf(os.path.join(root, a.region))
    # same run/figs naming as run_pipeline.sh: region-mc50.conf -> run-mc50
    tag = os.path.basename(a.region)[:-len(".conf")]
    tag = "" if tag == "region" else "-" + tag.removeprefix("region-")
    run = os.path.join(root, "run" + tag)
    figs = os.path.join(root, "figs" + tag)
    out = os.path.join(run, "s14")
    os.makedirs(out, exist_ok=True)
    year0 = float(rc.get("CATALOG_START", "2000-01-01")[:4])

    conf = E.read_conf(os.path.join(run, "s10", "hist-etas5pa.conf"))
    cat = E.read_delone2(os.path.join(run, "s10", "delone2.out"), conf)
    pr = E.neighbour_pairs(cat, conf["bi2"])
    fits = {"mk": E.read_mk(os.path.join(run, "s09", "hist-etas-mk.upda")),
            "5pa": E.read_5pa(os.path.join(run, "s10", "hist-etas5pa.upda"))}

    # stationary refit, started from the mk baselines
    st_json = os.path.join(out, "stationary.json")
    if os.path.exists(st_json) and not a.refit:
        with open(st_json) as fh:
            st_par = json.load(fh)["theta"]
        st_ll = E.loglik(cat, E.constant_fields(list(st_par.values()), cat.n),
                         pr, "fortran")["loglik"]
    else:
        f0 = fits["mk"]
        theta0 = [f0.mu.mean(), np.exp(np.log(f0.K).mean()), f0.c, 2 * f0.a[0],
                  f0.p[0], f0.d, f0.q[0]]
        st_par, st_ll = fit_stationary(cat, pr, theta0)
        with open(st_json, "w") as fh:
            json.dump(dict(theta=st_par, loglik=st_ll), fh, indent=2)
    fits = {"st": E.constant_fields(list(st_par.values()), cat.n), **fits}
    print("stationary refit:", {k: f"{v:.5g}" for k, v in st_par.items()}, f"logL={st_ll:.3f}")

    summary = {"n_target": int(cat.target().sum()), "stationary_refit": st_par, "models": {}}
    temporal = {}
    for key, f in fits.items():
        ll_f = E.loglik(cat, f, pr, "fortran")
        mass = E.region_mass(cat, f.at_events(cat.n), nseg=128)
        ll_x = E.loglik(cat, f, pr, "exact", mass=mass)
        tau = E.compensator_at_events(cat, f, mass)
        gaps = np.diff(np.concatenate([[0.0], tau]))
        ks_gap = stats.kstest(gaps, "expon")
        ks_uni = stats.kstest(tau / (ll_x["trig_int"] + ll_x["bg_int"]), "uniform")
        bg, trig = E.intensity_at_events(cat, f.at_events(cat.n), pr)
        lam = (bg + trig)[cat.target()]
        surv = float(np.sum(lam.min() / lam))
        temporal[key] = dict(tau=tau)
        n_par = {"st": 7}.get(key)
        summary["models"][key] = dict(
            loglik_fortran=ll_f["loglik"], loglik_exact=ll_x["loglik"],
            Lambda_T=ll_x["trig_int"] + ll_x["bg_int"],
            Lambda_T_background=ll_x["bg_int"],
            ks_gaps_D=ks_gap.statistic, ks_gaps_p=ks_gap.pvalue,
            ks_uniform_D=ks_uni.statistic, ks_uniform_p=ks_uni.pvalue,
            max_abs_N_minus_Lambda=float(np.max(np.abs(np.arange(1, len(tau) + 1) - tau))),
            thinning_expected_survivors=surv,
            aic=(2 * n_par - 2 * ll_f["loglik"]) if n_par else None)
        np.savetxt(os.path.join(out, f"tau-{key}.txt"), tau, fmt="%.6f")
        print(key, {k: (round(v, 4) if isinstance(v, float) else v)
                    for k, v in summary["models"][key].items()})

    fig_temporal(temporal, cat, year0, os.path.join(figs, "fig9-residuals-temporal.pdf"))

    vassign, nv = voronoi_assign(cat)
    gassign, ng, gshape = grid_assign(cat, 1.0)
    vres, gres = {}, {}
    for key, f in fits.items():
        v = cell_compensator(cat, f, vassign, nv, n_mc=a.n_mc)
        gres[key] = cell_compensator(cat, f, gassign, ng, n_mc=a.n_mc)
        vres[key] = v
        ks = stats.kstest(v, "gamma", args=(VORONOI_SHAPE, 0, 1 / VORONOI_SHAPE))
        u = stats.gamma.cdf(v, VORONOI_SHAPE, scale=1 / VORONOI_SHAPE)
        summary["models"][key].update(
            voronoi_sum=float(v.sum()), voronoi_ks_D=ks.statistic, voronoi_ks_p=ks.pvalue,
            voronoi_frac_below_1pct=float(np.mean(u < 0.01)),
            voronoi_frac_above_99pct=float(np.mean(u > 0.99)))
        np.savetxt(os.path.join(out, f"voronoi-{key}.txt"), v, fmt="%.6g")
        print(key, "voronoi sum", round(v.sum(), 1), "KS D", round(ks.statistic, 4),
              "p", f"{ks.pvalue:.3g}")

    coast = load_coast(os.path.join(root, "maps"))
    fig_voronoi(vres, cat, coast, os.path.join(figs, "fig10-residuals-voronoi.pdf"))

    dev = {}
    for k1, k2 in (("mk", "st"), ("5pa", "mk")):
        d = deviance_by_cell(cat, pr, fits[k1], fits[k2], gres[k1], gres[k2], gassign, ng)
        dev[(k1, k2)] = d
        exact = summary["models"][k1]["loglik_exact"] - summary["models"][k2]["loglik_exact"]
        summary[f"deviance_{k1}_vs_{k2}"] = dict(
            total=float(d.sum()), exact_loglik_difference=exact,
            cells_first_better=int((d > 0).sum()), cells_second_better=int((d < 0).sum()),
            cells=int(ng))
        np.savetxt(os.path.join(out, f"deviance-{k1}-vs-{k2}.txt"),
                   d.reshape(gshape[1], gshape[0]), fmt="%.4f")
        print(f"deviance {k1} vs {k2}: total {d.sum():.1f} (exact {exact:.1f}); "
              f"{(d > 0).sum()} cells favour {k1}, {(d < 0).sum()} favour {k2}")
    fig_deviance(dev, gshape, cat, coast, os.path.join(figs, "fig11-residuals-deviance.pdf"))
    with open(os.path.join(out, "residuals.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("wrote", out)


if __name__ == "__main__":
    main()
