"""Out-of-sample evaluation: fit on the training catalogue, score the test period.

The training run (run_pipeline.sh on a truncated catalogue) fixes the
tessellation and every field.  Nothing is refitted here.  Test-period events
enter the history as they occur; their parameters come from the frozen
fields, interpolated as the package's own forecast code (workpara5pa.f) and
our rasteriser (src/interp5pa.f) do:

    mu            linear in the value on the training Delaunay triangles
    K, alpha, p, q linear in the log

Anisotropy of test events follows the package's forecast convention (taken
from etas2aniso on the catalogue that contains them, which looks one day
ahead for each mainshock's own aftershocks); --isotropic-test replaces it
with circles, which is strictly prospective.

The log-likelihood is over [T1, T2] x A, the history before T1 entering only
through the triggering sum and its compensator:

    ln L = sum_{T1<=t_i<T2} ln lambda(t_i, x_i) - int_T1^T2 int_A lambda

Reference models, all fitted on the training period only:
    stationary      seven constants, scripts/refit_stationary.py on the training run
    poisson-mu      time-independent Poisson, shape of the 5pa background mu,
                    scaled to the training-period target count
    poisson-uniform time-independent Poisson, uniform over A, same count

Reported: ln L, information gain per event IG = (ln L_a - ln L_b) / N_test,
and the Rhoades et al. (2011) paired T-test on the per-event differences.

    python scripts/oos_eval.py --train run-mc50-train2020 --full run-mc50 --split 2021-01-01
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import etas_lik as E  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import matplotlib.dates  # noqa: E402

COLOURS = {"st": "#2a78d6", "mk": "#eb6834", "5pa": "#1baf7a", "poisson-mu": "#eda100"}
STYLES = {"st": ":", "mk": "--", "5pa": "-", "poisson-mu": "-."}
LABELS = {"st": "stationary ETAS", "mk": "hist-etas-mk", "5pa": "hist-etas5pa",
          "poisson-mu": "Poisson, mu-shaped", "poisson-uniform": "Poisson, uniform"}
INK, MUTED, GRID = "#1f1f1e", "#6b6a66", "#d9d8d4"


def catalogue_origin(hypo_path):
    """tsei2etas takes the first record of hypo.ts as day 0."""
    with open(hypo_path) as fh:
        v = fh.readline().split()
    y, mo, d, h, mi = map(int, v[:5])
    s = float(v[5])
    return dt.datetime(y, mo, d, h, mi) + dt.timedelta(seconds=s)


# ---------------------------------------------------------------- assembly
def combined_catalog(train, full, T1, T2, isotropic_test=False):
    """Training events as fitted, followed by the full catalogue's test events."""
    ntr = train.n
    old = full.t < T1
    if old.sum() != ntr or not np.allclose(full.t[old], train.t, atol=1e-6) \
            or not np.allclose(full.mag[old], train.mag, atol=1e-6):
        raise RuntimeError("training events do not match the full catalogue's pre-split events")
    # etas2aniso moves some mainshocks to their cluster centroid, and its
    # clustering sees the whole catalogue, so a few pre-split epicentres and
    # anisotropies differ between the two runs.  Training events keep the
    # values they were fitted with.
    moved = np.hypot(full.x[old] - train.x, full.y[old] - train.y)
    print(f"pre-split events relocated by etas2aniso between runs: {(moved > 1e-6).sum()}"
          f" (max {moved.max():.3f} deg)")
    new = (full.t >= T1) & (full.t < T2)
    cat = [np.concatenate([getattr(train, k), getattr(full, k)[new]])
           for k in ("x", "y", "mag", "t", "qa", "qb", "qc")]
    if isotropic_test:
        m = new.sum()
        cat[4][ntr:], cat[5][ntr:], cat[6][ntr:] = np.ones(m), np.zeros(m), np.ones(m)
    return E.Catalog(*cat, tx=train.tx, ty=train.ty, tz=T2, tstart=T1, xmg0=train.xmg0,
                     xmin=train.xmin, ymin=train.ymin, cosfac=train.cosfac)


def event_fields(train, f, cat):
    """Fields at every event of `cat`: vertex values for training events,
    interpolated on the training tessellation for the rest."""
    ntr = train.n
    tri = mtri.Triangulation(train.node_x, train.node_y, train.tri)
    xs, ys = cat.x[ntr:], cat.y[ntr:]

    def lin(v):
        out = mtri.LinearTriInterpolator(tri, v)(xs, ys)
        if np.ma.is_masked(out) and out.mask.any():
            raise RuntimeError("test event outside the training tessellation")
        return np.asarray(out)

    mu = np.concatenate([f.mu[:ntr], lin(f.mu)])
    rest = {k: np.concatenate([getattr(f, k)[:ntr], np.exp(lin(np.log(getattr(f, k))))])
            for k in ("K", "a", "p", "q")}
    return E.Fields(mu=mu, **rest, c=f.c, d=f.d)


def background_rate(train, f):
    """int_A mu, per day, on the training tessellation."""
    if f.const_mu:
        return float(f.mu[0]) * train.tx * train.ty
    return E.background_integral(train, f)


# ---------------------------------------------------------------- scoring
def score_etas(cat, fe, B, pr):
    """Per-test-event ln lambda, the cumulative compensator at each test event, and the total."""
    bg, trig = E.intensity_at_events(cat, fe, pr)
    tgt = np.nonzero(cat.target())[0]
    loglam = np.log((bg + trig)[tgt])
    S = E.region_mass(cat, fe, nseg=128)
    w = fe.K * S
    g0 = E.omori_int(cat.tstart - cat.t, fe.c, fe.p)
    comp = np.empty(len(tgt))
    for k, i in enumerate(tgt):
        j = slice(0, i)
        comp[k] = B * (cat.t[i] - cat.tstart) + np.sum(
            w[j] * (E.omori_int(cat.t[i] - cat.t[j], fe.c, fe.p[j]) - g0[j]))
    total = B * (cat.tz - cat.tstart) + np.sum(
        w * (E.omori_int(cat.tz - cat.t, fe.c, fe.p) - g0))
    return dict(loglam=loglam, comp=comp, Lambda=float(total))


def score_poisson(cat, rate_at, rate_total):
    tgt = cat.target()
    loglam = np.log(rate_at(cat.x[tgt], cat.y[tgt]))
    comp = rate_total * (cat.t[tgt] - cat.tstart)
    return dict(loglam=loglam, comp=comp, Lambda=float(rate_total * (cat.tz - cat.tstart)))


def lnL(s):
    return float(s["loglam"].sum() - s["Lambda"])


def t_test(a, b):
    """Rhoades et al. (2011): IG of a over b per event, with a 95% interval."""
    X = a["loglam"] - b["loglam"]
    N = len(X)
    ig = X.mean() - (a["Lambda"] - b["Lambda"]) / N
    s = X.std(ddof=1)
    half = stats.t.ppf(0.975, N - 1) * s / np.sqrt(N)
    return dict(IG=float(ig), lo=float(ig - half), hi=float(ig + half),
                T=float(ig / (s / np.sqrt(N))))


# ---------------------------------------------------------------- figure
def _gain(s, r):
    """Cumulative ln L of s minus r at each test event, then at T2."""
    ll = np.cumsum(s["loglam"] - r["loglam"]) - (s["comp"] - r["comp"])
    return np.concatenate([ll, [ll[-1] - (s["Lambda"] - r["Lambda"])
                                + (s["comp"][-1] - r["comp"][-1])]])


# the two large test-period mainshocks, both in the Molucca Sea region
M74 = dt.datetime(2026, 4, 1, 22, 48)
M78 = dt.datetime(2026, 6, 7, 23, 37)


def fig_cumulative(scores, cat, origin, out, width=6.5):
    """Authored at the width it is included at in the paper (a full text block),
    so the type sizes below are the printed sizes.  The period and event count
    belong in the caption, not in the figure."""
    fig, axes = plt.subplots(1, 2, figsize=(width, 0.4 * width))
    tgt = cat.target()
    days = np.concatenate([cat.t[tgt], [cat.tz]])
    when = [origin + dt.timedelta(days=float(d)) for d in days]
    panels = (("poisson-uniform", ("5pa", "mk", "st", "poisson-mu"), "(a)"),
              ("st", ("5pa", "mk"), "(b)"))
    for ax, (ref, keys, tag) in zip(axes, panels):
        fmt = "+.1f" if ref == "st" else "+.0f"
        for key in keys:
            ll = _gain(scores[key], scores[ref])
            ax.step(when, ll, where="post", color=COLOURS[key], ls=STYLES[key], lw=1.4,
                    label=f"{LABELS[key]}  ({ll[-1]:{fmt}})")
        ax.axhline(0, color=INK, lw=0.6)
        for t in (M74, M78):
            ax.axvline(t, color=MUTED, lw=0.6, ls=":")
        ax.set_ylabel("cumulative gain (nats)", fontsize=8, color=INK)
        ax.set_title(f"{tag} over {LABELS[ref]}", fontsize=8.5, color=INK, loc="left")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=7, length=2.5)
        ax.xaxis.set_major_locator(matplotlib.dates.YearLocator())
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
        ax.grid(color=GRID, lw=0.4)
        ax.set_axisbelow(True)
        ax.legend(fontsize=7, frameon=False, loc="upper left", handlelength=2.2)
    axes[1].annotate("M7.4, M7.8", xy=(M74, 0.09), xycoords=("data", "axes fraction"),
                     xytext=(-3, 0), textcoords="offset points",
                     ha="right", va="bottom", fontsize=7, color=MUTED)
    fig.tight_layout(pad=0.3, w_pad=1.5)
    fig.savefig(out)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", default="run-mc50-train2020")
    ap.add_argument("--full", default="run-mc50")
    ap.add_argument("--split", default="2021-01-01", help="first day of the test period")
    ap.add_argument("--hypo", default="catalog/hypo.ts")
    ap.add_argument("--isotropic-test", action="store_true")
    ap.add_argument("--figs", default="figs-mc50")
    ap.add_argument("--no-fig", action="store_true", help="scores only; leave fig12 alone")
    ap.add_argument("--out-name", default="oos", help="basename of the JSON in <train>/s14")
    a = ap.parse_args()
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tr_run, fu_run = os.path.join(root, a.train), os.path.join(root, a.full)

    origin = catalogue_origin(os.path.join(root, a.hypo))
    T1 = (dt.datetime.fromisoformat(a.split) - origin).total_seconds() / 86400.0

    ctr = E.read_conf(os.path.join(tr_run, "s10", "hist-etas5pa.conf"))
    cfu = E.read_conf(os.path.join(fu_run, "s10", "hist-etas5pa.conf"))
    train = E.read_delone2(os.path.join(tr_run, "s10", "delone2.out"), ctr)
    full = E.read_delone2(os.path.join(fu_run, "s10", "delone2.out"), cfu)
    T2 = cfu["tz"]
    cat = combined_catalog(train, full, T1, T2, a.isotropic_test)
    pr = E.neighbour_pairs(cat, ctr["bi2"])
    n_test = int(cat.target().sum())
    print(f"split {a.split} = day {T1:.3f}; train events {train.n}, test events {n_test}")

    fits = {"mk": E.read_mk(os.path.join(tr_run, "s09", "hist-etas-mk.upda")),
            "5pa": E.read_5pa(os.path.join(tr_run, "s10", "hist-etas5pa.upda"))}
    with open(os.path.join(tr_run, "s14", "stationary.json")) as fh:
        st = json.load(fh)["theta"]
    fits["st"] = E.constant_fields(list(st.values()), train.n)

    scores = {}
    for key, f in fits.items():
        fe = event_fields(train, f, cat) if not f.const_mu else \
            E.constant_fields(list(st.values()), cat.n)
        scores[key] = score_etas(cat, fe, background_rate(train, f), pr)

    # time-independent references, fitted on the training target period
    n_train_tgt = int(train.target().sum())
    rate_total = n_train_tgt / (T1 - ctr["tstart"])
    f5 = fits["5pa"]
    tri = mtri.Triangulation(train.node_x, train.node_y, train.tri)
    mu_interp = mtri.LinearTriInterpolator(tri, f5.mu)
    Bmu = background_rate(train, f5)
    scores["poisson-mu"] = score_poisson(
        cat, lambda x, y: np.asarray(mu_interp(x, y)) * rate_total / Bmu, rate_total)
    area = train.tx * train.ty
    scores["poisson-uniform"] = score_poisson(
        cat, lambda x, y: np.full(len(x), rate_total / area), rate_total)

    res = dict(split=a.split, T1=T1, T2=T2, n_train=train.n, n_train_target=n_train_tgt,
               n_test=n_test, isotropic_test=a.isotropic_test, models={}, tests={})
    for key, s in scores.items():
        res["models"][key] = dict(lnL=lnL(s), Lambda=s["Lambda"],
                                  IG_vs_uniform=(lnL(s) - lnL(scores["poisson-uniform"])) / n_test)
    for a_, b_ in (("5pa", "st"), ("mk", "st"), ("5pa", "mk"), ("st", "poisson-mu"),
                   ("5pa", "poisson-mu"), ("poisson-mu", "poisson-uniform")):
        res["tests"][f"{a_} vs {b_}"] = t_test(scores[a_], scores[b_])

    tag = "-isotropic" if a.isotropic_test else ""
    out = os.path.join(tr_run, "s14")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, f"{a.out_name}{tag}.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"{'model':18}{'ln L':>12}{'Lambda':>10}{'IG/event vs uniform':>22}")
    for k, m in res["models"].items():
        print(f"{k:18}{m['lnL']:12.2f}{m['Lambda']:10.1f}{m['IG_vs_uniform']:22.4f}")
    print(f"\n{'comparison':26}{'IG/event':>10}{'95% interval':>22}{'T':>8}")
    for k, t in res["tests"].items():
        print(f"{k:26}{t['IG']:10.4f}   [{t['lo']:8.4f}, {t['hi']:8.4f}]{t['T']:8.2f}")
    if not a.isotropic_test and not a.no_fig:
        fig_cumulative(scores, cat, origin, os.path.join(root, a.figs, "fig12-oos-gain.pdf"))


if __name__ == "__main__":
    main()
