#!/usr/bin/env python
"""
Completeness (Mc) and b-value analysis for the region catalogue.

Estimators implemented:
  MAXC  maximum-curvature: the mode of the incremental FMD (Wiemer & Wyss 2000)
  GFT   goodness-of-fit: lowest cutoff whose synthetic GR distribution
        reproduces the observed cumulative FMD to R% (Wiemer & Wyss 2000)
  MBS   b-value stability: lowest cutoff above which b stops drifting by more
        than one Shi & Bolt standard deviation (Cao & Gao 2002;
        Woessner & Wiemer 2005)

b is the Aki (1965) / Utsu (1965) maximum likelihood estimate with the Utsu
(1970) half-bin correction for magnitudes rounded to 0.1, and its uncertainty
is Shi & Bolt (1982).

Outputs figures + mc-results.json into the output directory.

Usage: python mc_analysis.py <catalog.csv> <outdir> [--target-mc 5.0]
           [--first-mc 6.0] [--bbox 118 132 -6 6]
"""
import csv, json, argparse, os, datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DM = 0.1  # magnitude bin width; USGS magnitudes are reported to 0.1


# ----------------------------------------------------------------- b-value
def b_aki(mags, mc, dm=DM):
    """Aki-Utsu MLE with Utsu's half-bin correction. Returns (b, sigma_b, n)."""
    m = mags[mags >= mc - dm / 2]
    n = m.size
    if n < 2:
        return np.nan, np.nan, n
    mean = m.mean()
    denom = mean - (mc - dm / 2)
    if denom <= 0:
        return np.nan, np.nan, n
    b = 1.0 / (np.log(10) * denom)
    # Shi & Bolt (1982)
    sig = 2.30 * b**2 * np.sqrt(((m - mean) ** 2).sum() / (n * (n - 1)))
    return b, sig, n


def a_value(mags, mc, b, dm=DM):
    n = (mags >= mc - dm / 2).sum()
    return np.log10(n) + b * mc if n > 0 else np.nan


# ------------------------------------------------------------- estimators
def maxc(mags, dm=DM):
    edges = np.arange(mags.min() - dm / 2, mags.max() + dm, dm)
    cnt, _ = np.histogram(mags, bins=edges)
    centres = edges[:-1] + dm / 2
    return round(float(centres[np.argmax(cnt)]), 1)


def gft(mags, dm=DM):
    """Return {90: Mc, 95: Mc} -- lowest cutoff reaching each R level."""
    out = {}
    cuts = np.round(np.arange(mags.min(), mags.max() - 0.5, dm), 1)
    for mco in cuts:
        sub = mags[mags >= mco - dm / 2]
        if sub.size < 50:
            continue
        b, _, n = b_aki(sub, mco, dm)
        if not np.isfinite(b):
            continue
        a = a_value(sub, mco, b, dm)
        bins = np.round(np.arange(mco, sub.max() + dm, dm), 1)
        obs = np.array([(sub >= x - dm / 2).sum() for x in bins], float)
        syn = 10 ** (a - b * bins)
        R = 100 * (1 - np.abs(obs - syn).sum() / obs.sum())
        for lvl in (95, 90):
            if lvl not in out and R >= lvl:
                out[lvl] = float(mco)
    return out


def mbs(mags, dm=DM, window=0.5):
    """b-value stability: |b(Mi) - mean b over [Mi, Mi+window]| <= sigma_b."""
    cuts = np.round(np.arange(mags.min(), mags.max() - 1.0, dm), 1)
    bs, sg = {}, {}
    for mco in cuts:
        b, s, n = b_aki(mags, mco, dm)
        if n >= 50 and np.isfinite(b):
            bs[mco], sg[mco] = b, s
    ks = sorted(bs)
    nw = int(round(window / dm))
    for i, mco in enumerate(ks):
        nxt = [bs[k] for k in ks[i:i + nw + 1]]
        if len(nxt) < nw + 1:
            break
        if abs(bs[mco] - np.mean(nxt)) <= sg[mco]:
            return float(mco)
    return None


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("catalog"); ap.add_argument("outdir")
    ap.add_argument("--target-mc", type=float, default=5.0)
    ap.add_argument("--first-mc", type=float, default=6.0,
                    help="threshold first tried, drawn for reference (0 = none)")
    # the USGS extract carries one stray California event (nc21364840) that
    # tsei2etas drops, so filter to the study box as the fits do
    ap.add_argument("--bbox", type=float, nargs=4, default=[118, 132, -6, 6],
                    metavar=("LON0", "LON1", "LAT0", "LAT1"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    lon0, lon1, lat0, lat1 = a.bbox
    rows = [r for r in csv.DictReader(open(a.catalog, encoding="utf-8"))
            if r.get("mag") and (r.get("type") or "earthquake") == "earthquake"
            and lon0 <= float(r["longitude"]) <= lon1
            and lat0 <= float(r["latitude"]) <= lat1]
    mags = np.array([float(r["mag"]) for r in rows])
    times = np.array([datetime.datetime.fromisoformat(r["time"].replace("Z", ""))
                      for r in rows])
    mtype = np.array([(r.get("magType") or "").strip().lower() for r in rows])
    order = np.argsort(times)
    mags, times, mtype = mags[order], times[order], mtype[order]

    res = {"n_events": int(mags.size),
           "mag_range": [float(mags.min()), float(mags.max())],
           "time_range": [str(times[0]), str(times[-1])],
           "target_mc": a.target_mc}

    # ---- estimators on the full catalogue
    res["maxc"] = maxc(mags)
    res["gft"] = gft(mags)
    res["mbs"] = mbs(mags)

    # ---- b at the estimated Mc and at the ETAS target threshold
    mc_pref = res["gft"].get(90, res["maxc"])
    for tag, mc in (("at_mc", mc_pref), ("at_target", a.target_mc)):
        b, s, n = b_aki(mags, mc)
        res[tag] = {"mc": float(mc), "b": float(b), "sigma_b": float(s), "n": int(n)}

    # ---- b as a function of cutoff
    cuts = np.round(np.arange(3.0, 7.0, DM), 1)
    bb, ss, nn = [], [], []
    for mco in cuts:
        b, s, n = b_aki(mags, mco)
        bb.append(b); ss.append(s); nn.append(n)
    bb, ss, nn = np.array(bb), np.array(ss), np.array(nn)

    # ---- annual Mc, to expose any change in network performance
    yrs, mc_y, b_y = [], [], []
    for y in range(times[0].year, times[-1].year + 1):
        sel = np.array([t.year == y for t in times])
        if sel.sum() < 100:
            continue
        sub = mags[sel]
        g = gft(sub).get(90)
        m_ = g if g is not None else maxc(sub)
        b_, _, _ = b_aki(sub, m_)
        yrs.append(y); mc_y.append(m_); b_y.append(b_)
    res["annual"] = {"year": yrs, "mc": mc_y, "b": b_y}
    res["annual_mc_max"] = float(np.nanmax(mc_y)) if mc_y else None

    # =================================================================
    # Figure 1 -- frequency-magnitude distribution
    # =================================================================
    edges = np.round(np.arange(mags.min() - DM / 2, mags.max() + DM, DM), 3)
    cnt, _ = np.histogram(mags, bins=edges)
    ctr = edges[:-1] + DM / 2
    cum = np.array([(mags >= x - DM / 2).sum() for x in ctr], float)

    # authored at its printed size (0.6 of the SN text width) so the type is
    # 7-8 pt on the page; no title, the caption carries it
    with plt.rc_context({"font.size": 7.5, "axes.labelsize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7}):
        fig, ax = plt.subplots(figsize=(3.1, 2.4))
        mcs = [v for v in (res["maxc"], res["gft"].get(90), res["mbs"]) if v is not None]
        ax.axvspan(min(mcs), max(mcs), color="#7b1fa2", alpha=.15, lw=0,
                   label=f"completeness $M_c$ {min(mcs):.1f}–{max(mcs):.1f}")
        ax.semilogy(ctr, cum, "o", ms=2.8, mfc="none", mec="#1f4e79",
                    mew=.8, label="cumulative $N(\\geq M)$")
        ax.semilogy(ctr, np.where(cnt > 0, cnt, np.nan), "s", ms=2.0,
                    color="#c0c0c0", mec="#909090", mew=.5, label="incremental")

        b_p = res["at_mc"]["b"]; a_p = a_value(mags, mc_pref, b_p)
        xs2 = np.linspace(mc_pref, mags.max(), 50)
        ax.semilogy(xs2, 10 ** (a_p - b_p * xs2), "--", color="#2e7d32", lw=1.1,
                    label=f"GR fit, $M\\geq{mc_pref:.1f}$ ($b={b_p:.2f}$)")
        b_t = res["at_target"]["b"]; a_t = a_value(mags, a.target_mc, b_t)
        xs = np.linspace(a.target_mc, mags.max(), 50)
        ax.semilogy(xs, 10 ** (a_t - b_t * xs), "-", color="#c00000", lw=1.3,
                    label=f"GR fit, $M\\geq{a.target_mc:.1f}$ ($b={b_t:.2f}$)")

        ax.axvline(a.target_mc, color="#c00000", lw=1.0)
        ax.text(a.target_mc - .05, 1.3, f"target M {a.target_mc:.1f}", color="#c00000",
                fontsize=6.5, rotation=90, va="bottom", ha="right")
        if a.first_mc:
            # in the legend: a rotated label here would sit on the squares
            ax.axvline(a.first_mc, color="#606060", lw=.9, ls="-.",
                       label=f"threshold first tried, M {a.first_mc:.1f}")

        ax.set_xlabel("magnitude $M$"); ax.set_ylabel("number of events")
        ax.legend(fontsize=6.3, frameon=False, loc="upper right",
                  handlelength=1.8, borderaxespad=.3, labelspacing=.35)
        ax.grid(alpha=.25, which="both", lw=.3)
        ax.set_xlim(3.0, mags.max() + .2)
        ax.set_ylim(.7, 6e5)
        fig.tight_layout(pad=.3); fig.savefig(f"{a.outdir}/fig-fmd.pdf"); plt.close(fig)

    # =================================================================
    # Figure 2 -- b-value vs cutoff (the MBS view)
    # =================================================================
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ok = np.isfinite(bb) & (np.array(nn) >= 50)
    ax.fill_between(cuts[ok], (bb - ss)[ok], (bb + ss)[ok],
                    color="#1f4e79", alpha=.18, lw=0)
    ax.plot(cuts[ok], bb[ok], "-", color="#1f4e79", lw=1.8)
    ax.axvline(a.target_mc, color="#c00000", lw=1.4)
    ax.text(a.target_mc + .04, ax.get_ylim()[1], " ETAS threshold",
            color="#c00000", fontsize=8, va="top")
    if res["mbs"]:
        ax.axvline(res["mbs"], color="#7b1fa2", ls=":", lw=1.4)
    ax.set_xlabel("cutoff magnitude $M_{co}$")
    ax.set_ylabel("$b(M_{co})$  $\\pm\\,1\\sigma$ (Shi & Bolt)")
    ax.set_title("b-value stability with cutoff", fontsize=10.5)
    ax.grid(alpha=.25, lw=.4)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/fig-bstab.pdf"); plt.close(fig)

    # =================================================================
    # Figure 3 -- annual Mc and b
    # =================================================================
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(yrs, mc_y, "o-", color="#1f4e79", ms=4, lw=1.5, label="annual $M_c$ (GFT 90%)")
    ax.axhline(a.target_mc, color="#c00000", lw=1.6, label=f"ETAS threshold {a.target_mc:.1f}")
    ax.fill_between([min(yrs), max(yrs)], 0, a.target_mc, color="#2e7d32",
                    alpha=.07, lw=0)
    ax.set_xlabel("year"); ax.set_ylabel("$M_c$")
    ax.set_ylim(3.5, 6.4)
    ax.set_title("Annual completeness magnitude against the ETAS threshold",
                 fontsize=10.5)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.grid(alpha=.25, lw=.4)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/fig-mctime.pdf"); plt.close(fig)

    # =================================================================
    # Figure 4 -- magnitude-type composition
    # =================================================================
    bands = [(3.0, 4.0), (4.0, 4.5), (4.5, 5.0), (5.0, 5.5),
             (5.5, 6.0), (6.0, 6.5), (6.5, 8.5)]
    lab, f_mb, f_mw = [], [], []
    comp = {}
    for lo, hi in bands:
        s = (mags >= lo) & (mags < hi)
        if s.sum() == 0:
            continue
        t = mtype[s]
        mb = np.char.startswith(t, "mb").sum()
        mw = np.char.startswith(t, "mw").sum()
        lab.append(f"{lo:.1f}–{hi:.1f}")
        f_mb.append(100 * mb / s.sum()); f_mw.append(100 * mw / s.sum())
        comp[f"{lo}-{hi}"] = {"n": int(s.sum()), "pct_mb": round(100*mb/s.sum(), 1),
                              "pct_mw": round(100*mw/s.sum(), 1)}
    res["magtype_composition"] = comp

    fig, ax = plt.subplots(figsize=(7.0, 3.8))
    x = np.arange(len(lab))
    ax.bar(x, f_mb, .62, label="$m_b$ (body wave)", color="#d99694")
    ax.bar(x, f_mw, .62, bottom=f_mb, label="$M_w$ (moment)", color="#1f4e79")
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8)
    ax.set_ylabel("per cent of events"); ax.set_xlabel("magnitude band")
    ax.set_title("Magnitude-type composition: the target set is effectively $M_w$",
                 fontsize=10.5)
    ax.legend(fontsize=8, frameon=False, loc="center left")
    ax.set_ylim(0, 108)
    fig.tight_layout(); fig.savefig(f"{a.outdir}/fig-magtype.pdf"); plt.close(fig)

    json.dump(res, open(f"{a.outdir}/mc-results.json", "w"), indent=2)

    # ---- console summary
    print(f"events            : {res['n_events']}  ({res['mag_range'][0]}–{res['mag_range'][1]})")
    print(f"MAXC              : {res['maxc']}")
    print(f"GFT               : {res['gft']}")
    print(f"MBS               : {res['mbs']}")
    print(f"b at Mc={mc_pref:<4}     : {res['at_mc']['b']:.3f} ± {res['at_mc']['sigma_b']:.3f}"
          f"  (n={res['at_mc']['n']})")
    print(f"b at target {a.target_mc}   : {res['at_target']['b']:.3f} ± "
          f"{res['at_target']['sigma_b']:.3f}  (n={res['at_target']['n']})")
    print(f"max annual Mc     : {res['annual_mc_max']}")
    print(f"figures + mc-results.json -> {a.outdir}")


if __name__ == "__main__":
    main()
