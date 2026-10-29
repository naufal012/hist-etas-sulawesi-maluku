#!/usr/bin/env python
"""
Collect ETAS results across target-magnitude thresholds and plot the trend.

Reads each run tree produced by run_pipeline.sh, takes the stationary
refit (s14/stationary.json, from refit_stationary.py -- not the package's
st-etas, whose compensator is defective) and the hist-etas-mk hierarchical
baselines, counts the
triggering pairs available at each threshold, and writes a comparison
figure plus a JSON/LaTeX summary.

The point of the figure is to test one prediction: if the low Omori
exponent at M>=6.0 is a sampling artefact rather than regional physics,
p should climb toward 1 as the threshold falls and triggered pairs
multiply. If p were flat, the diagnosis would be wrong.

Usage: python threshold_summary.py <outdir>
"""
import re, os, sys, json, csv, math, datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JAPAN = r"D:/Documents/Hierarchical ETAS/run"

RUNS = [  # (label, threshold, run directory)
    ("M>=6.0", 6.0, os.path.join(ROOT, "run")),
    ("M>=5.5", 5.5, os.path.join(ROOT, "run-mc55")),
    ("M>=5.0", 5.0, os.path.join(ROOT, "run-mc50")),
]
NAMES = ["mu", "K0", "c", "alpha", "p", "d", "q"]


def stationary(run):
    """Stationary refit from s14/stationary.json (scripts/refit_stationary.py).

    Returns (values, standard errors). The package's st-etas MLE is not used:
    its compensator integrates events at absolute lon/lat against [0,TX]x[0,TY].
    """
    p = os.path.join(run, "s14", "stationary.json")
    if not os.path.exists(p):
        return None, None
    d = json.load(open(p))
    key = {"K0": "K"}
    return ([d["theta"][key.get(k, k)] for k in NAMES],
            [d["se"][key.get(k, k)] for k in NAMES])


def stationary_st_etas(run):
    """Last 'mle =' block in st-etas.prt (defective; kept for reference only)."""
    p = os.path.join(run, "s05", "st-etas.prt")
    if not os.path.exists(p):
        return None
    txt = open(p, encoding="utf-8", errors="replace").read().replace("\r", "")
    blocks = re.findall(r"mle *=([^\n]*)\n([^\n]*)", txt)
    if not blocks:
        return None
    v = [float(x) for x in re.findall(r"[-+]?\d*\.\d+E[-+]\d+",
                                      blocks[-1][0] + " " + blocks[-1][1])]
    return v[:7] if len(v) >= 7 else None


def is_active(path, secs=180):
    """True if the file was written within `secs` -- i.e. a run still in progress.

    hist-etas-mk rewrites its .upda every time the simplex finds an improving
    solution, so the file exists and parses cleanly long before the search has
    converged. Reading it mid-run yields plausible but PROVISIONAL numbers.
    """
    try:
        return (datetime.datetime.now().timestamp() - os.path.getmtime(path)) < secs
    except OSError:
        return False


def hierarchical(run):
    """Baselines on lines 2-3 of hist-etas-mk.upda; header line 1 has w1 w2 ABIC m."""
    p = os.path.join(run, "s09", "hist-etas-mk.upda")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None, None
    if is_active(p):
        return "PROVISIONAL", None
    ls = open(p, encoding="utf-8", errors="replace").read().replace("\r", "").split("\n")
    v = [float(x) for x in re.findall(r"[-+]?\d*\.\d+E[-+]\d+", " ".join(ls[1:3]))]
    hdr = [float(x) for x in re.findall(r"[-+]?\d*\.\d+E[-+]\d+|\d+\.\d+|\b\d+\b", ls[0])]
    return (v[:7] if len(v) >= 7 else None), hdr


def nn_of(run):
    p = os.path.join(run, "s06", "delone2.out")
    if not os.path.exists(p):
        return None
    return int(open(p, encoding="utf-8").readline().split()[0])


def pair_counts(catalog, thresholds, dt_days=100, dkm=150):
    rows = [r for r in csv.DictReader(open(catalog, encoding="utf-8"))
            if r.get("mag") and (r.get("type") or "earthquake") == "earthquake"]
    out = {}
    for mm in thresholds:
        q = [r for r in rows if float(r["mag"]) >= mm]
        t = [datetime.datetime.fromisoformat(r["time"].replace("Z", "")) for r in q]
        o = np.argsort(t)
        t = [t[i] for i in o]
        x = [float(q[i]["longitude"]) for i in o]
        y = [float(q[i]["latitude"]) for i in o]
        m = [float(q[i]["mag"]) for i in o]
        c = tr = 0
        for i in range(len(t)):
            for j in range(i + 1, len(t)):
                d = (t[j] - t[i]).total_seconds() / 86400.0
                if d > dt_days:
                    break
                la = math.radians((y[i] + y[j]) / 2)
                if math.hypot((x[j] - x[i]) * 111.32 * math.cos(la),
                              (y[j] - y[i]) * 110.57) <= dkm:
                    c += 1
                    if m[j] < m[i]:
                        tr += 1
        out[mm] = {"n": len(t), "pairs": c, "triggered": tr}
    return out


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "threshold-summary")
    os.makedirs(outdir, exist_ok=True)

    def pack(lab, mm, run):
        st, st_se = stationary(run)
        pkg = stationary_st_etas(run)
        hi, hdr = hierarchical(run)
        status = "ok"
        if hi == "PROVISIONAL":
            hi, status = None, "provisional (run still in progress)"
        elif hi is None:
            status = "hierarchical stage not finished"
        return {"threshold": mm, "run": run, "nn": nn_of(run),
                "stationary": dict(zip(NAMES, st)) if st else None,
                "stationary_se": dict(zip(NAMES, st_se)) if st_se else None,
                "stationary_st_etas_defective": dict(zip(NAMES, pkg)) if pkg else None,
                "hierarchical": dict(zip(NAMES, hi)) if hi else None,
                "status": status, "mk_header": hdr}

    res = {lab: pack(lab, mm, run) for lab, mm, run in RUNS}
    res["Japan M>=6.0"] = pack("Japan M>=6.0", 6.0, JAPAN)

    cat = os.path.join(ROOT, "catalog", "usgs-mc.csv")
    res["pairs"] = pair_counts(cat, [6.0, 5.5, 5.0]) if os.path.exists(cat) else {}

    json.dump(res, open(os.path.join(outdir, "threshold-summary.json"), "w"), indent=2)

    # ---------------- figure: p and alpha against threshold ----------------
    th, ps, ph, as_, ah, trig, pse, ase = [], [], [], [], [], [], [], []
    for lab, mm, _ in RUNS:
        r = res[lab]
        if not r["stationary"]:
            continue
        th.append(mm)
        ps.append(r["stationary"]["p"]); as_.append(r["stationary"]["alpha"])
        pse.append(r["stationary_se"]["p"]); ase.append(r["stationary_se"]["alpha"])
        ph.append(r["hierarchical"]["p"] if r["hierarchical"] else np.nan)
        ah.append(r["hierarchical"]["alpha"] if r["hierarchical"] else np.nan)
        trig.append(res["pairs"].get(mm, {}).get("triggered", np.nan))

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.2))

    ax.errorbar(th, ps, yerr=pse, fmt="o-", color="#1f4e79", lw=1.9, ms=6, capsize=3,
                label="stationary (refit, ±1 s.e.)")
    ax.plot(th, ph, "s--", color="#c00000", lw=1.9, ms=6, label="hierarchical (hist-etas-mk)")
    if res["Japan M>=6.0"]["stationary"]:
        ax.errorbar([6.0], [res["Japan M>=6.0"]["stationary"]["p"]],
                    yerr=[res["Japan M>=6.0"]["stationary_se"]["p"]], fmt="o",
                    mfc="none", mec="#1f4e79", ecolor="#1f4e79", ms=10, mew=1.6, capsize=3)
        ax.plot([6.0], [res["Japan M>=6.0"]["hierarchical"]["p"]], "s",
                mfc="none", mec="#c00000", ms=10, mew=1.6)
        ax.annotate("Japan", (6.0, res["Japan M>=6.0"]["hierarchical"]["p"]),
                    textcoords="offset points", xytext=(9, -3), fontsize=8, color="#555")
        ax.set_ylim(top=1.08)
    ax.axhline(1.0, color="#888", ls=":", lw=1.2)
    ax.text(5.8, 1.008, "$p=1$", fontsize=8, color="#666")
    ax.set_xlabel("target magnitude threshold"); ax.set_ylabel("$p$")
    ax.invert_xaxis()
    ax.set_title("Omori exponent rises as the threshold falls", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    ax.grid(alpha=.25, lw=.4)

    ax2.errorbar(th, as_, yerr=ase, fmt="o-", color="#1f4e79", lw=1.9, ms=6, capsize=3,
                 label="stationary")
    ax2.plot(th, ah, "s--", color="#c00000", lw=1.9, ms=6, label="hierarchical")
    ax2.set_xlabel("target magnitude threshold"); ax2.set_ylabel(r"$\alpha$")
    ax2.invert_xaxis()
    # Neither alpha is monotonic: the refitted stationary one goes 0.77 -> 1.09
    # -> 0.95 (the monotonic 0.79 -> 0.92 -> 1.08 was st-etas's defective fit),
    # and the hierarchical baseline 1.44 -> 1.70 -> 0.96.
    ax2.set_title(r"$\alpha$: no monotonic trend in either fit", fontsize=10)
    ax2.legend(fontsize=8, frameon=False, loc="lower left")
    ax2.grid(alpha=.25, lw=.4)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fig-threshold-trend.pdf"))
    plt.close(fig)

    # ---------------- console table ----------------
    print(f"{'run':14}{'NN':>7}{'trig pairs':>12}{'p (stat)':>10}{'p (hier)':>10}"
          f"{'alpha (stat)':>14}{'alpha (hier)':>14}")
    for lab, mm, _ in RUNS:
        r = res[lab]
        if not r["stationary"]:
            print(f"{lab:14}{r['status']:>50}")
            continue
        h = r["hierarchical"]
        tp = res["pairs"].get(mm, {}).get("triggered", "")
        print(f"{lab:14}{r['nn'] or 0:7}{tp:>12}"
              f"{r['stationary']['p']:10.3f}"
              f"{(h['p'] if h else float('nan')):10.3f}"
              f"{r['stationary']['alpha']:14.3f}"
              f"{(h['alpha'] if h else float('nan')):14.3f}")
    r = res["Japan M>=6.0"]
    if r["stationary"]:
        h = r["hierarchical"]
        print(f"{'Japan 6.0':14}{r['nn'] or 0:7}{'':>12}"
              f"{r['stationary']['p']:10.3f}{(h['p'] if h else float('nan')):10.3f}"
              f"{r['stationary']['alpha']:14.3f}{(h['alpha'] if h else float('nan')):14.3f}")
    print(f"\nwrote {outdir}/fig-threshold-trend.pdf and threshold-summary.json")


if __name__ == "__main__":
    main()
