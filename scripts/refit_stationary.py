"""Refit the stationary ETAS model for one or more run directories.

The package's st-etas integrates the kernel with events at absolute lon/lat
against a box [0,TX]x[0,TY] (see README, stage 14), so its MLE is not a
stationary fit.  This refits the seven constants under the likelihood the
hierarchical programs use, all fields held constant, and writes
<run>/s14/stationary.json -- the same cache residuals.py reads.

Standard errors are from the observed information: a central-difference
Hessian of -log L in log-parameters, mapped back by the delta method.

    python scripts/refit_stationary.py run run-mc55 run-mc50 ../run
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import etas_lik as E  # noqa: E402
from residuals import fit_stationary  # noqa: E402

NAMES = ["mu", "K", "c", "alpha", "p", "d", "q"]


def read_st_mle(path):
    """The (broken) st-etas MLE, for the side-by-side table only."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    k = max(i for i, l in enumerate(lines) if l.strip().startswith("mle ="))
    return [float(v) for v in (lines[k].split("=")[1] + " " + lines[k + 1]).split()]


def hessian(fun, z, h=1e-3):
    n = len(z)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            e_i = np.eye(n)[i] * h
            e_j = np.eye(n)[j] * h
            H[i, j] = H[j, i] = (fun(z + e_i + e_j) - fun(z + e_i - e_j)
                                 - fun(z - e_i + e_j) + fun(z - e_i - e_j)) / (4 * h * h)
    return H


def refit(run):
    conf = E.read_conf(os.path.join(run, "s10", "hist-etas5pa.conf"))
    cat = E.read_delone2(os.path.join(run, "s10", "delone2.out"), conf)
    pr = E.neighbour_pairs(cat, conf["bi2"])
    mk_path = next(p for p in (os.path.join(run, "s09", "hist-etas-mk.upda"),
                               os.path.join(run, "s10", "hist-etas-mk.upda"))
                   if os.path.exists(p))
    mk = E.read_mk(mk_path)
    theta0 = [mk.mu.mean(), np.exp(np.log(mk.K).mean()), mk.c, 2 * mk.a[0],
              mk.p[0], mk.d, mk.q[0]]
    par, ll = fit_stationary(cat, pr, theta0)

    def nll(z):
        return -E.loglik(cat, E.constant_fields(np.exp(z), cat.n), pr, "fortran")["loglik"]

    z = np.log(list(par.values()))
    H = hessian(nll, z)
    cov = np.linalg.inv(H)
    se = np.sqrt(np.diag(cov)) * np.exp(z)          # delta method, log -> linear
    out = dict(theta=par, se=dict(zip(NAMES, se.tolist())), loglik=ll,
               aic=2 * 7 - 2 * ll, n_target=int(cat.target().sum()), n_events=cat.n,
               mk_upda=os.path.relpath(mk_path, run),
               hessian_min_eig=float(np.linalg.eigvalsh(H).min()))
    st_prt = os.path.join(run, "s05", "st-etas.prt")
    if os.path.exists(st_prt):
        out["st_etas_package_mle"] = dict(zip(NAMES, read_st_mle(st_prt)))
    os.makedirs(os.path.join(run, "s14"), exist_ok=True)
    with open(os.path.join(run, "s14", "stationary.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", help="run directories, relative to the project root")
    a = ap.parse_args()
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for r in a.runs:
        run = os.path.normpath(os.path.join(root, r))
        o = refit(run)
        t, s = o["theta"], o["se"]
        pkg = o.get("st_etas_package_mle", {})
        print(f"{r}: N={o['n_events']} target={o['n_target']}  logL={o['loglik']:.2f}  "
              f"minEig={o['hessian_min_eig']:.3g}")
        for k in NAMES:
            print(f"   {k:6s} {t[k]:.5g} ± {s[k]:.2g}" + (f"   (st-etas {pkg[k]:.5g})" if pkg else ""))


if __name__ == "__main__":
    main()
