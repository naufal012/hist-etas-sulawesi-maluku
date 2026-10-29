#!/usr/bin/env python
"""
Set the fitted HIST-ETAS fields against independent tectonic geometry.

Two questions, both asked of the M>=5.0 run:

1. Depth. HIST-ETAS is two-dimensional, so an intraslab event and a shallow
   one at the same epicentre can form a triggering pair. Slab2 (Hayes et al.
   2018) gives the slab-top depth under each epicentre; each target event is
   labelled upper-plate, interface or intraslab against the nearest slab
   surface, and the candidate triggered pairs (the same 150 km / 100 day /
   later-is-smaller definition as threshold_summary.py) are split by whether
   the two events share a label.

2. Surface structure. Distances from every grid pixel to the nearest mapped
   trench (USGS 2017 traces shipped with Slab2) or active fault (GEM Global
   Active Faults, harmonised) are used to ask whether the fields concentrate
   on structures, and how much of that is just where the events are.

Caveat that governs the reading: Slab2 was built partly from the same USGS
hypocentres, and mapped structures partly from seismicity, so agreement of
mu with them is expected by construction. It is a consistency check, not an
independent validation.

Inputs   maps/tectonics/slab2/*_dep_*.xyz, *_clp_*.csv
         maps/tectonics/trenches_usgs_2017_depths.csv
         maps/tectonics/gem_active_faults_harmonized.geojson
         maps/tectonics/gvp_holocene_box.csv  (GVP Volcanoes of the World,
                                         snapshot 2026-05-03, via the R
                                         package `volcanoes`)
         <run>/s03/work.etas, <run>/s13/interp5pa.pixel
Outputs  <run>/tect/tectonic-summary.json
         <run>/tect/slab-<code>.grid     lon lat depth(km, positive down), box only
         <run>/tect/structures.txt       seg kind lon lat, box only
         <run>/tect/events-class.txt     per target event
         <run>/tect/covariates.pixel     per grid pixel: distances to fault,
                                         subduction zone and volcano; slab
                                         depth, dip, strike (upper slab)
         <run>/tect/volcanoes.txt        Holocene volcanoes in the box

Usage: python tectonic_context.py [--run DIR] [--mc 5.0]
"""
import argparse, csv, glob, json, math, os
from collections import Counter

import numpy as np
from scipy.spatial import cKDTree

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TECT = os.path.join(ROOT, "maps", "tectonics")
BOX = (118.0, 132.0, -6.0, 6.0)
PAD = 1.0          # structures just outside the box still set distances inside it

# Depth labels. TOL is the half-thickness of the band around the slab top
# counted as "on the slab"; ZINT is the base of the seismogenic interface.
# Both are conventional choices, stated in the paper, not fitted.
TOL, ZINT = 20.0, 60.0
NEAR_KM = 50.0     # "near a structure"
DZ_FAR = 50.0      # a pair this far apart in depth is not plausibly one fault system
SLABS = ["sul", "hal", "sum", "cot", "phi"]   # png has no nodes inside the box


def km_xy(lon, lat, lat0=0.0):
    return (np.asarray(lon) * 111.32 * math.cos(math.radians(lat0)),
            np.asarray(lat) * 110.57)


def load_slab(code, what="dep"):
    """Slab2 grid as (lon, lat, values). Depth is flipped to positive down."""
    f = glob.glob(os.path.join(TECT, "slab2", f"{code}_slab2_{what}_*.xyz"))[0]
    a = np.loadtxt(f, delimiter=",")
    a[:, 0] = np.where(a[:, 0] > 180, a[:, 0] - 360, a[:, 0])
    ux, uy = np.unique(a[:, 0]), np.unique(a[:, 1])
    z = np.full((len(uy), len(ux)), np.nan)
    v = -a[:, 2] if what == "dep" else a[:, 2]
    z[np.searchsorted(uy, a[:, 1]), np.searchsorted(ux, a[:, 0])] = v
    return ux, uy, z


def nearest(grid, lon, lat):
    """Nearest-node value; used for strike, which must not be averaged
    across the 0/360 wrap."""
    ux, uy, z = grid
    i = np.clip(np.rint((lon - ux[0]) / (ux[1] - ux[0])).astype(int), 0, len(ux) - 1)
    j = np.clip(np.rint((lat - uy[0]) / (uy[1] - uy[0])).astype(int), 0, len(uy) - 1)
    out = z[j, i]
    off = (lon < ux[0]) | (lon > ux[-1]) | (lat < uy[0]) | (lat > uy[-1])
    out[off] = np.nan
    return out


def sample(grid, lon, lat):
    """Bilinear slab-top depth (km, positive down); NaN off the model."""
    ux, uy, z = grid
    i = np.searchsorted(ux, lon) - 1
    j = np.searchsorted(uy, lat) - 1
    out = np.full(len(lon), np.nan)
    ok = (i >= 0) & (i < len(ux) - 1) & (j >= 0) & (j < len(uy) - 1)
    i, j = i[ok], j[ok]
    tx = (lon[ok] - ux[i]) / (ux[i + 1] - ux[i])
    ty = (lat[ok] - uy[j]) / (uy[j + 1] - uy[j])
    out[ok] = ((1 - tx) * (1 - ty) * z[j, i] + tx * (1 - ty) * z[j, i + 1]
               + (1 - tx) * ty * z[j + 1, i] + tx * ty * z[j + 1, i + 1])
    return out


def classify(depth, slab_depths):
    """slab_depths: (n, nslab). Label against the nearest slab surface."""
    lab = np.empty(len(depth), dtype=object)
    zs = np.full(len(depth), np.nan)
    any_slab = np.isfinite(slab_depths).any(axis=1)
    diff = np.abs(slab_depths - depth[:, None])
    diff[~np.isfinite(diff)] = np.inf
    k = diff.argmin(axis=1)
    zs[any_slab] = slab_depths[any_slab, k[any_slab]]
    for n in range(len(depth)):
        z = depth[n]
        if not any_slab[n]:
            lab[n] = "no slab"
        elif z < zs[n] - TOL:
            lab[n] = "upper plate"
        elif z <= ZINT:
            lab[n] = "interface"
        else:
            lab[n] = "intraslab"
    nstack = np.isfinite(slab_depths).sum(axis=1)
    return lab, zs, nstack


def in_box(lon, lat, pad=0.0):
    return (BOX[0] - pad <= lon <= BOX[1] + pad) and (BOX[2] - pad <= lat <= BOX[3] + pad)


def structures():
    """List of (kind, Nx2 lon/lat) polylines touching the padded box."""
    segs = []
    rows = list(csv.DictReader(open(os.path.join(TECT, "trenches_usgs_2017_depths.csv"))))
    by = {}
    for r in rows:
        lon = float(r["lon"]); lon = lon - 360 if lon > 180 else lon
        by.setdefault((r["slab"], r["bound"]), []).append((lon, float(r["lat"])))
    for key, pts in by.items():
        pts = np.array(pts)
        if any(in_box(x, y, PAD) for x, y in pts):
            # the file lists each trench as consecutive points; break on jumps
            cut = np.where(np.hypot(*np.diff(pts, axis=0).T) > 0.5)[0] + 1
            for s in np.split(pts, cut):
                if len(s) > 1:
                    segs.append(("trench", s))
    g = json.load(open(os.path.join(TECT, "gem_active_faults_harmonized.geojson"),
                       encoding="utf-8"))
    for f in g["features"]:
        geom = f["geometry"]
        lines = [geom["coordinates"]] if geom["type"] == "LineString" else geom["coordinates"]
        kind = "subduction" if "Subduction" in (f["properties"].get("slip_type") or "") else "fault"
        for ln in lines:
            s = np.array([c[:2] for c in ln])
            if any(in_box(x, y, PAD) for x, y in s):
                segs.append((kind, s))
    return segs


def densify(s, step_km=2.0):
    out = [s[0]]
    for a, b in zip(s[:-1], s[1:]):
        d = math.hypot((b[0] - a[0]) * 111.32, (b[1] - a[1]) * 110.57)
        n = max(1, int(d / step_km))
        out += [a + (b - a) * t for t in np.linspace(0, 1, n + 1)[1:]]
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=os.path.join(ROOT, "run-mc50"))
    ap.add_argument("--mc", type=float, default=5.0)
    a = ap.parse_args()
    out = os.path.join(a.run, "tect")
    os.makedirs(out, exist_ok=True)
    res = {"tol_km": TOL, "zint_km": ZINT, "near_km": NEAR_KM, "dz_far_km": DZ_FAR}

    # ---- target events ------------------------------------------------
    w = np.loadtxt(os.path.join(a.run, "s03", "work.etas"), skiprows=1, usecols=(1, 2, 3, 4, 5))
    w = w[w[:, 2] >= a.mc - 0.05]
    lon, lat, mag, t, dep = w[:, 0], w[:, 1], w[:, 2], w[:, 3], -w[:, 4]
    res["n_target"] = int(len(w))

    # ---- slabs --------------------------------------------------------
    grids = {c: load_slab(c) for c in SLABS}
    sd = np.column_stack([sample(grids[c], lon, lat) for c in SLABS])
    lab, zs, nstack = classify(dep, sd)
    cnt = Counter(lab)
    res["classes"] = {k: int(v) for k, v in cnt.items()}
    res["stacked_epicentres"] = int((nstack >= 2).sum())
    res["depth_fixed_10_or_33_35"] = int(np.isin(np.round(dep, 1), [10.0, 33.0, 35.0]).sum())
    near_slab = {c: int(((np.abs(sd[:, i] - dep) <= TOL)).sum()) for i, c in enumerate(SLABS)}
    res["events_on_slab_by_model"] = near_slab
    with open(os.path.join(out, "events-class.txt"), "w") as fh:
        fh.write("lon lat mag day depth slabtop nslab class\n")
        for r in zip(lon, lat, mag, t, dep, zs, nstack, lab):
            fh.write("%.4f %.4f %.1f %.5f %.1f %.1f %d %s\n" % (*r[:7], r[7].replace(" ", "_")))

    for c in SLABS:
        ux, uy, z = grids[c]
        ix = (ux >= BOX[0] - PAD) & (ux <= BOX[1] + PAD)
        iy = (uy >= BOX[2] - PAD) & (uy <= BOX[3] + PAD)
        with open(os.path.join(out, f"slab-{c}.grid"), "w") as fh:
            fh.write("%d %d\n" % (ix.sum(), iy.sum()))
            for yy, row in zip(uy[iy], z[np.ix_(iy, ix)]):
                for xx, zz in zip(ux[ix], row):
                    fh.write("%.2f %.2f %s\n" % (xx, yy, "NA" if np.isnan(zz) else "%.1f" % zz))

    # ---- pairs ---------------------------------------------------------
    o = np.argsort(t)
    lon, lat, mag, t, dep, lab = lon[o], lat[o], mag[o], t[o], dep[o], lab[o]
    x, y = km_xy(lon, lat)
    pairs = []
    for i in range(len(t)):
        for j in range(i + 1, len(t)):
            if t[j] - t[i] > 100.0:
                break
            la = math.radians((lat[i] + lat[j]) / 2)
            if (math.hypot((lon[j] - lon[i]) * 111.32 * math.cos(la),
                           (lat[j] - lat[i]) * 110.57) <= 150.0 and mag[j] < mag[i]):
                pairs.append((i, j))
    pi, pj = np.array(pairs).T
    dz = np.abs(dep[pi] - dep[pj])
    mixed = lab[pi] != lab[pj]
    far = dz > DZ_FAR
    mol = ((lon[pi] >= 124) & (lon[pi] <= 128) & (lat[pi] >= -1) & (lat[pi] <= 5))
    res["pairs"] = {
        "n": int(len(pi)),
        "mixed_class": int(mixed.sum()),
        "dz_gt_far": int(far.sum()),
        "dz_gt_far_molucca": int((far & mol).sum()),
        "n_molucca": int(mol.sum()),
        "dz_median": float(np.median(dz)),
        "intraslab_parent_shallow_child": int(((lab[pi] == "intraslab") & (dep[pj] <= ZINT) & far).sum()),
        "shallow_parent_intraslab_child": int(((dep[pi] <= ZINT) & (lab[pj] == "intraslab") & far).sum()),
        "by_class_pair": {f"{p}->{c}": int(v) for (p, c), v in
                          Counter(zip(lab[pi], lab[pj])).most_common()},
    }

    # ---- surface structures -------------------------------------------
    segs = structures()
    with open(os.path.join(out, "structures.txt"), "w") as fh:
        fh.write("seg kind lon lat\n")
        for k, (kind, s) in enumerate(segs):
            for xx, yy in s:
                fh.write("%d %s %.4f %.4f\n" % (k, kind, xx, yy))
    res["structures"] = dict(Counter(k for k, _ in segs))
    P = np.vstack([densify(s) for _, s in segs])
    tree = cKDTree(np.column_stack(km_xy(P[:, 0], P[:, 1])))

    px = np.loadtxt(os.path.join(a.run, "s13", "interp5pa.pixel"))
    px = px[px[:, 2] > -998]
    dist = tree.query(np.column_stack(km_xy(px[:, 0], px[:, 1])))[0]
    near = dist <= NEAR_KM
    mu, K0, al, p = px[:, 2], px[:, 3], px[:, 4], px[:, 5]
    edist = tree.query(np.column_stack(km_xy(lon, lat)))[0]

    def top(v, q=0.9):
        return v >= np.quantile(v, q)
    # Pixels are equal in degrees, so cos(lat) weights convert to area.
    wA = np.cos(np.radians(px[:, 1]))
    res["fields_vs_structure"] = {
        "area_share_near": float(wA[near].sum() / wA.sum()),
        "events_share_near": float((edist <= NEAR_KM).mean()),
        "mu_integral_share_near": float((mu * wA)[near].sum() / (mu * wA).sum()),
        "mu_top10_share_near": float(near[top(mu)].mean()),
        "K0_top10_share_near": float(near[top(K0)].mean()),
        "alpha_top10_share_near": float(near[top(al)].mean()),
        "p_median_near": float(np.median(p[near])),
        "p_median_far": float(np.median(p[~near])),
        "K0_median_near": float(np.median(K0[near])),
        "K0_median_far": float(np.median(K0[~near])),
    }

    # Where the K0 and mu maxima sit, by slab context of the pixel.
    psd = np.column_stack([sample(grids[c], px[:, 0], px[:, 1]) for c in SLABS])
    def ctx(mask):
        c = Counter()
        for row in psd[mask]:
            have = [SLABS[i] for i in np.where(np.isfinite(row))[0]]
            c["+".join(have) or "none"] += 1
        return {k: round(v / mask.sum(), 3) for k, v in c.most_common()}
    res["top10_pixel_slab_context"] = {"mu": ctx(top(mu)), "K0": ctx(top(K0)),
                                       "all": ctx(np.ones(len(px), bool))}

    # Field value at each event's pixel, summarised by depth label.
    ptree = cKDTree(px[:, :2])
    k = ptree.query(np.column_stack([lon, lat]))[1]
    res["field_at_events_by_class"] = {
        c: {"n": int((lab == c).sum()),
            "K0_median": float(np.median(K0[k][lab == c])),
            "mu_median": float(np.median(mu[k][lab == c])),
            "p_median": float(np.median(p[k][lab == c]))}
        for c in cnt}

    # ---- covariate rasters ---------------------------------------------
    # Every pixel of the evaluation grid, inside the hull or not. Where slabs
    # are stacked the upper (shallowest) surface is reported, with its own
    # dip and strike, so the three slab panels describe one surface.
    allpx = np.loadtxt(os.path.join(a.run, "s13", "interp5pa.pixel"))
    gx, gy = allpx[:, 0], allpx[:, 1]
    inhull = allpx[:, 2] > -998
    gxy = np.column_stack(km_xy(gx, gy))
    kinds = np.array([k for k, s_ in segs for _ in densify(s_)])
    sub = np.isin(kinds, ["trench", "subduction"])
    d_fault = cKDTree(np.column_stack(km_xy(P[~sub, 0], P[~sub, 1]))).query(gxy)[0]
    d_sub = cKDTree(np.column_stack(km_xy(P[sub, 0], P[sub, 1]))).query(gxy)[0]
    vol = list(csv.DictReader(open(os.path.join(TECT, "gvp_holocene_box.csv"), encoding="utf-8")))
    V = np.array([[float(v["lon"]), float(v["lat"])] for v in vol])
    d_volc = cKDTree(np.column_stack(km_xy(V[:, 0], V[:, 1]))).query(gxy)[0]
    with open(os.path.join(out, "volcanoes.txt"), "w", encoding="utf-8") as fh:
        fh.write("lon lat name\n")
        for v in vol:
            fh.write("%s %s %s\n" % (v["lon"], v["lat"], v["name"].replace(" ", "_")))

    dep_all = np.column_stack([sample(grids[c], gx, gy) for c in SLABS])
    dip_all = np.column_stack([sample(load_slab(c, "dip"), gx, gy) for c in SLABS])
    str_all = np.column_stack([nearest(load_slab(c, "str"), gx, gy) for c in SLABS])
    has = np.isfinite(dep_all).any(axis=1)
    up = np.where(has, np.nanargmin(np.where(np.isfinite(dep_all), dep_all, np.inf), axis=1), 0)
    r = np.arange(len(gx))
    s_dep = np.where(has, dep_all[r, up], np.nan)
    s_dip = np.where(has, dip_all[r, up], np.nan)
    s_str = np.where(has, str_all[r, up], np.nan)

    # GSRM v2.1 (Kreemer et al. 2014): second invariant of the horizontal
    # strain-rate tensor, sqrt(exx^2 + eyy^2 + 2 exy^2), in 1e-9/yr, on
    # 0.25 x 0.2 deg cells. Only deforming cells are listed; a pixel with no
    # listed cell lies on a plate the model treats as rigid, and gets 0.
    gs = np.loadtxt(os.path.join(TECT, "gsrm", "GSRM_average_strain_v2.1.txt"), comments="#")
    gs = gs[(gs[:, 1] >= BOX[0] - PAD) & (gs[:, 1] <= BOX[1] + PAD)
            & (gs[:, 0] >= BOX[2] - PAD) & (gs[:, 0] <= BOX[3] + PAD)]
    e2 = np.sqrt(gs[:, 2] ** 2 + gs[:, 3] ** 2 + 2 * gs[:, 4] ** 2)
    dd, ii = cKDTree(np.column_stack([gs[:, 1] / 0.25, gs[:, 0] / 0.2])).query(
        np.column_stack([gx / 0.25, gy / 0.2]), p=np.inf)
    strain = np.where(dd <= 0.5 + 1e-6, e2[ii], 0.0)
    res["gsrm_rigid_pixel_share"] = float((strain == 0).mean())

    # CRUST1.0 (Laske et al. 2013), 1 x 1 deg cells from 89.5N, 179.5W with
    # longitude the inner loop. Crustal thickness as getCN1maps.f defines it:
    # bottom of water to Moho, i.e. ice + sediments + crystalline crust.
    cb = np.loadtxt(os.path.join(TECT, "crust1", "crust1.bnds")).reshape(180, 360, 9)
    row = np.clip(np.floor(90.0 - gy).astype(int), 0, 179)
    col = np.clip(np.floor(gx + 180.0).astype(int), 0, 359)
    crust = cb[row, col, 1] - cb[row, col, 8]

    np.savetxt(os.path.join(out, "covariates.pixel"),
               np.column_stack([gx, gy, inhull, d_fault, d_sub, d_volc, s_dep, s_dip, s_str,
                                strain, crust]),
               fmt=["%.4f", "%.4f", "%d", "%.1f", "%.1f", "%.1f", "%.1f", "%.1f", "%.1f",
                    "%.1f", "%.2f"],
               header="lon lat inhull d_fault d_sub d_volc slab_dep slab_dip slab_str"
                      " strain crust_thk",
               comments="")
    res["n_volcanoes"] = len(vol)

    # Rank correlation of each field with each covariate, inside the hull.
    # Strike is circular and left out.
    from scipy.stats import spearmanr
    F = {"log10_mu": np.log10(allpx[:, 2]), "log10_K0": np.log10(allpx[:, 3]),
         "alpha": allpx[:, 4], "p": allpx[:, 5]}
    C = {"d_fault": d_fault, "d_sub": d_sub, "d_volc": d_volc,
         "slab_dep": s_dep, "slab_dip": s_dip,
         "strain": strain, "crust_thk": crust}
    res["spearman_in_hull"] = {}
    for fk, fv in F.items():
        row = {}
        for ck, cv in C.items():
            m = inhull & np.isfinite(cv)
            row[ck] = round(float(spearmanr(fv[m], cv[m])[0]), 3)
        res["spearman_in_hull"][fk] = row
    res["slab_pixels_in_hull_share"] = float((inhull & has).sum() / inhull.sum())

    json.dump(res, open(os.path.join(out, "tectonic-summary.json"), "w"), indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
