#!/usr/bin/env python
"""
Generate every HIST-PPM stage config from region.conf + state.conf.

The pipeline learns three things only by running:
  TZ  -- total time span in days, on tsei2etas's own time axis
  NP  -- events surviving the magnitude threshold   (etas2aniso.out3 header)
  NN  -- target events after tessellation           (delone2.out header)
and, for the hierarchical stages, the st-etas MLE used as initial theta.

Those land in state.conf as the pipeline progresses; this script re-renders
the configs each time so NP/NN/TZ never drift out of sync between stages.

Stage 13 (conditional intensity and parameter rasterisation) adds three more
generated files -- node.conf, the *-forecast.conf origin files, and time.conf.
time.conf is the one to watch: the package's intensity program evaluates at
zzz = TZ + tday1, so an instant SNAP days after the start of the catalogue is
requested as tday1 = SNAP - TZ, which is normally negative here. That is the
whole reason stage 13 needs the causality patch in scripts/patch_intensity.py.
That drift is the single easiest way to get plausible-looking wrong output --
the package's own shipped Japan config carries the author's uncertainty about
it as a comment ("delone1 redued 310 to 308?").

Usage:  python gen_configs.py <region.conf> <state.conf> <outdir>
"""
import sys, os


def load(path):
    d = {}
    if not os.path.exists(path):
        return d
    for line in open(path, encoding="utf-8"):
        line = line.split("#")[0].strip()
        if line and "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def main():
    region, state, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    c = load(region)
    c.update(load(state))          # state overrides / supplies derived values

    xmin = float(c["XMIN"]); ymin = float(c["YMIN"])
    tx = float(c["TX"]); ty = float(c["TY"])
    mc = float(c["MC"]); clms = float(c["CLMS"]); xmg0 = float(c["XMG0"])
    zmin = float(c["ZMIN"]); bi2 = float(c["BI2"])
    tstart = float(c["TSTART"])
    tz = float(c.get("TZ", 0.0))
    np_ = int(float(c.get("NP", 0)))     # events after thresholding
    nn = int(float(c.get("NN", 0)))      # target events after tessellation
    theta = c.get("THETA_MK", c.get("THETA0")).split()
    delo_w1 = c.get("DELO_W1", "6.0d0")
    mc_utsu = mc - 0.05                  # Utsu correction for 0.1-rounded mags

    os.makedirs(outdir, exist_ok=True)
    w = lambda name, text: open(os.path.join(outdir, name), "w",
                                encoding="utf-8", newline="\n").write(text)

    # ---- stage 4 : anisotropic cluster identification -----------------
    w("etas2aniso.conf",
      f"./work.etas\n"
      f"{clms} {mc}           !clms cutm\n"
      f"1.                !one day span\n"
      f"read till above\n"
      f"0.0416666666667   !one hour span\n")

    # ---- stage 5 : stationary space-time ETAS -------------------------
    # NOTE: st-etas takes SEVEN fields on the region line (xmg1 is present);
    # the hist-etas configs take six. Do not copy one to the other.
    w("st-etas.conf",
      f"./etas2aniso.out3\n"
      f"17                        ! nfunct\n"
      f"{tx}  {ty}  {tz}  {np_} ! tx,ty,tz,nn\n"
      f"{xmin} {ymin}  {xmg0}  {zmin}  0.0  {tstart}  {bi2} "
      f"!xmin,ymin,xmg0,zmin,xmg1,tstar,bi2\n"
      f"7                        ! #parameters\n"
      f"{'   '.join(theta[:5])}\n"
      f"{'   '.join(theta[5:7])}\n"
      f"7 1  !ipr,igrd\n")

    # ---- stage 6 : Delaunay tessellation ------------------------------
    w("delone1.conf",
      f"   1.00E-15         EPS\n"
      f"   1000  7000       NEF0  NRG0\n"
      f" {xmin} {ymin}        xmin,ymin\n"
      f"  {tx}  {ty}          BXLX   BXLY\n"
      f"   {np_}             NP\n")

    # ---- stages 7 & 8 : b-value and background rate -------------------
    # NB: no trailing "." after xmin/ymin here. The package's Japan config reads
    # "128. 30. 5.95", where the dot is the decimal point of an integer-looking
    # literal -- not a separator. xmin is already formatted as "118.0", so
    # appending one yields "118.0." and the Fortran list-directed read crashes.
    for nm in ("delo2d-bvalues.conf", "delo2d-poisson.conf"):
        w(nm, f"{xmin} {ymin} {mc_utsu}\n{delo_w1}   !w1\n7  ! ipr\n")

    # ---- stage 9 : HIST-ETAS, varying mu and K0 -----------------------
    # init=0 (cold start): the .upda filename is read but never opened.
    w("hist-etas-mk.conf",
      f"delone2.out\n"
      f"{tx}  {ty}  {tz}  {nn}\n"
      f"{xmin}  {ymin}  {xmg0}   {zmin} {tstart}  {bi2}\n"
      f"0\n0\n1\n"
      f"hist-etas-mk.upda\n"
      f"0. 1. 1. \n"
      f"7\n"
      f"{' '.join(theta[:7])}\n"
      f"0  ipr(linear) 7 for printing the linear search\n"
      f"0 1.d0 1.d-0   nhesapp,dist,eps(simplex)\n")

    # ---- stage 10 : HIST-ETAS, five varying fields --------------------
    # Reads hist-etas-mk.upda for its starting field.
    w("hist-etas5pa.conf",
      f"delone2.out\n"
      f"{tx}      {ty}     {tz}   {nn}\n"
      f"{xmin}     {ymin}      {xmg0}     {zmin}  {tstart}     {bi2}\n"
      f"0\n0\n1\n"
      f"hist-etas-mk.upda\n"
      f"0. 1000.\n"
      f"10. 100. 1000.\n"
      f"0               !ipr  if ipr>0, print linear search procedure.\n"
      f"0.1d-3  0.1d-3  !tau1,tau2(davidn)\n"
      f"0.1d-3  0.1d-3  !eps1,eps2(davidn)\n"
      f"0 1.d0 0.5d-0    !nhesapp,dist,eps(simplex)\n")

    # ---- stage 13 : conditional intensity + parameter rasterisation ---
    # node.f writes one grid point per pixel CENTRE; interp5pa.f uses the
    # same convention so the two rasters align cell for cell.
    nx = int(float(c.get("NODE_NX", 200)))
    ny = int(float(c.get("NODE_NY", 200)))
    box = c.get("NODE_BOX", "").split()
    if len(box) == 4:
        glon1, glon2, glat1, glat2 = (float(v) for v in box)
    else:
        glon1, glon2 = xmin, xmin + tx
        glat1, glat2 = ymin, ymin + ty
    w("node.conf", f"""{glon1} {glon2} {glat1} {glat2}   !lon1 lon2 lat1 lat2
{nx} {ny}             !pixels in x and y
""")

    # Both intensity programs read tx, ty, tz then the region origin from
    # their own *-forecast.conf. The values are the same as the fitting
    # configs; they are separate files only because the package kept them so.
    for nm in ("hist-etas5pa-forecast.conf", "hist-etas-mk-forecast.conf"):
        w(nm, f"""{tx}   {ty}   {tz}
{xmin}   {ymin}
""")

    # Evaluation instant. SNAP is days from the start of the catalogue; the
    # program wants the offset from the end of the fitted span.
    snap = float(c.get("SNAP_DAYS", tz))
    w("time.conf", f"""{snap - tz:.5f}   ! zzz = TZ + this = day {snap:.5f}
""")

    print(f"configs written to {outdir}  (TZ={tz} NP={np_} NN={nn})")


if __name__ == "__main__":
    main()
