#!/usr/bin/env python
"""
Derive the stage-13 conditional-intensity program from the package original.

Source: HIST-PPM-V3/programs/Section13.3files5pa/intensityday2magnod5pa.f

The package program evaluates the HIST-ETAS conditional intensity on a pixel
grid at one instant.  It is written for the package's *forecast* setting: the
model is fitted to a catalogue spanning 0..TZ, a SECOND catalogue covering the
period after TZ is passed in through work.para, and the evaluation instant
zzz = TZ + tday1 always lies after every event in the fitted set.  Under that
assumption the fitted events need no causality test, and the program has none.

Here there is no second catalogue: the evaluation instant lies inside the
fitted span, so the test has to exist.  Four edits, all reversible, none of
which changes the arithmetic when the original's assumption does hold:

  1. CAUSALITY  -- skip fitted events with t_j >= t in the triggering sum.
     When t >= max(t_j) this skips nothing and the program is unchanged.
     Without it, events in the future of the evaluation instant contribute,
     which silently inflates lambda and can raise (t - t_j + c) to a negative
     power of a negative base.

  2. MEMORY     -- the fixed array bounds are cut to this region's size.
     The original declares common /jfun/ jxy(22345,22345): 2.0 GB of BSS for
     an array no statement in the file ever references.  It is dropped, and
     npa/nda/mdl are cut from 22345/23333 to values checked against
     delone2.out's header at run time by run_pipeline.sh.

  3. DEAD STORE -- xl(50000,60), 24 MB, likewise never referenced.

  4. DEBUG WRITE -- the "ramdai > 0.9999" diagnostic prints xoi(j) and xmin,
     neither of which is defined at that point (j is never assigned in the
     enclosing scope).  It is a latent out-of-bounds read that the Japan run
     never triggers because its rates stay below the threshold; ours exceed
     it near the 2026 M7.8.  Removed.

Usage:  python patch_intensity.py <package.f> <out.f> [--npa N] [--mdl N]
"""
import argparse
import os
import re
import sys


EDITS = []


def edit(name):
    def deco(fn):
        EDITS.append((name, fn))
        return fn
    return deco


@edit("causality guard on fitted events")
def _causality(src, o):
    # loop 120 walks l = 1..nn+n1; l <= nn is the fitted (delone2) set, whose
    # times live in zz(l).  The original tests only the work.para branch.
    old = ("      do 120 l=1,nn+n1\n"
           "      if(l.le.nn) then\n"
           "      xxi=xx(l)\n")
    new = ("      do 120 l=1,nn+n1\n"
           "      if(l.le.nn) then\n"
           "c     [patch] causality: a fitted event at or after the evaluation\n"
           "c     instant cannot have triggered anything observed at it.\n"
           "      if(zz(l).ge.zzz) go to 120\n"
           "      xxi=xx(l)\n")
    return src.replace(old, new, 1), src.count(old)


@edit("array bounds")
def _bounds(src, o):
    n = 0
    src, k = re.subn(r"parameter\(ldat=\d+, *npa=\d+, *mpa=npa\*5,mdl=\d+\)",
                     f"parameter(ldat={o.npa}, npa={o.npa}, mpa=npa*5,mdl={o.mdl})",
                     src)
    n += k
    src, k = re.subn(r"parameter\(npa=\d+,nda=\d+,mdl=\d+,nr=(\d+)\)",
                     rf"parameter(npa={o.npa},nda={o.npa},mdl={o.mdl},nr=\1)",
                     src)
    n += k
    return src, n


@edit("drop unreferenced common /jfun/ (2.0 GB)")
def _jfun(src, o):
    return re.subn(r"^      common /jfun/rmax\(nda\),jxy\(nda,nda\),jnn\(nda\)$",
                   "c     [patch] dropped: common /jfun/ is never referenced.",
                   src, flags=re.M)


@edit("drop unreferenced xl(50000,60) (24 MB)")
def _xl(src, o):
    return re.subn(r"^      dimension xl\(50000,60\)$",
                   "c     [patch] dropped: xl() is never referenced.",
                   src, flags=re.M)


@edit("remove debug write with undefined j / xmin")
def _debug(src, o):
    old = ("      if(ramdai.gt.0.9999) then\n"
           "      write(6,*) jj,l1,i,ramdai,rmdi,b4i\n"
           "      write(6,*) xxx,xoi(j)+xmin,yyy+ymin,zzz\n"
           "      endif\n")
    new = ("c     [patch] removed: printed xoi(j) with j undefined here.\n")
    return src.replace(old, new, 1), src.count(old)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source")
    ap.add_argument("out")
    ap.add_argument("--npa", type=int, default=8000,
                    help="max vertices/events (must exceed npext in delone2.out)")
    ap.add_argument("--mdl", type=int, default=16000,
                    help="max triangles (must exceed ndelo in delone2.out)")
    o = ap.parse_args()

    src = open(o.source, encoding="utf-8", errors="replace").read()
    total, failed = 0, []
    for name, fn in EDITS:
        src, k = fn(src, o)
        total += k
        if k == 0:
            failed.append(name)
        print(f"   [{k}] {name}")

    if failed:
        print("FATAL: these edits matched nothing -- the package source has "
              "changed and the patch is no longer valid:", file=sys.stderr)
        for f in failed:
            print(f"   - {f}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(os.path.abspath(o.out)), exist_ok=True)
    open(o.out, "w", encoding="utf-8", newline="\n").write(src)
    print(f"   -> {o.out}  ({total} edits, npa={o.npa} mdl={o.mdl})")


if __name__ == "__main__":
    main()
