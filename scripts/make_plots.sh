#!/usr/bin/env bash
# Generate the figures for this region.
#
# Uses the package's own R plotting scripts unmodified, but substitutes the
# world coastline for the Japan-only one (see scripts/drawmap.R).
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PKG="D:/Documents/Hierarchical ETAS/HIST-PPM-V3/programs"
RSCRIPT="/c/Program Files/R/R-4.5.1/bin/Rscript.exe"

REGION="$ROOT/region.conf"; RUN="$ROOT/run"; FIGS="$ROOT/figs"
while [ $# -gt 0 ]; do
  case "$1" in
    --region) REGION="$2"; shift 2;;
    --rundir) RUN="$2";    shift 2;;
    --figs)   FIGS="$2";   shift 2;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done
mkdir -p "$FIGS"

prep(){ # prep <dir> -- drop in the shared R helpers and map data
  mkdir -p "$1/maps"
  cp "$ROOT/scripts/drawmap.R" "$1/drawmap.R"
  cp "$PKG/Section7files/f2.R" "$1/f2.R"
  cp "$ROOT/maps/w1.pp.data" "$ROOT/maps/w1.br.data" "$1/maps/"
}

ok=0; fail=0
plot(){ # plot <dir> <script> <expected.pdf> <fig name>
  prep "$1"
  # The package's R scripts hardcode Japan's map origin (128, 30); delone2.out
  # holds coordinates relative to (XMIN, YMIN), so without this the values are
  # right but drawn at Japanese coordinates.
  python "$ROOT/scripts/localise_R.py" "$1/$2" "$REGION"
  ( cd "$1" && "$RSCRIPT" "$2" > "R-$2.log" 2>&1 )
  if [ -s "$1/$3" ]; then
    cp "$1/$3" "$FIGS/$4"; echo "   ok   $4"; ok=$((ok+1))
  else
    echo "   FAIL $4  (see $1/R-$2.log)" >&2; fail=$((fail+1))
  fi
}

# --- catalogue ---------------------------------------------------------
# Data exploration, before any model: epicentres, magnitude-time, cumulative
# count. Reads work.etas rather than the raw CSV so that what is plotted is
# exactly what tsei2etas handed to the rest of the pipeline.
MC=$(awk -F= '/^ *MC *=/{sub(/#.*/,"",$2); gsub(/ /,"",$2); print $2}' "$REGION" | tail -1)
TSTART=$(awk -F= '/^ *TSTART *=/{sub(/#.*/,"",$2); gsub(/ /,"",$2); print $2}' "$REGION" | tail -1)
CSTART=$(awk -F= '/^ *CATALOG_START *=/{sub(/#.*/,"",$2); gsub(/ /,"",$2); print $2}' "$REGION" | tail -1)
prep "$RUN/s03"
cp "$ROOT/scripts/mapkit.R" "$ROOT/scripts/catalog_plot.R" "$RUN/s03/"
if ( cd "$RUN/s03" && "$RSCRIPT" catalog_plot.R "$MC" "$CSTART" "$TSTART"        > R-catalog_plot.R.log 2>&1 ) && [ -s "$RUN/s03/catalog.pdf" ]; then
  cp "$RUN/s03/catalog.pdf" "$FIGS/fig0-catalog.pdf"; echo "   ok   fig0-catalog.pdf"; ok=$((ok+1))
else
  echo "   FAIL fig0-catalog.pdf (see $RUN/s03/R-catalog_plot.R.log)" >&2; fail=$((fail+1))
fi

# --- tessellation ------------------------------------------------------
cp "$PKG/Section11.1files/delone-plot.R" "$RUN/s06/" 2>/dev/null
cp "$RUN/s06/delone2.out" "$RUN/s06/" 2>/dev/null
plot "$RUN/s06" delone-plot.R delone-plot.pdf fig1-delaunay.pdf

# --- b-value field -----------------------------------------------------
cp "$PKG/Section7files/delo2d-bvalues.R" "$RUN/s07/"
plot "$RUN/s07" delo2d-bvalues.R delo2d-bvalues.pdf fig2-bvalues.pdf

# --- background rate ---------------------------------------------------
cp "$PKG/Section8files/delo2d-poisson.R" "$RUN/s08/"
plot "$RUN/s08" delo2d-poisson.R delo2d-poisson.pdf fig3-poisson.pdf

# --- HIST-ETAS mu, K0 --------------------------------------------------
mkdir -p "$RUN/s11.4"
cp "$PKG/Section11.4files/hist-etas-mk.R" "$RUN/s11.4/"
cp "$RUN/s09/hist-etas-mk.omap" "$RUN/s06/delone2.out" "$RUN/s11.4/" 2>/dev/null
plot "$RUN/s11.4" hist-etas-mk.R hist-etas-mk.pdf fig4-hist-etas-mk.pdf

# --- HIST-ETAS five fields ---------------------------------------------
mkdir -p "$RUN/s11.5"
cp "$PKG/Section11.5files/hist-etas5pa.R" "$RUN/s11.5/"
cp "$RUN/s10/hist-etas5pa.omap" "$RUN/s10/hist-etas5pa.upda" \
   "$RUN/s06/delone2.out" "$RUN/s11.5/" 2>/dev/null
plot "$RUN/s11.5" hist-etas5pa.R hist-etas5pa.pdf fig5-hist-etas5pa.pdf

echo "   figures: $ok ok, $fail failed -> $FIGS"
