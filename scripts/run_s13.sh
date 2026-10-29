#!/usr/bin/env bash
# =====================================================================
#  Stage 13 -- conditional intensity maps and rasterised parameter fields
#
#  Everything before stage 13 estimates the model. Stage 13 evaluates it:
#  it puts the fitted fields, and the space-time conditional intensity they
#  imply, onto a regular lon/lat grid so they can be drawn as maps rather
#  than as coloured dots at the tessellation vertices.
#
#  Two programs:
#    intensity5pa.exe  -- the package's intensityday2magnod5pa.f with the
#                         four edits documented in scripts/patch_intensity.py.
#                         Evaluates lambda(x, y, t) at one instant.
#    interp5pa.exe     -- src/interp5pa.f (ours). Evaluates mu, K, alpha, p, q
#                         on the same grid.
#
#  Both need stage 6 (delone2.out) and stage 10 (hist-etas5pa.upda).
#
#  Usage: bash scripts/run_s13.sh [--region F] [--rundir D] [--figs D]
#                                 [--snap DAY:TAG:LABEL] ... [--skip-plots]
# =====================================================================
set -euo pipefail
export PATH=/c/msys64/ucrt64/bin:$PATH

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PKG="D:/Documents/Hierarchical ETAS/HIST-PPM-V3/programs"
RSCRIPT="/c/Program Files/R/R-4.5.1/bin/Rscript.exe"
SRC13="$PKG/Section13.3files5pa/intensityday2magnod5pa.f"

REGION="$ROOT/region.conf"; RUN=""; FIGS=""; PLOTS=1
SNAPS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --region)     REGION="$2"; shift 2;;
    --rundir)     RUN="$2";    shift 2;;
    --figs)       FIGS="$2";   shift 2;;
    --snap)       SNAPS+=("$2"); shift 2;;
    --skip-plots) PLOTS=0; shift;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

tag="$(basename "$REGION" .conf)"
[ "$tag" = "region" ] && tag="" || tag="${tag#region-}"
STATE="$ROOT/state${tag:+-$tag}.conf"
[ -n "$RUN" ]  || RUN="$ROOT/run${tag:+-$tag}"
[ -n "$FIGS" ] || FIGS="$ROOT/figs${tag:+-$tag}"
S13="$RUN/s13"

say(){ printf '\n\033[1m== %s\033[0m\n' "$*"; }
cfg(){ awk -F= -v k="$1" '$1 ~ "^ *"k" *$" {sub(/#.*/,"",$2); gsub(/^ +| +$/,"",$2); print $2}' "$REGION" | tail -1; }
st(){  awk -F= -v k="$1" '$1 ~ "^ *"k" *$" {gsub(/^ +| +$/,"",$2); print $2}' "$STATE" | tail -1; }

TZ="$(st TZ)"
[ -n "$TZ" ] || { echo "FATAL: TZ missing from $STATE -- run stages 3..10 first" >&2; exit 1; }

for f in "$RUN/s06/delone2.out" "$RUN/s10/hist-etas5pa.upda"; do
  [ -s "$f" ] || { echo "FATAL: missing $f" >&2; exit 1; }
done

# Snapshots: --snap wins, else SNAPSHOTS from region.conf, else end of span.
if [ ${#SNAPS[@]} -eq 0 ]; then
  IFS=',' read -ra SNAPS <<< "$(cfg SNAPSHOTS)"
fi
if [ ${#SNAPS[@]} -eq 0 ] || [ -z "${SNAPS[0]:-}" ]; then
  SNAPS=("$TZ:end:end of target period")
fi

mkdir -p "$S13" "$FIGS"

say "stage 13  conditional intensity  ($(basename "$REGION"))"

# ---- configs ---------------------------------------------------------
python "$ROOT/scripts/gen_configs.py" "$REGION" "$STATE" "$S13" >/dev/null

# The array bounds of both programs are fixed at compile time; size them from
# the tessellation actually produced rather than trusting the package's Japan
# figures. A silent overflow here writes past the end of a COMMON block.
read -r NNH NPEXT NDELO _ <<<"$(head -1 "$RUN/s06/delone2.out")"
NPA=$(( NPEXT + NPEXT/4 + 500 ))
MDL=$(( NDELO + NDELO/4 + 500 ))
echo "   tessellation: nn=$NNH npext=$NPEXT ndelo=$NDELO  -> npa=$NPA mdl=$MDL"

# ---- inputs ----------------------------------------------------------
cp "$RUN/s06/delone2.out"        "$S13/"
cp "$RUN/s10/hist-etas5pa.upda"  "$S13/"
cp "$RUN/s10/hist-etas5pa.prt"   "$S13/"   # best ABIC for interp_plot.R's label
cp "$RUN/s10/hist-etas5pa.conf"  "$S13/"
cp "$ROOT/scripts/drawmap.R" "$ROOT/scripts/mapkit.R" \
   "$ROOT/scripts/intensity_plot.R" "$ROOT/scripts/interp_plot.R" "$S13/"
mkdir -p "$S13/maps"
cp "$ROOT/maps/w1.pp.data" "$ROOT/maps/w1.br.data" "$S13/maps/"

# work.para carries events that post-date the fitted span. There are none
# here -- the evaluation instants lie inside it -- but the program reads the
# file unconditionally, so it must exist and be empty.
: > "$S13/work.para"

# ---- build -----------------------------------------------------------
say "build"
python "$ROOT/scripts/patch_intensity.py" "$SRC13" "$S13/intensity5pa.f" \
       --npa "$NPA" --mdl "$MDL"
( cd "$S13" && gfortran -O2 -std=legacy -w -o intensity5pa.exe intensity5pa.f )

sed -e "s/parameter(npa=8000, mdl=16000)/parameter(npa=$NPA, mdl=$MDL)/" \
    "$ROOT/src/interp5pa.f" > "$S13/interp5pa.f"
( cd "$S13" && gfortran -O2 -std=legacy -w -o interp5pa.exe interp5pa.f )

cp "$PKG/Section13.3files5pa/node.f" "$S13/node.f"
( cd "$S13" && gfortran -O2 -std=legacy -w -o node.exe node.f && ./node.exe > node.prt )
echo "   grid points: $(wc -l < "$S13/node.dat")"

# ---- conditional intensity, one run per instant ----------------------
: > "$S13/snapshots.txt"
for s in "${SNAPS[@]}"; do
  day="${s%%:*}"; rest="${s#*:}"; nm="${rest%%:*}"; lab="${rest#*:}"
  [ "$lab" = "$nm" ] && lab="day $day"
  off=$(python -c "print('%.5f' % (float('$day') - float('$TZ')))")
  printf '%s   ! zzz = TZ + this = day %s\n' "$off" "$day" > "$S13/time.conf"
  say "  lambda at day $day  ($nm: $lab)"
  ( cd "$S13" && ./intensity5pa.exe > "intensity-$nm.prt" 2>&1 )
  mv "$S13/img1.data" "$S13/img1-$nm.data"
  printf '%s\t%s\t%s\n' "$nm" "$day" "$lab" >> "$S13/snapshots.txt"
  echo "   $(wc -l < "$S13/img1-$nm.data") pixels; log10 lambda range: \
$(awk '{if(NR==1||$4<a)a=$4; if(NR==1||$4>b)b=$4} END{printf "%.3f .. %.3f", a, b}' \
  "$S13/img1-$nm.data")"
done

# ---- parameter fields ------------------------------------------------
say "  parameter fields"
( cd "$S13" && ./interp5pa.exe > interp5pa.prt 2>&1 )
tail -1 "$S13/interp5pa.prt"

# ---- figures ---------------------------------------------------------
if [ "$PLOTS" -eq 1 ]; then
  say "figures"
  ( cd "$S13" && "$RSCRIPT" intensity_plot.R > R-intensity.log 2>&1 ) \
    && { cp "$S13/intensity.pdf" "$FIGS/fig6-intensity.pdf"; \
         echo "   ok   fig6-intensity.pdf"; sed -n '1,12p' "$S13/R-intensity.log"; } \
    || echo "   FAIL fig6-intensity.pdf (see $S13/R-intensity.log)" >&2
  ( cd "$S13" && "$RSCRIPT" interp_plot.R > R-interp.log 2>&1 ) \
    && { cp "$S13/interp5pa.pdf" "$FIGS/fig7-interp5pa.pdf"; \
         echo "   ok   fig7-interp5pa.pdf"; sed -n '1,12p' "$S13/R-interp.log"; } \
    || echo "   FAIL fig7-interp5pa.pdf (see $S13/R-interp.log)" >&2
fi

say "stage 13 done -- $S13"
