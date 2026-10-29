#!/usr/bin/env bash
# =====================================================================
#  HIST-ETAS pipeline driver -- parameterised by ../region.conf
#
#  Runs stages 3..10 in isolated directories, each fed from the ACTUAL
#  output of the previous stage. Three quantities are discoverable only
#  at run time (TZ, NP, NN) plus the st-etas MLE; each is extracted the
#  moment it exists, written to state.conf, and the stage configs are
#  re-rendered so nothing drifts out of sync.
#
#  Usage:  bash scripts/run_pipeline.sh [--from N] [--skip-plots]
#
#  Stage 13 (conditional intensity + rasterised parameter fields) lives in
#  scripts/run_s13.sh and is invoked from here; it can also be run alone
#  once stages 3..10 have completed.
# =====================================================================
set -euo pipefail
export PATH=/c/msys64/ucrt64/bin:$PATH

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT/bin"
RSCRIPT="/c/Program Files/R/R-4.5.1/bin/Rscript.exe"

FROM=0; PLOTS=1
REGION="$ROOT/region.conf"; RUN=""; HYPO=""; FIGS=""
while [ $# -gt 0 ]; do
  case "$1" in
    --from)       FROM="$2"; shift 2;;
    --skip-plots) PLOTS=0; shift;;
    --region)     REGION="$2"; shift 2;;
    --rundir)     RUN="$2"; shift 2;;
    --catalog)    HYPO="$2"; shift 2;;
    --figs)       FIGS="$2"; shift 2;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

# Each region variant gets its own state/run/figs tree, keyed off the region
# file's basename, so running several thresholds side by side cannot clobber.
tag="$(basename "$REGION" .conf)"          # e.g. "region-mc55" -> tag
[ "$tag" = "region" ] && tag="" || tag="${tag#region-}"
STATE="$ROOT/state${tag:+-$tag}.conf"
[ -n "$RUN" ]  || RUN="$ROOT/run${tag:+-$tag}"
[ -n "$FIGS" ] || FIGS="$ROOT/figs${tag:+-$tag}"
[ -n "$HYPO" ] || HYPO="$ROOT/catalog/hypo.ts"

# A function, not a string: ROOT contains a space ("Hierarchical ETAS"), and an
# unquoted string variable would word-split on it.
gen(){ python "$ROOT/scripts/gen_configs.py" "$REGION" "$STATE" "$1"; }

echo "region : $REGION"
echo "catalog: $HYPO"
echo "rundir : $RUN"

say(){ printf '\n\033[1m== %s\033[0m\n' "$*"; }
setstate(){ # setstate KEY VALUE  -- idempotent upsert into state.conf
  touch "$STATE"
  grep -v "^$1 *=" "$STATE" > "$STATE.tmp" 2>/dev/null || true
  echo "$1 = $2" >> "$STATE.tmp"; mv "$STATE.tmp" "$STATE"
}

mkdir -p "$RUN"/{s03,s04,s05,s06,s07,s08,s09,s10,s13}

# ---------------------------------------------------------------- 3 ---
if [ "$FROM" -le 3 ]; then
say "stage 3  tsei2etas  (catalogue -> etas format)"
cp "$BIN/tsei2etas.exe" "$RUN/s03/"
cp "$HYPO" "$RUN/s03/hypo.ts"
( cd "$RUN/s03" && ./tsei2etas.exe < hypo.ts > tsei2etas.prt 2>&1 )
# TZ = span in days on tsei2etas's OWN time axis (col 5), not a calendar
# difference -- its leap-year handling must not be second-guessed here.
TZ=$(awk 'NR>1 && NF>=5 {if($5+0>m) m=$5+0} END{printf "%.1f", (int(m)+1)}' "$RUN/s03/work.etas")
setstate TZ "$TZ"
echo "   events: $(($(wc -l < "$RUN/s03/work.etas")-1))   TZ=$TZ days"
fi
gen "$RUN/s04" >/dev/null

# ---------------------------------------------------------------- 4 ---
if [ "$FROM" -le 4 ]; then
say "stage 4  etas2aniso  (anisotropic clusters)"
cp "$BIN/etas2aniso.exe" "$RUN/s04/"
cp "$RUN/s03/work.etas" "$RUN/s04/"
( cd "$RUN/s04" && ./etas2aniso.exe > etas2aniso.prt 2>&1 )
NP=$(awk 'NR==1{print $1}' "$RUN/s04/etas2aniso.out3")
setstate NP "$NP"
echo "   events above Mc: NP=$NP"
fi
NP=$(awk -F= '/^NP/{gsub(/ /,"",$2);print $2}' "$STATE")
[ -n "$NP" ] && [ "$NP" -ge 30 ] || { echo "FATAL: NP=$NP too small to fit" >&2; exit 1; }
gen "$RUN/s05" >/dev/null; gen "$RUN/s06" >/dev/null

# ---------------------------------------------------------------- 5 ---
if [ "$FROM" -le 5 ]; then
say "stage 5  st-etas  (stationary ETAS, 7 params)"
cp "$BIN/st-etas.exe" "$RUN/s05/"
cp "$RUN/s04/etas2aniso.out3" "$RUN/s05/"
( cd "$RUN/s05" && ./st-etas.exe > st-etas.prt 2>&1 ) || true
# The MLE is printed as "mle =" then 5 values, continued by 2 on the next
# line. Feed it forward as the starting theta for the hierarchical stages.
MLE=$(awk '/mle =/{getline nx; sub(/.*mle *=/,""); print $0, nx}' "$RUN/s05/st-etas.prt" | tail -1)
if [ "$(echo $MLE | wc -w)" -eq 7 ]; then
  setstate THETA_MK "$MLE"
  echo "   MLE: $MLE"
  grep -E "log likelihood" "$RUN/s05/st-etas.prt" | tail -1
else
  echo "   WARNING: could not parse 7 MLE values; keeping THETA0 from region.conf" >&2
fi
fi
gen "$RUN/s06" >/dev/null

# ---------------------------------------------------------------- 6 ---
if [ "$FROM" -le 6 ]; then
say "stage 6  delone1 + delone2  (Delaunay tessellation)"
cp "$BIN/delone1.exe" "$BIN/delone2.exe" "$RUN/s06/"
cp "$RUN/s04/etas2aniso.out3" "$RUN/s06/"
( cd "$RUN/s06" && ./delone1.exe > delone1.prt 2>&1 && ./delone2.exe > delone2.prt 2>&1 )
read -r NN NPEXT NDELO _ <<<"$(head -1 "$RUN/s06/delone2.out")"
setstate NN "$NN"
echo "   target events NN=$NN  vertices=$NPEXT  triangles=$NDELO"
fi
for s in s07 s08 s09 s10; do gen "$RUN/$s" >/dev/null; done

# ------------------------------------------------------------- 7 & 8 ---
if [ "$FROM" -le 7 ]; then
say "stage 7  delo2d-bvalues  (location-dependent b)"
cp "$BIN/delo2d-bvalues.exe" "$RUN/s07/"; cp "$RUN/s06/delone2.out" "$RUN/s07/"
( cd "$RUN/s07" && ./delo2d-bvalues.exe > delo2d-bvalues.prt 2>&1 ) || true
echo "   ABIC: $(grep -oE 'abic *= *[-0-9.DE+]+' "$RUN/s07/delo2d-bvalues.prt" | tail -1)"
fi
if [ "$FROM" -le 8 ]; then
say "stage 8  delo2d-poisson  (background rate)"
cp "$BIN/delo2d-poisson.exe" "$RUN/s08/"; cp "$RUN/s06/delone2.out" "$RUN/s08/"
( cd "$RUN/s08" && ./delo2d-poisson.exe > delo2d-poisson.prt 2>&1 ) || true
echo "   ABIC: $(grep -oE 'abic *= *[-0-9.DE+]+' "$RUN/s08/delo2d-poisson.prt" | tail -1)"
fi

# ---------------------------------------------------------------- 9 ---
if [ "$FROM" -le 9 ]; then
say "stage 9  hist-etas-mk  (varying mu and K0)  -- slow"
cp "$BIN/hist-etas-mk.exe" "$RUN/s09/"; cp "$RUN/s06/delone2.out" "$RUN/s09/"
( cd "$RUN/s09" && ./hist-etas-mk.exe > hist-etas-mk.prt 2>&1 ) || true
if [ -s "$RUN/s09/hist-etas-mk.upda" ]; then
  echo "   w1 w2 ABIC m: $(head -1 "$RUN/s09/hist-etas-mk.upda")"
  echo "   simplex iterations: $(grep -c '####' "$RUN/s09/hist-etas-mk.prt" || true)"
else
  echo "   FAILED -- see $RUN/s09/hist-etas-mk.prt" >&2
fi
fi

# --------------------------------------------------------------- 10 ---
if [ "$FROM" -le 10 ]; then
say "stage 10  hist-etas5pa  (five varying fields)  -- slow"
cp "$BIN/hist-etas5pa.exe" "$RUN/s10/"; cp "$RUN/s06/delone2.out" "$RUN/s10/"
cp "$RUN/s09/hist-etas-mk.upda" "$RUN/s10/"
( cd "$RUN/s10" && ./hist-etas5pa.exe > hist-etas5pa.prt 2>&1 ) || true
if [ -s "$RUN/s10/hist-etas5pa.upda" ]; then
  echo "   header: $(head -2 "$RUN/s10/hist-etas5pa.upda" | tail -1)"
  echo "   simplex iterations: $(grep -c '####' "$RUN/s10/hist-etas5pa.prt" || true)"
else
  echo "   FAILED -- see $RUN/s10/hist-etas5pa.prt" >&2
fi
fi

# ------------------------------------------------------------- plots ---
if [ "$PLOTS" -eq 1 ]; then
say "figures (R)"
bash "$ROOT/scripts/make_plots.sh" --region "$REGION" --rundir "$RUN" --figs "$FIGS" || echo "   plots failed (non-fatal)" >&2
fi

# --------------------------------------------------------------- 13 ---
# Evaluation rather than estimation: put the fitted fields, and the
# conditional intensity they imply, on a regular grid. Needs s06 and s10.
if [ "$FROM" -le 13 ]; then
S13ARGS=(--region "$REGION" --rundir "$RUN" --figs "$FIGS")
[ "$PLOTS" -eq 1 ] || S13ARGS+=(--skip-plots)
bash "$ROOT/scripts/run_s13.sh" "${S13ARGS[@]}" || echo "   stage 13 failed (non-fatal)" >&2
fi

say "done -- outputs under $RUN"
