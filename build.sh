#!/usr/bin/env bash
# Build all HIST-PPM-V3 FORTRAN programs with MSYS2 UCRT64 gfortran.
# Each .f is compiled in place, so FILENAME.exe sits beside FILENAME.conf
# exactly as the manual (section 2.1) assumes.
export PATH=/c/msys64/ucrt64/bin:$PATH

FFLAGS="-O2 -std=legacy -fno-automatic -fallow-argument-mismatch -w"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ok=0; fail=0; failed=()

while IFS= read -r src; do
    dir=$(dirname "$src"); base=$(basename "$src" .f)
    if ( cd "$dir" && gfortran $FFLAGS -o "$base.exe" "$base.f" ) 2>"$ROOT/.build_err"; then
        ok=$((ok+1)); printf '  ok    %s\n' "${src#./}"
    else
        fail=$((fail+1)); failed+=("${src#./}")
        printf '  FAIL  %s\n' "${src#./}"
        sed 's/^/          /' "$ROOT/.build_err" | head -12
    fi
done < <(find "$ROOT/programs" -name '*.f' | sed "s|$ROOT/||" | sort)

rm -f "$ROOT/.build_err"
echo
echo "built $ok, failed $fail"
for f in "${failed[@]}"; do echo "  failed: $f"; done
