# HIST-PPM-V3 build on Windows

Source: https://bemlar.ism.ac.jp/ogata/HIST-PPM-V3/ (Ogata, Katsura, Tanemura,
Harte & Zhuang). Downloaded 2026-09-02, 407 files.

## Toolchain

MSYS2 UCRT64 + `gfortran` 16.2.0 (installed via
`pacman -S mingw-w64-ucrt-x86_64-gcc-fortran`; this also moved gcc 15.2.0 -> 16.2.0).
No WSL on this machine.

## Rebuilding

    ./build.sh

Compiles every `.f` in place, so `FILENAME.exe` lands next to `FILENAME.conf`,
as section 2.1 of the manual assumes. Flags:

    -O2 -std=legacy -fno-automatic -fallow-argument-mismatch -w

`-std=legacy` and `-fallow-argument-mismatch` are needed because the code is
F77 with argument-type mismatches that modern gfortran rejects by default.
`-fno-automatic` gives locals static storage (old F77 semantics) and keeps the
large local arrays off the stack.

38 of 39 sources build clean.

## Changes made to the upstream tree

1. **32 Linux ELF binaries renamed to `*.linux`.** The package ships compiled
   Linux executables at exactly the names the manual tells you to run
   (`./st-etas`, `./delone2`, ...). Left in place, bash finds the ELF file and
   fails. Renamed, `./st-etas` falls through to `st-etas.exe` and the manual's
   commands work verbatim.

2. **`programs/Section15files/bvalue2magsim.f` — two source fixes.**
   Original kept as `bvalue2magsim.f.orig`.

   Line 44 was missing a comma:

       read(2,*) npo, tx, ty amc, tpend     ->   read(2,*) npo, tx, ty, amc, tpend

   In fixed-form Fortran `ty amc` collapses to one variable `tyamc`, so the
   statement read 4 values from a 5-value record and `tpend` silently received
   the cutoff magnitude. The 5-field layout is confirmed by the writer,
   `Section14files/simNHPoisson.f:87` (`write(2,2) nev,tx1,ty,cm0,tend`,
   format `(i5,3f12.5,f14.5)`).

   Line 45 had a redundant paren group around the implied-DO io-list, which
   gfortran rejects:

       read(2,*) ((np(i),...,tp(i)), i=1,npo)  ->  read(2,*) (np(i),...,tp(i), i=1,npo)

## Known-bad upstream file

`programs/Section10files/t-etas5pa-detect.f` does not compile and cannot be
fixed: the file is not Fortran. It contains a captured `less` help screen
(overstrike backspace sequences, "SUMMARY OF LESS COMMANDS"). The local copy
md5-matches the server byte for byte (`e2d82a9e0e203879afd53d408f52b412`), so
this is corrupt at the source, not a bad download. It is a diagnostic helper,
not part of the main estimation pipeline. Ask ogata@ism.ac.jp for a good copy
if you need it.

## Verification

`Section5files/st-etas` was run against its shipped sample data and compared to
the author's reference output `st-etas.prt`:

  - MLE parameters: identical to all printed digits
  - AIC: identical (1695.1)
  - log-likelihood: agrees to 12 significant figures (840.532443492599 here vs
    840.532443490078 in the reference); our gradient is converged several
    orders tighter. The difference is optimizer float path on a newer
    compiler, not a porting error.

## Workflow note

Several programs write to `fort.2` (an unconnected unit auto-connected by the
runtime) rather than a named file — this is intended and documented in the
manual, sections 14.1, 15.1 and 16.1. Note that Section 15's `bvalue2magsim`
opens `simnhPoi.data`, which is Section 14's `fort.2` renamed; the manual's
"Reads:" line for section 15 omits this.
