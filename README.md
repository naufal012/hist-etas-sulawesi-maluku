# Hierarchical ETAS for Sulawesi–Maluku — reproduction package

Code, configuration, catalogue and fit results accompanying *Threshold
Adequacy, Starting Values and Residual Diagnostics of Hierarchical Space–Time
ETAS: An application to Sulawesi–Maluku* (N. L. Tasbihi, Department of
Statistics, Institut Teknologi Sepuluh Nopember), submitted.

## What this is

The analysis fits Ogata's HIST-ETAS model to 118–132°E, 6°S–6°N for
2000-01-01 to 2026-06-30, using the `HIST-PPM-V3` Fortran package. It runs as
nine stages, each in its own directory, each fed from the actual output of its
predecessor.

## What is and is not included

**Included.** The pipeline driver and config generator (`scripts/`), the
locally modified Fortran sources (`src/`), every stage configuration, the R
plotting scripts, the catalogue snapshot the paper's numbers come from, the
build notes, and the program output of the reported fits (see *Results*
below), so the paper's numbers can be checked without rerunning anything.

**Not included: `HIST-PPM-V3` itself.** The package is Ogata's and is
distributed from <https://bemlar.ism.ac.jp/ogata/HIST-PPM-V3/>. Redistributing
it here would be presumptuous about its licence, so modifications are supplied
as **diffs against the upstream tree** instead. Download the package, apply
the patches in `patches/`, then build.

## Reproducing

1. Download `HIST-PPM-V3` from the URL above (407 files as of 2026-09-02).
2. Apply `patches/*.patch`; see `BUILD-NOTES.md` for what each one does and
   why.
3. Build with `./build.sh` — MSYS2 UCRT64, `gfortran` 16.2.0, flags
   `-O2 -std=legacy -fno-automatic -fallow-argument-mismatch -w`.
4. Generate the stage configs and run the stages:
   `bash scripts/run_pipeline.sh --region region-mc50-ss.conf --catalog catalog/hypo-full.ts`
   `region-mc50-ss.conf` is the reported configuration. It differs from
   `region-mc50.conf` only in the starting values for `hist-etas-mk`: the
   `-ss` state file carries the stationary MLE refitted by
   `scripts/refit_stationary.py` (`s14/stationary.json`), because the
   package's `st-etas` integrates the kernel over a box that does not contain
   the epicentres and its estimate is not a valid start.
5. The threshold comparison in the paper's justification section uses
   `hist-etas-mk` at M>=6.0 (`run-ss/`) and M>=5.5 (`run-mc55-ss/`), and the
   package's JMA demo refitted the same way (`japan-control/`), each started
   from its own stationary MLE. `run/`, `run-mc50/` and `run-mc55/` are the
   earlier st-etas-started fits, kept for comparison.

6. *Not reported in the paper; kept for completeness.* A smoothing-sensitivity
   scan refits hist-etas5pa at the ABIC weights scaled by s = 4 ... 1/64:
   apply `patches/hist-etas5pa-fixw.patch` to a copy of the upstream
   `Section10files/hist-etas5pa.f` (it replaces the outer simplex by one
   inner solve), save it as `src/hist-etas5pa-fixw.f`, build it into `bin/`,
   then
   `python scripts/setup_smoothing_scan.py --run run-mc50-ss`, run
   `hist-etas5pa-fixw.exe > hist-etas5pa.prt` in each `run-mc50-ss/s15/s*`,
   and score them with `python scripts/smoothing_sensitivity.py`. The
   s = 1/32 solve diverges from the ABIC start; rerun it with
   `--tags 1_32 --from-tag 1_16`. Results: `run-mc50-ss/s15/summary.json`.

R 4.5.1 produced the figures.

## The catalogue

`catalog/usgs-mc.csv` is the USGS FDSN extract behind every reported number:
21 528 events, M 3.2-7.8, fetched in yearly slices.
`catalog/hypo-full.ts` is the same catalogue in the package's input format
and `catalog/work-mc50.etas` is the `tsei2etas` output the run used. The
sub-threshold events matter: the cluster shape matrices $S_j$ are estimated
from the aftershocks in the first day after each mainshock.

`catalog/usgs-raw.csv` / `hypo.ts` are an older M>=4.5 extract, kept only for
provenance. Do not start from them: they agree with the full catalogue above
M 4.5 but give different cluster anisotropy.

## Layout

    scripts/        pipeline driver, config generator, catalogue fetch,
                    Python likelihood port, residual and covariate analyses
    src/            new Fortran (interp5pa.f, the parameter-field rasteriser)
    patches/        diffs against the upstream HIST-PPM-V3 tree
    catalog/        catalogue snapshot + the work.etas actually used
    run-mc50-ss/    stage configs for the reported fit (M>=5.0)
    run-ss/, run-mc55-ss/   hist-etas-mk at M>=6.0 and M>=5.5
    japan-control/  hist-etas-mk config for the JMA demo, stationary start
    run*/           the earlier st-etas-started configs

## Results

Program output of the reported fits, as written by HIST-PPM-V3 (`.prt` is the
printed log, `.upda` the final estimates, `.omap` the fitted coefficients;
take ABIC from the `.prt`, where it is the last evaluation).

    run-mc50-ss/s09/   hist-etas-mk at M>=5.0 (Table 1), and delone2.out,
                       the Delaunay tessellation both hierarchical fits use
    run-mc50-ss/s10/   hist-etas5pa at M>=5.0 (Tables 1-2)
    run-mc50-ss/s13/   intensity grids at the three instants of Figure 3
                       (img1-*.data, intensity-*.prt) and the parameter-field
                       grid of Figure 4 (interp5pa.pixel, node.dat)
    run-mc50-ss/s14/   stationary refit (stationary.json) and the residuals of
                       Figures 6-7: transformed times (tau-*.txt), Voronoi
                       residuals (voronoi-*.txt), residuals.json
    run-ss/s09/, run-mc55-ss/s09/   hist-etas-mk at M>=6.0 and M>=5.5 (Table 3)
    japan-control/     hist-etas-mk on the package's JMA demo (Table 3)

## Licence

The author's own code and configurations are under the MIT licence
(`LICENSE`). HIST-PPM-V3 is not included and is not covered by it: the
`patches/` are differences against the authors' distribution and remain
subject to its terms. The USGS catalogue is a US government work in the
public domain.
