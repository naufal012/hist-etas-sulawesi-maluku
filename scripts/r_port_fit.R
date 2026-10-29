# End-to-end HIST-ETAS fit of the Sulawesi-Maluku catalogue at M>=5.0 with the
# R port (etas-jalilian, histetas()), as a test against the Fortran run-mc50.
#
# Unlike etas-jalilian's tools/validation-sulawesi/refit.R, which restarts the
# programs from their own .conf files, this starts from the catalogue: the
# stationary model is fitted by histetas() itself (the Fortran pipeline started
# hist-etas-mk from st-etas, whose compensator is defective), then mk, then 5pa.
# The anisotropy is estimated from the full sub-threshold catalogue, as the
# Fortran's etas2aniso does, and passed in.
#
#   Rscript scripts/r_port_fit.R            # ~30 min; progress on stderr
#
# Output: run-mc50-R/fit.rds, run-mc50-R/fit-log.txt
lib <- Sys.getenv("ETAS_LIB", "D:/Documents/Hierarchical ETAS/.Rlib-histetas")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(ETAS))
E <- asNamespace("ETAS")
stopifnot(exists(".hist_model", envir = E))

ROOT <- "D:/Documents/Hierarchical ETAS/sulawesi-maluku"
FRUN <- file.path(ROOT, "run-mc50")
OUT <- file.path(ROOT, "run-mc50-R")
dir.create(OUT, showWarnings = FALSE)
log <- function(...) {
  s <- sprintf(...)
  message(s)
  cat(s, "\n", file = file.path(OUT, "fit-log.txt"), append = TRUE, sep = "")
}
unlink(file.path(OUT, "fit-log.txt"))
log("ETAS from %s, %s", find.package("ETAS"), format(Sys.time()))

# ---- catalogue: the same USGS extract the Fortran pipeline read (run-mc50/s03/hypo.ts) ----
h <- read.table(file.path(ROOT, "catalog", "hypo-full.ts"),
                col.names = c("yr", "mo", "dy", "hr", "mi", "sec", "long", "lat",
                              "depth", "mag"))
quakes <- data.frame(
  date = sprintf("%04d-%02d-%02d", h$yr, h$mo, h$dy),
  time = sprintf("%02d:%02d:%05.2f", h$hr, h$mi, h$sec),
  long = h$long, lat = h$lat, mag = h$mag)
# work.etas measures time in days from the first event; TSTART = 730, TZ = 9676
t0 <- as.POSIXct(paste(quakes$date[1], quakes$time[1]), tz = "GMT")
day <- 86400
mkcat <- function(mth)
  suppressWarnings(catalog(quakes, time.begin = t0, study.start = t0 + 730 * day,
                           study.end = t0 + 9676 * day,
                           long.range = c(118, 132), lat.range = c(-6, 6),
                           mag.threshold = mth, flatmap = FALSE, roundoff = FALSE))
low <- mkcat(min(quakes$mag))
cat5 <- mkcat(5.0)
log("catalogue: %d events (M>=%.1f), %d at M>=5.0, %d target",
    nrow(low$revents), min(quakes$mag), nrow(cat5$revents),
    sum(cat5$revents[, "flag"] == 1))

# ---- anisotropy from the whole catalogue (etas2aniso: clms 5.5, cutm 5.0) ----
an <- E$.hist_aniso_catalog(low, clsm = 5.5, cutm = 5.0, window = 1)
keep <- low$revents[, "mm"] + low$mag.threshold >= 5.0
aniso <- an$all[keep, c("long", "lat", "dd", "ee", "ff")]
stopifnot(nrow(aniso) == nrow(cat5$revents),
          max(abs(an$all$mag[keep] - (cat5$revents[, "mm"] + 5))) < 1e-9)
log("anisotropy: %d clusters, %d anisotropic main shocks",
    length(unique(an$all$cluster[an$all$cluster > 0])), sum(aniso$ff != 0))
ref <- E$.hist_read_aniso(file.path(FRUN, "s04", "etas2aniso.out3"))
log("  vs Fortran etas2aniso.out3: %d rows, max |dd,ee,ff| diff %.2e",
    nrow(ref$events),
    max(abs(as.matrix(an$events[, c("dd", "ee", "ff")]) -
            as.matrix(ref$events[, c("dd", "ee", "ff")]))))

# ---- the fit ----------------------------------------------------------------
if (nzchar(Sys.getenv("DRYRUN"))) quit(save = "no")
t1 <- Sys.time()
fit <- histetas(cat5, fields = c("mu", "K", "alpha", "p", "q"), mode = "fortran",
                control = list(aniso = aniso, verbose = TRUE))
mins <- as.numeric(difftime(Sys.time(), t1, units = "mins"))
saveRDS(fit, file.path(OUT, "fit.rds"))

d2 <- E$.hist_read_delone2(file.path(FRUN, "s06", "delone2.out"))
cmp <- E$.hist_mesh_compare(fit$mesh, E$.hist_mesh_delone2(d2))
log("mesh: %d nodes, %d triangles; vs Fortran delone2.out: %d unmatched nodes, %d/%d unmatched triangles, max area diff %g",
    nrow(fit$mesh$nodes), nrow(fit$mesh$tri), cmp$nodes.unmatched,
    nrow(cmp$only.a), nrow(cmp$only.b), cmp$max.area.diff)
log("fit: %.1f min", mins)
log("stationary: %s  logL %.4f  (conv %d)",
    paste(signif(fit$stationary$theta, 6), collapse = " "),
    fit$stationary$loglik, fit$stationary$convergence)
log("mk:  theta %s", paste(signif(fit$mk$theta, 6), collapse = " "))
log("     w %s  ABIC %.6f  evals %d  simplex lines %d",
    paste(signif(fit$mk$w, 6), collapse = " "), fit$mk$abic,
    nrow(fit$mk$evals), nrow(fit$mk$simplex$trace))
log("5pa: w %s  ABIC %.6f  evals %d  simplex lines %d",
    paste(signif(fit$`5pa`$w, 6), collapse = " "), fit$`5pa`$abic,
    nrow(fit$`5pa`$evals), nrow(fit$`5pa`$simplex$trace))
log("5pa logL (fortran-mode quadrature) %.4f", fit$loglik)
