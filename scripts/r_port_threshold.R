# Threshold table (tab:trend) with the R port: stationary + hist-etas-mk at
# M>=6.0 and M>=5.5, from the same inputs the Fortran runs used
# (run/ = mc60 from catalog/hypo.ts, TZ 9675; run-mc55/ from hypo-full.ts, TZ 9676).
#
#   Rscript scripts/r_port_threshold.R 6.0
#   Rscript scripts/r_port_threshold.R 5.5
# Output: run-mc50-R/threshold-mc<..>.rds and .txt
lib <- Sys.getenv("ETAS_LIB", "D:/Documents/Hierarchical ETAS/.Rlib-histetas")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(ETAS))
E <- asNamespace("ETAS")
options(histetas.threads = 4)

ROOT <- "D:/Documents/Hierarchical ETAS/sulawesi-maluku"
mc <- as.numeric(commandArgs(TRUE)[1])
cfg <- if (mc == 6.0) list(hypo = "hypo.ts", tz = 9675, frun = "run", tag = "mc60") else
  if (mc == 5.5) list(hypo = "hypo-full.ts", tz = 9676, frun = "run-mc55", tag = "mc55") else
    stop("threshold must be 6.0 or 5.5")
OUT <- Sys.getenv("RPORT_OUT", file.path(ROOT, "run-mc50-R"))
dir.create(OUT, showWarnings = FALSE)

h <- read.table(file.path(ROOT, "catalog", cfg$hypo),
                col.names = c("yr", "mo", "dy", "hr", "mi", "sec", "long", "lat",
                              "depth", "mag"))
quakes <- data.frame(date = sprintf("%04d-%02d-%02d", h$yr, h$mo, h$dy),
                     time = sprintf("%02d:%02d:%05.2f", h$hr, h$mi, h$sec),
                     long = h$long, lat = h$lat, mag = h$mag)
t0 <- as.POSIXct(paste(quakes$date[1], quakes$time[1]), tz = "GMT")
mkcat <- function(mth)
  suppressWarnings(catalog(quakes, time.begin = t0, study.start = t0 + 730 * 86400,
                           study.end = t0 + cfg$tz * 86400,
                           long.range = c(118, 132), lat.range = c(-6, 6),
                           mag.threshold = mth, flatmap = FALSE, roundoff = FALSE))
low <- mkcat(min(quakes$mag))
catm <- mkcat(mc)
an <- E$.hist_aniso_catalog(low, clsm = mc + 0.5, cutm = mc, window = 1)
keep <- low$revents[, "mm"] + low$mag.threshold >= mc
aniso <- an$all[keep, c("long", "lat", "dd", "ee", "ff")]
ref <- E$.hist_read_aniso(file.path(ROOT, cfg$frun, "s04", "etas2aniso.out3"))

t1 <- Sys.time()
fit <- histetas(catm, fields = c("mu", "K"), mode = "fortran",
                control = list(aniso = aniso, verbose = TRUE))
mins <- as.numeric(difftime(Sys.time(), t1, units = "mins"))
saveRDS(fit, file.path(OUT, sprintf("threshold-%s.rds", cfg$tag)))
up <- E$.hist_read_upda(file.path(ROOT, cfg$frun, "s09", "hist-etas-mk.upda"))
out <- c(
  sprintf("M>=%.1f: %d events, %d target; aniso vs Fortran %d rows, %d anisotropic (Fortran %d), max diff %.2e",
          mc, nrow(catm$revents), sum(catm$revents[, "flag"] == 1), nrow(ref$events),
          sum(aniso$ff != 0), sum(ref$events$ff != 0),
          max(abs(as.matrix(an$events[, c("dd", "ee", "ff")]) -
                  as.matrix(ref$events[, c("dd", "ee", "ff")])))),
  sprintf("fit %.1f min", mins),
  sprintf("stationary theta %s logL %.4f", paste(signif(fit$stationary$theta, 6), collapse = " "),
          fit$stationary$loglik),
  sprintf("mk theta %s", paste(signif(fit$mk$theta, 6), collapse = " ")),
  sprintf("mk w %s ABIC %.4f evals %d lines %d", paste(signif(fit$mk$w, 6), collapse = " "),
          fit$mk$abic, nrow(fit$mk$evals), nrow(fit$mk$simplex$trace)),
  sprintf("Fortran mk theta %s", paste(signif(up$theta, 6), collapse = " ")),
  sprintf("Fortran mk w %s ABIC %.4f", paste(signif(up$weights, 6), collapse = " "), up$abic))
writeLines(out)
writeLines(out, file.path(OUT, sprintf("threshold-%s.txt", cfg$tag)))
