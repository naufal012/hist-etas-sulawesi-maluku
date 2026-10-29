# Jalilian's stationary-parameter ETAS (etas()) on the Sulawesi-Maluku M>=5.0
# catalogue, with his two spatial kernels:
#   mver = 1  inverse power  f = (q-1)/(pi s) (1 + r^2/s)^-q,  s = D exp(gamma (m - m0))
#   mver = 2  Gaussian       f = exp(-r^2/(2 s)) / (2 pi s)
# productivity A exp(alpha (m - m0)), normalised Omori (p-1)/c (1 + t/c)^-p,
# background nu * u(x, y) by stochastic declustering. Same catalogue, target
# window and units (degrees; cos 0 = 1) as the HIST-ETAS fits.
#   Rscript scripts/jalilian_mver.R 1|2   -> run-mc50-R/jalilian-mver<k>.rds/.txt
lib <- Sys.getenv("ETAS_LIB", "D:/Documents/Hierarchical ETAS/.Rlib-histetas-frozen")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(ETAS))
mver <- as.integer(commandArgs(TRUE)[1])
ROOT <- "D:/Documents/Hierarchical ETAS/sulawesi-maluku"
OUT <- file.path(ROOT, "run-mc50-R")
h <- read.table(file.path(ROOT, "catalog", "hypo-full.ts"),
                col.names = c("yr", "mo", "dy", "hr", "mi", "sec", "long", "lat",
                              "depth", "mag"))
quakes <- data.frame(date = sprintf("%04d-%02d-%02d", h$yr, h$mo, h$dy),
                     time = sprintf("%02d:%02d:%05.2f", h$hr, h$mi, h$sec),
                     long = h$long, lat = h$lat, mag = h$mag)
t0 <- as.POSIXct(paste(quakes$date[1], quakes$time[1]), tz = "GMT")
cat5 <- suppressWarnings(catalog(quakes, time.begin = t0, study.start = t0 + 730 * 86400,
                                 study.end = t0 + 9676 * 86400,
                                 long.range = c(118, 132), lat.range = c(-6, 6),
                                 mag.threshold = 5.0, flatmap = FALSE, roundoff = FALSE))
# optional 2nd argument "from2": start mver = 1 from the converged Gaussian fit
# (the package's default start, p = 1.3, ran into the p -> 1 edge of the
# normalised Omori density and its line search stalled), with q from the
# HIST-ETAS stationary fit (2.43), since the Gaussian has no q
param0 <- NULL
if (identical(commandArgs(TRUE)[2], "from2")) {
  g <- readRDS(file.path(OUT, "jalilian-mver2.rds"))$param
  param0 <- c(mu = g[["mu"]], A = g[["A"]], c = g[["c"]], alpha = g[["alpha"]],
              p = g[["p"]], D = g[["D"]], q = 2.43, gamma = g[["gamma"]])
  print(param0)
}
t1 <- Sys.time()
fit <- etas(cat5, param0 = param0, mver = mver, nthreads = 8, verbose = TRUE)
mins <- as.numeric(difftime(Sys.time(), t1, units = "mins"))
saveRDS(fit, file.path(OUT, sprintf("jalilian-mver%d.rds", mver)))
out <- c(sprintf("mver = %d: %.1f min, %d target events", mver, mins,
                 sum(cat5$revents[, "flag"] == 1)),
         paste("param:", paste(names(fit$param), signif(fit$param, 5), collapse = "  ")),
         paste("s.e.: ", paste(names(fit$param), signif(fit$asd[fit$itr, ], 3), collapse = "  ")),
         sprintf("opt$loglik %.3f  opt$aic %.3f  (declustering iterations %d)",
                 fit$opt$loglik, fit$opt$aic, fit$itr),
         sprintf("background share (sum of pb over targets) %.3f",
                 sum(fit$pb[cat5$revents[, "flag"] == 1]) / sum(cat5$revents[, "flag"] == 1)))
writeLines(out)
writeLines(out, file.path(OUT, sprintf("jalilian-mver%d.txt", mver)))
