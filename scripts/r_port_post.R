# After scripts/r_port_fit.R: writes the R-port fit into run-mc50-R/ in the
# Fortran pipeline's file formats, so the existing downstream scripts
# (residuals.py, tectonic_context.py, interp_plot.R, intensity_plot.R) produce
# the paper's numbers and figures from it without changes.
#
#   Rscript scripts/r_port_post.R
#
#   s09/hist-etas-mk.upda, s10/hist-etas5pa.upda   mk / 5pa layouts
#   s14/stationary.json                             histetas()'s stationary fit
#   s13/interp5pa.pixel, s13/img1-<tag>.data        predict() on the 210 x 180 grid
#   s13/interp5pa.pdf, s13/intensity.pdf            the package plotting scripts
#   r-port-summary.json                             numbers for the comparison draft
# Inputs that do not depend on the fit (settings, delone2.out, work.etas) are
# copied from run-mc50; r_port_fit.R checks that the port's mesh is identical.
lib <- Sys.getenv("ETAS_LIB", "D:/Documents/Hierarchical ETAS/.Rlib-histetas")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(ETAS))
E <- asNamespace("ETAS")

ROOT <- "D:/Documents/Hierarchical ETAS/sulawesi-maluku"
FRUN <- file.path(ROOT, "run-mc50")
RUN <- file.path(ROOT, "run-mc50-R")
fit <- readRDS(file.path(RUN, "fit.rds"))
for (s in c("s03", "s06", "s09", "s10", "s13", "s14")) dir.create(file.path(RUN, s), FALSE)
cp <- function(from, to) stopifnot(file.copy(file.path(FRUN, from), file.path(RUN, to),
                                             overwrite = TRUE))

# ---- .upda files ------------------------------------------------------------
fe <- function(x) sprintf("%20.12E", x)
wrap <- function(v, k = 4) {
  v <- fe(v)
  vapply(split(v, ceiling(seq_along(v) / k)), paste, "", collapse = "")
}
mk <- fit$mk
m <- length(mk$phi)
writeLines(c(sprintf("%s%s%12.2f%7d", fe(mk$w[1]), fe(mk$w[2]), mk$abic, m),
             wrap(mk$theta), wrap(as.vector(mk$phi))),
           file.path(RUN, "s09", "hist-etas-mk.upda"))
p5 <- fit$`5pa`
m5 <- length(p5$phi)
writeLines(c(paste(fe(p5$w[1:4]), collapse = ""),
             sprintf("%s%s%10d%s", fe(p5$w[5]), fe(p5$abic), m5, fe(1000)),
             wrap(p5$theta), wrap(as.vector(p5$phi))),
           file.path(RUN, "s10", "hist-etas5pa.upda"))
# round trip through the package's own reader
u <- E$.hist_read_upda(file.path(RUN, "s10", "hist-etas5pa.upda"))
stopifnot(max(abs(u$phi - p5$phi)) < 1e-10, max(abs(u$weights / p5$w - 1)) < 1e-10)
# the port's own tessellation, in delone2.out's layout. It differs from the
# Fortran's: the epicentres agree to 5e-6 deg (input rounding), enough to flip
# the diagonal of some nearly co-circular quadrilaterals.
ms <- fit$mesh
d <- fit$data
nd <- ms$nodes
ev <- seq_len(ms$npnt)
num <- function(k, v) if (k <= ms$npnt) v[k] else 0
tri <- ms$triangles
nb <- lapply(seq_len(ms$npext), function(i) {
  r <- which(tri[, 1] == i | tri[, 2] == i | tri[, 3] == i)
  sort(setdiff(unique(as.vector(tri[r, ])), i))
})
dl <- c(sprintf("%10d%10d%10d%12.5f%12.5f", ms$npnt, ms$npext, ms$ndelo, ms$tx, ms$ty),
        sprintf("%5d%12.8f%12.8f%8.2f%16.7f%9.4f%9.4f%9.4f", seq_len(ms$npext), nd$x, nd$y,
                c(d$mag, rep(0, ms$npext - ms$npnt)), c(d$t, rep(0, ms$npext - ms$npnt)),
                c(d$dd, rep(0, ms$npext - ms$npnt)), c(d$ee, rep(0, ms$npext - ms$npnt)),
                c(d$ff, rep(0, ms$npext - ms$npnt))),
        sprintf("%7d%8d%8d%8d%20.12E", seq_len(ms$ndelo), tri[, 1], tri[, 2], tri[, 3], ms$area),
        vapply(seq_along(nb), function(i) paste(sprintf("%7d", c(i, length(nb[[i]]), nb[[i]])),
                                               collapse = ""), ""),
        sprintf("%12d", max(lengths(nb))))
writeLines(dl, file.path(RUN, "s06", "delone2.out"))
file.copy(file.path(RUN, "s06", "delone2.out"), file.path(RUN, "s10"), overwrite = TRUE)
chk <- E$.hist_mesh_compare(E$.hist_mesh_delone2(E$.hist_read_delone2(file.path(RUN, "s06", "delone2.out"))), ms)
stopifnot(nrow(chk$only.a) == 0, nrow(chk$only.b) == 0, chk$max.area.diff < 1e-10)
cp("s10/hist-etas5pa.conf", "s10/hist-etas5pa.conf")
cp("s09/hist-etas-mk.conf", "s09/hist-etas-mk.conf")
cp("s03/work.etas", "s03/work.etas")

# ---- stationary -------------------------------------------------------------
st <- fit$stationary
jsonish <- function(x) paste0("{\n", paste(sprintf('    "%s": %.17g', names(x), x),
                                          collapse = ",\n"), "\n  }")
writeLines(c("{", paste0('  "theta": ', jsonish(st$theta), ","),
             sprintf('  "loglik": %.17g,', st$loglik),
             sprintf('  "aic": %.17g,', 14 - 2 * st$loglik),
             '  "source": "ETAS::histetas() stationary stage (R port), fortran mode"', "}"),
           file.path(RUN, "s14", "stationary.json"))

# ---- maps on the 210 x 180 grid (same pixel centres as node.f/interp5pa.f) --
g <- predict(fit, type = "fields", dimyx = c(180, 210))
pix <- expand.grid(lat = g$lat, long = g$long)[, 2:1]
F <- sapply(c("mu", "K", "alpha", "p", "q"), function(k) as.vector(t(g$fields[[k]])))
F[is.na(F)] <- -999
writeLines(sprintf("%10.4f %10.4f %14.6E %14.6E %14.6E %14.6E %14.6E", pix$long, pix$lat,
                   F[, 1], F[, 2], F[, 3], F[, 4], F[, 5]),
           file.path(RUN, "s13", "interp5pa.pixel"))
snaps <- read.table(file.path(FRUN, "s13", "snapshots.txt"), sep = "\t",
                    col.names = c("tag", "day", "label"), stringsAsFactors = FALSE)
TZ <- 9676
rates <- list()
for (k in seq_len(nrow(snaps))) {
  li <- predict(fit, type = "intensity", time = snaps$day[k], dimyx = c(180, 210))
  lam <- as.vector(t(li$intensity))
  writeLines(sprintf("%10.5f %9.4f %9.4f %9.5f", snaps$day[k] - TZ, pix$long, pix$lat,
                     log10(lam)), file.path(RUN, "s13", sprintf("img1-%s.data", snaps$tag[k])))
  a <- (14 / 210) * (12 / 180)
  i <- which.max(lam)
  rates[[snaps$tag[k]]] <- list(day = snaps$day[k], regional_rate = sum(lam) * a,
                                peak = lam[i], peak_long = pix$long[i], peak_lat = pix$lat[i])
}
for (f in c("snapshots.txt", "hist-etas5pa-forecast.conf", "interp_plot.R",
            "intensity_plot.R", "drawmap.R", "mapkit.R", "node.conf"))
  cp(file.path("s13", f), file.path("s13", f))
file.copy(file.path(FRUN, "s13", "maps"), file.path(RUN, "s13"), recursive = TRUE)
file.copy(file.path(RUN, "s06", "delone2.out"), file.path(RUN, "s13"), overwrite = TRUE)
file.copy(file.path(RUN, "s10", "hist-etas5pa.upda"), file.path(RUN, "s13"), overwrite = TRUE)
# interp_plot.R takes the ABIC label from the .prt's abic= lines
writeLines(sprintf("  abic=  %.6f  -l= 0", p5$evals[, "abic"]),
           file.path(RUN, "s13", "hist-etas5pa.prt"))

# ---- Table 2: fields at the event vertices -----------------------------------
n <- fit$mesh$npnt
fv <- sapply(c("mu", "K", "alpha", "p", "q"),
             function(k) p5$theta[[k]] * exp(p5$phi[seq_len(n), k]))
tab <- lapply(colnames(fv), function(k) {
  v <- fv[, k]
  list(min = min(v), median = median(v), max = max(v), ratio = max(v) / min(v),
       weight = unname(p5$w[[k]]))
})
names(tab) <- colnames(fv)

esc <- function(x) if (is.list(x)) {
  if (is.null(names(x))) paste0("[", paste(vapply(x, esc, ""), collapse = ", "), "]") else
    paste0("{", paste(sprintf('"%s": %s', names(x), vapply(x, esc, "")), collapse = ", "), "}")
} else if (is.character(x)) sprintf('"%s"', x) else if (length(x) > 1)
  esc(as.list(x)) else sprintf("%.10g", x)
cnt <- function(f) list(evals = nrow(f$evals), simplex_lines = nrow(f$simplex$trace),
                        abic = f$abic, w = as.list(f$w), theta = as.list(f$theta))
summ <- list(stationary = list(theta = as.list(st$theta), loglik_fortran = st$loglik,
                               counts = as.list(st$counts)),
             mk = cnt(mk), `5pa` = cnt(p5), loglik_5pa_fortran = fit$loglik,
             fields_at_vertices = tab, intensity = rates,
             mesh = list(nodes = nrow(fit$mesh$nodes), triangles = nrow(fit$mesh$tri)))
writeLines(esc(summ), file.path(RUN, "r-port-summary.json"))

# ---- figures -----------------------------------------------------------------
owd <- setwd(file.path(RUN, "s13"))
for (s in c("interp_plot.R", "intensity_plot.R")) {
  out <- system2(file.path(R.home("bin"), "Rscript"), s, stdout = TRUE, stderr = TRUE)
  writeLines(out, paste0("R-", sub("_plot.R", "", s), ".log"))
  cat(out, sep = "\n")
}
setwd(owd)
