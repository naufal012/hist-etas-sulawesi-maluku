# Stage 13, figure 2: the five hist-etas5pa parameter fields rasterised onto
# the lon/lat grid by src/interp5pa.f.
#
# The package plots these as coloured dots at the tessellation vertices
# (Section11.5files/hist-etas5pa.R). Evaluating the piecewise-linear field on
# a grid instead shows what the model actually asserts between the events, and
# makes the smoothing that ABIC chose visible.
#
# Input   interp5pa.pixel   lon lat mu K alpha p q  (-999 outside the hull)
#         delone2.out       vertices, for the epicentre overlay
#         hist-etas5pa-forecast.conf
# Output  interp5pa.pdf

source("drawmap.R")
source("mapkit.R")

fc <- scan("hist-etas5pa-forecast.conf", quiet = TRUE, comment.char = "!")
lon0 <- fc[4]; lat0 <- fc[5]

num <- scan("delone2.out", nlines = 1, quiet = TRUE)
eqs <- matrix(scan("delone2.out", skip = 1, nlines = num[1], quiet = TRUE),
              ncol = 8, byrow = TRUE)
eq <- list(lon = eqs[, 2] + lon0, lat = eqs[, 3] + lat0)

d <- matrix(scan("interp5pa.pixel", quiet = TRUE), ncol = 7, byrow = TRUE)
d[d <= -998] <- NA

fields <- list(
  list(col = 3, log = TRUE,  ttl = expression(mu ~ "  (" * deg^-2 ~ d^-1 * ")")),
  list(col = 4, log = TRUE,  ttl = expression(K[0])),
  list(col = 5, log = FALSE, ttl = expression(alpha)),
  list(col = 6, log = FALSE, ttl = expression(p)),
  list(col = 7, log = FALSE, ttl = expression(q))
)

# Three columns of (map, key) pairs over two rows, authored at the width the
# report includes it at. The alternative shapes either scale down far enough
# to make the axis text illegible or cost most of a page.
W <- 6.9
wid <- rep(c(4, 1.15), 3)
pw <- W * 4 / sum(wid) - 3.0 * 10 / 72
rowh <- panel_height(pw, diff(range(d[, 1])), diff(range(d[, 2])),
                     mar.lines = 4.2)

open_pdf("interp5pa.pdf", width = W, height = 2 * rowh)
layout(matrix(seq_len(12), nrow = 2, byrow = TRUE), widths = wid)

for (k in seq_along(fields)) {
  f <- fields[[k]]
  v <- d[, f$col]
  if (f$log) {
    v <- log10(v)
    pal <- hcl.colors(64, "Viridis", rev = TRUE)
  } else {
    pal <- hcl.colors(64, "Spectral", rev = TRUE)
  }
  fin <- v[is.finite(v)]
  # Trim the extreme 1% at each end: the field is unconstrained in the empty
  # corners of the hull and a single runaway triangle would flatten the scale.
  zlim <- as.numeric(quantile(fin, c(0.01, 0.99), names = FALSE))
  if (diff(zlim) <= 0) zlim <- range(fin)
  g <- as_grid(d[, 1], d[, 2], pmin(pmax(v, zlim[1]), zlim[2]))
  panel_map(g, zlim, pal, title = "", old = eq,
            corner = paste0("(", letters[k], ")"),
            xlab = "longitude (E)", ylab = "latitude (N)")
  if (f$log) {
    tk <- log_ticks(zlim[1], zlim[2], max.n = 5)
    panel_key(zlim, pal, at = tk$at, labels = tk$labels)
  } else {
    panel_key(zlim, pal, at = pretty(zlim, 5))
  }
  mtext(f$ttl, side = 3, line = 0.25, cex = 0.75, adj = 0)
}
# The sixth cell pair carries the reference parameters the fields multiply,
# so the maps can be read without going back to hist-etas5pa.upda.
up1 <- scan("hist-etas5pa.upda", nlines = 2, quiet = TRUE)
up2 <- scan("hist-etas5pa.upda", skip = 2, nlines = 2, quiet = TRUE)
par(mar = c(2.4, 0.4, 1.8, 0.4))
plot.new(); plot.window(c(0, 1), c(0, 1))
txt <- c(
  sprintf("mu[0]  == %s", format(up2[1], digits = 3, scientific = TRUE)),
  sprintf("K[0]   == %s", format(up2[2], digits = 3, scientific = TRUE)),
  sprintf("c      == %s", format(up2[3], digits = 3, scientific = TRUE)),
  sprintf("alpha  == %s", format(up2[4], digits = 3)),
  sprintf("p      == %s", format(up2[5], digits = 3)),
  sprintf("d      == %s", format(up2[6], digits = 3)),
  sprintf("q      == %s", format(up2[7], digits = 3)))
text(0.02, 0.93, "reference parameters", adj = c(0, 1), cex = 0.8, font = 2)
for (i in seq_along(txt))
  text(0.06, 0.80 - 0.085 * (i - 1),
       parse(text = gsub("e([-+]?[0-9]+)", "%*% 10^{\\1}",
                         gsub(" +", " ", txt[i]))),
       adj = c(0, 1), cex = 0.78)
# The ABIC on line 2 of hist-etas5pa.upda is the simplex's LAST evaluation,
# while the coefficients and weights saved there are from its BEST one, so
# take the best ABIC from the program's own output (run_s13.sh copies it in).
abic <- up1[6]
if (file.exists("hist-etas5pa.prt")) {
  ab <- grep("^ +abic=", readLines("hist-etas5pa.prt"), value = TRUE)
  abic <- min(as.numeric(sub(" *-l=.*", "", sub(".*abic= *", "", ab))))
} else warning("hist-etas5pa.prt not found: ABIC label is the .upda's last evaluation")
text(0.02, 0.14, sprintf("ABIC = %.1f", abic), adj = c(0, 1), cex = 0.78)
text(0.02, 0.05, sprintf("%d vertices", as.integer(up1[7]) / 5),
     adj = c(0, 1), cex = 0.78)
plot.new()
dev.off()

cat(sprintf("interp5pa.pdf: %.2f x %.2f in. field ranges inside the hull:\n",
            W, 2 * rowh))
nm <- c("mu", "K0", "alpha", "p", "q")
for (i in seq_along(fields)) {
  v <- d[, fields[[i]]$col]; v <- v[is.finite(v)]
  cat(sprintf("   %-6s %12.5g .. %12.5g   median %12.5g   ratio %8.3g\n",
              nm[i], min(v), max(v), median(v), max(v) / min(v)))
}
