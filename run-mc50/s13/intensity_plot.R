# Stage 13, figure 1: the HIST-ETAS conditional intensity on a lon/lat grid,
# at each requested instant, on one shared colour scale.
#
# Input   img1-<tag>.data   from intensity5pa.exe:
#                             col1 = t - TZ (days), col2 = lon, col3 = lat,
#                             col4 = log10 lambda  [events deg^-2 day^-1]
#         snapshots.txt     one "<tag> <day> <label>" line per snapshot,
#                           written by run_s13.sh
#         delone2.out       target events, for the epicentre overlay
#         hist-etas5pa-forecast.conf   region origin
# Output  intensity.pdf

source("drawmap.R")
source("mapkit.R")

snaps <- read.table("snapshots.txt", header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE)
names(snaps) <- c("tag", "day", "label")

fc <- scan("hist-etas5pa-forecast.conf", quiet = TRUE, comment.char = "!")
lon0 <- fc[4]; lat0 <- fc[5]

num <- scan("delone2.out", nlines = 1, quiet = TRUE)
eqs <- matrix(scan("delone2.out", skip = 1, nlines = num[1], quiet = TRUE),
              ncol = 8, byrow = TRUE)
eq.lon <- eqs[, 2] + lon0
eq.lat <- eqs[, 3] + lat0
eq.t   <- eqs[, 5]
eq.mag <- eqs[, 4]

# Events still inside this window before the instant are drawn as circles:
# beyond it the Omori term has decayed far enough that the event is no longer
# what the reader is being shown.
RECENT <- 90

grids <- lapply(snaps$tag, function(tg) {
  d <- matrix(scan(sprintf("img1-%s.data", tg), quiet = TRUE),
              ncol = 4, byrow = TRUE)
  as_grid(d[, 2], d[, 3], d[, 4])
})

# One scale for every panel, so the panels are comparable. Clip the bottom of
# the range: the lower tail is empty ocean far from any vertex and would
# otherwise absorb most of the colour resolution.
all.v <- unlist(lapply(grids, function(g) as.vector(g$z)))
all.v <- all.v[is.finite(all.v)]
zlim <- c(floor(quantile(all.v, 0.02, names = FALSE)), ceiling(max(all.v)))
# Optional fixed range, "Rscript intensity_plot.R -5 2", so figures drawn from
# different fits (retrospective vs prospective) share one colour scale.
args <- as.numeric(commandArgs(trailingOnly = TRUE))
if (length(args) == 2) zlim <- args

pal <- hcl.colors(64, "YlOrRd", rev = TRUE)
tk <- log_ticks(zlim[1], zlim[2])

np <- nrow(snaps)
# Authored at the width it is included at, so nothing is scaled down: a 9.4 in
# figure squeezed into a 6.5 in text block puts the axis labels near 4.7 pt.
W <- 6.9
PS <- 10
wid <- c(rep(4, np), 1.3)
pw <- W * 4 / sum(wid) - 3.0 * PS / 72          # plotting width of one panel
H <- panel_height(pw, diff(range(grids[[1]]$x)), diff(range(grids[[1]]$y)),
                  mar.lines = 4.2, ps = PS)

open_pdf("intensity.pdf", width = W, height = H, pointsize = PS)
layout(matrix(seq_len(np + 1), nrow = 1), widths = wid)
for (i in seq_len(np)) {
  # Only events that precede the instant may be drawn: the intensity shown is
  # conditional on exactly that history.
  before <- eq.t < snaps$day[i]
  recent <- before & eq.t > snaps$day[i] - RECENT
  panel_map(grids[[i]], zlim, pal,
            title = snaps$label[i],
            corner = paste0("(", letters[i], ")"),
            old = list(lon = eq.lon[before], lat = eq.lat[before]),
            new = list(lon = eq.lon[recent], lat = eq.lat[recent],
                       cex = 0.20 + 0.30 * pmax(0, eq.mag[recent] - 4.9)),
            axes.y = (i == 1),
            xlab = "longitude (E)",
            ylab = if (i == 1) "latitude (N)" else "")
}
panel_key(zlim, pal, at = tk$at, labels = tk$labels)
mtext(expression(lambda), side = 3, line = 0.25, cex = 0.75, adj = 0)
mtext(expression(deg^-2 ~ d^-1), side = 1, line = 1.2, cex = 0.6, adj = 0)
dev.off()

cat(sprintf("intensity.pdf: %d panels, %.2f x %.2f in, log10 lambda in [%.2f, %.2f]\n",
            np, W, H, min(all.v), max(all.v)))
for (i in seq_len(np)) {
  z <- grids[[i]]$z
  cat(sprintf("   %-5s day %9.3f  max log10 lambda %7.3f  at (%.2f, %.2f)  %d events in prior %d d\n",
              snaps$tag[i], snaps$day[i], max(z[is.finite(z)]),
              grids[[i]]$x[which(z == max(z), arr.ind = TRUE)[1, 1]],
              grids[[i]]$y[which(z == max(z), arr.ind = TRUE)[1, 2]],
              sum(eq.t < snaps$day[i] & eq.t > snaps$day[i] - RECENT), RECENT))
}
