# Data-exploration figure: what the catalogue looks like before any model is
# fitted to it.
#
# (a) epicentres, sized by magnitude and coloured by depth -- the depth
#     colouring is what separates the shallow arc seismicity from the
#     Molucca Sea and Sangihe slabs that project onto the same map;
# (b) magnitude against time, with the target threshold drawn in;
# (c) cumulative count of target events, whose slope is the thing the
#     background rate mu has to reproduce.
#
# Input   work.etas   from tsei2etas: id, lon, lat, mag, t(days), -depth, y, m, d
#         region.conf via --args
# Output  catalog.pdf

source("drawmap.R")
source("mapkit.R")

args <- commandArgs(trailingOnly = TRUE)
mc <- if (length(args) >= 1) as.numeric(args[1]) else 5.0
t0 <- if (length(args) >= 2) args[2] else "2000-01-01"
tstart <- if (length(args) >= 3) as.numeric(args[3]) else 730

x <- read.table("work.etas", skip = 1, header = FALSE)
names(x)[1:6] <- c("id", "lon", "lat", "mag", "t", "negdep")
x$dep <- -x$negdep
origin <- as.Date(t0)
x$date <- origin + x$t

tgt <- x$mag >= mc & x$t >= tstart

dep.br <- c(-1, 30, 70, 150, 300, 1000)
dep.pal <- hcl.colors(5, "Zissou 1")
dep.k <- cut(x$dep, dep.br, labels = FALSE)

W <- 6.9
open_pdf("catalog.pdf", width = W, height = 3.25)
layout(matrix(c(1, 2, 1, 3), nrow = 2, byrow = TRUE), widths = c(5.0, 4.4))

# ---- (a) epicentres --------------------------------------------------
par(mar = c(2.6, 2.8, 1.8, 0.6), mgp = c(1.5, 0.35, 0), tcl = -0.22)
plot.new()
plot.window(xlim = range(x$lon), ylim = range(x$lat), xaxs = "i", yaxs = "i",
            asp = 1)
drawmap.world(col = grey(0.35), lwd = 0.5)
ord <- order(x$mag)
points(x$lon[ord], x$lat[ord], pch = 21, lwd = 0.25,
       cex = 0.10 + 0.24 * (x$mag[ord] - min(x$mag)),
       col = grey(0.25), bg = adjustcolor(dep.pal[dep.k[ord]], alpha.f = 0.65))
axis(1, cex.axis = 0.72); axis(2, cex.axis = 0.72, las = 1); box(lwd = 0.8)
mtext("longitude (E)", side = 1, line = 1.35, cex = 0.62)
mtext("latitude (N)", side = 2, line = 1.6, cex = 0.62)
mtext("(a)", side = 3, line = 0.35, adj = 0, cex = 0.72, font = 2)
title(main = sprintf("%d events, M %.1f-%.1f", nrow(x), min(x$mag), max(x$mag)),
      cex.main = 0.85, font.main = 1, line = 0.4)
legend("bottomleft", bty = "n", cex = 0.6, pt.cex = 0.9, pch = 21,
       col = grey(0.25), pt.lwd = 0.3,
       pt.bg = adjustcolor(dep.pal, alpha.f = 0.75),
       legend = c("0-30 km", "30-70", "70-150", "150-300", "> 300"),
       bg = "white", inset = c(0.01, 0.01), y.intersp = 0.9)
legend("topright", bty = "n", cex = 0.6, pch = 21, col = grey(0.25),
       pt.bg = grey(0.85), pt.lwd = 0.3,
       pt.cex = 0.10 + 0.24 * (c(4, 6, 8) - min(x$mag)),
       legend = c("M 4", "M 6", "M 8"), inset = c(0.01, 0.01),
       y.intersp = 1.6)

# ---- (b) magnitude against time --------------------------------------
par(mar = c(1.4, 2.8, 1.8, 0.6))
plot(x$date, x$mag, type = "h", lwd = 0.25, col = grey(0.55),
     xlab = "", ylab = "", axes = FALSE, ylim = c(min(x$mag), max(x$mag) + 0.2))
abline(h = mc, col = "firebrick", lwd = 0.9, lty = 2)
abline(v = origin + tstart, col = "steelblue", lwd = 0.9, lty = 3)
big <- x$mag >= 7.0
points(x$date[big], x$mag[big], pch = 21, cex = 0.6, lwd = 0.4,
       col = "black", bg = "white")
axis.Date(1, x$date, cex.axis = 0.72); axis(2, cex.axis = 0.72, las = 1)
box(lwd = 0.8)
mtext("magnitude", side = 2, line = 1.6, cex = 0.62)
mtext("(b)", side = 3, line = 0.35, adj = 0, cex = 0.72, font = 2)
text(origin + tstart, max(x$mag) + 0.1, " target period starts",
     adj = 0, cex = 0.6, col = "steelblue")
text(min(x$date), mc, sprintf(" M %.1f target threshold", mc),
     adj = c(0, -0.5), cex = 0.6, col = "firebrick")

# ---- (c) cumulative target count -------------------------------------
par(mar = c(2.6, 2.8, 1.4, 0.6))
plot(x$date[tgt], seq_len(sum(tgt)), type = "l", lwd = 1.1, col = "black",
     xlab = "", ylab = "", axes = FALSE)
axis.Date(1, x$date, cex.axis = 0.72); axis(2, cex.axis = 0.72, las = 1)
box(lwd = 0.8)
mtext(sprintf("cumulative events, M >= %.1f", mc), side = 2, line = 1.6,
      cex = 0.62)
mtext("(c)", side = 3, line = 0.25, adj = 0, cex = 0.72, font = 2)
# A homogeneous Poisson process would trace the chord between the endpoints;
# the departure from it is the clustering the ETAS triggering term absorbs.
n <- sum(tgt)
lines(range(x$date[tgt]), c(1, n), col = "firebrick", lwd = 0.8, lty = 2)
legend("topleft", bty = "n", cex = 0.62, lty = c(1, 2), lwd = c(1.1, 0.8),
       col = c("black", "firebrick"),
       legend = c("observed", "constant rate"))
dev.off()

cat(sprintf("catalog.pdf: %d events total, %d target (M >= %.1f, t >= %g d)\n",
            nrow(x), n, mc, tstart))
cat(sprintf("   span %s .. %s;  depth %.0f .. %.0f km\n",
            min(x$date), max(x$date), min(x$dep), max(x$dep)))
