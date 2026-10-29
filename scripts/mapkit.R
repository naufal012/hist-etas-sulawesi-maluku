# Shared map/panel helpers for the stage-13 figures.
#
# The package's own figures are built on filled2.contour (f2.R), which calls
# layout() internally and therefore cannot be combined with mfrow to place
# several maps on one page. These helpers split the two jobs -- draw a raster
# panel, draw a colour key -- so a page can be composed with a single layout()
# and every panel can share one scale.
#
# Expects drawmap.R (world coastline) to have been sourced already.

# cairo_pdf embeds the fonts it uses; the default pdf() device only
# references the base-14 set, so a plotmath mu or lambda arrives in the
# report as whatever Symbol substitute the reader's viewer happens to pick.
open_pdf <- function(file, width, height, pointsize = 10) {
  if (isTRUE(capabilities("cairo"))) {
    cairo_pdf(file, width = width, height = height, pointsize = pointsize,
              onefile = TRUE)
  } else {
    pdf(file, width = width, height = height, pointsize = pointsize)
  }
}

# Height in inches that makes a panel of `pw` inches of plotting width come
# out at the given degree extent with square degrees, plus its margins.
panel_height <- function(pw, dlon, dlat, mar.lines = 4.6, ps = 10) {
  pw * dlat / dlon + mar.lines * ps / 72
}

# Turn columns (lon, lat, v) into a list(x, y, z) grid.
#
# The coordinates arrive through a Fortran f10.4 field, so the spacing is
# regular only to four decimals and image(useRaster=TRUE) rejects it. The
# cell INDEX is recovered from the values as written, then the axes are
# rebuilt as exact arithmetic sequences.
as_grid <- function(lon, lat, v) {
  ux <- sort(unique(lon))
  uy <- sort(unique(lat))
  nx <- length(ux); ny <- length(uy)
  z <- matrix(NA_real_, nrow = nx, ncol = ny)
  z[cbind(match(lon, ux), match(lat, uy))] <- v
  list(x = seq(ux[1], ux[nx], length.out = nx),
       y = seq(uy[1], uy[ny], length.out = ny),
       z = z)
}

# One map panel: raster, coastline, optional epicentres, box and axes.
#
# `old` and `new` are both epicentre sets: `old` is drawn as faint dots for
# spatial context, `new` as open circles sized by magnitude. In the intensity
# figure `new` holds the events recent enough to still be driving the
# triggering term, which is what the reader needs to connect map to history.
panel_map <- function(g, zlim, pal, title = "", old = NULL, new = NULL,
                      xlab = "", ylab = "", axes.x = TRUE, axes.y = TRUE,
                      coast.col = grey(0.1), coast.lwd = 0.5,
                      corner = "", mar = c(2.4, 2.6, 1.8, 0.4)) {
  par(mar = mar, mgp = c(1.5, 0.35, 0), tcl = -0.22)
  plot.new()
  plot.window(xlim = range(g$x), ylim = range(g$y), xaxs = "i", yaxs = "i")
  image(g$x, g$y, g$z, col = pal, zlim = zlim, add = TRUE, useRaster = TRUE)
  drawmap.world(col = coast.col, lwd = coast.lwd)
  if (!is.null(old) && length(old$lon))
    points(old$lon, old$lat, pch = 16, cex = 0.12, col = rgb(0, 0, 0, 0.30))
  if (!is.null(new) && length(new$lon))
    points(new$lon, new$lat, pch = 1, cex = new$cex, col = "black", lwd = 0.7)
  if (axes.x) axis(1, cex.axis = 0.72) else axis(1, labels = FALSE)
  if (axes.y) axis(2, cex.axis = 0.72, las = 1) else axis(2, labels = FALSE)
  box(lwd = 0.8)
  if (nzchar(title))
    title(main = title, cex.main = 0.85, font.main = 1, line = 0.4)
  if (nzchar(corner))
    mtext(corner, side = 3, line = 0.35, adj = 0, cex = 0.72, font = 2)
  if (nzchar(xlab)) mtext(xlab, side = 1, line = 1.35, cex = 0.62)
  if (nzchar(ylab)) mtext(ylab, side = 2, line = 1.55, cex = 0.62)
}

# Vertical colour key. `at`/`labels` in the same units as zlim.
panel_key <- function(zlim, pal, at = pretty(zlim), labels = at, title = "",
                      mar = c(2.4, 0.3, 1.8, 2.9)) {
  par(mar = mar, mgp = c(1.5, 0.35, 0), tcl = -0.22)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = zlim, xaxs = "i", yaxs = "i")
  br <- seq(zlim[1], zlim[2], length.out = length(pal) + 1)
  rect(0, br[-length(br)], 1, br[-1], col = pal, border = NA)
  box(lwd = 0.8)
  keep <- at >= zlim[1] & at <= zlim[2]
  axis(4, at = at[keep], labels = labels[keep], las = 1, cex.axis = 0.7)
  if (nzchar(title)) mtext(title, side = 3, line = 0.3, cex = 0.7, adj = 0)
}

# Tick positions/labels for a base-10 log scale over [lo, hi] in log10 units.
# Thinned to at most `max.n` decades so the key stays readable.
log_ticks <- function(lo, hi, max.n = 6) {
  e <- seq(ceiling(lo), floor(hi))
  if (length(e) == 0) e <- round((lo + hi) / 2)
  step <- max(1, ceiling(length(e) / max.n))
  e <- rev(rev(e)[seq(1, length(e), by = step)])
  list(at = e, labels = parse(text = paste0("10^{", e, "}")))
}
