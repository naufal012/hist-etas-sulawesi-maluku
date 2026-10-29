# Tectonic covariates on the stage-13 grid, one panel each:
#   (a) distance to the nearest active fault      (GEM, non-subduction)
#   (b) distance to the nearest subduction zone   (GEM subduction + USGS trenches)
#   (c) distance to the nearest Holocene volcano  (GVP)
#   (d) slab-top depth, (e) dip, (f) strike       (Slab2, upper slab where stacked)
#
# Run in <run>/tect after scripts/tectonic_context.py, with drawmap.R,
# mapkit.R and maps/ alongside.
#
# Input   covariates.pixel, structures.txt, volcanoes.txt, events-class.txt
# Output  covariates.pdf

source("drawmap.R")
source("mapkit.R")

cv <- read.table("covariates.pixel", header = TRUE)
st <- read.table("structures.txt", header = TRUE, stringsAsFactors = FALSE)
vo <- read.table("volcanoes.txt", header = TRUE, stringsAsFactors = FALSE)
ev <- read.table("events-class.txt", header = TRUE, stringsAsFactors = FALSE)
eq <- list(lon = ev$lon, lat = ev$lat)

segs <- split(st, st$seg)
is_sub <- vapply(segs, function(s) s$kind[1] != "fault", logical(1))

draw_lines <- function(which) {
  for (s in segs[which]) {
    lines(s$lon, s$lat, col = grey(1, 0.85), lwd = 1.6)
    lines(s$lon, s$lat, col = "black", lwd = 0.6)
  }
}
draw_volc <- function()
  points(vo$lon, vo$lat, pch = 24, cex = 0.55, lwd = 0.5,
         bg = "white", col = "black")

# Strike wraps at 0/360, so it gets a cyclic palette: constant chroma and
# luminance, hue going once round.
cyclic <- hcl(h = seq(0, 360, length.out = 65)[-65] + 20, c = 55, l = 62)

panels <- list(
  list(col = "d_fault",  ttl = "km", pal = hcl.colors(64, "Viridis"),
       lim = c(0, 300), extra = function() draw_lines(!is_sub)),
  list(col = "d_sub",    ttl = "km", pal = hcl.colors(64, "Viridis"),
       lim = c(0, 300), extra = function() draw_lines(is_sub)),
  list(col = "d_volc",   ttl = "km", pal = hcl.colors(64, "Viridis"),
       lim = c(0, 600), extra = draw_volc),
  list(col = "slab_dep", ttl = "km", pal = hcl.colors(64, "Viridis", rev = TRUE),
       lim = c(0, 650), extra = NULL),
  list(col = "slab_dip", ttl = "deg", pal = hcl.colors(64, "Viridis", rev = TRUE),
       lim = c(0, 70), extra = NULL),
  list(col = "slab_str", ttl = "deg", pal = cyclic,
       lim = c(0, 360), extra = NULL, at = seq(0, 360, 90))
)
names_txt <- c("distance to fault", "distance to subduction zone",
               "distance to volcano", "slab depth", "slab dip", "slab strike")

W <- 6.9
wid <- rep(c(4, 1.15), 3)
pw <- W * 4 / sum(wid) - 3.0 * 10 / 72
# The top row shares its longitude axis with the bottom row, so it carries
# ticks only and a shallower bottom margin; that is height the page needs.
top.mar <- 0.5
rowh <- panel_height(pw, diff(range(cv$lon)), diff(range(cv$lat)),
                     mar.lines = 4.2)
toph <- rowh - (2.4 - top.mar) * 10 / 72

open_pdf("covariates.pdf", width = W, height = toph + rowh)
layout(matrix(seq_len(12), nrow = 2, byrow = TRUE), widths = wid,
       heights = c(toph, rowh))

for (k in seq_along(panels)) {
  f <- panels[[k]]
  v <- pmin(pmax(cv[[f$col]], f$lim[1]), f$lim[2])
  g <- as_grid(cv$lon, cv$lat, v)
  top <- k <= 3
  panel_map(g, f$lim, f$pal, old = eq,
            corner = paste0("(", letters[k], ")"),
            xlab = if (top) "" else "longitude (E)", ylab = "latitude (N)",
            axes.x = !top,
            mar = c(if (top) top.mar else 2.4, 2.6, 1.8, 0.4))
  if (!is.null(f$extra)) { f$extra(); box(lwd = 0.8) }
  # Placed in data units, clear of the bold corner letter.
  mtext(names_txt[k], side = 3, line = 0.35, adj = 0, cex = 0.66,
        at = min(g$x) + 0.14 * diff(range(g$x)))
  panel_key(f$lim, f$pal, at = if (is.null(f$at)) pretty(f$lim, 4) else f$at,
            mar = c(if (top) top.mar else 2.4, 0.3, 1.8, 2.9))
  mtext(f$ttl, side = 3, line = 0.25, cex = 0.62, adj = 0)
}
dev.off()
cat(sprintf("covariates.pdf: %.2f x %.2f in\n", W, toph + rowh))
