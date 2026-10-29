# Coastline drawing for a non-Japan region.
#
# The package's own drawmap.R reads jp.pp.dat / jp.br.dat, which cover Japan
# only and draw nothing over Indonesia. This replacement uses the world
# coastline that ships with the package (w1.pp.data + w1.br.data, found in
# Section13.4filesPlots/maps) and is clipped by the plot's own xlim/ylim.
#
# w1.pp.data is an Nx2 matrix of (lat, lon) -- latitude FIRST -- and
# w1.br.data holds the vertex count of each separate coastline segment, so
# the segments are recovered by a cumulative sum. Drawing all segments as one
# polyline would string unrelated islands together with spurious lines.
#
# drawmap.japan is defined as an alias so the package's unmodified plotting
# scripts (delo2d-bvalues.R, hist-etas-mk.R, ...) work as they are.

pb <- matrix(scan('./maps/w1.pp.data', quiet = TRUE), ncol = 2, byrow = TRUE)
pr <- scan('./maps/w1.br.data', quiet = TRUE)

pr1 <- cumsum(pr)
pr0 <- c(1, pr1[-length(pr1)] + 1)

drawmap.world <- function(...) {
  for (i in seq_along(pr1)) {
    seg <- pr0[i]:pr1[i]
    points(pb[seg, 2], pb[seg, 1], type = 'l', ...)
  }
}

# The upstream scripts call drawmap.japan(inner=TRUE/FALSE, ...); the world
# coastline has no inner/outer distinction, so absorb and ignore that arg.
drawmap.japan <- function(inner = FALSE, ...) drawmap.world(...)
