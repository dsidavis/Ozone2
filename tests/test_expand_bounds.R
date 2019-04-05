bounds = list(dos = c(5, 2500), a = c(-0.2, -1e-5), k = c(1e5, 0.2))
ans = OzoneExposure:::expand_bounds(bounds, 5)

stopifnot(nrow(ans) == 5^3)
 
