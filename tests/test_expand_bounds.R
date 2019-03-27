bounds = list(dos = c(5, 2500), a = c(-0.2, 0), k = c(0, 0.2))
ans = expand_bounds(bounds, 5)

stopifnot(nrow(ans) == 5^3)
 
