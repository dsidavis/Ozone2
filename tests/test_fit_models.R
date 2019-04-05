if(FALSE){
library(OzoneExposure)
    # requires the data to run
tars = list.files("~/Downloads", pattern = "*bdat$", full.names = TRUE)
d = OzoneExposure:::readSAS(tars)
dd = mungeSAS(d)

# Fit McDonnell
fit_McDonnell(dd, n_opt = 2)

# fit Schelegle
fit_Schelegle(dd, bounds = list(dos = c(900,901), a = c(-0.01,-0.011), k = c(0.002,0.0021), sigma = c(5,5.5)), n_inter = 2)

fit_Schelegle(dd, bounds = list(dos = c(900,901), a = c(-0.01,-0.011), k = c(0.002,0.0021), sigma = c(5,5.5)), n_inter = 2, cores = 2L)

}
