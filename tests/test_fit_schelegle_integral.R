if(FALSE){


    library(OzoneExposure)
    # requires the data to run
tars = list.files("~/Downloads", pattern = "*bdat$", full.names = TRUE)
d = OzoneExposure:::readSAS(tars)
dd = mungeSAS(d)
    system.time({ff = fit_Schelegle(dd, model = stanmodels$schelegle2, cores = 3, n_interval = 10)})
    system.time({f2 = fit_Schelegle(dd, cores = 3, n_interval = 10)})
    
#     ff = fit_Schelegle_integral(dd, verbose = TRUE)

#     ff = rstan::optimizing(stanmodels$schelegle2, c(dd, sigma_U = 0.1), verbose = TRUE)
#     ff = rstan::sampling(stanmodels$schelegle2, c(dd, sigma_U = 0.1), verbose = TRUE)

# # fit Schelegle
# fit_Schelegle(dd, bounds = list(dos = c(900,901), a = c(-0.01,-0.011), k = c(0.002,0.0021), sigma = c(5,5.5)), n_inter = 2)

# fit_Schelegle(dd, bounds = list(dos = c(900,901), a = c(-0.01,-0.011), k = c(0.002,0.0021), sigma = c(5,5.5)), n_inter = 2, cores = 2L)

}
