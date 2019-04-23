# Testing the Stan version of Ed's model

library(OzoneExposure)

time = c(1, 50, 60, 110, 120, 170, 215, 265, 275, 325, 335, 385, 399)
ve = c(38.3, rep(c(38.3, 15.0792453), 6))
o3 = c(0, 0.0793933333, 0.0793933333, 0.07981, 0.07981,
       0.08029, 0.08029,
       0.07958, 0.07958,
       0.0803966666666, 0.0803966666666,
       0.07977333333333, 0.07977333333333) 

expected_ans = c(0, 0, 0, 0, 0, 0, 0, -1.72087088354959e-94, -1.90017077526133e-49, 
-4.43485912124668, -4.19394040921599, -6.22098178707323, -5.47539345520911
)

data_list = list(max_timepts = length(time),
         max_n_dFEV1 = length(time),
         n_obs = 1L,
         n_ind = 1L,
         n_dFEV1 = as.array(length(expected_ans)),
         n_timepts = as.array(length(time)),
         ind = as.array(1L),
         age = 23,
         BMI = 23,
         BSA = 5L, # Not used by schelegle 
         Ve = t(as.matrix(ve)),
         Cm = t(as.matrix(o3)),
         Cs = t(as.matrix(o3)), # not used by schelegle
         Time = t(as.matrix(time)),
         dFEV1_measure_idx = t(as.matrix(seq(length(time)))),
         dFEV1 = t(as.matrix(-1 * expected_ans + rnorm(length(expected_ans), 0, 0.1))),
         k = log(log(2)/35),
         a = log(0.0246633),
         dos = log(1400)
         )

# rstan::stan_rdump(names(data_list), file = "test_schelegle.rdump", envir = as.environment(data_list))
mod = rstan::stan_model("../src/stan_files/schelegle2.stan")
rstan::sampling(mod, data = data_list, verbose = TRUE)#, init = list(sigma = 0.1, U_a = as.array(1e-5), U_dos = as.array(1e-5), U_a = as.array(1e-5)), chains = 1)
