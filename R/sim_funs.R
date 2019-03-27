# Simulation of data with normal error for testing

sim_data = function(O3, Ve, CD, t, #inputs
                    Dos, K, A, n, # parameters
                    sigma) #error 

{
    uos = UOS(O3, Ve, CD, t, Dos)
    x = FEV(uos, K, A)

    rnorm(n, x, sigma)       
}


sim_experiment = function(O3 = rep(rep(c(0.123, 0), each = 2), 3),
                          Ve = rep(c(30, 13), length.out = 12),
                          t_stop = c(50, 63, 113, 126, 176, 200,
                                     260, 273, 323, 336, 386, 399),
                          dos = 1100, k = 0.02, a = -0.02,
                          sdO3 = 0,
                          sdVe = 2,
                          sdFEV1 = 0.5)
{
    O3 = rnorm(length(O3), mean = O3, sd = sdO3)
    Ve = rnorm(length(Ve), mean = Ve, sd = sdVe)
    ans = experimentFEV1(O3, Ve, t_stop, dos, k, a)
    data.frame(FEV1 = rnorm(length(ans), mean = ans, sd = sdFEV1),
               O3 = O3,
               Ve = Ve,
               t_stop = t_stop)
    
}
