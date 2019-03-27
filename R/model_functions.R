# Structural functions for each model

################################################################################
# From Schelegle et al.

UOS = function(O3, Ve, t, # inputs
                DR = O3 * Ve * 1.96,
                FrDos = DR / Dos,
                CumFrDos = cumsum(rep(FrDos, length(t))),
                Dos) #parameter
    # Vectorized means for calc UOS over a number of t timepoints
    # Reflects the order of calculations and intermediate variables in
    # the Excel spreadsheet, not the Matlab code
{
    
    DR / (1 + exp(-20 * (t - (t / CumFrDos))))
}
 
deltaX = function(UOS, n_t = length(UOS), fev_base, #input
               K # parameters
               )
    # Takes vectors of UOS and calculates the change in X,
    # a scaled dFEV1, i.e., dFEV1 = X * A
    # Reflects the Excel spreadsheet calculations, not Matlab
{
    r = (1 - exp(-K))
    ceoss = UOS / K
    x = numeric(n_t)

    # Every person starts with a baseline
    x[1] = fev_base
    
    if(n_t > 1)
        for(i in 2:n_t)
            # sum as we go
            x[i] = x[i-1] + ((ceoss[i - 1] - x[i - 1]) * r)
    
    x
}

experimentFEV1 = function(O3, Ve, t_stop, Dos, K, A)
    # Calculates the dFEV1 over an entire experiment (~400 min long)
    # O3, Ve are vectors with the measurement for the time interval
    # t_stop are the stop points (in min) for each associated interval, with t_stop[1] == 1
    # This avoids having to expand the measurements into a min by min measurement
{
    dFEV1 = numeric(length(t_stop) - 1)
    x_previous = FrDos_previous = 0 # start at 0 delta FEV1
    
    if(t_stop[1] != 0)
        t_stop = c(0, t_stop)
    
    for(i in seq(length(t_stop) - 1)) {
        t = (t_stop[i]+1):(t_stop[i+1])
        dr = O3[i] * Ve[i] * 1.96
        FrDos = dr / Dos
        CumFrDos = FrDos_previous + cumsum(rep(FrDos, length(t)))
        uos = UOS(O3[i], Ve[i], t, DR = dr, CumFrDos = CumFrDos, Dos = Dos)
        tmp = deltaX(uos, fev_base = x_previous, K = K)
        FrDos_previous = CumFrDos[length(CumFrDos)]
        # browser()
        x_previous = dFEV1[i] = tmp[length(tmp)]
    }
    
    dFEV1 * A
}

################################################################################
# Not Used
if(FALSE){
cuml_integral = function(x, t)
    # Uses the trapzoid - integral approx. over a vector of x and t vals
{
    cumsum((c(0,x[-length(x)]) + x) / 2 * (diff(c(0,t))))
}


UOS = function(O3, Ve, t, # inputs
               DR = O3 * Ve * 1.96,
               CD = cumsum(rep(DR, length(t))),
               Dos) #parameter
    # Vectorized means for calc UOS over a number of t timepoints
{
    DR / (1 + exp(-20 * (t - (Dos * (t / CD)))))
}

FEV1 = function(dX, A){
    dX * A 
}

}
