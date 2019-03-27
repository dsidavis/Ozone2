

get_XB5 = function(Cm, Cs, Vs, Ve,
                   BSA, Time,
                   B5, B6)
    # Direct translation of Stan code, from SAS code
{
    N = length(Cm)
    XB5 = numeric(N)
    
    for(j in 1:N){
        Vm = (Ve[j] ^ B6) / (BSA ^ B6)
        Ta = ifelse(j > 1, Time[j-1], 0.0)
        Tb = Time[j]
        TD = (Tb - Ta)
        XB5_previous = ifelse(j > 1, XB5[j-1], 0)
        
        XB5[j] = XB5_previous * (exp(-B5 * TD)) +
            (Cm[j] * Vm * (B5^-1)) * (1 - exp(-B5 * TD)) +
            (Cm[j] * Vs * (B5^-2)) * (((1 - B5 * Ta) * exp(-B5 * TD)) - (1 -B5 * Tb)) +
            (Cs[j] * Vm * (B5^-2)) * (((1 - B5 * Ta) * exp(-B5 * TD)) - (1 -B5 * Tb)) +
            (Cs[j] * Vs * (B5^-3)) * (((-2 + (2 * B5 * Ta) -
                                        (B5^2 * Ta^2)) * exp(-B5 * TD)) -
                                      (-2 + (2 * B5 * Tb) - (B5^2 * Tb^2)))
    }
    
    return(XB5)
}

get_pop_median = function(XB5, age_c, BMI_c, 
                          B1, B2, B3,
                          B4, B8, B9)
{
	N = length(XB5)
	XB5G = numeric(N)
    
	F1 = B1 + B2 * age_c + B8 * BMI_c
	T2 = 1 + B4

    Median = numeric(N)
	XB5G = numeric(N)
    XB5G[XB5 >= B9] = XB5[XB5 >= B9] - B9
    
	T1 = 1 + B4 * exp(-B3 * XB5G)


    Median = F1 * (1/T1 - 1/T2)
	
	return(Median)
}


