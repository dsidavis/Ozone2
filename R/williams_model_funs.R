################################################################################
# From McDonnell et al, 2013

get_XB5 = function(Cm, Cs, Vs, BSA, Time)
{
    
    # Cm = Ozone concentration
    # Cs = O3_slope
    # Vs = 0
    # Vm = (Ve[J] ^ B6) / (BSA[J]^B6

    XB5 = rep(0, Time)
    
    for(J in 1:Time) {
        Vm = (Ve[J] ^ B6) / (BSA[J]^B6)
        Ta = T[J - 1]
        Tb = T[J]
        TD = (Tb - Ta)
        
        XB5[J] = XB5[J - 1] * ( exp( -B5 * TD ) ) +
            ( Cm * Vm * (B5^-1) ) * ( 1 - exp(-B5 * TD) ) +
            ( Cm * Vs * (B5^-2) ) * (( ( 1 - B5 * Ta ) * exp(-B5 * TD) ) - (1-B5*Tb)) +
            ( Cs * Vm * (B5^-2) ) * (( ( 1 - B5 * Ta ) * exp(-B5 * TD) ) - (1-B5*Tb)) +
            ( Cs * Vs * (B5^-3) ) * (( (-2 + (2 * B5 * Ta) -
                                        (B5^2 * Ta^2)) * exp(-B5 * TD)) -
                                     (-2 + (2 * B5 * Tb) - (B5^2 * Tb^2)))
    }
    return(XB5)
}

predict_dFEV1 = function(XB5, B1, B2, B3, B4, B8, B9)
{
    XB5G = (XB5 - B9) * (B9 <= XB5)

    F1 = B1 + B2 * (age-23.8) + B8 * (BMI - 23.1)
    T1 = 1 + B4 * exp(-B3 * XB5G)
    T2 = 1 + B4
    Median = F1 * (1/T1 - 1/T2)
    return(Median)
}

################################################################################
# Not used
if(FALSE){
FEV2 = function(Ui, M)
    # Non-linear form of the model
    # The E term is error, modeled elsewhere
{
    exp(Ui) * M
}


M = function(b1, b2, b3, b4, X, age, bmi)
    # Calculates the M values with covars
{
    (b1 + b2 * (age - 23.8) + b8 * (bmi - 23.1)) / (1 + b4 * exp(-b3 * X)) -
        (b1 + b2 * (age -23.8) + b8 * (bmi - 23.1)) / (1 + b4)
}

              
}
