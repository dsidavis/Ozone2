# Functions to fit a model and tune parameters by optimization


get_FEV1_prediction = function(df, Dos, K, A)

{
    # dd = by(df, list(df$person, df$protocolNum), function(x){
        # browser()
        experimentFEV1(df$O3, df$VE, df$endTime, Dos = Dos, K = K, A = A)
    # })

    # do.call(rbind, dd)
}

get_error = function(pred, obs, fun = MSE)
{
    fun(pred, obs)
}

SSE = function(x, y)
{
    (x - y) ^ 2
}

loglik = function(x, y, sigma)
{
    if(sigma <= 0) return(NA)
    # finish this later
    -sum(dnorm(y, x, sigma, log = TRUE))
}


fit_FEV1 = function(pars, df, cost_fun = c("sse", "loglik"))
{
    pred = get_FEV1_prediction(df, Dos = pars[1], K = pars[2], A = pars[3])
    if(cost_fun == "loglik" & length(pars) != 4)
        stop("To fit with loglik, sigma must be the 4th element in pars")

    i = !is.na(df$dFEV1)
    # Just trying out the switch() fun
    switch(cost_fun,
           sse = sum(SSE(pred[i], df$dFEV1[i]), na.rm = TRUE),
           loglik = sum(loglik(pred[i], df$dFEV1[i], pars[4])))
}


fit_FEV1b = function(pars, d, cost = "sse")

{
    sum(sapply(d, function(x, pars)
        fit_FEV1(pars, df = x, cost_fun = cost), pars = pars), na.rm = TRUE)
}

