fit_McDonnell = function(data, sigma_u = 0.1,
                        model = OzoneExposure::stanmodels$mcdonnell,
                        n_optim = 1L, cores = 1L, ...)
{
    data = c(data, sigma_U = sigma_u)
    if(cores > 1){
        cl = makeCluster(cores, "FORK")
        on.exit(stopCluster(cl))
        ans = parSapply(cl, seq(n_optim), function(i)
            rstan::optimizing(model, data, ...)$par)
    } else 
        ans = replicate(n_optim, rstan::optimizing(model, data, ...)$par)
    
    as.data.frame(t(ans))
            
}

fit_Schelegle = function(data,
                         bounds = list(dos = c(5, 2500), a = c(-0.2, -1e-5),
                                       k = c(1e-6, 0.2), sigma = c(1e-5,5)),
                         n_interval = 50L,
                         model = OzoneExposure::stanmodels$schelegle,
                         cores = 1L, ...)
{
    grid = expand_bounds(bounds, n_interval)
       
    if(cores > 1){
        cl = makeCluster(cores, "FORK")
        on.exit(stopCluster(cl))
        ans = parApply(cl, grid, 1, function(x) 
            fit_eds(model, data, x, ...))
    } else
        ans = apply(grid, 1, function(x)
            fit_eds(model, data, x, ...))
    
    return(cbind(grid, aic = ans))
}

fit_eds = function(mod, data, pars, ...)
{
    data = c(data, pars)
    fit = rstan::sampling(mod, data, algorithm = "Fixed_param",
                          iter = 1, chains = 1, ...)
    rstan::extract(fit, pars = "aic")$aic
}


        
expand_bounds = function(b, n)
{
    if(!all(sapply(b, length) == 2))
        stop("Please specify the bounds as the upper and lower bounds for each parameter")
    if(is.null(names(b)) | !all(names(b) %in% c("dos","a","k","sigma")))
        stop("Please provide bounds with the names 'dos','a','k','sigma'")
    if(any(b$k <= 0))
        stop("The rate, 'k', cannot be less than or equal to zero")
    tmp = lapply(b, function(x)
        seq(x[1], x[2], length.out = n))
    expand.grid(tmp)
}
