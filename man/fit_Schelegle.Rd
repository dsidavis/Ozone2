\name{fit_Schelegle}
\alias{fit_Schelegle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit the model from Schelegle et al.
}
\description{
This function 'fits' the model from Schelegle et al papers. Since the
number of parameters is relatively small and the response surface has
many local optimum, this function conducts a grid search.
}
\usage{
fit_Schelegle(data, bounds = list(dos = c(5, 2500), a = c(-0.2, 0), k = c(0, 0.2), sigma = c(0.01, 5)), n_interval = 50L, model = OzoneExposure::stanmodels$schelegle, cores = 1L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
list, result of \code{mungeSAS}
}
  \item{bounds}{
list, with the upper and lower bounds for each parameter ("dos", "a",
"k", and "sigma"). There must only two numbers for each parameter, and
each parameter must be named.
}
  \item{n_interval}{
integer, the number of intervals to create between each parameter's
upper and lower bound. 
}
  \item{model}{
a compiled Stan model
}
  \item{cores}{
the number of CPU cores to use. On OSX and Linux, any value > 1 will
result in parallel computations. Currently, parallel ops are not
available under Windows.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A data.frame with
  \item{dos}{}
  \item{k}{}
  \item{a}{}
  \item{sigma}{}
  \item{aic}{the associated AIC value for each parameter combination}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Matt Espe
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, bounds = list(dos = c(5, 2500), a = c(-0.2, 0), 
    k = c(0, 0.2), sigma = c(0.01, 5)), n_interval = 50L, model = OzoneExposure::stanmodels$schelegle, 
    cores = 1L) 
{
    grid = expand_bounds(bounds, n_interval)
    if (cores > 1) {
        cl = makeCluster(cores, "FORK")
        on.exit(stopCluster(cl))
        ans = parApply(cl, grid, 1, function(x) fit_eds(model, 
            data, x))
    }
    else ans = apply(grid, 1, function(x) fit_eds(model, data, 
        x))
    return(cbind(grid, aic = ans))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
