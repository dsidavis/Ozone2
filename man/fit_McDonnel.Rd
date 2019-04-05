\name{fit_McDonnell}
\alias{fit_McDonnell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit the model from McDonnell et al.
}
\description{
This will fit the model from McDonnell using optimization. The default
arguments will fit the model a single time. However, it might be
desirable to fit the model multiple times, and evaluate which of the
optimum is the global optimum. By default, Stan initializes the
optimizer at random values.
}
\usage{
fit_McDonnell(data, sigma_u = 0.1, model = OzoneExposure::stanmodels$mcdonnell, n_optim = 1L, cores = 1L, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
	a \code{list} produced by \code{mungeSAS}
}
  \item{sigma_u}{
	numeric positive, the standard deviation of the individual random effects.
}
  \item{model}{
	a compiled Stan model
}
  \item{n_optim}{
	integer, number of times to run the optimizer.
}
  \item{cores}{
	integer, number of CPU cores to use. Only functional on OSX and
	Linux (can be fixed later).
}
  \item{\dots}{
	additional args passed to \code{rstan::optimizing}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
function (data, sigma_u = 0.1, model = OzoneExposure::stanmodels$mcdonnell, 
    n_optim = 1L, cores = 1L, ...) 
{
    data = c(data, sigma_U = sigma_u)
    if (cores > 1) {
        cl = makeCluster(cores, "FORK")
        on.exit(stopCluster(cl))
        parSapply(cl, seq(n_optim), function(i) rstan::optimizing(model, 
            data, ...)$par)
    }
    else replicate(n_optim, rstan::optimizing(model, data, ...)$par)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
