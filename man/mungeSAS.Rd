\name{mungeSAS}
\alias{mungeSAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format ozone exposure data from SAS for Stan models
}
\description{
This function takes a \code{data.frame} from \code{readSAS}, and formats
it to be used by a Stan model.
}
\usage{
mungeSAS(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
	\code{data.frame} from \code{readSAS} (or similarly formatted).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list, with the following:
%%  If it is a LIST, use

\item{max_timepts}{integer, maximum number of measurement timepoints for
  all observations}
\item{max_n_dFEV1}{integer, maximum number of measured dFEV1 overall measurements}
\item{n_obs}{integer, the number of observations}
\item{n_ind}{integer, the number of individuals}
\item{n_dFEV1}{integer vector, the number of dFEV1 measurement points for each
  observation, i.e., how often dFEV1 was measured per observation}
\item{n_timepts}{integer vector, the number of other
  measurements for each observation}
\item{ind}{integer vector, with one value per observation specifying the
  individual ID for that observation}
\item{age}{numeric vector, with one value per observation specifying the
  individual's age at that observation. There is one per observation, as
  individuals age between experiments}
\item{BMI}{numeric vector, with one value per observation}
\item{BSA}{numeric vector, with one value per observation}
\item{Ve}{matrix of dims [max_timepts, n_obs], with each column
  corresponding the values for a single experimental observation. Values
  are padded with 0 where there number of observed points is less than max_timepts}
\item{Cm}{matrix, O3 measurements in similar configuration as Ve}
\item{Cs}{matrix, O3 slope measurements, similar configuration as Ve}
\item{Time}{matrix, integer values of the time, in minutes, when each
  Ve, Cm, etc. was measured. Extra spaces are padded with zeros.}
\item{dFEV1_measure_idx}{matrix, the index of which times above dFEV1
  was measured during}
\item{dFEV1}{matrix, measured delta FEV1 at each time point}
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
function (df) 
{
    tmp = split(df, paste(df$ID, df$STUDY, df$Lab, df$EXPOSURE))
    t_vars = lapply(tmp, extract_t_vars)
    ind_vars = as.data.frame(do.call(rbind, lapply(tmp, extract_ind_vars)))
    dFEV1 = lapply(tmp, extract_dFEV1)
    t_vars = collapse_results(t_vars)
    n_timepts = lapply(t_vars, function(x) apply(x, 2, function(x) sum(!is.na(x))))
    n_timepts = apply(do.call(rbind, n_timepts), 2, unique)
    t_vars = lapply(t_vars, trim_vars)
    t_vars = lapply(t_vars, function(x) {
        x[is.na(x)] = 0
        x
    })
    dFEV1 = collapse_results(dFEV1)
    dFEV1$DELFEV1[is.na(dFEV1$DELFEV1)] = 0
    n_obs = length(tmp)
    n_ind = length(unique(ind_vars$ID))
    n_dFEV1 = sapply(tmp, nrow)
    list(max_timepts = max(n_timepts), max_n_dFEV1 = max(n_dFEV1), 
        n_obs = n_obs, n_ind = n_ind, n_dFEV1 = n_dFEV1, n_timepts = n_timepts, 
        ind = to_id(ind_vars$ID), age = as.numeric(ind_vars$AGE), 
        BMI = as.numeric(ind_vars$BMI), BSA = as.numeric(apply(t_vars$BSA, 
            2, function(x) unique(x[x != 0]))), Ve = t(as.matrix(t_vars$Ve)), 
        Cm = t(as.matrix(t_vars$O3_mean)), Cs = t(as.matrix(t_vars$O3_slope)), 
        Time = t(as.matrix(t_vars$T)), dFEV1_measure_idx = t(dFEV1$TIME_ID), 
        dFEV1 = t(dFEV1$DELFEV1))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
