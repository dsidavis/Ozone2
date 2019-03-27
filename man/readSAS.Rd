\name{readSAS}
\alias{readSAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Read a collection of SAS formatted Ozone experiment files into R.
}
\description{
  \code{readSAS} is a convienence function to read experimental data
  from O3 exposure experiments into R.
}
\usage{
readSAS(files)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
	character vector, a collection of the file paths to the SAS files to
	be read.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{\code{data.frame}}{
	A \code{data.frame} with the files all together.
  }
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
function (files) 
{
    d = lapply(files, function(x) as.data.frame(haven::read_sas(x)))
    do.call(rbind, d)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
