\name{stanmodels}
\alias{stanmodels}
\docType{data}
\title{
Compiled Stan Models for Human Ozone Exposure
}
\description{
This object contains two different models for the impact of ozone
exposure to FEV1. Each model attempts to predict the change in FEV1
(delta-FEV1) over time. 
}
\usage{data("stanmodels")}
\format{
}
\details{
There are two models available. The first, mcdonnel, is based on the
papers by McDonnel et al. The second, schelegle, is based on the papers
by Schelegle et al.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(stanmodels)
## maybe str(stanmodels) ; plot(stanmodels) ...
}
\keyword{datasets}
