\name{ComputeSummaryStatistics}
\alias{ComputeSummaryStatistics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes summary likelihood statistics}
\description{
  The function computes summary likelihood statistics for an object of type
  \code{\link{SimStudyObject-class}} that is returned from the function \code{\link{SimulationStudy}}.
}
\usage{
ComputeSummaryStatistics(simObj, altHyp = "unrelated")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simObj}{Object of type \code{\link{SimStudyObject-class}}}
  \item{altHyp}{The alternative hypothesis to be used when doing likelihood
    ratio statistics}
}
%\details{
%   If necessary, more details than the description above 
%}
\value{
  \item{posterior}{Three dimensional array, where
    element (i,j,k) is the average posterior probability for simulation
    sub study i,
    true model j and alternative model k. The sub studies correspond to
    the different number of markers specified in \code{\link{SimulationStudy}}
  }
  \item{posterior.sd}{Standard deviation of the average posterior
    values. Same data structure as for posterior.}
  \item{likrat}{Average likelihood ratio between true models and \code{altHyp}. Same data structure as for posterior.}
  \item{likrat.sd}{Standard deviation of the average likelihood ratio. Same data structure as for posterior.}
  \item{classrat}{Average number of times the likelihood ratio is positive. Same data structure as for posterior.}
  \item{classrat.sd}{Standard deviation of the \code{classrat} statistics. Same data structure as for posterior.}
}
%\references{ put references to the literature/web site here  }
\author{Øivind Skare \email{oivind.skare@medisin.uio.no}}
%\author{ who you are }
%\note{ further notes 
%
% Make other sections like Warning with \section{Warning }{....} 
%}
\seealso{\code{\link{SimStudyObject-class}}, \code{\link{SimulationStudy}}}
%\examples{
%## See example in \code{\link{SimulationStudy}}
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ kwd1 }
%\keyword{ kwd2 }% __ONLY ONE__ keyword per line
