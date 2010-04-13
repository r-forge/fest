\name{SimStudyObject-class}
\docType{class}
\alias{SimStudyObject-class}
%\alias{SimStudyObject}
\alias{initialize,SimStudyObject-method}

%\alias{initialize,SimStudyObject-method}

\title{Class "SimStudyObject"}
\description{An object returned by the 'SimulationStudy' function. Used
  as input to \code{\link{ComputeSummaryStatistics}}.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("SimStudyObject", ...)}.
  describe objects here 
}
\section{Slots}{
  \describe{
    \item{\code{posterior}:}{Two dimensional list structure, where
      posterior[[i]][[j]] is a matrix of posterior probabilites for data
      set i and true model j. Element (k,l) of the matrix contains the
      posterior probability of alternative model k for simulation j.}
    \item{\code{logLik}:}{Same data structure as for posterior
      containing the log likelihood values.}
    \item{\code{nsim}:}{Number of simulations.}
    \item{\code{nmarker}:}{Number of markers.}
    \item{\code{maf}:}{Minor allele frequency}
    \item{\code{freqThreshold}:}{A frequency threshold that have been
      used in selecting markers. Only markers with minor allele
      frequency larger than this threshold have been considered for selection.}
    \item{\code{model}:}{Object of class \code{\link{Model-class}}.}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "SimStudyObject")}: ... }
  }
}
%\references{Reference here?}
\author{Øivind Skare \email{oivind.skare@medisin.uio.no}}
%\note{ further notes}
%
% Make other sections like Warning with \section{Warning }{....}
%
\seealso{
  \code{\link{SimulationStudy}}
  %  or \code{\linkS4class{CLASSNAME}} for links to other classes
}
%\examples{
%## See example in \code{\link{SimulationStudy}}
%}
\keyword{classes}
