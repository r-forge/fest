%\name{Model-class}
\name{Model-class}
\docType{class}
\alias{Model-class}
%\alias{Model}

\title{Class "Model"}
\description{Contains true and alternative pairwise family relations. An object
  of class \code{Model} is used as input to
  \code{\link{SimulationStudy}}.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Model", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{true}:}{True family relations.}
    \item{\code{alternative}:}{Alternative family relations. Specfied for
    each true family relation.}
  }  
}
\section{Methods}{
  No methods defined with class "Model" in the signature.
}
\details{
  Models are set by the function \code{\link{SetModels}}.
}
%\references{
%  Reference to future paper..
%}
\author{Oivind Skare \email{oivind.skare@medisin.uio.no}}
%\note{ further notes }

%Make other sections like Warning with \section{Warning }{....} 

\seealso{
  \code{\link{SetModels}}, \code{\link{SimulationStudy}}
  %	or \code{\linkS4class{CLASSNAME}} for links to other classes
}
\examples{
  models1 <- SetModels(trueModels=paste("HS-", 1:6, sep=""),
                       altModels=c("true", "unrelated"))
  modelList <- paste("HS-", 1:5, sep="")
  models2 <- SetModels(trueModels=modelList, altModels=c(modelList, "unrelated"))
}
\keyword{classes}
