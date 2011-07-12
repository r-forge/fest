\name{Pedigree-class}
\docType{class}
\alias{Pedigree-class}
\alias{initialize,Pedigree-method}

\title{Class "Pedigree"}
\description{Contains information about pedigree structure.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Pedigree", ...)}.
  Internal class used in FEST.
}
\section{Slots}{
  \describe{
    \item{\code{fam}:}{Integer vector with family numbers.}
    \item{\code{subj}:}{Integer vector with subject numbers.}
    \item{\code{fatherID}:}{Integer vector with ID of father.}
    \item{\code{motherID}:}{Integer vector with ID of mother.}
    \item{\code{gender}:}{Integer vector specifying gender (1=father, 2=mother).}
    \item{\code{affected}:}{Integer vector with affected status
      (0=unaffected, 1=affected).}
    \item{\code{typed}:}{Integer vector with indices of typed
      persons.}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "Pedigree")}: ... }
  }
}
%\references{ put references to the literature/web site here}
\author{Øivind Skare \email{oivind.skare@medisin.uio.no}}

%\note{further notes}

%Make other sections like Warning with \section{Warning }{....} 

\seealso{
  \code{\link{PlotPedigree}}
}
\examples{
showClass("Pedigree")
}
\keyword{classes}
