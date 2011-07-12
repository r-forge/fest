\name{PlotPedigree}
\alias{PlotPedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot of a pedigree}
\description{
  Given a family relation or a pedigree file, the function plots the
  pedigree.
}
\usage{
PlotPedigree(reltype, plotfile, devtype = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reltype}{Character string specifying the family relation. Either
  'reltype' or 'pedfile' must be specified. 'reltype' must be one of the
  family relations listed in \code{\link{SetModels}}.}
  \item{plotfile}{The plot is written to this file if 'devtype' is
    different from NULL.}
  \item{devtype}{Device type. Possible values: 'pdf', 'postscript',
    NULL. If NULL, the plot is not written to file 'plotfile'.}
}
%\details{
%   If necessary, more details than the description above 
%}
%\value{
%  Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
%}
%\references{ put references to the literature/web site here  }
\author{Øivind Skare \email{oivind.skare@medisin.uio.no}}
%\note{ further notes 
%
% Make other sections like Warning with \section{Warning }{....} 
%}
%\seealso{ objects to See Also as \code{\link{help}},  }
\examples{
PlotPedigree("S-1-4")
PlotPedigree("HS-3-4")
PlotPedigree("PC-3")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
