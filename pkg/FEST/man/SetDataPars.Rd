\name{SetDataPars}
\alias{SetDataPars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set input parameters for function RealStudy}
\description{
  Set input parameters for function RealStudy.
}
\usage{
SetDataPars(path, chrdirs = NULL, prefixInputFiles = NULL, suffixPed = ".pre",
 format = c("qtdt", "linkage"), famList = NULL, individualsTyped = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{The data directory.}
  \item{chrdirs}{Sub directories relative to the data directory: the
    data directories for the different chromosomes. A sub directory
    may contain data on several chromosomes. Default is NULL, i.e. no sub
    directories, the data files are in the directory given by 'path'.}
  \item{prefixInputFiles}{Prefix for int input files. If prefixInputFiles=NULL, then
  it is assumed that there are unique input files corresponding to the default suffixes:  .dat
  (linkage format), .dat, .freq, .map (qtdt format).} 
  \item{suffixPed}{Suffix for the pedigree files.}
  \item{format}{'qtdt' or 'linkage'}
  \item{famList}{Vector: familes included in the analysis}
  \item{individualsTyped}{Matrix: Number of rows equal number of
    families included in the analysis. Row i: the two individuals which
    relationship is to be tested}
}
%\details{
%  If necessary, more details than the description above 
%}
\value{
  \item{chrdirs}{Data directories for the different chromosomes}
  \item{suffixPed}{Suffix for the pedigree files.}
  \item{format}{'qtdt' or 'linkage'}
  \item{fam}{Vector: familes included in the analysis}
  \item{individualsTyped}{Matrix: Number of rows equal number of
    families included in the analysis. Row i: the two individuals which
    relationship is to be tested}
}
%\references{ put references to the literature/web site here  }
\author{Øivind Skare \email{oivind.skare@medisin.uio.no}}
%\note{ further notes 
%
%  Make other sections like Warning with \section{Warning }{....} 
%}
\seealso{\code{\link{RealStudy}}}
%\examples{
%## See \code{\link{RealStudy}}
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data}
