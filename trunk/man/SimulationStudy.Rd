\name{SimulationStudy}
\alias{SimulationStudy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation study of family relationships}
\description{
  Assume a given set of true and alternative family relationships. Use
  simulations and exact likelihood computations to compute likelihood
  and posterior values. These computations are done for different
  number of markers.
}
\usage{
SimulationStudy(models, chr=c(1:22), nmarker = c(22 * c(1, 10, 100,
1000), 5e+05), nsim = c(1000, 1000, 1000, 1000, 400), maf = numeric(),
frequencyData = NULL, freqThreshold = c(rep(0.1, 4), 0),
saveMerlinFiles=FALSE, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{Object of type \code{Model} that specifies true and alternative family
    relations.}
  \item{chr}{Either a vector or a list. If chr is a vector it contains the chromosomes included in
    the simulation study. If it is list it should have same length as the
    nmarker vector. Each element of the list is a vector containing
    chromosomes included in the simulation study.}
  \item{nmarker}{Vector of number of markers. A simulation study is done
  for each number of markers.}
  \item{nsim}{Number of simulations. A vector with same length as the
    nmarker vector.}
  \item{maf}{Minor allele frequency. Same for all SNPs, must be
    specified if \code{frequencyData} is not specified.}
  \item{frequencyData}{A list with frequency information for each chromosome. See
    \code{\link{affy}} for a description of the format.}
  \item{freqThreshold}{Selects a sub set of the SNPs in
    \code{frequencyData}: only SNPs with minor allele frequency >
    \code{freqThreshold} are retained. Vector with same length as the
    nmarker vector.}
  \item{saveMerlinFiles}{If TRUE the files used as input to the
    likelihood computations in merlin are saved. Default value is FALSE.}
  \item{verbose}{If TRUE, information about simulations are
    output to screen. Default TRUE.}
}
%\details{
%   If necessary, more details than the description above 
%}
\value{
  An object of type \code{\link{SimStudyObject-class}}.
}
\references{
  \url{http://folk.uio.no/thoree/FEST}
  
  Øivind Skare, Nuala Sheehan, and Thore Egeland
  Identification of distant family relationships
  Bioinformatics Advance Access published on July 6, 2009.}
\author{Øivind Skare \email{oivind.skare@medisin.uio.no}}
%\note{ further notes 
%
% Make other sections like Warning with \section{Warning }{....} 
%}
\seealso{\code{\link{SetModels}}}
\examples{
  set.seed(17)
  models <- SetModels(trueModels=paste("HS-", 1:6, sep=""),
                      altModels=c("true", "unrelated"))
  nsim <- rep(10, 2)
  nmarker <- 22*c(1, 10)
  chr <- c(1:22)
  simObj1 <- SimulationStudy(models, chr=chr, nmarker=nmarker, nsim=nsim, maf=0.5)
  stat1 <- ComputeSummaryStatistics(simObj1)
  ## Average posterior results for no of markers=22:
  print(round(stat1$posterior[1,,],4)) # rows: true models, columns: alternative models
###        HS-1   HS-2   HS-3   HS-4   HS-5 HS-6 unrelated
### HS-1 0.5284     NA     NA     NA     NA   NA    0.4716
### HS-2     NA 0.4821     NA     NA     NA   NA    0.5179
### HS-3     NA     NA 0.4945     NA     NA   NA    0.5055
### HS-4     NA     NA     NA 0.5004     NA   NA    0.4996
### HS-5     NA     NA     NA     NA 0.5001   NA    0.4999
### HS-6     NA     NA     NA     NA     NA  0.5    0.5000
  ## No of markers=220:
  print(round(stat1$posterior[2,,],4))

  ## Simulation study using Affymetrix frequency data
  ## The complete 500K data file may be downloaded from
  ## 'http://folk.uio.no/thoree/FEST/affy.RData'
  ## load("affy.RData")

  ## A small subset of the Affymetrix frequency data. 100 markers on
  ## each chromsome
  data(affy.subset)
  simObj2 <- SimulationStudy(models, chr=chr, nmarker=nmarker, nsim=nsim,
                             frequencyData=affy.subset, freqThreshold=c(0.1,0.1))
  stat2 <- ComputeSummaryStatistics(simObj2)
  ## Average posterior results (Affymetrix):
  ## No of markers=22:
  print(round(stat1$posterior[1,,],4)) # rows: true models, columns: alternative models
  ## No of markers=220:
  print(round(stat1$posterior[2,,],4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{iteration}
\keyword{htest}% __ONLY ONE__ keyword per line
