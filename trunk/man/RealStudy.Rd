\name{RealStudy}
\alias{RealStudy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analysis of alternative hypothesized family relations for observed marker data}
\description{
  Based on observed marker data, likelihood and posterior values are
  computed for alternative hypothesized
  family relations.
}
\usage{
RealStudy(altModels, dataPars, saveMerlinFiles, limitCentiMorgan=0,
          freqThreshold=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{altModels}{Vector of strings that gives alternative family
    relations. See \code{\link{SetModels}} for a list of the
    family relations that can be used.}
  \item{dataPars}{List: Output from function \code{\link{SetDataPars}}}
  \item{saveMerlinFiles}{If TRUE the files used as input to the
    likelihood computations in merlin are saved. Default value is
    FALSE.}
  \item{limitCentiMorgan}{Markers are thinned such
    that the distance between two consecutive markers are larger than this
    limit.}
  \item{freqThreshold}{Includes only markers with minor allele frequency larger
    than this threshold.}
}
%\details{
%  If necessary, more details than the description above 
%}
\value{
  \item{logLiks}{log likelihood values}
  \item{posterior}{posterior values}
}
%\references{Reference to future paper.. and web site..}
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
%\seealso{ objects to See Also as \code{\link{help}},  }
\examples{
## The example can not be run, user must supply data
\dontrun{pathDataDir <- "../Data/RealData/"} ## path to data directory
\dontrun{chrdirs <- paste("Chr", 1:22, sep="")}
\dontrun{dataPars <- SetDataPars(pathDataDir, chrdirs=chrdirs,
prefixInputFiles="", format="linkage", famList=32,
individualsTyped=rbind(c(32,33), c(31,33), c(30,33), c(26,33)))}
\dontrun{realObj <- RealStudy(altModels=c("S-1", "HS-1", "S-2", "S-3","unrelated"),
 dataPars)}

## A simple test example: Two persons assumed to be half-sibs, 2 markers
## Make first data files (qtdt format)
ped <- rbind(c(1, 1, 0, 0, 2, 0, 0, 0, 0),
             c(1, 2, 0, 0, 1, 0, 0, 0, 0),
             c(1, 3, 0, 0, 2, 0, 0, 0, 0),
             c(1, 4, 2, 1, 1, 1, 2, 1, 1),
             c(1, 5, 2, 3, 1, 1, 1, 2, 2))
dat <- rbind("M  locus1", "M  locus2")
freq <- rbind("M locus1", "F 0.4 0.6", "M locus2", "F 0.7 0.3")
map <- rbind("CHROMOSOME MARKER POSITION",
             "1          locus1     123.4",
             "1          locus2     136.2")
write.table(ped, file="test1.ped", col.names=FALSE, row.names=FALSE)
write.table(dat, file="test1.dat", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(freq, file="test1.freq", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(map, file="test1.map", col.names=FALSE, row.names=FALSE, quote=FALSE)

## Analysis of this data set
mypath <- "."
chrdirs <- NULL
suffixPed <- ".ped"
format <- "qtdt"
famList <- 1
individualsTyped <- cbind(4,5)
prefixInputFiles <- "test1"
dataPars <- SetDataPars(path=mypath, chrdirs=chrdirs,
                        suffixPed=suffixPed,
                        prefixInputFiles=prefixInputFiles, format=format,
                        famList=famList, individualsTyped=individualsTyped)
realObj <- RealStudy(altModels=c("HS-1", "HS-2", "HS-3",
                                 "S-1", "unrelated"), dataPars)
## realObj$posterior
##         HS-1      HS-2      HS-3        S-1 unrelated
##4-5 0.1440720 0.2432792 0.2669072 0.07134156 0.2743999

## Take away the second locus and compare three family relations that
## should give equal likelihood for 1 marker (and more generally for
## unlinked markers)

dat <- rbind("M  locus1", "S2  locus2") # S2: skips the locus
write.table(dat, file="test1.dat", col.names=FALSE, row.names=FALSE, quote=FALSE)
realObj <- RealStudy(altModels=c("HS-1", "S-1-2", "PC-2",
                                 "unrelated"), dataPars)
## realObj$logLiks
##        HS-1     S-1-2      PC-2 unrelated
##4-5 -2.448768 -2.448768 -2.448768 -2.566551
## realObj$posterior
##         HS-1     S-1-2      PC-2 unrelated
##4-5 0.2571429 0.2571429 0.2571429 0.2285714
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{iteration}
\keyword{htest}
