\name{SimulationStudyLD}
\alias{SimulationStudyLD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation study of family relationships}
\description{
  Assume a given set of true and alternative family relationships. Use
  simulations and exact likelihood computations to compute likelihood
  and posterior values. These computations are done for different
  number of markers.
}
\usage{
SimulationStudyLD(models,
                  nmarker=c(10,100,1000, Inf),
                  nsim=c(1000,1000,1000,400),
                  freqThreshold=c(rep(0.1, 4), 0.0),
                  saveMerlinFiles=FALSE,
                  verbose=TRUE,
                  frequencyFile=NULL,
                  hapdatFile=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{Object of type \code{Model} that specifies true and alternative family
    relations.}
  \item{nmarker}{Vector of number of markers. A simulation study is done
    for each number of markers.}
  \item{nsim}{Number of simulations. A vector with same length as the
    nmarker vector.}
  \item{freqThreshold}{Selects a sub set of the SNPs in
    \code{frequencyData}: only SNPs with minor allele frequency >
    \code{freqThreshold} are retained. Vector with same length as the
    nmarker vector.}
  \item{saveMerlinFiles}{If TRUE the files used as input to the
    likelihood computations in merlin are saved. Default value is FALSE.}
  \item{verbose}{If TRUE, information about simulations are
    output to screen. Default TRUE.}
  \item{frequencyFile}{Name of file to load. The file contains a
    data frame object containing frequency information
    for each marker for the chromosome. See
    \code{\link{InitializeLD}} for a description of the data frame object.
%    \describe{
%      \item{\code{SNP}}{
%	Name of SNP
%      }
%      \item{\code{cM}}{
%	Distance in centi Morgan
%      }
%      \item{\code{A}}{
%	Frequency of allele A
%      }
%      \item{\code{C}}{
%	Frequency of allele C
%}
%}
  }
  \item{hapdatFile}{Name of file to load. The file contains a list of
    haplotype data object needed for simulating haplotypes with
    haplosim. Each haplotype data object corresponds to a chunk of the
    chromosome. Inside each chunk, markers are assumed to be in LD,
    while markers from different chunks are assumed to be in LE.
    See \code{\link{InitializeLD}} for a description of the data object.
%    \describe{
%      item{\code{freqs}}{
%	Allele frequencies
%      }
%      item{\code{cor}}{
%	Correlation matrix (LD coefficients)
%      }
%      item{\code{div}}{
%	Locus-specific diversity measure
%      }
%      item{\code{cov}}{
%	Covariance matrix for the normal distribution
%      }
%    }
  }
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
\seealso{\code{\link{InitializeLD}}}
\examples{
\dontrun{set.seed(17)}
  \dontrun{models <- SetModels(trueModels=c("S-2", "HS-2", "unrelated"),
                      altModels=c("S-2", "HS-2", "unrelated"))}
  \dontrun{nmarkerLD <- c(1,10,100,1000)}
  ## Initialize LD simulation (Note: may be time consuming)
  \dontrun{initFiles <- InitializeLD(chunkSize=1000, ldPath="/home/skare/Common/Genetikk/LinkedMarkerProject/JennyData/", ldFile="ld_chr22_CEU.txt")}

  \dontrun{simObjLD <- SimulationStudyLD(models, chr=chr, nmarker=nmarker,
  nsim=nsim, frequencyFile=initFiles$file.affy, hapdatFile=initFiles$file.hapdat)}

  \dontrun{statLD <- ComputeSummaryStatistics(simObjLD)}
  ## Average posterior results for no of markers=1000:
  \dontrun{print(round(statLD$posterior[4,,],4))} # rows: true models, columns: alternative models
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{iteration}
\keyword{htest}% __ONLY ONE__ keyword per line
