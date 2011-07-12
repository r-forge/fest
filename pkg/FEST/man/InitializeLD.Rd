\name{InitializeLD}
\alias{InitializeLD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initialize LD information}
\description{
  Initialize LD information....
}
\usage{
InitializeLD(chunkSize=1000, ldPath="", ldFile="ld_chr22_CEU.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chunkSize}{Number of markers in each chunk.}
  \item{ldPath}{Path to linkage disequilbrium file (hapmat format).}
  \item{ldFile}{Linkage disequilbrium file.}
}
%\details{
%   If necessary, more details than the description above 
%}
\value{
  A list containing names of files containing frequency and haplotype
  data information.
  The corresponding frequency data object is a data frame with columns:
%  \describe{
    \item{\code{SNP}}{
      Name of SNP
    }
    \item{\code{cM}}{
      Distance in centi Morgan
    }
    \item{\code{A}}{
      Frequency of allele A
    }
    \item{\code{C}}{
      Frequency of allele C
    }
%  }
  The haplotype file contains a list of
  haplotype data objects needed for simulating haplotypes with
  haplosim. Each haplotype data object corresponds to a chunk of the
  chromosome. Inside each chunk, markers are assumed to be in LD,
  while markers from different chunks are assumed to be in LE.
  Each haplotype object is a list containing:
%  \describe{
    item{\code{freqs}}{
      Allele frequencies
    }
    item{\code{cor}}{
      Correlation matrix (LD coefficients)
    }
    item{\code{div}}{
      Locus-specific diversity measure
    }
    item{\code{cov}}{
      Covariance matrix for the normal distribution
    }
%  }
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
\seealso{\code{\link{SimulationStudyLD}}}
\examples{
    \dontrun{tidInit <- system.time(initFiles <- InitializeLD(chunkSize=1000,
                                                     ldPath="/home/skare/Common/Genetikk/LinkedMarkerProject/JennyData/"))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
