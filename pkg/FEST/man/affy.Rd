\name{affy}
\alias{affy}
%\alias{affymetrix}
\docType{data}
\title{Affymetrix 500K frequency data}
\description{
   Allele frequency data from the Affymetrix 500K chip derived from the
   HapMap data. It contains frequency information
   for all 22 autosomal chromosomes.
   See also \code{\link{affy.subset}}.
   %It was downloaded    from the Chiamo home page
   %\url{http://www.stats.ox.ac.uk/\~marchini/software/gwas/chiamo.html}.
 }
\usage{## download file from 'http://folk.uio.no/thoree/FEST/affy.RData'}
%\usage{download.file("http://folk.uio.no/thoree/FEST/affy.RData", "affy.RData")}
\format{
List of length 22 where element i (affy[[i]]) is a data frame with the following columns:
\describe{
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
}}}

%\details{
%   If necessary, more details than the description above 
%}
\source{
  \url{http://www.stats.ox.ac.uk/~marchini/software/gwas/chiamo.html}
}
%\references{
%   possibly secondary sources and usages 
%}
\examples{
## download file from 'http://folk.uio.no/thoree/FEST/affy.RData'
\dontrun{load("affy.RData")}
## Average frequencies for chromosome 1:
\dontrun{f <- apply(affy[[1]][,-c(1,2)],2,mean)}
\dontrun{print(f)}

## Minor allele frequency for chromosome 1
\dontrun{maf1 <- pmin(affy[[1]][,3], affy[[1]][,4])}
## Empirical cumulative distribution function
\dontrun{plot(ecdf(maf1), xlim=c(0,0.5))}
## Some statistics
\dontrun{stat <- c(length(maf1), mean(maf1>0.01), mean(maf1>0.02), mean(maf1>0.05))}
}
\keyword{datasets}
