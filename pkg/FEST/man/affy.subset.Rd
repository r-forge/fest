\name{affy.subset}
\alias{affy.subset}
%\alias{affymetrix}
\docType{data}
\title{Subset of Affymetrix 500K frequency data}
\description{
   A subset of the allele frequency data from the Affymetrix 500K chip derived from the
   HapMap data. It contains frequency information of 100 markers
   for all 22 autosomal chromosomes.
   %It was downloaded    from the Chiamo home page
   %\url{http://www.stats.ox.ac.uk/\~marchini/software/gwas/chiamo.html}.
 }
\usage{data(affy.subset)}
\format{
  See \code{\link{affy}}.}

%\details{
%   If necessary, more details than the description above 
%}
%\source{
%  \url{http://www.stats.ox.ac.uk/~marchini/software/gwas/chiamo.html}
%}
%\references{
%   possibly secondary sources and usages 
%}
\examples{
data(affy.subset)
## Average frequencies for chromosome 1:
f <- apply(affy.subset[[1]][,-c(1,2)],2,mean)
print(f)

## Minor allele frequency for chromosome 1
maf1 <- pmin(affy.subset[[1]][,3], affy.subset[[1]][,4])
## Empirical cumulative distribution function
plot(ecdf(maf1), xlim=c(0,0.5))
## Some statistics
stat <- c(length(maf1), mean(maf1>0.01), mean(maf1>0.02), mean(maf1>0.05))
}
\keyword{datasets}
