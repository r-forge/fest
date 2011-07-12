\name{ThinMerlinInputFiles}
\Rdversion{1.1}
\alias{ThinMerlinInputFiles}
\title{Thinning of Merlin Input files}
\description{
  Merlin input files are thinned based on thresholds for minor allele freqency and minimal distance (in cM) between markers.

}
\usage{
ThinMerlinInputFiles(mapfile, datfile, freqfile, pedfile, nNotMarker = 5,
                     limitCentiMorgan = 0, freqThreshold = 0, suffix = "_thinned")
}

\arguments{
  \item{mapfile}{
File of the type filename.map which contains the genetic position of each marker (chromosome and position in cM).

}
  \item{datfile}{
File of the type filename.dat which contains the marker identifiers (e.g. rsIDs).
}
  \item{freqfile}{
File of the type filename.freq which contains the allele frequencies of each marker.

}
  \item{pedfile}{
File of the type filename.ped which contains family information and genotype data for all individuals.

}
  \item{nNotMarker}{
    Number of columns in pedfile that precedes the marker
    columns. Default value equal 5. %% COMMENT: This is ad-hoc and
    %% a preliminary solution, I should
    %% find a way to automatically find this number from the information
    %% stored in the datfile.
}
  \item{limitCentiMorgan}{
The data set is thinned such that the distance between two consecutive markers is larger than this limit. Default value is 0.
}
  \item{freqThreshold}{
Markers with a minor allele frequency smaller or equal to this limit will be removed. Default value is 0. 
}
 \item{suffix}{
Suffix that will be added to the filenames of the thinned files.
}
}
\author{
Øivind Skare \email{oivind.skare@medisin.uio.no}
}

\seealso{
  \code{\link{SortMerlinInputFiles}}
}
\examples{
## The example can not be run, user must supply data.
## Replace "file" with the name of your files.

## nMarkerThinned <- ThinMerlinInputFiles("file.map", "file.dat", "file.freq", "file.ped",
##                                         limitCentiMorgan=0.1, freqThreshold=0.1,
##                                         suffix="_thinned")
## file.rename("file.map_thinned", "file2.map") 
## file.rename("file.dat_thinned", "file2.dat") 
## file.rename("file.freq_thinned", "file2.freq") 
## file.rename("file.ped_thinned", "file2.ped")  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{data}% __ONLY ONE__ keyword per line
