\name{SetModels}
\alias{SetModels}
\title{Specify true and alternative family relations}
\description{
  Specify true and alternative family relations.
}
\usage{
SetModels(trueModels, altModels = c("true", "unrelated"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trueModels}{a character array of the true family relations. See details.}
  \item{altModels}{a character array of alternative hypothised family relations. See details.}
}
\details{
  The true and alternative pairwise family relations are among the following:
  \describe{
    \item{\code{HS-k-l}:}{Generalised half-sib
      relationships where k=1, 2, ... and l=1, 2,... 'HS-k-l'
      specifies a relation where two persons sharing one common
      ancestor; one person is k generation apart and the other is
      l generation apart. Examples: 'HS-1-1' is a
      half-sib relation, 'HS-2-2' specifies two persons
      sharing a grandparent.}
    \item{\code{HS-k}:}{Same as HS-k-k}
    \item{\code{S-k-l}:}{Generalised sibling relationships where
      k=1, 2 ... and l=1, 2,... 'S-k-l'
      specifies a relation where two persons are sharing one common
      ancestor; one person is k+1 generation apart and the other is
      l+1 generation apart. Examples: 'S-1-1'
      specifies sibs, 'S-2-2' first cousins, 'S-1-2' an
      uncle-niece relation.}
    \item{\code{S-k}:}{Same as S-k-k}
    \item{\code{PC-k}:}{Generalised sibling relationships where
      k=1, 2 ... 'PC-k' specifies a relation where the
      second person is the ancestor of first person, k generations
      apart. Examples: 'PC-1' specifies
      a parent-child relation, 'PC-2' a
      grandparent-grandchild relation.}
    \item{\code{unrelated:}}{Two unrelated persons}
  }
  In addition \code{altModels} may also contain the following
  specifications:
  \describe{
    \item{\code{true}:}{Same models as specified in \code{trueModels}}
    \item{\code{upper}:}{A generation-less relation than the true
      model specified in \code{trueModels}}
    \item{\code{lower}:}{A generation-more relation than the true
      model specified in \code{trueModels}}
  }
  The available family relations are constrained by default setting in
  merling for maximal allowed pedigree complexity allowed by merlin
  \url{http://www.sph.umich.edu/csg/abecasis/Merlin/index.html}. For the
  moment we allow HS-k-l and S-k-l where k+l<= 25, and
  PC-k where k <= 25.
}
\value{
  An object of \code{class} \code{\linkS4class{Model}}.
}
%\references{Reference to future paper.. and web site..}
\author{Øivind Skare \email{oivind.skare@medisin.uio.no}}
%\note{ further notes 
%
% Make other sections like Warning with \section{Warning }{....} 
%}
\seealso{
  \code{\linkS4class{Model}}
}

\examples{
## Example 1
  modelList <- paste("HS-", 1:5, sep="")
  models <- SetModels(trueModels=modelList, altModels=c(modelList, "unrelated"))
## Example 2
  models <- SetModels(trueModels=paste("HS-", 1:6, sep=""),
                      altModels=c("true", "unrelated"))
}
\keyword{models}
