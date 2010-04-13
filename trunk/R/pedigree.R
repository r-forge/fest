setMethod("initialize", "Pedigree", function(.Object, ...) {
  .args <- list(...)
  if ((nargs() == 2 && hasArg(gender)) ||
      (nargs() == 3 && hasArg(gender) && hasArg(typed))) {
    nSubj <- length(.args$gender)
    .Object@fam <- as.integer(rep(1, nSubj)) # default family number: 1
    .Object@subj <- as.integer(1:nSubj) # subjects numbered 1, 2, ..., nSubj
    .Object@gender <- as.integer(.args$gender)
    if (!hasArg(typed))
      .Object@typed <- integer()
    else
      .Object@typed <- as.integer(.args$typed)
    .Object@affected <- as.integer(rep(0, nSubj)) # affected default = 0
    .Object@affected[.Object@typed] <- as.integer(1)
    .Object@fatherID <- as.integer(rep(0, nSubj))
    .Object@motherID <- as.integer(rep(0, nSubj))
    validObject(.Object)
  }
  else {
    callNextMethod()
  }
  .Object
})


setMethod("print", "Pedigree", function(x) {
  df <- cbind(fam=x@fam, subj=x@subj, fatherID=x@fatherID, motherID=x@motherID,
              gender=x@gender, affected=x@affected)
  print(df)
})


## aD-hoc when setMethod("initialize", ...) does not work
InitializePed <- function(.Object)
{
  nSubj <- length(.Object@gender)
  .Object@fam <- as.integer(rep(1, nSubj)) # default family number: 1
  .Object@subj <- as.integer(1:nSubj) # subjects numbered 1, 2, ..., nSubj
  .Object@affected <- as.integer(rep(0, nSubj)) # affected default = 0
  .Object@affected[.Object@typed] <- as.integer(1)
  .Object@fatherID <- as.integer(rep(0, nSubj))
  .Object@motherID <- as.integer(rep(0, nSubj))
  .Object
}

##setGeneric("writePed", function(ped, pedfile) {standardGeneric("writePed")})
##setMethod("writePed", signature(ped="pedigree", pedfile="character"), function(ped, pedfile) {
WritePed <- function(ped, pedfile=NULL) {
  pedTable <- cbind(ped@fam, ped@subj, ped@fatherID, ped@motherID, ped@gender, ped@affected)
  
  if (!is.null(pedfile))
    write.table(pedTable, file=pedfile, quote=FALSE, col.names=FALSE, row.names=FALSE)
  invisible(pedTable)
}
          
CheckClass <- function() {
  gender <- as.integer(c(2,1,2,1,1))
  nSubj <- length(gender)
  ped <- new("Pedigree", fam=as.integer(rep(1, nSubj)), subj=as.integer(1:nSubj),
             fatherID=as.integer(rep(0, nSubj)), motherID=as.integer(rep(0, nSubj)),
             gender=gender, affected=as.integer(c(0,0,0,1,1)))
  validObject(ped)
  ped <- new("Pedigree", gender=gender, typed=c(nSubj-1, nSubj))
  validObject(ped)
  print(ped)
  WritePed(ped, "test.dat")
  ped
}

MakePedigree <- function(mod, pedfile=NULL) {
  if (!(mod %in% availableModels))
    stop("Family relation ", mod, " not available for analysis\n")

  if (substring(mod, 1, 2) == "HS") {
    sp <- unlist(strsplit(mod, "-"))
    gen <- as.integer(sp[-1])
    if (length(gen) == 1) {
      gen <- rep(gen, 2)
    }
    ped <- MakePedigreeGeneralisedHalfSib(gen[1], gen[2])
  }
  else if (substring(mod, 1, 1) == "S") {
    sp <- unlist(strsplit(mod, "-"))
    gen <- as.integer(sp[-1])
    if (length(gen) == 1) {
      gen <- rep(gen, 2)
    }
    ped <- MakePedigreeUncleNiece(gen[1], gen[2])
  }
  else if (substring(mod, 1, 2) == "PC") {
    sp <- unlist(strsplit(mod, "-"))
    gen <- as.integer(sp[-1])
    ped <- MakePedigreeParentChild(gen)
  }
   else if (mod == "unrelated") {
    ped <- MakePedigreeNotRelated(gender=c(1,2))
  }

  if (!is.null(pedfile)) {
    WritePed(ped, pedfile)
  }
  
  ped
}
          
MakePedigreeNotRelated <- function(gender) {
  nSubj <- length(gender)
  ped <- new("Pedigree", gender=as.integer(gender), typed=as.integer(c(nSubj-1, nSubj)))
  
  invisible(ped)
}

## Make pedigree of not related individuals without marker data
##   gender: Gender of individuals
##   pedfile: the pedigree is written to this file

MakePedigreeHalfSib <- function(extraGenerations, genderVersion=2) {
  if (genderVersion == 1)
    gender <- c(2,1,2, rep(c(2,1,1,2), extraGenerations-1), 1,1) # assumed pedigree pattern
  else
    gender <- c(2,1,2, rep(c(2,1,1,2), extraGenerations-1), 1,2) # assumed pedigree pattern

  nSubj <- length(gender)
  ped <- new("Pedigree", gender=as.integer(gender), typed=as.integer(c(nSubj-1, nSubj)))
  ##  ped <- InitializePed(ped)
  
  fatherID <- c(0, 0, 0) ## ID's of first generation
  motherID <- c(0, 0, 0)
  gens <- seq(from=0, by=1, len=extraGenerations-1)
  father <- c(2,2)
  mother <- c(1,3)
  for (i in gens) {    
    motherID <- c(motherID, c(0, mother, 0))
    mother <- 4*(i+1) + c(0,3)

    fatherID <- c(fatherID, c(0, father, 0))
    father <- 4*(i+1) + c(1,2)
  }
  motherID <- c(motherID, mother)
  fatherID <- c(fatherID, father)

  ped@motherID <- as.integer(motherID)
  ped@fatherID <- as.integer(fatherID)
  

  invisible(ped)
}

MakePedigreeSibling <- function(extraGenerations) {
  gender <- c(1,2, rep(c(1,2,1,2), extraGenerations-1), 1,2) # assumed pedigree pattern

  nSubj <- length(gender)
  ped <- new("Pedigree", gender=as.integer(gender), typed=as.integer(c(nSubj-1, nSubj)))
  ##  ped <- InitializePed(ped)

  ## determine mother and father id's
  fatherID <- c(0, 0)
  motherID <- c(0, 0)
  gens <- seq(from=0, by=1, len=extraGenerations-1)
  father <- c(1,1)
  mother <- c(2,2)
  for (i in gens) {
    motherID <- c(motherID, c(0, mother, 0))
    mother <- 4*(i+1) + c(0,2)

    fatherID <- c(fatherID, c(0, father, 0))
    father <- 4*(i+1) + c(-1,1)
  }
  motherID <- c(motherID, mother)
  fatherID <- c(fatherID, father)

  ped@motherID <- as.integer(motherID)
  ped@fatherID <- as.integer(fatherID)

  invisible(ped)
}


MakePedigreeParentChild <- function(extraGenerations) {
  gender <- c(rep(c(1, 2), extraGenerations), 2) # assumed pedigree pattern

  nSubj <- length(gender)
  ped <- new("Pedigree", gender=as.integer(gender), typed=as.integer(c(1, nSubj)))

  ## determine mother and father id's
  fatherID <- c(0, 0)
  motherID <- c(0, 0)
  gens <- seq(from=0, by=1, len=extraGenerations-1)
  father <- 1
  mother <- 2
  for (i in gens) {
    motherID <- c(motherID, c(mother, 0))
    mother <- 2*(i+1) + 2

    fatherID <- c(fatherID, c(father, 0))
    father <- 2*(i+1) + 1
  }
  motherID <- c(motherID, mother)
  fatherID <- c(fatherID, father)

  ped@motherID <- as.integer(motherID)
  ped@fatherID <- as.integer(fatherID)

  invisible(ped)
}

MakePedigreeUncleNiece <- function(gen1, gen2) {
  ## determine mother and father id's
  fatherID <- c(0, 0)
  motherID <- c(0, 0)

  gens <- seq(from=0, by=1, len=gen1-1)
  father <- 1
  mother <- 2
  for (i in gens) {
    motherID <- c(motherID, c(0, mother))
    mother <- 2*(i+1) + 2

    fatherID <- c(fatherID, c(0, father))
    father <- 2*(i+1) + 1
  }
  motherID <- c(motherID, mother)
  fatherID <- c(fatherID, father)

  gender1 <- c(rep(c(1, 2), gen1), 1) # assumed pedigree pattern
  nSubj1 <- length(motherID)
  
  gens <- seq(from=0, by=1, len=gen2-1)
  father <- 1
  mother <- 2
  for (i in gens) {
    motherID <- c(motherID, c(mother, 0))
    mother <- 2*i + nSubj1 + 2

    fatherID <- c(fatherID, c(father, 0))
    father <- 2*i + nSubj1 + 1
  }
  motherID <- c(motherID, mother)
  fatherID <- c(fatherID, father)
  
  gender <- c(gender1, rep(c(1, 2), gen2-1), 2) # assumed pedigree pattern

  nSubj <- length(gender)
  ped <- new("Pedigree", gender=as.integer(gender), typed=as.integer(c(nSubj1, nSubj)))

  ped@motherID <- as.integer(motherID)
  ped@fatherID <- as.integer(fatherID)

  invisible(ped)
}

MakePedigreeGeneralisedHalfSib <- function(gen1, gen2) {
  ## determine mother and father id's
  fatherID <- c(0, 0, 0)
  motherID <- c(0, 0, 0)

  gens <- seq(from=0, by=1, len=gen1-1)
  father <- 2
  mother <- 1
  for (i in gens) {
    motherID <- c(motherID, c(0, mother))
    mother <- 2*(i+1) + 2

    fatherID <- c(fatherID, c(0, father))
    father <- 2*(i+1) + 3
  }
  motherID <- c(motherID, mother)
  fatherID <- c(fatherID, father)

  gender1 <- c(c(2,1,2), rep(c(2, 1), gen1-1), 1) # assumed pedigree pattern
  nSubj1 <- length(motherID)
  
  gens <- seq(from=0, by=1, len=gen2-1)
  father <- 2
  mother <- 3
  for (i in gens) {
    motherID <- c(motherID, c(mother, 0))
    mother <- 2*i + nSubj1 + 2

    fatherID <- c(fatherID, c(father, 0))
    father <- 2*i + nSubj1 + 1
  }
  motherID <- c(motherID, mother)
  fatherID <- c(fatherID, father)
  
  gender <- c(gender1, rep(c(1, 2), gen2-1), 2) # assumed pedigree pattern

  nSubj <- length(gender)
  ped <- new("Pedigree", gender=as.integer(gender), typed=as.integer(c(nSubj1, nSubj)))

  ped@motherID <- as.integer(motherID)
  ped@fatherID <- as.integer(fatherID)

  invisible(ped)
}

InitialPedigree <- function(relation, nmarker, prefix.ped) {
  if (!(relation %in% availableModels))
    stop("Family relation ", relation, " not available for analysis\n")

  ##  prefix <- relation
  inped <- paste(prefix.tmpfiles, prefix.ped, ".ped", sep="")
  ped <- MakePedigree(relation, pedfile=inped)
  typed <- ped@typed
  
  cmd <- paste("perl", paste(get("FEST.perlpath", envir=topenv(), inherits = FALSE), "makepedigree.pl", sep=""),
               inped, prefix.ped, nmarker, paste(typed, collapse=" "))
  system(cmd)
  file.remove(inped)

  length(ped@fam)
}
