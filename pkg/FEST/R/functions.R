####################################
####################################
## PUBLIC FUNCTIONS (NAMESPACE)
####################################
####################################

SetModels <- function(trueModels,
                      altModels=c("true", "unrelated")) {
  available <- trueModels %in% availableModels
  if (!all(available))
    stop("True Models not available for analysis: ", trueModels[!available], "\n")
  if (is.character(altModels)) {
    available <- altModels %in% availableAltModels
    if (!all(available))
      stop("Alternative models not available for analysis: ", altModels[!available], "\n")
    
    nTrue <- length(trueModels)
    altModelsList <- vector("list", nTrue)
    nAlt <- length(altModels)
    for (i in 1:nTrue) {
      reltype <- unlist(strsplit(trueModels[i], "-"))[1]
      if (reltype == "HS" ||
          reltype == "S" ||
          reltype == "PC") {
        extraGen <- as.integer(unlist(strsplit(trueModels[i], "-"))[-1])
##        if (length(extraGen) == 1
      }
      else
        extraGen <- NULL
      for (j in 1:nAlt) {
        if (altModels[j] == "true")
          altModelsList[[i]] <- c(altModelsList[[i]], trueModels[i])
        else if (altModels[j] == "lower") {
          if (all(extraGen > 1)) {
            if (length(extraGen) == 1)
              altModelsList[[i]] <- c(altModelsList[[i]], paste(reltype, "-", extraGen-1, sep=""))
            else
              altModelsList[[i]] <- c(altModelsList[[i]], paste(reltype, "-", extraGen[1]-1, "-", extraGen[2]-1, sep=""))
          }
        }
        else if (altModels[j] == "upper") {
          ## Must check if upper model is available
          if (length(extraGen) == 1)
            altMod <- paste(reltype, "-", extraGen+1, sep="")
          else
            altMod <- paste(reltype, "-", extraGen[1]+1, "-", extraGen[2]+1, sep="")
          available <- altMod %in% availableModels
          if (available) {
            altModelsList[[i]] <- c(altModelsList[[i]], altMod)
          }
        }
        else
          altModelsList[[i]] <- c(altModelsList[[i]], altModels[j])
      }
    }
  }
  else if (is.list(altModels)) {
    n <- length(altModels)
    for (i in 1:n) {
      available <- altModels[[i]] %in% availableModels
      if (!all(available))
        stop("Alternative models not available for analysis: ", altModels[[i]][!available], "\n")
    }
    altModelsList <- altModels
  }
  else {
    stop("altModels should either be a character vector or a list")
  }
  model <- new("Model", true=trueModels,
               alternative=altModelsList)
  model
}

SimulationStudyLD <- function(models,
                              ##                              chr=22,
                              nmarker=c(10,100,1000, Inf),
                              nsim=c(1000,1000,1000,400),
                              frequencyData=NULL,
                              freqThreshold=c(rep(0.1, 4), 0.0),
                              limitCentiMorgan=0,
                              saveMerlinFiles=FALSE,
                              verbose=TRUE) {
  
  op <- options(stringAsFactors=FALSE)
  on.exit(options(op))

  chr <- 22 ## Temporary restricted to chromosome 22
  
  if (length(nmarker) != length(nsim)) {
    stop("Input variables nmarker and nsim have unequal lengths")
  }
  if (min(nsim) <= 0) {
    stop("Input variable nsim must be >= 1")
  }
  if (length(nmarker) != length(nsim)) {
    stop("Input variables nmarker and nsim have unequal lengths")
  }

  if (is.vector(chr)) {
    nchr <- rep(length(chr), length(nmarker))
  }
  else if (is.list(chr)) {
    nchr <- lapply(chr, length)
    if (length(nchr) != length(nsim)) {
      stop("List chr and and vector nsim have unequal lengths")
    }
    p <- length(nmarker)

    chr.list <- vector("list", p)
    for (i in 1:p) {
      chr.list[[i]] <- chr
    }
    chr <- chr.list
  }
  else {
    stop("Input variable chr should be a vector or a list")
  }
  if (any(nchr > nmarker)) {
    stop("Number of markers less than number of chromosomes")
  }

  if (length(nmarker) != length(freqThreshold)) {
    stop("Input variables nmarker and freqThreshold have unequal lengths")
  }

  if (is.null(frequencyData)) {
    affy22.sub <- InitializeLD(chunkSize=1000)
    ## check hapdat$freq med affy22.sub
  }
  else {
    affy22.sub <- frequencyData
  }

###Thinning of frequencyData
  
  affy22 <- vector("list", 22)
  affy22[[22]] <- affy22.sub
    
  neach <- nmarker/nchr
  simres <- SimulationStudyWithFrequencyData(models, chr=chr, neach=neach,
                                             nsim=nsim, frequencyData=affy22, #frequencyData,
                                             freqThreshold=freqThreshold,
                                             limitCentiMorgan=limitCentiMorgan,
                                             saveMerlinFiles=saveMerlinFiles,
                                             verbose=verbose, LD=TRUE)
  
  lnliks <- simres$lnliks
  nmarker <- simres$nmarker
  simObject <- new("SimStudyObject", logLik=lnliks, nsim=as.integer(nsim),
                   nmarker=as.integer(nmarker), model=models, maf=numeric(),
                   freqThreshold=freqThreshold)
  simObject
}


SimulationStudy <- function(models,
                            chr=c(1:22),
                            nmarker=c(22*c(1,10,100,1000),500000),
                            nsim=c(1000,1000,1000,1000,400),
                            maf=numeric(),
                            frequencyData=NULL,
                            freqThreshold=c(rep(0.1, 4), 0.0),
                            ##                            limitCentiMorgan=0,
                            saveMerlinFiles=FALSE,
                            verbose=TRUE) {

  op <- options(stringAsFactors=FALSE)
  on.exit(options(op))

  if (is.null(frequencyData) && length(maf)==0) {
    stop("Either affyData (affymetrix data object) or maf (minor allele frequency) should be specified\n")
  }
  if (length(nmarker) != length(nsim)) {
    stop("Input variables nmarker and nsim have unequal lengths")
  }
  if (min(nsim) <= 0) {
    stop("Input variable nsim must be >= 1")
  }
  if (length(nmarker) != length(nsim)) {
    stop("Input variables nmarker and nsim have unequal lengths")
  }

  if (is.vector(chr)) {
    p <- length(nmarker)
    nchr <- rep(length(chr), p)
    chr.list <- vector("list", p)
    for (i in 1:p) {
      chr.list[[i]] <- chr
    }
    chr <- chr.list
  }
  else if (is.list(chr)) {
    nchr <- lapply(chr, length)
    if (length(nchr) != length(nsim)) {
      stop("List chr and and vector nsim have unequal lengths")
    }
  }
  else {
    stop("Input variable chr should be a vector or a list")
  }
  if (any(nchr > nmarker)) {
    stop("Number of markers less than number of chromosomes")
  }
##  if (LD) {
##    if (is.null(frequencyData)) {
##      stop("If LD then a frequency data object should be specified")
##    }
##  }
 
  if (!is.null(frequencyData)) {
    if (length(nmarker) != length(freqThreshold)) {
      stop("Input variables nmarker and freqThreshold have unequal lengths")
    }
    
    neach <- nmarker/nchr
    limitCentiMorgan <- 0 ## Not yet in use
    simres <- SimulationStudyWithFrequencyData(models, chr=chr, neach=neach,
                                               nsim=nsim, frequencyData=frequencyData,
                                               freqThreshold=freqThreshold,
                                               limitCentiMorgan=limitCentiMorgan,
                                               saveMerlinFiles=saveMerlinFiles,
                                               verbose=verbose, LD=FALSE)
    lnliks <- simres$lnliks
    nmarker <- simres$nmarker
  }
  else {
    if (maf <= 0.0 || maf >= 1.0) {
      stop("Minor allele frequency (maf) should be > 0 and < 1")
    }

   lnliks <- SimulationStudyWithVariableMafs(models, nchr=nchr, nmarker=nmarker, nsim=nsim,
                                              mafs=maf, saveMerlinFiles=saveMerlinFiles,
                                              verbose=verbose)
  }
  
  simObject <- new("SimStudyObject", logLik=lnliks, nsim=as.integer(nsim),
                   nmarker=as.integer(nmarker), model=models, maf=maf,
                   freqThreshold=freqThreshold)
  simObject
}


SetDataPars <- function(path,
                        chrdirs=NULL,
                        prefixInputFiles=NULL,
                        suffixPed=".pre",
                        format=c("qtdt","linkage"),
                        famList=NULL,
                        individualsTyped=NULL) {

  if (!is.null(chrdirs)) {
    chrdirs <- file.path(path, chrdirs)
  }
  else {
    chrdirs <- path
  }

  if (is.null(prefixInputFiles)) {
    cat("Message: Input parameter 'prefixInputFiles' not specified.\n")
    cat("         Assumes that there are input files with suffixes: .dat (linkage format)\n")
    cat("                                                           .dat, .freq, .map (qtdt format).\n")
  }
    
  if (!is.null(individualsTyped) & !is.matrix(individualsTyped)) {
    stop("Input parameter 'individualsTyped' should be a matrix.")
  }
  
  if (!is.null(famList)) {
    if (length(famList) == 1 && !is.null(individualsTyped)) {
      famList <- rep(famList, nrow(individualsTyped))
    }
    if (length(famList) != nrow(individualsTyped))
      stop("Length of famList differs from nrow(individualsTyped)")
  }
  
  list(chrdirs=chrdirs, prefixInputFiles=prefixInputFiles, suffixPed=suffixPed,
       format=format, fam=famList, individualsTyped=individualsTyped)
}



RealStudy <- function(altModels,
                      dataPars,
                      saveMerlinFiles=FALSE,
                      limitCentiMorgan=0,
                      freqThreshold=0) {
  
  op <- options(stringAsFactors=FALSE)
  on.exit(options(op))
  
  available <- altModels %in% availableModels
  if (!all(available)) {
    stop("Following alternative models not available for analysis: ",
         altModels[!available], "\n")
  }
  
  subjList <- dataPars$individualsTyped
  famList <- dataPars$fam
  nCases <- nrow(subjList)
  nAlt <- length(altModels)
  nDir <- length(dataPars$chrdirs)


  ## TEMPORARY: # markers to be thinned
  nmarkers <- NULL
  
  
  lnliks <- matrix(0, nrow=nCases, ncol=nAlt)
  for (i in 1:nCases) {
    for (j in 1:nAlt) {
      for (k in 1:nDir) {
        ll <- RunRealAnalyses(altModel=altModels[j], datadir=dataPars$chrdirs[k],
                              famList=famList[i], subjList=subjList[i,],
                              prefixInput=dataPars$prefixInputFiles,
                              suffixPed=dataPars$suffixPed, format=dataPars$format,
                              saveMerlinFiles=saveMerlinFiles,
                              limitCentiMorgan=limitCentiMorgan,
                              freqThreshold=freqThreshold,
                              nmarker=nmarkers[i])
        lnliks[i,j] <- lnliks[i,j] + ll$lnlik
      }
    }
  }

  colnames(lnliks) <- altModels
  rownames(lnliks) <- apply(subjList, 1, paste, collapse="-")
  
  posterior <- t(apply(lnliks, 1, Posterior))
  res <- list(logLiks=lnliks, posterior=posterior, nmarker=ll$nmarker)
  res
}


ComputeSummaryStatistics <- function(simObj, #colNames, rowNames,
                                     altHyp="unrelated") {
  nDataSets <- length(simObj@nmarker)
  nTrue <- length(simObj@logLik[[1]])

  altModels <- sort(AlternativeModels(simObj@model@alternative))
  nAltAll <- length(altModels)
  
  likratMean <- array(NA, dim=c(nDataSets, nTrue, nAltAll))
  likratMean.sd <- array(NA, dim=c(nDataSets, nTrue, nAltAll))
  likratMedian <- array(NA, dim=c(nDataSets, nTrue, nAltAll))
  classRate <- array(NA, dim=c(nDataSets, nTrue, nAltAll))
  classRate.sd <- array(NA, dim=c(nDataSets, nTrue, nAltAll))
  posteriorMean <- array(NA, dim=c(nDataSets, nTrue, nAltAll))
  posteriorMean.sd <- array(NA, dim=c(nDataSets, nTrue, nAltAll))

  dimNames <- list(simObj@nmarker,
                   simObj@model@true,
                   altModels)
  dimnames(likratMean) <- dimNames
  dimnames(likratMean.sd) <- dimNames
  dimnames(likratMedian) <- dimNames
  dimnames(classRate) <- dimNames
  dimnames(classRate.sd) <- dimNames
  dimnames(posteriorMean) <- dimNames
  dimnames(posteriorMean.sd) <- dimNames

  for (i in 1:nDataSets) {
    for (j in 1:nTrue) {
      ll <- simObj@logLik[[i]][[j]]
      po <- simObj@posterior[[i]][[j]]
      nSim <- ncol(ll)

      alts <- rownames(po)

      ind.alt <- which(alts %in% altHyp)
      if (length(ind.alt) > 0) {
        lrat <- ll[-ind.alt,,drop=FALSE] - ll[ind.alt,]
        alts.sub <- alts[-ind.alt]
      
        lmean <- apply(lrat, 1, mean)
        lsd <- apply(lrat, 1, sd)/sqrt(nSim)   
        likratMean[i,j,alts.sub] <- lmean
        likratMean.sd[i,j,alts.sub] <- lsd
      

        lmedian <- apply(lrat, 1, median)
        likratMedian[i,j,alts.sub] <- lmedian

        crate <- apply(lrat>0, 1, mean)
        crate.sd <- apply(lrat>0, 1, sd)/sqrt(nSim)
        classRate[i,j,alts.sub] <- crate
        classRate.sd[i,j,alts.sub] <- crate.sd
      }

      pmean <- apply(po, 1, mean)
      psd <- apply(po, 1, sd)/sqrt(nSim)
      posteriorMean[i,j,alts] <- pmean
      posteriorMean.sd[i,j,alts] <- psd
    }
  }


  stat <- list(likrat=likratMean, likrat.sd=likratMean.sd, median=likratMedian,
               classrat=classRate, classrat.sd=classRate,
               posterior=posteriorMean, posterior.sd=posteriorMean.sd)
  stat
}


####################################
####################################
## INTERNAL FUNCTIONS
####################################
####################################

setMethod("initialize", "SimStudyObject", function(.Object, ...) {
  .args <- list(...)

  .Object <- callNextMethod();

  if (hasArg(logLik)) {
    .Object@posterior <- ComputePosterior(.args$logLik)
  }
  
  .Object
})

## Simulates family-relations and computes likelihood statistics using merlin
SimulationStudyWithFixedMaf <- function(model,
                                        nchr,
                                        nmarker=c(22, 22*10, 22*100, 22*1000, 500000),
                                        nsim=c(1000,1000,1000,1000,400),
                                        maf=0.5,
                                        saveMerlinFiles=FALSE,
                                        verbose=TRUE) {
  startSeed <- 107
  nMarkerNumbers <- length(nmarker)
  simRes <- vector("list", nMarkerNumbers)
  ##  cat("Minor allele frequency ", maf, "\n")
  for (j in 1:nMarkerNumbers) {
    if (verbose) {
      cat("Number of markers ", format(as.character(nmarker[j]), width=6, justify="left"), " - ", date(), "\n")
    }
    simRes[[j]] <- RunAnalysesAll(model=model, nmarker=nmarker[j], nchr=nchr[j], nsim=nsim[j], maf=maf,
                                  seed=startSeed, saveMerlinFiles=saveMerlinFiles)
  }
  if (verbose) {
    cat("Simulation study finished  - ", date(), "\n")
  }

  simRes
}

## Simulates family-relations and computes likelihood statistics using merlin
SimulationStudyWithVariableMafs <- function(model,
                                            nchr,
                                            nmarker=22*c(1,10,100),
                                            nsim=c(1000,1000,1000),
                                            mafs=c(0.1,0.25,0.5),
                                            saveMerlinFiles=FALSE,
                                            verbose=TRUE) {
  nMaf <- length(mafs)
  simRes <- NULL
  for (i in 1:nMaf) {
    maf <- mafs[i]
    if (verbose) {
      cat("Minor allele frequency ", maf, "\n")
    }
    simRes <- c(simRes,
                SimulationStudyWithFixedMaf(model=model, nchr=nchr, nmarker=nmarker, nsim=nsim, maf=mafs[i],
                                            saveMerlinFiles=saveMerlinFiles, verbose=verbose))
  }

  simRes
}


## Simulation study for the affymetrix data
SimulationStudyWithFrequencyData <- function(model, chr,
                                             neach=c(1, 10, 100, 1000, Inf, Inf),
                                             frequencyData=NULL,
                                             freqThreshold=c(0.1, 0.1, 0.1, 0.1, 0.1, 0.0),
                                             limitCentiMorgan=0,
                                             nsim=c(1000,1000,1000,1000,400,400),
                                             saveMerlinFiles=FALSE,
                                             verbose=TRUE,
                                             LD=FALSE) {
  
  nMarkerNumbers <- length(neach)
  lnliks <- vector("list", nMarkerNumbers)
  nmarker <- rep(NA, nMarkerNumbers)
  
  startSeed <- 107
  if (verbose) {
    cat("Run analysis for subsets of frequency data.\n")
  }
  for (i in 1:nMarkerNumbers) {
    if (verbose) {
      cat("Number of SNPs for each chromosome: ",
          format(as.character(neach[i]), width=5, justify="left"), date(), "\n")
    }
    
    subFreq <- SelectSNPs(frequencyData, neach=neach[i], chr=chr[[i]], threshold=freqThreshold[i],
                          limitCentiMorgan=limitCentiMorgan)

    if (LD) {
      SelectHapdat(subFreq)
    }
    WriteInputFiles(subFreq, chr=chr[[i]])
    nmarker <- sum(unlist(lapply(subFreq, nrow)))
    nchr <- length(chr[[i]])
    lnliks[[i]] <- RunAnalysesAll(model=model, nmarker=nmarker, nchr=nchr,
                                  nsim=nsim[i], seed=startSeed,
                                  useInputFiles=TRUE, saveMerlinFiles=saveMerlinFiles,
                                  LD=LD)
  }
  
  if (verbose) {
    cat("Simulation study finished                 ", date(), "\n")
  }

  
  simres <- list(lnliks=lnliks, nmarker=nmarker)
  simres
}

##SaveToFile <- function(saveFile, simRes, model, nDataSets, nDone, nmarker, nsim,
##                       maf=numeric(), freqThreshold=numeric()) {
##  simObject.save <- new("SimStudyObject", logLik=simRes, nsim=as.integer(nsim),
##                        nmarker=as.integer(nmarker), model=model, maf=maf,
##                        freqThreshold=freqThreshold)
##  save(nDone, simObject.save, file=saveFile)
##}

CheckMerlinRun <- function(merlinObj) {
  ind.warning <- grep("warning", merlinObj, ignore.case = TRUE)
  if (length(ind.warning)>0) {
    for (i in 0:3) {
      ind.warning <- union(ind.warning, ind.warning+i)
    }
    ind.warning <- sort(ind.warning)
    warning("MERLIN WARNING: ", paste(merlinObj[ind.warning], collapse="\n"))
  }
  ind.error <- grep("error", merlinObj, ignore.case = TRUE)
  ind.error.info <- grep("--error", merlinObj, ignore.case = TRUE)
  ind.error <- setdiff(ind.error, ind.error.info)
  if (length(ind.error)>0) {
    for (i in 0:3) {
      ind.error <- union(ind.error, ind.error+i)
    }
    ind.error <- sort(ind.error)
    stop("MERLIN ERROR: ", paste(merlinObj[ind.error], collapse="\n"))
  }
  ind.badinheritance <- grep("bad inheritance",  merlinObj, ignore.case = TRUE)
  bad.inheritance <- length(ind.badinheritance) > 0
  invisible(bad.inheritance)
}

## Remove temporary simulated merlin files (merlin*.ped/freq/dat/map)
CleanMerlinFiles <- function(prefixes, suffixes=c("ped", "freq", "dat", "map")) {
  for (prefix in prefixes) {
    for (suffix in suffixes) {
      file.remove(list.files(pattern=paste("^", prefix, ".*\\.", suffix, "$", sep="")))
    }
  }
##  file.remove(list.files(pattern=paste("^", prefix.ped, ".*\\.freq$", sep="")))
##  file.remove(list.files(pattern=paste("^", prefix.ped, ".*\\.dat$", sep="")))
##  file.remove(list.files(pattern=paste("^", prefix.ped, ".*\\.map$", sep="")))
##  system(paste("rm -f ", prefix.ped, "*.ped", sep=""))
##  system(paste("rm -f ", prefix.ped, "*.dat", sep=""))
##  system(paste("rm -f ", prefix.ped, "*.freq", sep=""))
##  system(paste("rm -f ", prefix.ped, "*.map", sep=""))
}

SimulatePedigreesWithLD <- function(trueModel, nsim, seed, nSimInEachPedigree=5,
                                    nmarker=NULL, altModels=NULL, nrows=NULL) {
  prefix.ped <- paste(prefix.tmpfiles, "simped_", sep="")
  for (i in 1:nsim) {
    simpedfile <- paste(prefix.ped, i, ".ped", sep="")
    SimulatePedigreeWithLD(trueModel=trueModel, simpedfile=simpedfile, nmarker=nmarker)
  }
    
  MergeSimulatedPedigrees(prefixSimFiles=prefix.ped, trueModel=trueModel,
                          nSimInEachPedigree=nSimInEachPedigree, nmarker=nmarker,
                          altModels=altModels, nrows=nrows)
}

SimulatePedigreeWithLD <- function(trueModel, simpedfile, nmarker) {
  
  map <- read.table(paste(prefixMerlin.input, ".map", sep=""), header=TRUE)
  ped <- read.table(paste(prefix.tmpfiles, trueModel, ".ped", sep=""), header=FALSE)[,1:6]
  typed <- ped[,6]
  ped <- ped[,1:5]
  ## NB: freq is already used to make the hapdat object
  
  nFounders <- sum(ped[,3] == 0)

  hap <- SimulateHaplotypesLD(nmarker, 2*nFounders)

  markers <- ConditionalSimulation(hap, map$POS, ped)
  markers[!typed,] <- 0 # Set markers of untyped individuals = 0 (required by MergeSimulatedPedigrees)
  affected <- rep(1, nrow(ped))
  ped.sim <- cbind(ped, markers, affected)

  write.table(ped.sim, file=simpedfile, col.names=FALSE, row.names=FALSE)
}

InitializeLD <- function(chunkSize=1000, ldFile="ld_chr22_CEU.txt") {
  cat("Select Hapmap LD subset: Only SNPs in Affymetrix 500K are retained\n")
  cmd <- paste("perl hapmapLD_affy.pl", ldFile, "../Data/Affy500K_Allele_Frequency_Files/chr22.freq.txt")
  system(cmd)
  cat("Divide LD file into chunks of size ", chunkSize , "\n")
  ldFile.affy <- paste(ldFile, ".affy", sep="")
  cmd <- paste("perl hapmapLD_divide.pl", ldFile.affy, chunkSize)
  res <- system(cmd, intern=TRUE)
  nchunk <- as.integer(unlist(strsplit(res, " "))[1])

  ##  datadir <- "."
  ##  files <- list.files(datadir, pattern=glob2rx("ld*.chunk*"))
  files <- paste(ldFile.affy, ".chunk", 1:nchunk, sep="")

  ##  nchunk <- length(files)
  cat("Number of chunks: ", nchunk, "\n")

  ##  download.file("http://folk.uio.no/thoree/FEST/affy.RData", "affy.RData")
  affy <- NULL # define affy such that R CMD check do not report warnings.
  load("affy.RData") ## get affy
  affy22 <- affy[[22]] #SNP         cM     A     C

  hapdat <- vector("list", nchunk)
  ind.affy <- NULL
  for (i in 1:nchunk) {
    ldObj <- read.table(files[i], header=TRUE)

    snps <- unique(ldObj$snp1)

    ind <- which(affy22[,"SNP"] %in% snps)
    affy22.sub <- affy22[ind,]
    ind.affy <- c(ind.affy, ind)

    snps.affy <- affy22[ind,"SNP"]

    indLD <- (ldObj[,"snp1"] %in% snps.affy) & (ldObj[,"snp2"] %in% snps.affy)
    ldObj.sub <- ldObj[indLD,]
    
    nsnps <- length(snps.affy)
    r2mat <- matrix(0,nrow=nsnps, ncol=nsnps)
    rownames(r2mat) <- snps.affy
    colnames(r2mat) <- snps.affy

    nLD <- nrow(ldObj.sub)
    snp1 <- ldObj.sub[,"snp1"]
    snp2 <- ldObj.sub[,"snp2"]
    r2 <- ldObj.sub[,"R2"]
    for (j in 1:nLD) {
      r2mat[snp1[j],snp2[j]] <- r2[j]
      r2mat[snp2[j],snp1[j]] <- r2[j]
    }
    diag(r2mat) <- 1

    cat("Make Haplodata object for chunk ", i, "\n")
    freqs <- affy22.sub[,"A"]
    names(freqs) <- snps.affy
    hapdat[[i]] <- MakeHaplodataObject(freqs, r2mat)
  }
  affy22.hapmap <- affy22[ind.affy, ]
  
  assign("hapdat", hapdat, envir=topenv())

  file.hapdat <- "hapdat.Rdata"
  cat("Write HapMap LD correlation data to ", file.hapdat, "\n")
  save(hapdat, file=file.hapdat)
  affy22.hapmap
}

MakeHaplodataObject <- function(freqs, cor) {
  P <- freqs
  Q <- qnorm(P)
  nloci <- length(freqs)
  null.mat <- matrix(0, nrow = nloci, ncol = nloci)
  vmat <- .C("covariance", as.integer(nloci), as.double(cor), 
             as.double(P), as.double(Q), rlt = as.double(null.mat), 
             PACKAGE = "hapsim")$rlt
  V <- matrix(vmat, nrow = nloci, ncol = nloci)
  if (!checkpd(V)) 
    V <- makepd(V)
  return(list(freqs = P, cor = cor, cov = V, div = NULL))
}


SelectHapdat <- function(subFreq) {
  hapdat <- NULL # Temporary, to avoid warnings in 'R CMD check'
  hapdat.sub <- hapdat
  for (i in 1:length(hapdat.sub)) {
    nm <- names(hapdat.sub[[i]]$freqs)
    
    ind <- which(nm %in% subFreq$SNP)

    hapdat.sub[[i]]$freqs <- hapdat.sub[[i]]$freqs[ind]
    hapdat.sub[[i]]$cor <- hapdat.sub[[i]]$cor[ind,ind]
    hapdat.sub[[i]]$cov <- hapdat.sub[[i]]$cov[ind,ind]
  }

  assign("hapdat.sub", hapdat.sub, envir=topenv())
}

#
SimulateHaplotypesLD <- function(nmarker, nhap) {
  u <- NULL
  hapdat.sub <- NULL # Temporary, to avoid warnings in 'R CMD check'
  for (i in 1:length(hapdat.sub)) {
    hsim <- haplosim.fix(nhap, hapdat.sub[[i]], force.polym=FALSE, summary=FALSE)
    u <- cbind(u, hsim$data+1)
  }
  
  u
}

## Simulated pedigrees conditional on founder haplotypes
## hap: founder haplotypes
## mapPos: genetic map position in centiMorgans
## ped: pedigree
ConditionalSimulation <- function(hap, mapPos, ped) {
  combine <- function(child, father, mother) {
    if (!finished[father])
      combine(father, ped[father, 3], ped[father, 4])
    if (!finished[mother])
      combine(mother, ped[mother, 3], ped[mother, 4])

    ## Simulate recombination events for mother and fathers genes
    ran <- range(mapPos)
    lambda <- diff(ran)/100
    pIndex <- -1
    for (person in c(mother, father)) {
      n.recomb <- rpois(1, lambda=lambda)
      if (n.recomb > 0) {
        recomb.pos <- runif(n.recomb, min=ran[1], max=ran[2])
        index <- findInterval(recomb.pos, mapPos) + 1
        
        u <- rep(0, nmarker)
        u[index] <- 1
        u <- (cumsum(u) %% 2) + 1
      }
      else {
        u <- rep(1, nmarker)
      }

      for (j in 1:nmarker) { # could this be vectorized?
        markers[child,2*j+pIndex] <- markers[person, 2*j-2 + u[j]]
      }
      pIndex <- pIndex + 1
    }
    markers
  }

  nmarker <- ncol(hap)
  nmarker2 <- 2*nmarker
  nsubj <- nrow(ped)
  
  founders <- ped[,3] == 0
  markers <- matrix(NA, nrow=nsubj, ncol=nmarker2)
  markers[founders,] <- hap

  finished <- rep(FALSE, nsubj)
  finished[founders] <- TRUE
  for (i in 1:nrow(ped)) {
    if (!finished[i])
      markers <- combine(i, ped[i,3], ped[i, 4])
  }

  markers
}

## Simulate nsim pedigrees (and corresponding freq,dat,map files)
## based on pedigree structure given in the 'inped' file
SimulatePedigrees <- function(trueModel, nsim, seed, nSimInEachPedigree=5,
                              nmarker=NULL, altModels=NULL, nrows=NULL) { 
  ## Remove previous simulated merlin data files
  prefix.ped <- "merlin-"

  ## remove temporary simulated files
  on.exit(CleanMerlinFiles(prefix.ped))

  ## Simulate nsim pedigree files using MERLIN
  merlinSimulateCommand <- paste("merlin ",
                                 paste(paste(c("-d ", "-m ", "-f "), prefixMerlin.input,
                                             c(".dat", ".map", ".freq"), sep=""), collapse=" "),
                                 " -p ", paste(prefix.tmpfiles, trueModel, ".ped", sep=""),
                                 " --simulate -r ", seed," --reruns:", nsim," --save",sep="")
  merlinObj.tmp <- system(merlinSimulateCommand, intern=TRUE)
  CheckMerlinRun(merlinObj.tmp)

  MergeSimulatedPedigrees(prefix.ped, trueModel, nSimInEachPedigree, nmarker, altModels, nrows)
}

MergeSimulatedPedigrees <- function(prefixSimFiles, trueModel, nSimInEachPedigree, nmarker, altModels, nrows) {
  ## Merge simulated pedigree files by calling perl script 'merge.pl'.
  ## Makes pedigree files for 
  if (is.null(nrows)) {
    nrows.str <- system(paste("sed -n \"$=\" ", paste(prefix.tmpfiles, trueModel, ".ped", sep="")), intern=TRUE)
    nrows <- as.integer(nrows.str)
  }
  nAlt <- length(altModels)
  filenames.altModels <- paste(prefix.tmpfiles, altModels, sep="")
  cmd <- paste("perl", paste(get("FEST.perlpath", envir=topenv(), inherits = FALSE), "merge.pl", sep=""),
               prefixSimFiles, nrows, nmarker, nSimInEachPedigree, #prefixMerlin.merged,
               paste(filenames.altModels, collapse =" "),
               paste(paste(prefixMerlin.mergedAlt, 1:nAlt, sep=""), collapse=" "), sep=" ")
  system(cmd)
}

RunAnalysesAll <- function(model, nmarker, nchr, nsim, seed,
                           maf=NULL, useInputFiles=FALSE, saveMerlinFiles=FALSE, LD=FALSE) {
  InitialMerlinInputFiles <- function(nchr, neach, n.alleles=2, freq, chrLength=190) {
    cmd <- paste("perl", paste(get("FEST.perlpath", envir=topenv(), inherits = FALSE), "makeMerlinInputs.pl", sep=""),
                 prefixMerlin.input, chrLength, neach, nchr, n.alleles, paste(freq, collapse=" "))
    system(cmd)
  }

  if (!saveMerlinFiles) {
    on.exit(CleanMerlinFiles(prefixMerlin.input, c("dat", "freq", "map")))
  }
  
  if (!useInputFiles && is.null(maf)) {
    stop("Either affyData (affymetrix data object) or maf (minor allele frequency) should be specified\n")
  }
  
  if (!useInputFiles) {
    neach <- as.integer(nmarker/nchr)
    nmarker <- nchr * neach

    ## makes initial input files: for.dat, for.freq, for.map (preceeded with prefix.tmpfiles)
    InitialMerlinInputFiles(nchr, neach, 2, c(maf, 1-maf))
  }

  nTrue <- length(model@true)

  lnliks <- vector("list", nTrue)
  for (i in 1:nTrue) {
    lnliks[[i]] <- RunAnalyses(trueModel=model@true[i], altModels=model@alternative[[i]],
                               nsim=nsim, seed=seed, nchr=nchr, nmarker=nmarker, saveMerlinFiles=saveMerlinFiles, LD=LD)
  }

  lnliks
}

## 1) Simulate pedigrees 
## 2) Compute likelihoods for simulated pedigrees for 
##   1) a given pedigree structure (given by the 'inped' file)
##   2) subjects are not related
RunAnalyses <- function(trueModel, altModels, nsim=4, seed=107, nchr=2,
                        nmarker=4, nSimInEachPedigree=5, saveMerlinFiles=FALSE, LD=FALSE) {
  ## Simulate pedigrees and corresponding data files by call to Merlin
  ## Compute likelihoods by call to Merlin
  lnlik.allSplits <- NULL

  trueAndAltmodels <- unique(c(trueModel, altModels))

  prefix.trueAndAltmodels <- paste(prefix.tmpfiles, trueAndAltmodels, sep="")
  
  ## remove temporary files
  on.exit(CleanMerlinFiles(prefix.trueAndAltmodels, "ped"))
  if (!saveMerlinFiles) {
    on.exit(CleanMerlinFiles(prefixMerlin.mergedAlt, "ped"), add=TRUE)
  }
  
  ## Make initial pedigree files for true and alternative family relations
  nrows <- NULL
  for (i in 1:length(trueAndAltmodels)) {
    mod <- trueAndAltmodels[i]
    prefix <- prefix.trueAndAltmodels[i]
    nrows <- c(nrows,
               InitialPedigree(mod, nmarker=nmarker, prefix))
  }
  names(nrows) <- trueAndAltmodels
  
  if (nsim < nSimInEachPedigree) {
    nSplits <- 1
    nSimInEachPedigree <- nsim
  }
  else
    nSplits <- as.integer(nsim/nSimInEachPedigree)

  nAlt <- length(altModels)
  for (i in 1:nSplits) {
    ## Simulate pedigree files and merge these in (file names preceeded with prefix.tmpfiles)
    ## msimAlt1.ped, msimAlt2.ped, ... (nAlt files)
    if (LD) {
      SimulatePedigreesWithLD(trueModel=trueModel, nsim=nSimInEachPedigree, seed=seed,
                              nSimInEachPedigree=nSimInEachPedigree, nmarker=nmarker, altModels=altModels,
                              nrows=nrows[trueModel])
    }
    else {
      SimulatePedigrees(trueModel=trueModel, nsim=nSimInEachPedigree, seed=seed,
                        nSimInEachPedigree=nSimInEachPedigree, nmarker=nmarker, altModels=altModels,
                        nrows=nrows[trueModel])
    }
    
    lnlik <- NULL
    
    for (j in 1:nAlt) {
      merlinCommand <- paste("merlin ",
                             paste(paste(c("-d ", "-m ", "-f "), prefixMerlin.input,
                                         c(".dat", ".map", ".freq"), sep=""), collapse=" "),
                             " -p ", prefixMerlin.mergedAlt, j, ".ped --quiet --lik --per --mega 1024", sep="")
      res.all <- system(merlinCommand, intern=TRUE, wait=TRUE)
      bad.inheritance <- CheckMerlinRun(res.all)

      ## Extract likelihoods from Merlin output and compute likelihood ratios
      if (bad.inheritance) {
        lnlik <- rbind(lnlik, rep(-Inf, nSimInEachPedigree)) ## 
      }
      else {
        lnlik <- rbind(lnlik, ExtractLikelihood(res.all, 1, nSimInEachPedigree, nchr)[1,]) ## assumes 1 family
      }
     }
    seed <- seed + 1
    lnlik.allSplits <- cbind(lnlik.allSplits, lnlik)
  }
  rownames(lnlik.allSplits) <- altModels
  
  ## Return likelihoods
  lnlik.allSplits
}


ExtractLikelihood <- function(res.all, n.fam, nsim, nchr) {
  ind <- grep("lnLikelihood =", res.all)
  lnlik <- as.numeric(substring(res.all[ind], first=18))
  lnlik.array <- array(lnlik, dim=c(n.fam, nsim, nchr))

  lnlik.sim <- apply(lnlik.array, c(1,2), sum)
  lnlik.sim
}



AlternativeModels <- function(altModels) {
  allMod <- NULL
  n <- length(altModels)
  for (i in 1:n) {
    allMod <- c(allMod, altModels[[i]])
  }
  unique(allMod)
}

MakeDimNames <- function(nmarker, altModels, mafs=NULL) {
  colNames <-  altModels
  rowNames <- NULL
  if (length(mafs)>1) {
    for (maf in mafs)
      rowNames <- c(rowNames,
                    paste(paste("(Maf ", format(round(maf,2),nsmall=2), ")", sep=""), nmarker))
  }
  else {
    rowNames <- nmarker
  }
  list(row=rowNames, col=colNames)
}

Posterior <- function(logLiks) {
  max.l <- max(logLiks)  
  exp(logLiks - max.l)/sum(exp(logLiks - max.l))
}


ComputePosterior <- function(logLiks) {
  nDataSets <- length(logLiks)
  nTrue <- length(logLiks[[1]])

  posterior <- vector("list", nDataSets)
  for (i in 1:nDataSets) {
    posterior[[i]] <- vector("list", nTrue)
    for (j in 1:nTrue) {
      x <- logLiks[[i]][[j]]
      if (!is.null(x)) {
        nAlt <- nrow(x)
        nSim <- ncol(x)

        posterior[[i]][[j]] <- matrix(NA, nrow=nAlt, ncol=nSim)
        for (l in 1:nSim) {
          max.l <- max(x[,l])
          posterior[[i]][[j]][,l] <- exp(x[,l]-max.l)/sum(exp(x[,l]-max.l))
        }
        rownames(posterior[[i]][[j]]) <- rownames(x)
      }
      else {
        posterior[[i]][[j]] <- NULL
      }
    }
  }
  posterior
}


###MakeLatexTable <- function(posterior, latexSuffix="") {
###  latex.table(round(posterior,3), paste("posterior",latexSuffix, sep=""),
###              rowlabel="", caption = "Posterior probabilities are shown. The first column shows the markers used.")
###}


##########################################################################
##########################################################################
##
## Real data analysis: Comparing possible family relationships
##
##########################################################################
##########################################################################

## pedfile0: gives hypothetized family relation
## datadirs: dat file for different chromosomes (linkage format)
## Must run separately for each chromosome

RunRealAnalyses <- function(altModel, datadir=file.path("../Data/KajaData/Merlin/Kr22/"),
                            famList, subjList, typed=FALSE, suffixPed=".pre",
                            prefixInput=NULL,
                            format=c("qtdt","linkage"), saveMerlinFiles=FALSE,
                            limitCentiMorgan=NULL, freqThreshold=NULL, nmarker.thin=NULL) {

  ComputeLikelihoods <- function(pedfile, datfile, n.fam, format,
                                 freqfile=NULL, mapfile=NULL) {
    if (format == "linkage") {
      merlinCommand.rel <- paste("merlin  -d ", datfile, " -p ", pedfile, " --quiet --lik --per", sep="")
    }
    else {
      merlinCommand.rel <- paste("merlin  -d ", datfile, " -p ", pedfile,
                                 " -f ", freqfile, " -m ", mapfile,
                                 " --quiet --lik --per", sep="")
    }
    res.all.rel <- system(merlinCommand.rel, intern=TRUE, wait=TRUE)
    bad.inheritance <- CheckMerlinRun(res.all.rel)
    
    ## Extract likelihoods from Merlin output and compute likelihood ratios
    if (bad.inheritance) {
      lnlik.rel <- -Inf
    }
    else {
      cmd <- paste("perl", paste(get("FEST.perlpath", envir=topenv(), inherits = FALSE), "numberOfChr.pl", sep=""),
                    mapfile)
      nchr <- as.integer(system(cmd, intern=TRUE))

      lnlik.rel <- ExtractLikelihood(res.all.rel, n.fam, 1, nchr) ## assumes 1 family
    }
    lnlik.rel
  }

  format <- match.arg(format)

  ##  cat(pedfile0, "", datadir, "\n")
  nFam <- length(famList)

  if (format == "linkage") {
    if (!is.null(prefixInput)) {
      datfile <- file.path(datadir, paste(prefixInput, ".dat", sep=""))
      pedfile <- file.path(datadir, paste(prefixInput, suffixPed, sep=""))
    }
    else {
      datfile <- list.files(datadir, pattern=glob2rx("*.dat"),
                            full.names=TRUE)
      pedfile <- list.files(datadir, pattern=glob2rx(paste("*", suffixPed, sep="")),
                            full.names=TRUE)
    }
  }
  else {
    if (!is.null(prefixInput)) {
      datfile <- file.path(datadir, paste(prefixInput, ".dat", sep=""))
      pedfile <- file.path(datadir, paste(prefixInput, suffixPed, sep=""))
      freqfile <- file.path(datadir, paste(prefixInput, ".freq", sep=""))
      mapfile <- file.path(datadir, paste(prefixInput, ".map", sep=""))
    }
    else {
      datfile <- list.files(datadir, pattern=glob2rx("*.dat"),
                            full.names=TRUE)
      pedfile <- list.files(datadir, pattern=glob2rx(paste("*", suffixPed, sep="")),
                            full.names=TRUE)    
      freqfile <- list.files(datadir, pattern=glob2rx("*.freq"),
                             full.names=TRUE)
      mapfile <- list.files(datadir, pattern=glob2rx("*.map"),
                            full.names=TRUE)
    }
    if  (!file.exists(freqfile)) {
      stop(freqfile, " does not exist")
    }
    if  (!file.exists(mapfile)) {
      stop(mapfile, " does not exist")
    }
 }
  if  (!file.exists(datfile)) {
    stop(datfile, " does not exist")
  }
  if  (!file.exists(pedfile)) {
    stop(pedfile, " does not exist")
  }
  
  pedAlt <- WritePed(MakePedigree(altModel))
  

  ##  ped <- read.table(pedfile, header=FALSE)
  n.col <- length(scan(pedfile, nlines=1, quiet = TRUE))
  ped <-matrix(scan(pedfile, quiet = TRUE), ncol=n.col, byrow=TRUE)

  nColPed <- ncol(ped)
  if (ncol(ped) %% 2 == 1) {
    ## Trait missing
    nNotMarkerFirst <- 5
    nTrait <- 0
  }
  else {
    nNotMarkerFirst <- 5
    nTrait <- 1
  }
  nMarker2AndTrait <- ncol(ped) - nNotMarkerFirst
  nSubjAlt <- nrow(pedAlt)

  ## MERGE pedfiles
  ## For each family: The subjects typed in pedAlt should get data from
  ##                  the corresponding two subjects of ped

  pedSpecified <- NULL
  for (i in 1:nFam) {
    f <- famList[i]
    typed <- which(pedAlt[,6] == 1)

    if (nSubjAlt>2) {
      pedAltMinusTyped <- cbind(pedAlt[-typed,1:nNotMarkerFirst, drop=FALSE],
                                matrix(0, nrow=nSubjAlt-length(typed), ncol=nMarker2AndTrait))
      colnames(pedAltMinusTyped) <- names(ped)
    }
    else
      pedAltMinusTyped <- NULL
    
    pedAltTyped <- cbind(pedAlt[typed,1:nNotMarkerFirst, drop=FALSE])

    ind.typed.ped <- which(is.element(ped[,2], subjList) & ped[,1] == f)
    pedAltTypedWithMarkers <- as.matrix(cbind(pedAltTyped, ped[ind.typed.ped,(nNotMarkerFirst+1):nColPed]))
    pedf <- rbind(pedAltMinusTyped, pedAltTypedWithMarkers)
    pedf[,1] <- f
    ##    pedf[,2] <- 1:nSubjAlt
    pedSpecified <- rbind(pedSpecified, pedf)
  }

  
  pedfileSpecified <- paste(prefix.tmpfiles, altModel, ".ped", sep="")
  write.table(pedSpecified, col.names=FALSE, row.names=FALSE, file=pedfileSpecified)

  if (!saveMerlinFiles) {
    on.exit(file.remove(pedfileSpecified))
  }
  else {
    cat("Save merlin pedfile to ", pedfileSpecified, "\n")
  }

  nMarker <- (nMarker2AndTrait-nTrait)/2
  
  if (limitCentiMorgan > 0 || freqThreshold > 0 || !is.null(nmarker.thin)) {
    suffix <- "_thinned"
    nMarker <- ThinMerlinInputFiles(mapfile=mapfile, datfile=datfile, freqfile=freqfile, pedfile=pedfileSpecified,
                                    nNotMarker=nNotMarkerFirst, limitCentiMorgan=limitCentiMorgan, freqThreshold=freqThreshold, suffix=suffix)
    
    pedfileSpecified <- paste(pedfileSpecified,  suffix, sep="")
    datfile <- paste(datfile,  suffix, sep="")
    freqfile <- paste(freqfile,  suffix, sep="")
    mapfile <- paste(mapfile,  suffix, sep="")

##    nMarker <- as.integer(system(paste("perl -pe '}{$_=$.'", datfile), intern=TRUE))
    nMarker <- length(count.fields(datfile))
    ##    cat("Thinned number of markers from ", nr1, " to ", nr2, "\n")
  }
   
  lnlik <- ComputeLikelihoods(pedfile=pedfileSpecified, datfile=datfile,
                              mapfile=mapfile, freqfile=freqfile,
                              n.fam=nFam, format=format)

  
  list(lnlik=lnlik, nmarker=nMarker)
}

ThinMerlinInputFiles <- function(mapfile, datfile, freqfile, pedfile, nNotMarker=5,  limitCentiMorgan=0, freqThreshold=0, suffix="_thinned") {
  nmarker.thin <- NULL
  if (is.null(nmarker.thin)) {
    ##    nmarker.thin <- as.integer(system(paste("perl -pe '}{$_=$.'", datfile), intern=TRUE))
    nmarker.thin <- length(count.fields(datfile))
  }
  cmd <- paste("perl", paste(get("FEST.perlpath", envir=topenv(), inherits = FALSE), "thinning.pl", sep=""),
               pedfile, datfile, mapfile, freqfile, nNotMarker, limitCentiMorgan, freqThreshold, nmarker.thin, suffix)
  system(cmd)
  
  datfile.thinned <- paste(datfile, suffix, sep="")
  ##  nMarker <- as.integer(system(paste("perl -pe '}{$_=$.'", datfile.thinned), intern=TRUE))
  nMarker <- length(count.fields(datfile.thinned))
  invisible(nMarker)
}

#######################################################################
#######################################################################
### Functions for reading Affymetrix frequency files
### and genetic maps,and for selecting subset of data
### and writing input merlin files
#######################################################################
#######################################################################

## Read the 500 K Affymetrix map and frequency data
## Make for.dat, for.map, for.freq (preceeded with prefixMerlin.input)
ReadAffy <- function(chr=c(1:22), convertToCM=TRUE, extrapolate=FALSE) {
  genMap <- ReadGeneticMaps(chr)

  if (extrapolate) {
    warning("extrapolate option is not implemented (should include call to approxExtrap)")
    extrapolate <- FALSE
  }
  path <- file.path("../Data/Affy500K_Allele_Frequency_Files/")
  current <- setwd(path)
  files <- paste("chr", chr, ".freq.txt", sep="") ## TEMPORARY
  n <- length(files)
  
  affy <- vector("list", n)
  nMarker.sum <- 0
  nSnp0.1.sum <- 0
  length.sum <- 0
  mat.print <- NULL
  for (i in 1:n) {
    f <- files[i]
    a <- read.table(f, header=FALSE)

    ind.na <- a[,5] == -1 # missing
    ind01 <- a[,5] == 0 | a[,5] == 1 # allele frequencies: 0,1
    a <- a[!ind.na & !ind01,]
    
    ind <- a[,3] == "A" | a[,3] == "T"

    a5 <- a[,5]
    a5[!ind] <- 1 - a[!ind,5]
    a6 <- 1 - a5

    if (convertToCM) { # convert to centimorgan
      ##      ind <- which(min(a[,2]) < min(genMap[[i]][,1]) | max(a[,2] > max(genMap[[i]][,1])))      
      ##      if (extrapolate)
      ##        a[,2] <- approxExtrap(genMap[[i]][,1], genMap[[i]][,2], a[,2], ties="ordered")$y
      ##      else {
      a[,2] <- approx(genMap[[i]][,1], genMap[[i]][,2], a[,2], ties="ordered")$y # NA for those outside SNP's outside interpolation range
      ##      }

      nOutsideRange <- sum(is.na(a[,2]))
      if (nOutsideRange > 0) {
        print(paste(nOutsideRange, "SNPs outside range for chromosome", i))
      }
    }
    
    affy[[i]] <- na.omit(cbind(a[,1:2], A=a5, C=a6))
    names(affy[[i]]) <- c("SNP", "cM", "A", "C")

    nSnp0.1 <- sum(pmin(a5,a6) > 0.1)
    cat("Cromosome ", i, " Number of SNP's ", nrow(affy[[i]]), " cM-range ", range(affy[[i]][,2]), " length ", diff(range(affy[[i]][,2])), "\n")
    mat.print <- rbind(mat.print, c(nrow(affy[[i]]), nSnp0.1, diff(range(affy[[i]][,2]))))

    nMarker.sum <- nMarker.sum + nrow(affy[[i]])
    nSnp0.1.sum <- nSnp0.1.sum + nSnp0.1
    length.sum <- length.sum + diff(range(affy[[i]][,2]))
  }
  setwd(current)

##  mat.print <- rbind(mat.print, c(nMarker.sum, nSnp0.1.sum, length.sum))
##  mat.print <- rbind(mat.print, c(nMarker.sum/n, nSnp0.1.sum/n, length.sum/n))
##  colnames(mat.print) <- c("Number of SNPs(MAF $>$ 0.0)", "Number of SNPs (MAF $>$ 0.1)", "Length (cM)")
##  rownames(mat.print) <- c(1:n, "Total", "Average")
##  latex.table(mat.print, "affySummary", caption = "Summary of 500K Affymetrix data.")

  cat("Total number of markers: ", nMarker.sum, "\n")
  cat("Total chromosome length: ", length.sum, "\n")
  cat("Average chromosome length: ", length.sum/n, "\n")
  

  invisible(affy)
}


## Writes dat, freq and map input file to Merlin
## freqObj: list where each element correspond to a chromosome
## freqObj[[i]] is a data frame with columns: SNP, cM, A, C
WriteInputFiles <- function(freqObj, chr) {
  n <- length(freqObj)
### WRITE DAT FILE
  ##  cat("Write dat file..\n")
  datfile <- paste(prefixMerlin.input, ".dat", sep="")
  appendFile <- FALSE
  for (i in 1:n) {
    nr <- nrow(freqObj[[i]])
    markernames <- as.character(freqObj[[i]][,1])
    tab <- data.frame(trait=rep("M", nr), marker=markernames)
    write.table(tab, file=datfile, col.names=FALSE, row.names=FALSE, quote=FALSE, append=appendFile)
    appendFile <- TRUE
  }
  
  write("A locus1", file=datfile, append=TRUE)

### WRITE FREQ FILE
  ##  cat("Write freq file..\n")
  freqfile <- paste(prefixMerlin.input, ".freq", sep="")
  appendFile <- FALSE
  for (i in 1:n) {
    nr <- nrow(freqObj[[i]])
    markernames <- as.character(freqObj[[i]][,1])
    tabm <- cbind(trait=rep("M", nr), markernames, rep("", nr))

    tabf <- cbind(rep("F", nr), freqObj[[i]][,3], freqObj[[i]][,4])

    names(tabf) <- c("trait", "freq1", "freq2")
    names(tabm) <- c("trait", "freq1", "freq2")
    tab <- rbind(tabm, tabf)
  
    ind.par <- (1:nr)*2
    ind.odd <- (1:nr)*2-1
    
    ind <- rep(NA, 2*nr)
    ind[ind.odd] <- 1:nr
    ind[ind.par] <- (nr+1):(2*nr)

    tab.print <- data.frame(tab[ind,])
    
    
    write.table(tab.print, file=freqfile, col.names=FALSE, row.names=FALSE, quote=FALSE, append=appendFile)
    appendFile <- TRUE
  }


### WRITE MAP FILE
  ##  cat("Write map file..\n")
  mapfile <- paste(prefixMerlin.input, ".map", sep="")
  write("CHR  MARKER    POS", file=mapfile)
  for (i in 1:n) {
    nr <- nrow(freqObj[[i]])
    write.table(cbind(rep(chr[i], nr), freqObj[[i]][,1:2]), file=mapfile, col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)
  }
  invisible(freqObj)
}


ReadGeneticMaps <- function(chr=c(1:22)) {
  filenames <- paste("../Data/GeneticMaps/smooth_map_b36/chr", chr, ".sm.map2", sep="")
  par(mfrow=c(5,5))
  n <- length(filenames)
  genMap <- vector("list", n)
  for (i in 1:n) {
    f <- filenames[i]
    x <- read.table(f, header=TRUE)
    plot(x[,2], x[,4], type="l", xlab="BP", ylab="CM", main=paste("Chromosome", i))
    genMap[[i]] <- x[,c(2,4)]
  }
  genMap
}


SelectSNPs <- function(freqObj, neach, chr=c(1:22), threshold=0.1, limitCentiMorgan=0, printInfo=FALSE) {
  nchr <- length(chr)
  freqObjSub <- vector("list", nchr)
  for (i in 1:nchr) {
    a <- freqObj[[chr[i]]]
    ind <- which(a[,3] >= threshold & a[,3] <= 1-threshold)
    n <- length(ind)
    ##    n.intervals <- neach
    if (!is.finite(neach)) {
      ind.pos <- 1:n # all
    }
    else if (neach > 1) {
      ind.pos <- 0:(neach-1)
      ind.pos <- as.integer(((n-1)/(neach-1)) * ind.pos + 1)
    }
    else {
      ind.pos <- as.integer((n+1)/2)
    }
    freqObjSub[[i]] <- a[ind[ind.pos],]
    if (printInfo) {
      cat("Summary for chromosome ", i, ": ")
      print(summary(apply(freqObjSub[[i]][,3:4], 1, min)))
      ##    print(dim(freqObjSub[[i]]))
    }
  }
  invisible(freqObjSub)
}

SummaryAffymetrix <- function(affy) {
  n <- NULL
  maf.min <- NULL
  
  for (i in 1:22) {
    maf <- pmin(affy[[i]][,3], affy[[i]][,4])
    n <- rbind(n, c(length(maf), sum(maf>0), sum(maf>0.001), sum(maf>0.01), sum(maf>0.02), sum(maf>0.05)))
    maf.min <- min(maf.min, maf)
  }
  n.sum <- apply(n, 2, sum)
  
  n.sum
}

SortMerlinInputFiles <- function(mapfile, datfile, freqfile, pedfile, nNotMarker=5, prefix="sorted_", chr=NULL, excludeSNP=NULL) {
  map <- read.table(mapfile, header=TRUE)
  key <- map[,1]*1e8 + map[,3]
  ind <- order(key)

  if (!is.null(chr)) {
    indchr <- which(map[ind, 1] %in% chr)
    ind <- ind[indchr]
  }
  if (!is.null(excludeSNP)) {
    indsnp <- which(!(map[ind,2] %in% excludeSNP))
    ind <- ind[indsnp]
    
  }

  write.table(map[ind,], file=paste(prefix, mapfile, sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)

  dat <- read.table(datfile, header=FALSE)
  write.table(dat[ind,], file=paste(prefix, datfile, sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)

  nRec <- length(scan(freqfile, nlines=2, quiet = TRUE, what=character()))
  nFreq <- nRec - 3
  freqfile.sorted <- paste(prefix, freqfile, sep="")
  file.copy(freqfile, freqfile.sorted, overwrite=TRUE)
  com <- paste("perl -pi~ -e \"s/M/M ", paste(rep("2.01", nFreq-1), collapse=" "), "/g\" ", freqfile.sorted, sep="") ##makes freq easy to read
  system(com)
  freq <- read.table(freqfile.sorted, header=FALSE)

  n <- length(ind)
  ind.freq <- rep(NA, len=2*n)
  ind.freq[seq(1,2*n, by=2)] <- 2*ind -1
  ind.freq[seq(2,2*n, by=2)] <- 2*ind
  write.table(freq[ind.freq,], file=freqfile.sorted, quote=FALSE, row.names=FALSE, col.names=FALSE)
  com <- paste("perl -pi~ -e \"s/M ", paste(rep("2.01", nFreq-1), collapse=" "), "/M/g\" ",freqfile.sorted, sep="")
  system(com) ##fixes freq file back

  ind.ped <- c(1:nNotMarker, ind.freq + nNotMarker)
  ##  nr <- as.integer(system(paste("perl -pe '}{$_=$.'", pedfile), intern=TRUE))
  ##  nr <- count.fields(pedfile)
  nc <- length(scan(pedfile, nlines=1))
  ##  ped <- matrix(scan(pedfile, quiet=TRUE), byrow=TRUE, nrow=nr)
  ped <- matrix(scan(pedfile, quiet=TRUE), byrow=TRUE, ncol=nc)
  write.matrix(ped[,ind.ped], file=paste(prefix, pedfile, sep="")) #, quote=FALSE, row.names=FALSE, col.names=FALSE)
}
