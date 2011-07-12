setClass("Pedigree",
         representation(fam="integer", subj="integer", fatherID="integer", motherID="integer", gender="integer", affected="integer", typed="integer"),
         prototype=list(fam=integer(), subj=integer(), fatherID=integer(), motherID=integer(), gender=integer(), affected=integer(), typed=integer()),
         validity= function(object) {
##           f <- names(getSlots("Pedigree"))
##           n <- length(f)
##           for (i in 2:n)
##             if (f[i] != "typed")
##               if(length(eval(parse(text=paste("object@",f[1], sep="")))) !=
##                  length(eval(parse(text=paste("object@",f[i], sep="")))))
##                 stop("specified ", f[1], " and ", f[i], " of different lengths")
           TRUE
         })

setClass("Model",
         representation(true="character", alternative="list"))


setClass("SimStudyObject",
         representation(posterior="list", logLik="list", nsim="integer", nmarker="integer",
                        maf="numeric", freqThreshold="numeric",
                        model="Model"),
         validity=function(object) TRUE)
