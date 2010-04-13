## Halfsib-k-l complexity: k + l - 1
## Sibling-k-l complexity: k + l - 1 (should have been k + l + 2 by the 2*n-f formula...)
## Parentchild-k complexity: k - 1
maxBit <- 24
maxExtraGen.halfsib <- as.integer((maxBit+1)/2)
maxExtraGen.sib <-  as.integer((maxBit+1)/2)
maxExtraGen.parentchild <- maxBit
availableModels <- c("unrelated",
                     paste("HS-", 1:maxExtraGen.halfsib, sep=""),
                     paste("S-", 1:maxExtraGen.sib, sep=""),
                     paste("PC-", 1:maxExtraGen.parentchild, sep=""))
for (i in 1:maxBit) {
  gen <- 1:(maxBit+1-i)
  availableModels <- c(availableModels,
                       paste("HS-", i, "-", gen, sep=""))
}
for (i in 1:maxBit) {
  gen <- 1:(maxBit+1-i)
  availableModels <- c(availableModels,
                       paste("S-", i, "-", gen, sep=""))
}


availableAltModels <- c(availableModels, "true", "lower", "upper")
FEST.perlpath <- ""
prefix.tmpfiles <- "__TMP__"
prefixMerlin.merged <- paste(prefix.tmpfiles, "msim", sep="")
prefixMerlin.mergedAlt <- paste(prefix.tmpfiles, "msimAlt", sep="")
prefixMerlin.input <- paste(prefix.tmpfiles, "for", sep="")

.onLoad <- function(libname, pkgname)
{
  ##  require(quantreg)
  ##
  
  FEST.path <- .find.package("FEST")
  local.perlpath <- file.path("/exec/")
  FEST.perlpath <- paste(FEST.path, local.perlpath, sep="")
  ##  pkgEnv <- topenv()
  
  assign("FEST.perlpath", FEST.perlpath, envir=topenv())
  ##  assignInNamespace("FEST.perlpath", FEST.perlpath,  ns = "FEST")

  ##  assignInNamespace("haplosim", "haplosim.fix", pos="package:hapsim")
  ##  reassignInPackage("haplosim", pkgName="hapsim", haplosim.fix)
  haplosim <- haplosim.fix # hapsim does not have a name space
  
  ## Check if merlin is installed
  merlinOutput <- system("merlin", intern=TRUE, ignore.stderr=TRUE)
  if (length(merlinOutput) == 0) {
    warning("merlin is not installed! Must be installed for FEST to be operative.")
  }
  perlOutput <- system("perl --version", intern=TRUE, ignore.stderr=TRUE)
  if (length(perlOutput) == 0) {
    warning("perl is not installed! Must be installed for FEST to be operative.")
  }
    
}
