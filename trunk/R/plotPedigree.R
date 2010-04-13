##PlotPedigree <- function(reltype, pedfile, plotfile, devtype=NULL) {
PlotPedigree <- function(reltype, plotfile, devtype=NULL) {
  pedfile <- paste(prefix.tmpfiles, reltype, sep="")
  on.exit(file.remove(pedfile))
  MakePedigree(reltype, pedfile)
  
  ##p1 <- scan("pedigreeFil", what=list(0,0,0,0,0,"",""))
  p1 <- scan(pedfile, what=list(0,0,0,0,0,0), quiet = TRUE)
  p2 <- as.data.frame(p1)
  ##names(p2) <- c("id","fid","mid","sex","aff","GABRB1","D4S1645")
  names(p2) <- c("fam", "id","fid","mid","sex","aff")
  n <- nrow(p2)
  ##  status <- c(rep(1, n-2), 0, 0)
  status <- p2$aff == 0

  p3 <- pedigree.mine(p2$id, p2$fid, p2$mid, p2$sex-1, p2$aff, status)
  id <- p3$id
  n <- length(id)
  id <- vector("character", n)
  id[status==FALSE] <- c("A", "B")
  if (!is.null(devtype)) {
    if (devtype == "pdf")
      pdf(plotfile)
    else if (devtype == "postscript")
      postscript(plotfile, horizontal=FALSE, onefile=FALSE) ## encapsulated postscript
  }
  par(xpd=TRUE)
  plot.pedigree.mine(p3, id=id, cex=2, symbolsize=2)
  if (!is.null(devtype))
    dev.off()
}

pedigree.mine <-   function (id, dadid, momid, sex, affected, status, relations) 
{
  n <- length(id)
  if (length(momid) != n) 
    stop("Mismatched lengths, id and momid")
  if (length(dadid) != n) 
    stop("Mismatched lengths, id and momid")
  if (length(sex) != n) 
    stop("Mismatched lengths, id and sex")
  if (any(is.na(id))) 
    stop("Missing value for the id variable")
  if (!is.numeric(id) && any(as.character(id) == "")) 
    stop("Empty string is not allowed as the id variable")
  if (is.factor(sex)) 
    sex <- as.character(sex)
  codes <- c("male", "female", "unknown", "terminated")
  if (is.character(sex)) 
    sex <- charmatch(casefold(sex, upper = F), codes, nomatch = 3)
  if (min(sex) == 0) 
    sex <- sex + 1
  sex <- ifelse(sex < 1 | sex > 4, 3, sex)
  if (all(sex > 2)) 
    stop("Invalid values for 'sex'")
  else if (mean(sex == 3) > 0.25) 
    warning("More than 25% of the gender values are 'unknown'")
  sex <- factor(sex, 1:4, labels = codes)
  father <- match(dadid, id, nomatch = 0)
  if (any(sex[father] != "male")) {
    who <- unique((id[father])[sex[father] != "male"])
    stop(paste("Id not male, but is a father:", paste(who, 
                                                      collapse = " ")))
  }
  mother <- match(momid, id, nomatch = 0)
  if (any(sex[mother] != "female")) {
    who <- unique((id[mother])[sex[mother] != "female"])
    stop(paste("Id not female, but is a mother:", paste(who, 
                                                        collapse = " ")))
  }
  depth <- kindepth(id, momid, dadid, align = T)
  temp <- list(id = id, momid = momid, dadid = dadid, sex = sex, 
               depth = depth)
  if (!missing(affected)) {
    if (is.matrix(affected)) {
      if (nrow(affected) != n) 
        stop("Wrong number of rows in affected")
      if (is.logical(affected)) 
        affected <- 1 * affected
    }
    else {
      if (length(affected) != n) 
        stop("Wrong length for affected")
      if (is.logical(affected)) 
        affected <- as.numeric(affected)
      if (is.factor(affected)) 
        affected <- as.numeric(affected) - 1
    }
    ##        if (max(affected) > min(affected)) 
    affected <- affected - min(affected)
    if (!all(affected == 0 | affected == 1 | affected == 
             2)) 
      stop("Invalid code for affected status")
    temp$affected <- affected
  }
  if (!missing(status)) {
    if (length(status) != n) 
      stop("Wrong length for status")
    if (any(status != 0 & status != 1)) 
      stop("Invalid status code")
    temp$status <- status
  }
  if (!missing(relations)) {
    if (is.matrix(relations)) {
      if (ncol(relations) != 3) 
        stop("Relations matrix must have 3 columns")
      id1 <- relations[, 1]
      id2 <- relations[, 2]
      code <- relations[, 3]
    }
    else if (is.list(relations)) {
      id1 <- relations$id1
      id2 <- relations$id2
      code <- relations$code
      if (is.null(id1) || is.null(id2) || is.null(code)) 
        stop("Relations list must have id1, id2, and code")
      if (length(id1) != length(id2) || (length(id1) != 
                  length(code))) 
        stop("Id1, id2 and code in the relations list are different lengths")
    }
    else stop("Relations argument must be a matrix or a list")
    if (any(code != floor(code)) || min(code) < 1 || max(code) > 
        4) 
      stop("Invalid relationship code")
    temp1 <- match(id1, id, nomatch = 0)
    temp2 <- match(id2, id, nomatch = 0)
    if (any(temp1 == 0 | temp2 == 0)) 
      stop("Subjects in relationships that are not in the pedigree")
    if (any(temp1 == temp2)) {
      who <- temp1[temp1 == temp2]
      stop(paste("Subject", id[who], "is their own spouse or twin"))
    }
    if (any(code < 3)) {
      twins <- (code < 3)
      if (any(momid[temp1[twins]] != momid[temp2[twins]])) 
        stop("Twins with non-identical parents")
      if (any(dadid[temp1[twins]] != dadid[temp2[twins]])) 
        stop("Twins with non-identical parents")
    }
    if (any(code == 1)) {
      mztwins <- (code == 1)
      if (any(sex[mztwins] != sex[mztwins])) 
        stop("MZ Twins with different genders")
    }
    temp$relation <- cbind(temp1, temp2, code)
  }
  temp <- c(temp, list(hints = autohint(temp)))
  oldClass(temp) <- "pedigree"
  temp
}

plot.pedigree.mine <- function (x, id = x$id, sex = x$sex, status = x$status, affected = x$affected, 
                                cex = 1, col = rep(1, length(x$id)), symbolsize = 1, branch = 0.6, 
                                packed = T, align = packed, width = 8, density = c(-1, 50, 
                                                                         70, 90), mar = c(4.1, 1, 4.1, 1), angle = c(90, 70, 50, 
                                                                                                             0), keep.par = F, ...) 
{
  maxlev <- max(x$depth) + 1
  n <- length(x$depth)
  sex <- as.numeric(sex)
  sex[is.na(sex)] <- 3
  if (any(sex < 1 | sex > 4)) 
    stop("Invalid sex code")
  if (length(sex) != n) 
    stop("Wrong length for sex")
  if (is.null(status)) 
    status <- rep(0, n)
  else {
    if (!all(status == 0 | status == 1 | status == 2)) 
      stop("Invalid status code")
    if (length(status) != n) 
      stop("Wrong length for status")
  }
  if (!is.null(id)) {
    if (length(id) != n) 
      stop("Wrong length for id")
  }
  if (is.null(affected)) {
    affected <- matrix(0, nrow = n)
  }
  else {
    if (is.matrix(affected)) {
      if (nrow(affected) != n) 
        stop("Wrong number of rows in affected")
      if (is.logical(affected)) 
        affected <- 1 * affected
    }
    else {
      if (length(affected) != n) 
        stop("Wrong length for affected")
      if (is.logical(affected)) 
        affected <- as.numeric(affected)
      if (is.factor(affected)) 
        affected <- as.numeric(affected) - 1
    }
    if (max(affected) > min(affected)) 
      affected <- matrix(affected - min(affected), nrow = n)
    else affected <- matrix(affected, nrow = n)
    if (!all(affected == 0 | affected == 1 | affected == 
             2)) 
      stop("Invalid code for affected status")
  }
  if (!missing(status)) {
    if (length(status) != n) 
      stop("Wrong length for affected")
    if (any(status != 0 & status != 1)) 
      stop("Invalid status code")
  }
  symbol <- c(0, 1, 5, 2)[sex]
  points.sym <- function(x, y, symbol, status, affected, size, 
                         col, aspect, angle, density, adj) {
    circle <- function(cx, cy, r, code = 0, col = 1, angle, 
                       density, ...) {
      z <- (0:360 * pi)/180
      pin <- par()$pin
      usr <- par()$usr
      adj <- (pin[2]/pin[1])/((usr[4] - usr[3])/(usr[2] -  usr[1]))
      x <- sin(z) * r + cx
      y <- cos(z) * r * 1/adj + cy
      if (sum(code, na.rm = T) == 0) 
        polygon(x, y, border = T, density = 0, col = col, 
                ...)
      else {
        if (length(code) == 1) 
          polygon(x, y, border = T, col = col, density = density[1], 
                  angle = angle[1], ...)
        if (length(code) == 2) {
          polygon(x, y, border = T, density = 0, ...)
          z <- (0:180 * pi)/180
          x <- sin(z) * r + cx
          y <- cos(z) * r * 1/adj + cy
          polygon(x, y, border = T, col = col, density = density[2] * 
                  code[2], angle = angle[2], ...)
          z <- (180:360 * pi)/180
          x <- sin(z) * r + cx
          y <- cos(z) * r * 1/adj + cy
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
        }
        if (length(code) == 3) {
          polygon(x, y, border = T, density = 0, ...)
          z <- (0:90 * pi)/180
          x <- c(cx, sin(z) * r + cx)
          y <- c(cy, cos(z) * r * 1/adj + cy)
          polygon(x, y, border = T, col = col, density = code[3] * 
                  density[3], angle = angle[3], ...)
          z <- (180:270 * pi)/180
          x <- c(cx, sin(z) * r + cx)
          y <- c(cy, cos(z) * r * 1/adj + cy)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          z <- (270:360 * pi)/180
          x <- c(cx, sin(z) * r + cx)
          y <- c(cy, cos(z) * r * 1/adj + cy)
          polygon(x, y, border = T, col = col, density = code[2] * 
                  density[2], angle = angle[2], ...)
        }
        if (length(code) == 4) {
          polygon(x, y, border = T, density = 0, ...)
          z <- (0:90 * pi)/180
          x <- c(cx, sin(z) * r + cx)
          y <- c(cy, cos(z) * r * 1/adj + cy)
          polygon(x, y, border = T, col = col, density = code[3] * 
                  density[3], angle = angle[3], ...)
          z <- (90:180 * pi)/180
          x <- c(cx, sin(z) * r + cx)
          y <- c(cy, cos(z) * r * 1/adj + cy)
          polygon(x, y, border = T, col = col, density = code[4] * 
                  density[4], angle = angle[4], ...)
          z <- (180:270 * pi)/180
          x <- c(cx, sin(z) * r + cx)
          y <- c(cy, cos(z) * r * 1/adj + cy)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          z <- (270:360 * pi)/180
          x <- c(cx, sin(z) * r + cx)
          y <- c(cy, cos(z) * r * 1/adj + cy)
          polygon(x, y, border = T, col = col, density = code[2] * 
                  density[2], angle = angle[2], ...)
        }
        if (length(code) > 4) 
          stop("Can only plot up to 4 levels of codes")
      }
      invisible()
    }
    square <- function(cx, cy, r, code = 0, col = 1, angle, 
                       density, ...) {
      pin <- par()$pin
      usr <- par()$usr
      adj <- (pin[2]/pin[1])/((usr[4] - usr[3])/(usr[2] -  usr[1]))
      x <- cx + c(-r, -r, r, r)
      y <- cy + (1/adj) * c(-r, r, r, -r)
      if (sum(code, na.rm = T) == 0) 
        polygon(x, y, border = T, density = 0, col = col, 
                ...)
      else {
        if (length(code) == 1) 
          polygon(x, y, border = T, col = col, density = density[1], 
                  angle = angle[1], ...)
        if (length(code) == 2) {
          polygon(x, y, border = T, density = 0, ...)
          x <- cx + c(-r, -r, 0, 0)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          x <- cx + c(0, 0, r, r)
          polygon(x, y, border = T, col = col, density = density[2] * 
                  code[2], angle = angle[2], ...)
        }
        if (length(code) == 3) {
          polygon(x, y, border = T, density = 0, ...)
          x <- cx + c(-r, -r, 0, 0)
          y <- cy + (1/adj) * c(-r, 0, 0, -r)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          x <- cx + c(-r, -r, 0, 0)
          y <- cy + (1/adj) * c(0, r, r, 0)
          polygon(x, y, border = T, col = col, density = code[2] * 
                  density[2], angle = angle[2], ...)
          x <- cx + c(0, 0, r, r)
          y <- cy + (1/adj) * c(0, r, r, 0)
          polygon(x, y, border = T, col = col, density = density[3] * 
                  code[3], angle[3], ...)
        }
        if (length(code) == 4) {
          polygon(x, y, border = T, density = 0, ...)
          x <- cx + c(-r, -r, 0, 0)
          y <- cy + (1/adj) * c(-r, 0, 0, -r)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          x <- cx + c(-r, -r, 0, 0)
          y <- cy + (1/adj) * c(0, r, r, 0)
          polygon(x, y, border = T, col = col, density = code[2] * 
                  density[2], angle = angle[2], ...)
          x <- cx + c(0, 0, r, r)
          y <- cy + (1/adj) * c(0, r, r, 0)
          polygon(x, y, border = T, col = col, density = code[3] * 
                  density[3], angle = angle[3], ...)
          x <- cx + c(0, 0, r, r)
          y <- cy + (1/adj) * c(-r, 0, 0, -r)
          polygon(x, y, border = T, col = col, density = code[4] * 
                  density[4], angle[4], ...)
        }
        if (length(code) > 4) 
          stop("Can only plot up to 4 levels of codes")
      }
      invisible()
    }
    diamond <- function(cx, cy, r, code = 0, col = 1, angle, 
                        density, ...) {
      pin <- par()$pin
      usr <- par()$usr
      adj <- (pin[2]/pin[1])/((usr[4] - usr[3])/(usr[2] -    usr[1]))
      x <- cx + c(-r, 0, r, 0)
      y <- cy + (1/adj) * c(0, r, 0, -r)
      if (sum(code, na.rm = T) == 0) 
        polygon(x, y, border = T, density = 0, col = col, 
                ...)
      else {
        if (length(code) == 1) 
          polygon(x, y, border = T, col = col, density = density[1], 
                  angle = angle[1], ...)
        if (length(code) == 2) {
          polygon(x, y, border = T, density = 0, ...)
          x <- cx + c(-r, 0, 0)
          y <- cy + (1/adj) * c(0, r, -r)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          x <- cx + c(0, 0, r)
          y <- cy + (1/adj) * c(r, -r, 0)
          polygon(x, y, border = T, col = col, density = density[2] * 
                  code[2], angle = angle[2], ...)
        }
        if (length(code) == 3) {
          polygon(x, y, border = T, density = 0, ...)
          x <- cx + c(-r, 0, 0)
          y <- cy + (1/adj) * c(0, 0, -r)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          x <- cx + c(-r, 0, 0)
          y <- cy + (1/adj) * c(0, r, 0)
          polygon(x, y, border = T, col = col, density = code[2] * 
                  density[2], angle = angle[2], ...)
          x <- cx + c(0, 0, r)
          y <- cy + (1/adj) * c(0, r, 0)
          polygon(x, y, border = T, col = col, density = density[3] * 
                  code[3], angle = angle[3], ...)
        }
        if (length(code) == 4) {
          polygon(x, y, border = T, density = 0, ...)
          x <- cx + c(-r, 0, 0)
          y <- cy + (1/adj) * c(0, 0, -r)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          x <- cx + c(-r, 0, 0)
          y <- cy + (1/adj) * c(0, r, 0)
          polygon(x, y, border = T, col = col, density = code[2] * 
                  density[2], angle = angle[2], ...)
          x <- cx + c(0, 0, r)
          y <- cy + (1/adj) * c(0, r, 0)
          polygon(x, y, border = T, col = col, density = code[3] * 
                  density[3], angle = angle[3], ...)
          x <- cx + c(0, 0, r)
          y <- cy + (1/adj) * c(-r, 0, 0)
          polygon(x, y, border = T, col = col, density = code[4] * 
                  density[4], angle = angle[4], ...)
        }
        if (length(code) > 4) 
          stop("Can only plot up to 4 levels of codes")
      }
      invisible()
    }
    triangle <- function(cx, cy, r, code = 0, col = 1, angle, 
                         density, ...) {
      pin <- par()$pin
      usr <- par()$usr
      adj <- (pin[2]/pin[1])/((usr[4] - usr[3])/(usr[2] - usr[1]))
      r <- r * 1.25
      a <- 3 * r/sqrt(3)
      b <- r/2
      x <- cx + c((-1/2) * a, (1/2) * a, 0)
      y <- cy + (1/adj) * c(-b, -b, r)
      if (sum(code, na.rm = T) == 0) 
        polygon(x, y, border = T, density = 0, col = col, 
                ...)
      else {
        if (length(code) == 1) 
          polygon(x, y, border = T, col = col, density = density[1], 
                  angle = angle[1], ...)
        if (length(code) == 2) {
          polygon(x, y, border = T, density = 0, ...)
          x <- cx + c((-1/2) * a, 0, 0)
          y <- cy + (1/adj) * c(-b, -b, r)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          x <- cx + c(0, (1/2) * a, 0)
          y <- cy + (1/adj) * c(-b, -b, r)
          polygon(x, y, border = T, col = col, density = density[2] * 
                  code[2], angle = angle[2], ...)
        }
        if (length(code) == 3) {
          polygon(x, y, border = T, density = 0, ...)
          midx <- (r * (0.5) * a)/(b + r)
          x <- cx + c((-1/2) * a, -midx, 0, 0)
          y <- cy + (1/adj) * c(-b, 0, 0, -b)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          x <- cx + c(-midx, 0, 0)
          y <- cy + (1/adj) * c(0, r, 0)
          polygon(x, y, border = T, col = col, density = code[2] * 
                  density[2], angle = angle[2], ...)
          x <- cx + c(0, 0, midx)
          y <- cy + (1/adj) * c(0, r, 0)
          polygon(x, y, border = T, col = col, density = density[3] * 
                  code[3], angle = angle[3], ...)
        }
        if (length(code) == 4) {
          polygon(x, y, border = T, density = 0, ...)
          midx <- (r * (0.5) * a)/(b + r)
          x <- cx + c((-1/2) * a, -midx, 0, 0)
          y <- cy + (1/adj) * c(-b, 0, 0, -b)
          polygon(x, y, border = T, col = col, density = code[1] * 
                  density[1], angle = angle[1], ...)
          x <- cx + c(-midx, 0, 0)
          y <- cy + (1/adj) * c(0, r, 0)
          polygon(x, y, border = T, col = col, density = code[2] * 
                  density[2], angle = angle[2], ...)
          x <- cx + c(0, 0, midx)
          y <- cy + (1/adj) * c(0, r, 0)
          polygon(x, y, border = T, col = col, density = code[3] * 
                  density[3], angle = angle[3], ...)
          x <- cx + c(0, midx, (1/2) * a, 0)
          y <- cy + (1/adj) * c(0, 0, -b, -b)
          polygon(x, y, border = T, col = col, density = code[4] * 
                  density[4], angle = angle[4], ...)
        }
        if (length(code) > 4) 
          stop("Can only plot up to 4 levels of codes")
      }
      invisible()
    }
    for (i in seq(along=x)[!is.na(match(symbol, c(0, 15)))]) {
      square(x[i], y[i], size, code = (affected[i, , drop = F]), 
             col = col[i], angle = angle, density = density)
    }
    for (i in seq(along=x)[!is.na(match(symbol, c(1, 16)))]) {
      circle(x[i], y[i], size, code = (affected[i, , drop = F]), 
             col = col[i], angle = angle, density = density)
    }
    for (i in seq(along=x)[!is.na(match(symbol, c(5, 18)))]) {
      diamond(x[i], y[i], size, code = (affected[i, , drop = F]), 
              col = col[i], angle = angle, density = density)
    }
    for (i in seq(along=x)[!is.na(match(symbol, c(2, 17)))]) {
      triangle(x[i], y[i], size, code = (affected[i, , 
                                                  drop = F]), col = col[i], angle = angle, density = density)
    }
    who <- (status == 1)
    if (any(who)) {
      deltax <- size * cos(pi/4) + aspect
      deltay <- (1/adj) * (size * sin(pi/4) + aspect)
      segments(x[who] - deltax, y[who] - deltay, x[who] + 
               deltax, y[who] + deltay)
    }
  }
  plist <- align.pedigree(x, packed = packed, width = width, 
                          align = align)
  who <- (plist$nid > 0)
  xp <- plist$pos[who]
  yp <- -(row(plist$nid))[who]
  np <- plist$nid[who]

  oldpar <- par(mar=c(0,0,0,0))
  if (!keep.par)
    on.exit(par(oldpar))
  ##    textoff <- -(1/adj) * (radius) - 0.5 * cex * (par()$cin[2]/7)
  radius0 <-  0.08 * symbolsize
  xlim <- range(xp) + 2.5*c(-radius0, radius0)
  ylim <- range(yp) + 2.5*c(-radius0, radius0)
  plot(xp, yp, axes = F, type = "n", xlab = "", ylab = "", xlim=xlim, ylim=ylim,
       ...)
  pin <- par()$pin
  usr <- par()$usr
  adj <- (pin[2]/pin[1])/((usr[4] - usr[3])/(usr[2] - usr[1]))
  symbolsize <- symbolsize * adj
  radius <- 0.08 * symbolsize
  deltaInch <- diff(ylim)/pin[2]
  textoff <- (radius/adj) + 0.5 * cex * par()$cin[2] * deltaInch
  textoff <- 1.1 * textoff
  delta <- 0.1
  aspect <- 0.05 * symbolsize
  points.sym(xp, yp, symbol[np], status[np], affected[np, , 
                                                      drop = F], radius, col[np], aspect, angle = angle, density = density, 
             adj = adj)
  if (!is.null(id)) {
    text(xp, yp - textoff, id[np], cex = cex, col = col[np])
  }
  for (i in 1:maxlev) {
    if (any(plist$spouse[i, ] > 0)) {
      temp <- (1:ncol(plist$spouse))[plist$spouse[i, ] > 
                                     0]
      segments(plist$pos[i, temp] + radius, rep(-i, length(temp)), 
               plist$pos[i, temp + 1] - radius, rep(-i, length(temp)))
    }
  }
  for (i in 2:maxlev) {
    zed <- unique(plist$fam[i, ])
    zed <- zed[zed > 0]
    for (fam in zed) {
      xx <- plist$pos[i - 1, fam + 0:1]
      parentx <- mean(xx)
      if (!is.null(plist$twins)) {
        tw.left <- (plist$twins[i, ] > 0 & plist$twins[i, 
                                                       ] < 4) & plist$fam[i, ] == fam
        mz.left <- (plist$twins[i, ] == 1) & plist$fam[i, 
                                 ] == fam
        un.left <- (plist$twins[i, ] == 3) & plist$fam[i, 
                                 ] == fam
        tw.right <- rep(F, length(tw.left))
        famlst <- plist$fam[i, ]
        twloc <- (1:length(tw.right))[tw.left]
        famloc <- (1:length(tw.right))[famlst == fam]
        flag <- NULL
        for (i in twloc) flag <- c(flag, min(famloc[i < 
                                                    famloc]))
        tw.right[flag] <- T
        mz.right <- rep(F, length(mz.left))
        mzloc <- (1:length(mz.right))[mz.left]
        flag <- NULL
        for (i in mzloc) flag <- c(flag, min(famloc[i < 
                                                    famloc]))
        mz.right[flag] <- T
        un.right <- rep(F, length(un.left))
        unloc <- (1:length(un.right))[un.left]
        flag <- NULL
        for (i in unloc) flag <- c(flag, min(famloc[i < 
                                                    famloc]))
        un.right[flag] <- T
      }
      if (is.null(plist$twins)) {
        twn <- length(plist$fam[i, ] == fam)
        tw.left <- rep(F, twn)
        tw.right <- rep(F, twn)
        mz.left <- rep(F, twn)
        mz.right <- rep(F, twn)
        un.left <- rep(F, twn)
        un.right <- rep(F, twn)
      }
      tw <- tw.left | tw.right
      mz <- mz.left | mz.right
      un <- un.left | un.right
      who <- (plist$fam[i, ] == fam) & !tw
      xx <- plist$pos[i, who]
      yy <- rep(-i, length = sum(who))
      ww <- plist$nid[i, who]
      xx.l <- plist$pos[i, tw.left]
      yy.l <- rep(-i, length = sum(tw.left))
      ww.l <- plist$nid[i, tw.left]
      xx.r <- plist$pos[i, tw.right]
      yy.r <- rep(-i, length = sum(tw.right))
      ww.r <- plist$nid[i, tw.right]
      segments(xx, yy + (1/adj) * radius, xx, yy + 3 * 
               delta, col = col[ww])
      xx.lm <- plist$pos[i, mz.left]
      yy.lm <- rep(-i, length = sum(mz.left))
      xx.rm <- plist$pos[i, mz.right]
      yy.rm <- rep(-i, length = sum(mz.right))
      xx.lu <- plist$pos[i, un.left]
      yy.lu <- rep(-i, length = sum(un.left))
      xx.ru <- plist$pos[i, un.right]
      yy.ru <- rep(-i, length = sum(un.right))
      who <- (plist$fam[i, ] == fam)
      xx <- plist$pos[i, who]
      yy <- rep(-i, length = sum(who))
      ww <- plist$nid[i, who]
      xx2 <- plist$pos[i, ]
      xx2 <- xx2[who]
      if (sum(tw) >= 2) {
        tw.who <- tw[who]
        n <- length(tw.who)
        n2 <- sum(tw.who)
        flagL <- sum(tw.who[1:2]) == 2
        flagR <- sum(tw.who[c(n, n - 1)]) == 2
        xx2.save <- xx2
        if (flagL) 
          xx2[tw.who][1:2] <- rep(mean(xx2.save[tw.who][1:2]), 
                                  2)
        if (flagR) 
          xx2[tw.who][c(n2, n2 - 1)] <- rep(mean(xx2.save[tw.who][c(n2, 
                                                                    n2 - 1)]), 2)
      }
      segments(min(xx2), 3 * delta - i, max(xx2), 3 * delta - 
               i)
      x1 <- mean(range(xx))
      y1 <- 3 * delta - i
      if (branch == 0) 
        segments(x1, y1, parentx, -(i - 1))
      else {
        y2 <- 1 - i
        x2 <- parentx
        ydelta <- ((y2 - y1) * branch)/2
        segments(c(x1, x1, x2), c(y1, y1 + ydelta, y2 - 
                                  ydelta), c(x1, x2, x2), c(y1 + ydelta, y2 - 
                                                            ydelta, y2))
      }
    }
  }
  arcconnect <- function(x, y) {
    xx <- seq(x[1], x[2], length = 15)
    yy <- seq(y[1], y[2], length = 15) + 0.5 - (seq(-7, 7))^2/98
    lines(xx, yy, lty = 2)
  }
  xx <- table(np)
  xx <- xx[xx > 1]
  if (length(xx) > 0) {
    multi <- as.numeric(names(xx))
    for (i in 1:length(xx)) {
      x2 <- xp[np == multi[i]]
      y2 <- yp[np == multi[i]]
      nn <- xx[i]
      for (j in 1:(nn - 1)) arcconnect(x2[j + 0:1], y2[j + 
                                                       0:1])
    }
  }
  ckall <- x$id[is.na(match(x$id, x$id[plist$nid]))]
  if (length(ckall > 0)) 
    cat("Did not plot the following people:", ckall, "\n")
  tmp <- plist$nid[plist$nid != 0]
  xp2 <- xp[order(tmp)]
  yp2 <- yp[order(tmp)]
  invisible(list(plist = plist, object = x, x = xp2, y = yp2, 
                 textoff = textoff, symbolsize = symbolsize))
}
