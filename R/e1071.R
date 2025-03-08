"bclust" <-
  function (x, centers = 2, iter.base = 10, minsize = 0,
            dist.method = "euclidian", hclust.method = "average",
            base.method = "kmeans", base.centers = 20,
            verbose = TRUE, final.kmeans = FALSE, docmdscale=FALSE,
            resample=TRUE, weights=NULL, maxcluster=base.centers, ...)
  {
    if (!require("mva")) stop("Could not load required package mva")
    if (!require("class"))
      stop("Could not load required package class from bundle VR")

    xr <- nrow(x)
    xc <- ncol(x)
    CLUSFUN <- get(base.method)

    object <- list(allcenters =
                     matrix(0, ncol = xc, nrow = iter.base * base.centers),
                   allcluster = NULL,
                   hclust = NULL,
                   members = NULL,
                   cluster = NULL,
                   centers = NULL,
                   iter.base = iter.base,
                   base.centers = base.centers,
                   prcomp = NULL,
                   datamean = apply(x, 2, mean),
                   colnames = colnames(x),
                   dist.method = dist.method,
                   hclust.method = hclust.method,
                   maxcluster = maxcluster)

    class(object) <- "bclust"

    optSEM <- getOption("show.error.messages")
    if(is.null(optSEM)) optSEM <- TRUE
    on.exit(options(show.error.messages = optSEM))

    if (verbose) cat("Committee Member:")
    for (n in 1:iter.base) {
      if (verbose){
        cat(" ", n, sep = "")
      }
      if(resample){
        x1 <- x[sample(xr, replace = TRUE, prob=weights), ]
      }
      else{
        x1 <- x
      }

      for(m in 1:20){
        if(verbose) cat("(",m,")",sep="")
        options(show.error.messages = FALSE)
        tryres <- try(CLUSFUN(x1, centers = base.centers, ...))
        if(!inherits(tryres, "try-error")) break
      }
      options(show.error.messages = optSEM)
      if(m==20)
        stop("Could not find valid cluster solution in 20 replications\n")

      object$allcenters[((n - 1) * base.centers + 1):(n * base.centers),] <-
        tryres$centers
    }
    object$allcenters <- object$allcenters[complete.cases(object$allcenters),]
    object$allcluster <- knn1(object$allcenters, x,
                              factor(1:nrow(object$allcenters)))

    if(minsize > 0){
      object <- prune.bclust(object, x, minsize=minsize)
    }

    if (verbose)
      cat("\nComputing Hierarchical Clustering\n")
    object <- hclust.bclust(object, x = x, centers = centers,
                            final.kmeans = final.kmeans,
                            docmdscale=docmdscale)
    object
  }

"centers.bclust" <- function (object, k)
{
  centers <- matrix(0, nrow = k, ncol = ncol(object$allcenters))
  for (m in 1:k) {
    centers[m, ] <-
      apply(object$allcenters[object$members[,k-1] == m, , drop = FALSE], 2, mean)
  }
  centers
}

"clusters.bclust" <- function (object, k, x=NULL)
{
  if(missing(x))
    allcluster <- object$allcluster
  else
    allcluster <- knn1(object$allcenters, x,
                       factor(1:nrow(object$allcenters)))

  return(object$members[allcluster, k - 1])
}


"hclust.bclust" <-
  function (object, x, centers, dist.method = object$dist.method,
            hclust.method = object$hclust.method, final.kmeans = FALSE,
            docmdscale = FALSE, maxcluster=object$maxcluster)
  {
    if (!require("mva")) stop("Could not load required package mva")

    d <- dist(object$allcenters, method = dist.method)
    if(hclust.method=="diana"){
      if (!require("cluster")) stop("Could not load required package cluster")
      object$hclust <- as.hclust(diana(d, diss=TRUE))
    }
    else
      object$hclust <- hclust(d, method = hclust.method)

    if(docmdscale){
      object$cmdscale <- cmdscale(d)
    }

    object$members <- cutree(object$hclust, 2:maxcluster)
    object$cluster <- clusters.bclust(object, centers)
    object$centers <- centers.bclust(object, centers)
    if (final.kmeans) {
      kmeansres <- kmeans(x, centers = object$centers)
      object$centers <- kmeansres$centers
      object$cluster <- kmeansres$cluster
    }
    object
  }



"plot.bclust" <-
  function (x, maxcluster=x$maxcluster,
            main = deparse(substitute(x)), ...)
  {
    if (!require("mva")) stop("Could not load required package mva")

    opar <- par(c("mar", "oma"))
    on.exit(par(opar))

    par(oma = c(0, 0, 3, 0))
    layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE))
    par(mar = c(0, 4, 4, 1))
    plot(x$hclust, labels = FALSE, hang = -1)
    x1 <- 1:maxcluster
    x2 <- 2:maxcluster
    y <- rev(x$hclust$height)[x1]
    z <- abs(diff(y))
    par(mar = c(4, 4, 1, 2))
    plot(x1, ((y - min(y))/(max(y) - min(y))), ty = "l", xlab = "",
         ylab = "", ylim = c(0, 1))
    lines(x2, z/sum(z), col = "grey")
    text(x2, z/sum(z), labels = as.character(x2))
    #    lx2 <- length(x2)
    #    abline(h=qexp(.95, rate = length(x2)), lty=3, col="grey")
    #    abline(h=qexp(.95^(1/lx2), rate = length(x2)), lty=3, col="grey")
    mtext(main, outer = TRUE, cex = 1.5)
    layout(1)
  }

"boxplot.bclust" <-
  function (x, n = nrow(x$centers), bycluster = TRUE,
            main = deparse(substitute(x)), oneplot=TRUE,
            which=1:n, ...)
  {
    if (!require("mva")) stop("Could not load required package mva")
    N <- length(which)

    opar <- par(c("mfrow", "oma", "mgp","xpd"))
    on.exit(par(opar))
    par(xpd=NA)

    memb <- x$members[, (n - 1)]
    tmemb <- table(memb)
    cendf <- as.data.frame(x$allcenters)
    ylim <- range(x$allcenters)

    if (bycluster) {
      if(oneplot){
        if (N <= 3) {
          par(mfrow = c(N, 1))
        }
        else {
          par(mfrow = c(ceiling(N/2), 2))
        }
      }
      tcluster <- table(clusters.bclust(x, n))
      for (k in which) {
        boxplot(cendf[memb == k, ], col = "grey",
                names = rep("",ncol(cendf)),
                ylim = ylim, ...)
        if (x$datamean) {
          lines(x$datamean, col = "red")
        }
        if(!is.null(x$colnames)){
          text(1:length(x$colnames)+0.2,
               par("usr")[3],
               adj=1,srt=35,
               paste(x$colnames, "  "))
        }

        title(main = paste("Cluster ", k, ": ", tmemb[k],
                           " centers, ", tcluster[k], " data points", sep = ""))
      }
    }
    else {
      a <- ceiling(sqrt(ncol(cendf)))
      if(oneplot){
        par(mfrow = c(a, ceiling(ncol(cendf)/a)))
      }
      memb <- as.factor(memb)
      for (k in 1:ncol(cendf)) {
        boxplot(cendf[, k] ~ memb, col = "grey", ylim = ylim, ...)
        title(main = x$colnames[k])
        abline(h = x$datamean[k], col = "red")
      }
    }
  }

### prune centers that contain not at least minsize data points
prune.bclust <- function(object, x, minsize=1, dohclust=FALSE, ...){

  ok <- FALSE
  while(!all(ok)){
    object$allcluster <- knn1(object$allcenters, x,
                              factor(1:nrow(object$allcenters)))

    ok <- table(object$allcluster) >= minsize
    object$allcenters <- object$allcenters[ok, ]
  }
  if(dohclust){
    object <- hclust.bclust(object, x, nrow(object$centers), ...)
  }
  object
}

bincombinations <- function(p) {

  retval <- matrix(0, nrow=2^p, ncol=p)

  for(n in 1:p){
    retval[,n] <- rep(c(rep(0, (2^p/2^n)), rep(1, (2^p/2^n))),
                      length=2^p)
  }
  retval
}



cmeans <- function (x, centers, iter.max = 100, verbose = FALSE,
                    dist = "euclidean", method = "cmeans",
                    m=2, rate.par = NULL)
{
  xrows <- dim(x)[1]
  xcols <- dim(x)[2]
  xold <- x
  perm <- sample(xrows)
  x <- x[perm, ]
  ## initial values are given
  if (is.matrix(centers))
    ncenters <- dim(centers)[1]
  else {
    ## take centers random vectors as initial values
    ncenters <- centers
    centers <- x[rank(runif(xrows))[1:ncenters], ]+0.001
  }

  dist <- pmatch(dist, c("euclidean", "manhattan"))
  if (is.na(dist))
    stop("invalid distance")
  if (dist == -1)
    stop("ambiguous distance")

  method <- pmatch(method, c("cmeans", "ufcl"))
  if (is.na(method))
    stop("invalid clustering method")
  if (method == -1)
    stop("ambiguous clustering method")

  if (method == 2) {
    if (missing(rate.par)) {
      rate.par <- 0.3
    }
  }

  initcenters <- centers
  ## dist <- matrix(0, xrows, ncenters)
  ## necessary for empty clusters
  pos <- as.factor(1:ncenters)
  rownames(centers) <- pos
  iter <- integer(1)

  if (method == 1) {
    retval <- .C("cmeans",
                 xrows = as.integer(xrows),
                 xcols = as.integer(xcols),
                 x = as.double(x),
                 ncenters = as.integer(ncenters),
                 centers = as.double(centers),
                 iter.max = as.integer(iter.max),
                 iter = as.integer(iter),
                 verbose = as.integer(verbose),
                 dist = as.integer(dist-1),
                 U = double(xrows*ncenters),
                 m = as.double(m),
                 ermin = double(1),
                 PACKAGE = "e1071")
  }
  else if (method == 2) {
    retval <- .C("ufcl",
                 xrows = as.integer(xrows),
                 xcols = as.integer(xcols),
                 x = as.double(x),
                 ncenters = as.integer(ncenters),
                 centers = as.double(centers),
                 iter.max = as.integer(iter.max),
                 iter = as.integer(iter),
                 verbose = as.integer(verbose),
                 dist = as.integer(dist-1),
                 U = double(xrows*ncenters),
                 m = as.double(m),
                 rate.par = as.double(rate.par),
                 ermin = double(1),
                 PACKAGE = "e1071")
  }

  centers <- matrix(retval$centers, ncol = xcols, dimnames = dimnames(initcenters))

  U <- retval$U
  U <- matrix(U, ncol=ncenters)
  #  clusterU <- max.col(U)
  clusterU <- apply(U,1,which.max)
  clusterU <- clusterU[order(perm)]
  U <- U[order(perm),]

  clustersize <- as.integer(table(clusterU))

  retval <- list(centers = centers, size = clustersize, cluster = clusterU,
                 iter= retval$iter - 1 , membership=U,
                 withinerror = retval$ermin, call = match.call())

  class(retval) <- c("fclust")
  return(retval)
}



#predict.fsegmentation <- function(clobj, x){

#  xrows<-dim(x)[1]
#  xcols<-dim(x)[2]
#  ncenters <- clobj$ncenters
#  cluster <- integer(xrows)
#  clustersize <- integer(ncenters)
#  f <- clobj$m


#  if(dim(clobj$centers)[2] != xcols){
#    stop("Number of variables in cluster object and x are not the same!")
#  }


#  retval <- .C("fuzzy_assign",
#               xrows = as.integer(xrows),
#               xcols = as.integer(xcols),
#               x = as.double(x),
#               ncenters = as.integer(ncenters),
#               centers = as.double(clobj$centers),
#               dist = as.integer(clobj$dist-1),
#              U = double(xrows*ncenters),
#               f = as.double(f))



#  U <- retval$U
#  U <- matrix(U, ncol=ncenters)
#  clusterU <- apply(U,1,which.max)
#  clustersize <- as.integer(table(clusterU))


#   clobj$iter <- NULL
#  clobj$cluster <- clusterU
#  clobj$size <- retval$clustersize
#  clobj$membership <- U

#  return(clobj)
#}
countpattern <- function(x, matching=FALSE)
{
  nvar <- dim(x)[2]
  n <- dim(x)[1]

  ## build matrix of all possible binary vectors
  b <- matrix(0, 2^nvar, nvar)
  for (i in 1:nvar)
    b[, nvar+1-i] <- rep(rep(c(0,1),c(2^(i-1),2^(i-1))),2^(nvar-i))

  namespat <- b[,1]
  for (i in 2:nvar)
    namespat <- paste(namespat, b[,i], sep="")

  xpat <- x[,1]
  for (i in 2:nvar)
    xpat <- 2*xpat+x[,i]
  xpat <- xpat+1

  pat <- tabulate(xpat, nb=2^nvar)
  names(pat) <- namespat

  if (matching)
    return(list(pat=pat, matching=xpat))
  else
    return(pat)
}
cshell <- function (x, centers, iter.max = 100, verbose = FALSE,
                    dist = "euclidean", method = "cshell",
                    m=2, radius= NULL)
{
  xrows <- dim(x)[1]
  xcols <- dim(x)[2]
  xold <- x
  perm <- sample(xrows)
  x <- x[perm, ]
  ## initial values are given
  if (is.matrix(centers))
    ncenters <- dim(centers)[1]
  else {
    ## take centers random vectors as initial values
    ncenters <- centers
    centers <- x[rank(runif(xrows))[1:ncenters], ]+0.001
  }

  ##initialize radius
  if (missing(radius))
    radius <- rep(0.2,ncenters)
  else
    radius <- as.double(radius)

  dist <- pmatch(dist, c("euclidean", "manhattan"))
  if (is.na(dist))
    stop("invalid distance")
  if (dist == -1)
    stop("ambiguous distance")

  method <- pmatch(method, c("cshell"))
  if (is.na(method))
    stop("invalid clustering method")
  if (method == -1)
    stop("ambiguous clustering method")

  initcenters <- centers
  ## dist <- matrix(0, xrows, ncenters)
  ## necessary for empty clusters
  pos <- as.factor(1:ncenters)
  rownames(centers) <- pos
  iter <- integer(1)

  flag <- integer(1)


  retval <- .C("cshell",
               xrows = as.integer(xrows),
               xcols = as.integer(xcols),
               x = as.double(x),
               ncenters = as.integer(ncenters),
               centers = as.double(centers),
               iter.max = as.integer(iter.max),
               iter = as.integer(iter),
               verbose = as.integer(verbose),
               dist = as.integer(dist-1),
               U = double(xrows*ncenters),
               UANT = double(xrows*ncenters),
               m = as.double(m),
               ermin = double(1),
               radius = as.double(radius),
               flag = as.integer(flag),
               PACKAGE = "e1071")

  centers <- matrix(retval$centers, ncol = xcols, dimnames = dimnames(initcenters))


  radius <- as.double(retval$radius)
  U <- retval$U
  U <- matrix(U, ncol=ncenters)
  UANT <- retval$UANT
  UANT <- matrix(UANT, ncol=ncenters)


  iter <- retval$iter
  flag <- as.integer(retval$flag)


  ##Optimization part
  while (((flag == 1) || (flag==4)) && (iter<=iter.max)){

    flag <- 3


    system <- function (spar=c(centers,radius), x, U, m, i){
      k <- dim(x)[1]
      d <- dim(x)[2]
      nparam<-length(spar)

      v<-spar[1:(nparam-1)]
      r<-spar[nparam]

      ##distance matrix x_k - v_i
      distmat <- t(t(x)-v)

      ##norm from x_k - v_i
      normdist <- distmat[,1]^2
      for (j in 2:d)
        normdist<-normdist+distmat[,j]^2
      normdist <- sqrt(normdist)

      ##equation 5
      op <- sum( (U[,i]^m) * (normdist-r) )^2
      ##equation 4
      equationmatrix <- ((U[,i]^m) * (1-r/normdist))*distmat
      op<- op+apply(equationmatrix, 2, sum)^2

    }


    for (i in 1:ncenters){
      spar <- c(centers[i,],radius[i])
      npar <- length(spar)

      optimres <- optim(spar ,system, method="CG", x=x, U=U, m=m, i=i)
      centers[i,] <- optimres$par[1:(npar-1)]
      radius[i] <- optimres$par[npar]
    }


    retval <- .C("cshell",
                 xrows = as.integer(xrows),
                 xcols = as.integer(xcols),
                 x = as.double(x),
                 ncenters = as.integer(ncenters),
                 centers = as.double(centers),
                 iter.max = as.integer(iter.max),
                 iter = as.integer(iter-1),
                 verbose = as.integer(verbose),
                 dist = as.integer(dist-1),
                 U = as.double(U),
                 UANT = as.double(UANT),
                 m = as.double(m),
                 ermin = double(1),
                 radius = as.double(radius),
                 flag = as.integer(flag),
                 PACKAGE = "e1071")

    flag<-retval$flag
    if (retval$flag!=2)
      flag<-1


    centers <- matrix(retval$centers, ncol = xcols, dimnames = dimnames(initcenters))

    radius <- as.double(retval$radius)
    U <- retval$U
    U <- matrix(U, ncol=ncenters)
    UANT <- retval$UANT
    UANT <- matrix(UANT, ncol=ncenters)



    iter <- retval$iter

  }

  centers <- matrix(retval$centers, ncol = xcols, dimnames = dimnames(initcenters))

  U <- retval$U
  U <- matrix(U, ncol=ncenters)


  clusterU <- apply(U,1,which.max)
  clusterU <- clusterU[order(perm)]
  U <- U[order(perm),]

  clustersize <- as.integer(table(clusterU))
  radius <- as.double(retval$radius)

  retval <- list(centers = centers, radius=radius,
                 size = clustersize, cluster = clusterU,
                 iter = retval$iter - 1, membership=U,
                 withinerror = retval$ermin,
                 call = match.call())

  class(retval) <- c("cshell", "fclust")
  return(retval)
}


#predict.cshell <- function( clobj, x){

#  xrows<-dim(x)[1]
#  xcols<-dim(x)[2]
#  ncenters <- clobj$ncenters
#  cluster <- integer(xrows)
#  clustersize <- integer(ncenters)
#  f <- clobj$m
#  radius <- clobj$radius

#  if(dim(clobj$centers)[2] != xcols){
#    stop("Number of variables in cluster object and x are not the same!")
#  }


#  retval <- .C("cshell_assign",
#               xrows = as.integer(xrows),
#               xcols = as.integer(xcols),
#               x = as.double(x),
#               ncenters = as.integer(ncenters),
#               centers = as.double(clobj$centers),
#               dist = as.integer(clobj$dist-1),
#               U = double(xrows*ncenters),
#               f = as.double(f),
#               radius = as.double(radius))



#  U <- retval$U
#  U <- matrix(U, ncol=ncenters)

#  clusterU <- apply(U,1,which.max)
#  clustersize <- as.integer(table(clusterU))


#  clobj$iter <- NULL
#  clobj$cluster <- clusterU
#  clobj$size <- retval$clustersize
#  clobj$membership <- U

#  return(clobj)
#}










rdiscrete <- function (n, probs, values = 1:length(probs), ...)
{
  sample(values, size=n, replace=TRUE, prob=probs)
}



ddiscrete <- function (x, probs, values = 1:length(probs))
{

  if (length(probs) != length(values))
    stop("ddiscrete: probs and values must have the same length.")
  if (sum(probs < 0) > 0)
    stop("ddiscrete: probs must not contain negative values.")
  if (!is.array(x) && !is.vector(x) && !is.factor(x))
    stop("ddiscrete: x must be an array or a vector or a factor.")

  p <- probs/sum(probs)

  y <- as.vector(x)
  l <- length(y)
  z <- rep(0,l)

  for (i in 1:l)
    if (any(values == y[i]))
      z[i] <- p[values == y[i]]

  z <- as.numeric(z)
  if (is.array(x))
    dim(z) <- dim(x)

  return(z)
}


pdiscrete <- function (q, probs, values = 1:length(probs))
{

  if (length(probs) != length(values))
    stop("pdiscrete: probs and values must have the same length.")
  if (sum(probs < 0) > 0)
    stop("pdiscrete: probs must not contain negative values.")
  if (!is.array(q) & !is.vector(q))
    stop("pdiscrete: q must be an array or a vector")

  p <- probs/sum(probs)

  y <- as.vector(q)
  l <- length(y)
  z <- rep(0,l)

  for (i in 1:l)
    z[i] <- sum(p[values <= y[i]])

  z <- as.numeric(z)
  if (is.array(q))
    dim(z) <- dim(q)

  return(z)
}

qdiscrete <- function (p, probs, values = 1:length(probs))
{

  if (length(probs) != length(values))
    stop("qdiscrete: probs and values must have the same length.")
  if (sum(probs < 0) > 0)
    stop("qdiscrete: probs must not contain negative values.")
  if (!is.array(p) & !is.vector(p))
    stop("qdiscrete: p must be an array or a vector")

  probs <- cumsum(probs)/sum(probs)

  y <- as.vector(p)
  l <- length(y)
  z <- rep(0,l)

  for (i in 1:l)
    z[i] <- length(values) - sum(y[i] <= probs) + 1

  z <- as.numeric(z)
  z <- values[z]
  if (is.array(p))
    dim(z) <- dim(p)

  return(z)
}



element <- function(x, i) {

  if(!is.array(x))
    stop("x is not an array")

  ni <- length(i)
  dx <- dim(x)

  if(length(i)!=length(dx))
    stop("Wrong number of subscripts")

  if(ni==1){
    return(x[i])
  }
  else{
    m1 <- c(i[1], i[2:ni]-1)
    m2 <- c(1,cumprod(dx)[1:(ni-1)])
    return(x[sum(m1*m2)])
  }
}

fclustIndex <- function ( y, x, index= "all" )
{

  clres <- y
  ###########################################################################
  ################SESSION 1: MEASURES#########################################
  ###########################################################################

  gath.geva <- function (clres,x)#for m=2
  {
    xrows <- dim(clres$me)[1]
    xcols <- dim(clres$ce)[2]
    ncenters <- dim(clres$centers)[1]
    scatter <- array(0.0, c(xcols, xcols, ncenters))
    scatternew <- array(0.0, c(xcols, xcols, ncenters))
    fhv <-as.double(0)
    apd <-as.double(0)
    pd <-  as.double(0)
    control <- as.double(0)

    for (i in 1:ncenters){
      paronomastis <- as.double(0)
      paronomastis2 <- as.double(0)

      for (j in 1:xrows){
        paronomastis <- paronomastis+clres$me[j,i]

        diff <- x[j,]-clres$ce[i,]
        scatternew[,,i] <- clres$me[j,i]*(t(t(diff))%*%t(diff))
        scatter[,,i] <- scatter[,,i]+scatternew[,,i]
      }#xrows
      scatter[,,i] <- scatter[,,i]/paronomastis


      for (j in 1:xrows){
        diff <- x[j,]-clres$ce[i,]
        control <- (t(diff)%*%solve(scatter[,,i]))%*%t(t(diff))
        if (control<1.0)
          paronomastis2 <- paronomastis2+clres$me[j,i]
        ##   else
        ##     cat("...")
      }#xrows
      fhv <- fhv+sqrt(det(scatter[,,i]))

      apd <- apd+paronomastis2/sqrt(det(scatter[,,i]))

      pd <- pd+paronomastis2

    }#ncenters
    pd <- pd/fhv
    apd <- apd/ncenters

    retval <- list(fuzzy.hypervolume=fhv,average.partition.density=apd, partition.density=pd)
    return(retval)
  }

  xie.beni <- function(clres){#for all m
    xrows <- dim(clres$me)[1]
    minimum<--1
    error <- clres$within#sd
    ncenters <- dim(clres$centers)[1]
    for (i in 1:(ncenters-1)){
      for (j in (i+1):ncenters){
        diff<- clres$ce[i,]-clres$ce[j,]
        diffdist <- t(diff)%*%t(t(diff))
        if (minimum==-1)
          minimum <- diffdist
        if (diffdist<minimum)
          minimum <- diffdist
      }}
    xiebeni <- error/(xrows*minimum)
    return(xiebeni)
  }

  fukuyama.sugeno <- function(clres){#for all m
    xrows <- dim(clres$me)[1]
    ncenters <-dim(clres$centers)[1]
    error <- clres$within#sd
    k2<-as.double(0)

    meancenters <- apply(clres$ce,2,mean)
    for (i in 1:ncenters){
      paronomastis3 <- as.double(0)
      for (j in 1:xrows){
        paronomastis3 <- paronomastis3+(clres$me[j,i]^2)}

      diff <- clres$ce[i,]-meancenters
      diffdist <- t(diff)%*%t(t(diff))
      k2 <- k2 +paronomastis3*diffdist
    }#ncenters

    fukuyamasugeno<-error-k2
    return(fukuyamasugeno)
  }

  partition.coefficient <- function(clres){
    xrows <- dim(clres$me)[1]

    partitioncoefficient <- sum(apply(clres$me^2,1,sum))/xrows
    return(partitioncoefficient)
  }

  partition.entropy <- function(clres){
    xrows <- dim(clres$me)[1]
    ncenters <- dim(clres$centers)[1]
    partitionentropy <- 0.0
    for (i in 1:xrows){
      for (k in 1:ncenters){
        if (clres$me[i,k]!=0.0)
          partitionentropy<- partitionentropy+(clres$me[i,k]*log(clres$me[i,k]))
      }}
    partitionentropy<-partitionentropy/((-1)*xrows)
    ##partitionentropy <- sum(apply(clres$me*log(clres$me),1,sum))/((-1)*xrows)
    return(partitionentropy)
  }

  separation.index <- function(clres, x)
  {
    xrows <- dim(clres$me)[1]
    xcols <- dim(x)[2]
    ncenters <- dim(clres$centers)[1]
    maxcluster <- double(ncenters)
    minimum <- -1.0
    ##hardpartition <- matrix(0,xrows,ncenters)

    ##for (i in 1:xrows)
    ## hardpartition[i,clres$cl[i]] <- 1
    for (i in 1:ncenters){
      maxcluster[i] <- max(dist(matrix(x[clres$cl==i],ncol=xcols)))
    }
    maxdia <- maxcluster[rev(order(maxcluster))[1]]
    for (i in 1:(ncenters-1)){
      for (j in (i+1):(ncenters)){
        for (m in 1:xrows){
          if (clres$cl[m]==i){
            for (l in 1:xrows){
              if (clres$cl[l]==j){
                diff <- x[m,]-x[l,]
                diffdist <- sqrt(t(diff)%*%t(t(diff)))
                fraction <- diffdist/maxdia
                if (minimum==-1)
                  minimum <- fraction
                if (fraction<minimum)
                  minimum <- fraction
              }}}}}}
    return(minimum)
  }


  proportion.exponent <- function(clres)
  {
    k <- dim(clres$centers)[2]
    xrows <- dim(clres$me)[1]

    bexp <- as.integer(1)
    for (j in 1:xrows){
      greatint <- as.integer(1/max(clres$me[j,]))
      aexp <- as.integer(0)
      for (l in 1:greatint){
        aexp <- aexp + (-1)^(l+1)*(gamma(k+1)/(gamma(l+1)*gamma(k-l+1)))*(1-l*max(clres$me[j,]))^(k-1)
        ##if (aexp==0.0){
        ##print("aexp=0.0")
        ##print(j)
        ##}
      }
      bexp <- bexp * aexp

    }
    proportionexponent <- -log(bexp)
    return(proportionexponent)
  }

  ###########################################################################
  ################SESSION 2: MAIN PROGRAM####################################
  ###########################################################################

  index <- pmatch(index, c("gath.geva", "xie.beni", "fukuyama.sugeno",
                           "partition.coefficient", "partition.entropy",
                           "proportion.exponent", "separation.index", "all"))


  if (is.na(index))
    stop("invalid clustering index")
  if (index == -1)
    stop("ambiguous index")

  vecallindex <- numeric(9)

  if (any(index==1) || (index==8)){
    gd <- gath.geva(clres, x)
    vecallindex[1] <- gd$fuzzy
    vecallindex[2] <- gd$average
    vecallindex[3] <- gd$partition}
  if (any(index==2) || (index==8))
    vecallindex[4] <- xie.beni(clres)
  if (any(index==3) || (index==8))
    vecallindex[5] <- fukuyama.sugeno(clres)
  if (any(index==4) || (index==8))
    vecallindex[6] <- partition.coefficient(clres)
  if (any(index==5) || (index==8))
    vecallindex[7] <- partition.entropy(clres)
  if (any(index==6) || (index==8))
    vecallindex[8] <- proportion.exponent(clres)
  if (any(index==7) || (index==8)){
    require(mva)
    vecallindex[9] <- separation.index(clres,x)}

  names(vecallindex) <- c("fhv", "apd", "pd", "xb", "fs", "pc", "pe",
                          "pre", "si")


  if (index < 8)
    vecallindex <- vecallindex[index]

  return(vecallindex)
}


hamming.distance <- function(x,y){

  z<-NULL

  if(is.vector(x) && is.vector(y)){
    z <- sum(as.logical(x) != as.logical(y))
  }
  else{
    z <- matrix(0,nrow=nrow(x),ncol=nrow(x))
    for(k in 1:(nrow(x)-1)){
      for(l in (k+1):nrow(x)){
        z[k,l] <- hamming.distance(x[k,], x[l,])
        z[l,k] <- z[k,l]
      }
    }
    dimnames(z) <- list(dimnames(x)[[1]], dimnames(x)[[1]])
  }
  z
}



















hamming.window <- function (n)
{
  if (n == 1)
    c <- 1
  else
  {
    n <- n-1
    c <- 0.54 - 0.46*cos(2*pi*(0:n)/n)
  }
  return(c)
}
hanning.window <- function (n)
{
  if (n == 1)
    c <- 1
  else
  {
    n <- n-1
    c <- 0.5 - 0.5*cos(2*pi*(0:n)/n)
  }
  return(c)
}
ica <- function(X, lrate, epochs=100, ncomp=dim(X)[2],
                fun="negative")
{
  if (!is.matrix(X))
  {
    if (is.data.frame(X))
      X <- as.matrix(X)
    else
      stop("ica: X must be a matrix or a data frame")
  }
  if (!is.numeric(X))
    stop("ica: X contains non numeric elements")

  m <- dim(X)[1]
  n <- dim(X)[2]

  Winit <- matrix(rnorm(n*ncomp), ncomp, n)
  W <- Winit

  if (!is.function(fun))
  {
    funlist <- c("negative kurtosis", "positive kurtosis",
                 "4th moment")
    p <- pmatch(fun, funlist)
    if (is.na(p))
      stop("ica: invalid fun")
    funname <- funlist[p]
    if (p == 1) fun <- tanh
    else if (p == 2) fun <- function(x) {x - tanh(x)}
    else if (p == 3) fun <- function(x) {sign(x)*x^2}
  }
  else funname <- as.character(substitute(fun))

  for (i in 1:epochs)
    for (j in 1:m)
    {
      x <- X[j,, drop=FALSE]
      y <- W%*%t(x)
      gy <- fun(y)
      W <- W + lrate*gy%*%(x-t(gy)%*%W)
    }
  colnames(W) <- NULL
  pr <- X%*%t(W)
  retval <- list(weights = W, projection = pr, epochs = epochs,
                 fun = funname, lrate = lrate, initweights = Winit)
  class(retval) <- "ica"
  return(retval)
}


print.ica <- function(x, ...)
{
  cat(x$epochs, "Trainingssteps with a learning rate of", x$lrate, "\n")
  cat("Function used:", x$fun,"\n\n")
  cat("Weightmatrix\n")
  print(x$weights, ...)
}

plot.ica <- function(x, ...) pairs(x$pr, ...)

impute <- function(x, what=c("median", "mean")){

  what <- match.arg(what)

  if(what == "median"){
    retval <-
      apply(x, 2,
            function(z) {z[is.na(z)] <- median(z, na.rm=TRUE); z})
  }
  else if(what == "mean"){
    retval <-
      apply(x, 2,
            function(z) {z[is.na(z)] <- mean(z, na.rm=TRUE); z})
  }
  else{
    stop("`what' invalid")
  }
  retval
}
interpolate <- function(x, a, adims=lapply(dimnames(a), as.numeric),
                        method="linear"){

  if(is.vector(x)) x<- matrix(x, ncol=length(x))

  if(!is.array(a))
    stop("a is not an array")

  ad <- length(dim(a))

  method <- pmatch(method, c("linear", "constant"))
  if (is.na(method))
    stop("invalid interpolation method")

  if(any(unlist(lapply(adims, diff))<0))
    stop("dimensions of a not ordered")

  retval <- rep(0, nrow(x))
  bincombi <- bincombinations(ad)

  convexcoeff <- function(x, y) {
    ok <- y>0
    x[ok] <- y[ok]-x[ok]
    x
  }

  for(n in 1:nrow(x)){

    ## the `leftmost' corner of the enclosing hypercube
    leftidx <- rep(0, ad)
    xabstand <- rep(0, ad)
    aabstand <- rep(0, ad)

    for(k in 1:ad){
      if(x[n,k] < min(adims[[k]]) || x[n,k] > max(adims[[k]]))
        stop("No extrapolation allowed")
      else{
        leftidx[k] <- max(seq(adims[[k]])[adims[[k]] <= x[n,k]])
        ## if at the right border, go one step to the left
        if(leftidx[k] == length(adims[[k]]))
          leftidx[k] <- leftidx[k] - 1

        xabstand[k] <- x[n,k] - adims[[k]][leftidx[k]]
        aabstand[k] <- adims[[k]][leftidx[k]+1] -
          adims[[k]][leftidx[k]]
      }
    }

    coefs <- list()
    if(method==1){
      for(k in 1:(2^ad)){
        retval[n] <- retval[n] +
          element(a, leftidx+bincombi[k,]) *
          prod((aabstand-
                  convexcoeff(xabstand,
                              aabstand*bincombi[k,]))/aabstand)
      }
    }
    else if(method==2){
      retval[n] <- element(a, leftidx)
    }
  }

  names(retval) <- rownames(x)
  retval
}

kurtosis <- function (x, na.rm = FALSE)
{
  if (na.rm)
    x <- x[!is.na(x)]
  sum((x-mean(x))^4)/(length(x)*var(x)^2) - 3
}


lca <- function(x, k, niter=100, matchdata=FALSE, verbose=FALSE)
{

  ## if x is a data matrix -> create patterns
  if (is.matrix(x))
  {
    if (matchdata)
    {
      x <- countpattern(x, matching=TRUE)
      xmat <- x$matching
      x <- x$pat
    }
    else
      x <- countpattern(x, matching=FALSE)
  }
  else   ## if no data ist given, matchdata must be FALSE
    matchdata <- FALSE

  n <- sum(x)
  npat <- length(x)
  nvar <- round(log(npat)/log(2))

  ## build matrix of all possible binary vectors
  b <- matrix(0, 2^nvar, nvar)
  for (i in 1:nvar)
    b[, nvar+1-i] <- rep(rep(c(0,1),c(2^(i-1),2^(i-1))),2^(nvar-i))

  ## initialize probabilities
  classprob <- runif(k)
  classprob <- classprob/sum(classprob)
  names(classprob) <- 1:k
  p <- matrix(runif(nvar*k), k)


  pas <- matrix(0, k, npat)
  classsize <- numeric(k)

  for (i in 1:niter)
  {
    for (j in 1:k)
    {
      ## P(pattern|class)
      mp <- t(b)*p[j,]+(1-t(b))*(1-p[j,])
      pas[j,] <- drop(exp(rep(1,nvar)%*%log(mp))) # column product
    }
    ##  P(pattern|class)*P(class)
    pas <- t(t(pas)*classprob)

    ## P(class|pattern)
    sump <- drop(rep(1,k)%*%pas)  # column sums
    pas <- t(t(pas)/sump)

    spas <- t(t(pas)*x)
    classsize <- drop(spas%*%rep(1,npat))  # row sums
    classprob <- classsize/n
    p <- pas%*%(x*b)/classsize
    if (verbose)
      cat("Iteration:", i, "\n")
  }

  for (j in 1:k)
  {
    mp <- t(b)*p[j,]+(1-t(b))*(1-p[j,])
    pas[j,] <- drop(exp(rep(1,nvar)%*%log(mp)))*classprob[j]
    # column product
  }

  ## LogLikelihood
  pmust <-  drop(rep(1,k)%*%pas)  # column sums
  ll <- sum(x*log(pmust))

  ## Likelihoodquotient
  xg0 <- x[x>0]
  ll0 <- sum(xg0*log(xg0/n))
  lq <- 2*(ll0-ll)

  ## bic
  bic <- -2*ll+log(n)*(k*(nvar+1)-1)
  bicsat <- -2*ll0+log(n)*(2^nvar-1)

  ## chisq
  ch <- sum((x-n*pmust)^2/(n*pmust))

  ## P(class|pattern)
  sump <- drop(rep(1,k)%*%pas)  # column sums
  pas <- t(t(pas)/sump)

  mat <- max.col(t(pas))
  if (matchdata)
    mat <- mat[xmat]

  colnames(p) <- 1:nvar
  rownames(p) <- 1:k

  lcaresult <- list(classprob=classprob, p=p, matching=mat,
                    logl=ll, loglsat=ll0,
                    chisq=ch, lhquot=lq, bic=bic, bicsat=bicsat, n=n,
                    np=(k*(nvar+1)-1), matchdata=matchdata)

  class(lcaresult) <- "lca"
  return(lcaresult)
}


print.lca <- function(x, ...)
{
  cat("LCA-Result\n")
  cat("----------\n\n")

  cat("Datapoints:", x$n, "\n")
  cat("Classes:   ", length(x$classprob), "\n")
  cat("Probability of classes\n")
  print(round(x$classprob,3))

  cat("Itemprobabilities\n")
  print(round(x$p,2))
}

summary.lca <- function(object, ...)
{
  nvar <- ncol(object$p)
  object$npsat <- 2^nvar-1
  object$df <- 2^nvar-1-object$np
  object$pvallhquot <- 1-pchisq(object$lhquot,object$df)
  object$pvalchisq <- 1-pchisq(object$chisq,object$df)
  object$k <- length(object$classprob)

  ## remove unnecessary list elements
  object$classprob <- NULL
  object$p <- NULL
  object$matching <- NULL

  class(object) <- "summary.lca"
  return(object)
}


print.summary.lca <- function(x, ...)
{
  cat("LCA-Result\n")
  cat("----------\n\n")

  cat("Datapoints:", x$n, "\n")
  cat("Classes:   ", x$k, "\n")

  cat("\nGoodness of fit statistics:\n\n")
  cat("Number of parameters, estimated model:", x$np, "\n")
  cat("Number of parameters, saturated model:", x$npsat, "\n")
  cat("Log-Likelihood, estimated model:      ", x$logl, "\n")
  cat("Log-Likelihood, saturated model:      ", x$loglsat, "\n")

  cat("\nInformation Criteria:\n\n")
  cat("BIC, estimated model:", x$bic, "\n")
  cat("BIC, saturated model:", x$bicsat, "\n")

  cat("\nTestStatistics:\n\n")
  cat("Likelihood ratio:  ", x$lhquot,
      "  p-val:", x$pvallhquot, "\n")
  cat("Pearson Chi^2:     ", x$chisq,
      "  p-val:", x$pvalchisq, "\n")
  cat("Degress of freedom:", x$df, "\n")
}


bootstrap.lca <- function(l, nsamples=10, lcaiter=30, verbose=FALSE)
{

  n <- l$n
  classprob <- l$classprob
  nclass <- length(l$classprob)
  prob <- l$p
  nvar <- ncol(l$p)
  npat <- 2^nvar

  ## build matrix of all possible binary vectors
  b <- matrix(0, npat, nvar)
  for (i in 1:nvar)
    b[, nvar+1-i] <- rep(rep(c(0,1),c(2^(i-1),2^(i-1))),2^(nvar-i))

  ll <- lq <- ll0 <- ch <- numeric(nsamples)

  for (i in 1:nsamples)
  {
    ## generate data
    cm <- sample(1:nclass, size=n, replace=TRUE, prob=classprob)
    x <- matrix(runif(n*nvar), nrow=n)
    x <- (x<prob[cm,])+0
    x <- countpattern(x, matching=FALSE)

    if (verbose)
      cat ("Start of Bootstrap Sample", i, "\n")
    lc <- lca(x, nclass, niter=lcaiter, verbose=verbose)

    ll[i] <- lc$logl
    ll0[i] <- lc$loglsat
    lq[i] <- lc$lhquot
    ch[i] <- lc$chisq

    if (verbose)
      cat("LL:", ll[i], " LLsat:", ll0[i],
          " Ratio:", lq[i], " Chisq:", ch[i], "\n")
  }


  lratm <- mean(lq)
  lrats <- sd(lq)
  chisqm <- mean(ch)
  chisqs <- sd(ch)

  zl <- (l$lhquot-lratm)/lrats
  zc <- (l$ch-chisqm)/chisqs
  pzl <- 1-pnorm(zl)
  pzc <- 1-pnorm(zc)
  pl <- sum(l$lhquot<lq)/length(lq)
  pc <- sum(l$ch<ch)/length(ch)

  lcaboot <- list(logl=ll, loglsat=ll0, lratio=lq,
                  lratiomean=lratm, lratiosd=lrats,
                  lratioorg=l$lhquot, zratio=zl,
                  pvalzratio=pzl, pvalratio=pl,
                  chisq=ch, chisqmean=chisqm, chisqsd=chisqs,
                  chisqorg=l$ch, zchisq=zc,
                  pvalzchisq=pzc, pvalchisq=pc,
                  nsamples=nsamples, lcaiter=lcaiter)
  class(lcaboot) <- "bootstrap.lca"
  return(lcaboot)
}

print.bootstrap.lca <- function(x, ...)
{
  cat("Bootstrap of LCA\n")
  cat("----------------\n\n")

  cat ("Number of Bootstrap Samples:    ", x$nsamples, "\n")
  cat ("Number of LCA Iterations/Sample:", x$lcaiter, "\n")

  cat("Likelihood Ratio\n\n")
  cat("Mean:", x$lratiomean, "\n")
  cat("SDev:", x$lratiosd, "\n")
  cat("Value in Data Set:", x$lratioorg, "\n")
  cat("Z-Statistics:     ", x$zratio, "\n")
  cat("P(Z>X):           ", x$pvalzratio, "\n")
  cat("P-Val:            ", x$pvalratio, "\n\n")

  cat("Pearson's Chisquare\n\n")
  cat("Mean:", x$chisqmean, "\n")
  cat("SDev:", x$chisqsd, "\n")
  cat("Value in Data Set:", x$chisqorg, "\n")
  cat("Z-Statistics:     ", x$zchisq, "\n")
  cat("P(Z>X):           ", x$pvalzchisq, "\n")
  cat("P-Val:            ", x$pvalchisq, "\n\n")
}

predict.lca <- function(object, x, ...)
{
  if (object$matchdata)
    stop("predict.lca: only possible, if lca has been called with matchdata=FALSE")
  else
  {
    x <- countpattern(x, matching=TRUE)
    return(object$matching[x$matching])
  }
}
classAgreement <-
  function (tab, match.names=FALSE)
  {
    n <- sum(tab)
    ni <- apply(tab, 1, sum)
    nj <- apply(tab, 2, sum)

    ## patch for matching factors
    if (match.names && !is.null(dimnames(tab))) {
      lev <- intersect (colnames (tab), rownames(tab))
      p0 <- sum(diag(tab[lev,lev]))/n
      pc <- sum(ni[lev] * nj[lev])/n^2
    } else { # cutoff larger dimension
      m <- min(length(ni), length(nj))
      p0 <- sum(diag(tab[1:m, 1:m]))/n
      pc <- sum(ni[1:m] * nj[1:m])/n^2
    }
    n2 <- choose(n, 2)
    rand <- 1 + (sum(tab^2) - (sum(ni^2) + sum(nj^2))/2)/n2
    nis2 <- sum(choose(ni[ni > 1], 2))
    njs2 <- sum(choose(nj[nj > 1], 2))
    crand <- (sum(choose(tab[tab > 1], 2)) - (nis2 * njs2)/n2)/((nis2 +
                                                                   njs2)/2 - (nis2 * njs2)/n2)
    list(diag = p0, kappa = (p0 - pc)/(1 - pc), rand = rand,
         crand = crand)
  }

matchClasses <- function(tab, method = "rowmax", iter=1, maxexact=9,
                         verbose=TRUE){

  methods <- c("rowmax", "greedy", "exact")
  method <- pmatch(method, methods)

  rmax <- apply(tab,1,which.max)
  myseq <- 1:ncol(tab)
  cn <- colnames(tab)
  rn <- rownames(tab)
  if(is.null(cn)){
    cn <- myseq
  }
  if(is.null(rn)){
    rn <- myseq
  }

  if(method==1){
    retval <- rmax
  }
  if(method==2 | method==3){
    if(ncol(tab)!=nrow(tab)){
      stop("Unique matching only for square tables.")
    }


    dimnames(tab) <- list(myseq, myseq)
    cmax <- apply(tab,2,which.max)
    retval <- rep(NA, ncol(tab))
    names(retval) <- colnames(tab)

    baseok <- cmax[rmax]==myseq
    for(k in myseq[baseok]){
      therow <- (tab[k,])[-rmax[k]]
      thecol <- (tab[, rmax[k]])[-k]
      if(max(outer(therow, thecol, "+")) < tab[k, rmax[k]]){
        retval[k] <- rmax[k]
      }
      else{
        baseok[k] <- FALSE
      }
    }



    if(verbose){
      cat("Direct agreement:", sum(baseok),
          "of", ncol(tab), "pairs\n")
    }

    if(!all(baseok)){
      if(method==3){
        if(sum(!baseok)>maxexact){
          method <- 2
          warning(paste("Would need permutation of", sum(!baseok),
                        "numbers, resetting to greedy search\n"))
        }
        else{
          iter <- gamma(ncol(tab)-sum(baseok)+1)
          if(verbose){
            cat("Iterations for permutation matching:", iter, "\n")
          }
          perm <- permutations(ncol(tab)-sum(baseok))
        }
      }

      ## rest for permute matching
      if(any(baseok)){
        rest <- myseq[-retval[baseok]]
      }
      else{
        rest <- myseq
      }

      for(l in 1:iter){
        newretval <- retval
        if(method == 2){
          ok <- baseok
          while(sum(!ok)>1){
            rest <- myseq[!ok]
            k <- sample(rest, 1)
            if(any(ok)){
              rmax <- tab[k, -newretval[ok]]
            }
            else{
              rmax <- tab[k,]
            }
            newretval[k] <- as.numeric(names(rmax)[which.max(rmax)])
            ok[k] <- TRUE
          }
          newretval[!ok] <- myseq[-newretval[ok]]
        }
        else{
          newretval[!baseok] <- rest[perm[l,]]
        }

        if(l>1){
          agree <- sum(diag(tab[,newretval]))/sum(tab)
          if(agree>oldagree){
            retval <- newretval
            oldagree <- agree
          }
        }
        else{
          retval <- newretval
          agree <- oldagree <- sum(diag(tab[,newretval]))/sum(tab)
        }
      }
    }
  }

  if(verbose){
    cat("Cases in matched pairs:",
        round(100*sum(diag(tab[,retval]))/sum(tab), 2), "%\n")
  }

  if(any(as.character(myseq)!=cn)){
    retval <- cn[retval]
  }
  names(retval) <- rn

  retval
}

compareMatchedClasses <- function(x, y,
                                  method="rowmax", iter=1, maxexact=9,
                                  verbose=FALSE)
{
  if(missing(y)){
    retval <- list(diag=matrix(NA, nrow=ncol(x), ncol=ncol(x)),
                   kappa=matrix(NA, nrow=ncol(x), ncol=ncol(x)),
                   rand=matrix(NA, nrow=ncol(x), ncol=ncol(x)),
                   crand=matrix(NA, nrow=ncol(x), ncol=ncol(x)))
    for(k in 1:(ncol(x)-1)){
      for(l in (k+1):ncol(x)){
        tab <- table(x[,k], x[,l])
        m <- matchClasses(tab, method=method, iter=iter,
                          verbose=verbose, maxexact=maxexact)
        a <- classAgreement(tab[,m])
        retval$diag[k,l] <- a$diag
        retval$kappa[k,l] <- a$kappa
        retval$rand[k,l] <- a$rand
        retval$crand[k,l] <- a$crand
      }
    }
  }
  else{
    x <- as.matrix(x)
    y <- as.matrix(y)
    retval <- list(diag=matrix(NA, nrow=ncol(x), ncol=ncol(y)),
                   kappa=matrix(NA, nrow=ncol(x), ncol=ncol(y)),
                   rand=matrix(NA, nrow=ncol(x), ncol=ncol(y)),
                   crand=matrix(NA, nrow=ncol(x), ncol=ncol(y)))
    for(k in 1:ncol(x)){
      for(l in 1:ncol(y)){
        tab <- table(x[,k], y[,l])
        m <- matchClasses(tab, method=method, iter=iter,
                          verbose=verbose, maxexact=maxexact)
        a <- classAgreement(tab[,m])
        retval$diag[k,l] <- a$diag
        retval$kappa[k,l] <- a$kappa
        retval$rand[k,l] <- a$rand
        retval$crand[k,l] <- a$crand
      }
    }
  }

  retval
}


permutations <- function(n) {

  if(n ==1)
    return(matrix(1))
  else if(n<2)
    stop("n must be a positive integer")

  z <- matrix(1)
  for (i in 2:n) {
    x <- cbind(z, i)
    a <- c(1:i, 1:(i - 1))
    z <- matrix(0, ncol=ncol(x), nrow=i*nrow(x))
    z[1:nrow(x),] <- x
    for (j in 2:i-1) {
      z[j*nrow(x)+1:nrow(x),] <- x[, a[1:i+j]]
    }
  }
  dimnames(z) <- NULL
  z
}
moment <- function(x, order = 1, center = FALSE, absolute = FALSE,
                   na.rm = FALSE) {
  if (na.rm)
    x <- x[!is.na(x)]
  if (center)
    x <- x - mean(x)
  if (absolute)
    x <- abs(x)
  sum(x ^ order) / length(x)
}
plot.stft <- function (x, col = gray (63:0/63), ...)
{
  x <- x$values
  image(x=1:dim(x)[1], y=1:dim(x)[2], z=x, col=col, ...)
}
print.fclust <- function (x, ...)
{
  clobj <- x
  if (!is.null(clobj$iter))
    cat("\n                            Clustering on Training Set\n\n\n")
  else
    cat("\n                              Clustering on Test Set\n\n\n")

  cat("Number  of  Clusters: ", dim(clobj$centers)[1], "\n")
  cat("Sizes   of  Clusters: ", clobj$size, "\n\n")
  #    cat("Centers of  Clusters: ", "\n")
  #matrix(clobj$centers, ncol=dim(clobj$ce)[2]), "\n\n")
  #    clobj$centers

  #    cat("Learning Parameters:",clobj$call,"\n\n")

  #    if (clobj$iter < clobj$call$iter.max)
  #      cat("Algorithm converged after", clobj$iter, "iterations.\n")
  #    else
  #      cat("Algorithm did not converge after", clobj$iter, "iterations.\n")
}


rbridge <- function(end=1, frequency=1000) {

  z <- rwiener(end=end, frequency=frequency)
  ts(z - time(z)*as.vector(z)[frequency],
     start=1/frequency, frequency=frequency)
}

read.octave <- function (file, quiet=FALSE) {

  nr <- 0
  nc <- 0

  if(!quiet)
    cat("Header: ")

  head <- scan(file=file,what=character(),nlines=4, sep=":", quiet=quiet)
  if(length(head) != 8){
    stop("Header seem to be corrupt")
  }
  for(k in 1:4){
    if(head[2*k-1] == "# rows"){
      nr <- as.integer(head[2*k])
    }
    else if(head[2*k-1] == "# columns"){
      nc <- as.integer(head[2*k])
    }
  }

  if(!quiet)
    cat("Data  : ")

  z <- scan(file=file,skip=4,quiet=quiet)
  if(length(z) != nc*nr){
    stop("Wrong number of data elements")
  }

  if((nr>1) && (nc>1)){
    if(!quiet)
      cat(paste("Matrix:", nr, "rows,", nc, "columns\n"))

    z<-matrix(z, nrow=nr, ncol=nc, byrow=TRUE)
  }
  else if(!quiet){
    cat("Vector:", nr*nc, "elements\n")
  }
  z
}

rectangle.window <- function (n)
  rep (1, n)
rwiener <- function(end=1, frequency=1000) {

  z<-cumsum(rnorm(end*frequency)/sqrt(frequency))
  ts(z, start=1/frequency, frequency=frequency)
}

scaclust <- function (x, centers, iter.max = 100, verbose = FALSE,
                      method = "ad", theta= NULL)
{
  xrows <- dim(x)[1]
  xcols <- dim(x)[2]
  xold <- x
  perm <- sample(xrows)
  x <- x[perm, ]
  ## initial values are given
  if (is.matrix(centers))
    ncenters <- dim(centers)[1]
  else {
    ## take centers random vectors as initial values
    ncenters <- centers
    centers <- x[rank(runif(xrows))[1:ncenters], ]+0.001
  }
  ##  method <- pmatch(method, c("ad", "mtv", "sand","ml"))
  method <- pmatch(method, c("ad", "mtv", "sand","mlm"))
  if (is.na(method))
    stop("invalid clustering method")
  if (method == -1)
    stop("ambiguous clustering method")

  if (method == 1) {
    beta <- 1/xcols
    taf <- 0 }
  if (method == 2) {
    beta <- 0.5
    taf <- xcols/2 }
  if (method == 3) {
    beta <- 1/xcols
    taf <- 1 }
  if (method == 4){
    beta <- 0.0
    taf <- -1 }


  ##initialize theta
  ## if (method != 4){
  if (missing(theta))
    theta <- rep(1.0,ncenters)
  else
    theta <- as.double(theta)

  ##}


  initcenters <- centers
  ## dist <- matrix(0, xrows, ncenters)
  ## necessary for empty clusters
  pos <- as.factor(1:ncenters)
  rownames(centers) <- pos
  iter <- integer(1)

  ##  if ((method == 1) || (method == 2) || (method == 3)){
  retval <- .C("common",
               xrows = as.integer(xrows),
               xcols = as.integer(xcols),
               x = as.double(x),
               ncenters = as.integer(ncenters),
               centers = as.double(centers),
               itermax = as.integer(iter.max),
               iter = as.integer(iter),
               verbose = as.integer(verbose),
               U = double(xrows*ncenters),
               beta = as.double(beta),
               taf = as.double(taf),
               theta = as.double(theta),
               ermin = double(1),
               PACKAGE = "e1071")

  centers <- matrix(retval$centers, ncol = xcols, dimnames = dimnames(initcenters))

  U <- retval$U
  U <- matrix(U, ncol=ncenters)
  ##  clusterU <- max.col(U)
  clusterU <- apply(U,1,which.max)
  clusterU <- clusterU[order(perm)]
  U <- U[order(perm),]

  clustersize <- as.integer(table(clusterU))

  retval <- list(centers = centers, size = clustersize,
                 cluster = clusterU, iter = retval$iter - 1, membership=U,
                 withinerror = retval$ermin, call = match.call())

  class(retval) <- c("fclust")
  return(retval)
}



allShortestPaths <- function(x){

  x <- as.matrix(x)
  x[is.na(x)] <- .Machine$double.xmax
  x[is.infinite(x) & x>0] <- .Machine$double.xmax
  if(ncol(x) != nrow(x))
    stop("x is not a square matrix")
  n <- ncol(x)

  z <- .C("e1071_floyd",
          as.integer(n),
          double(n^2),
          as.double(x),
          integer(n^2),
          PACKAGE = "e1071")
  z <- list(length = matrix(z[[2]], n),
            middlePoints = matrix(z[[4]]+1, n))
  z$length[z$length == .Machine$double.xmax] <- NA
  z
}

extractPath <- function(obj, start, end){

  z <- integer(0)
  path <- function(i, j){

    k <- obj$middlePoints[i, j]
    if (k != 0)
    {
      path(i,k);
      z <<- c(z, k)
      path(k,j);
    }
  }
  path(start,end)
  c(start, z, end)
}

sigmoid <- function(x) 1/(1 + exp(-x))

dsigmoid <- function(x) sigmoid(x) * (1 - sigmoid(x))

d2sigmoid <- function(x) dsigmoid(x) * (1 - 2 * sigmoid(x))
skewness <- function (x, na.rm = FALSE)
{
  if (na.rm)
    x <- x[!is.na(x)]
  sum((x-mean(x))^3)/(length(x)*sd(x)^3)
}

read.matrix.csr <- function (file, fac = TRUE, ncol = NULL) {
  library(methods)
  if(!require(SparseM))
    stop("Need `SparseM' package!")

  con <- file(file)
  open(con)
  y <- vector()
  x <- new("matrix.csr")
  x@ia <- 1
  i <- 1
  maxcol <- 1
  while (isOpen(con) & length(buf <- readLines (con, 1)) > 0) {
    s <- strsplit(buf, "[ ]+",extended=TRUE)[[1]]

    ## y
    if (length(grep(":", s[1])) == 0) {
      y[i] <- if (fac) s[1] else as.numeric(s[1])
      s <- s[-1]
    }

    ## x-values
    if (length(s)) {
      tmp <- do.call("rbind", strsplit(s, ":"))
      x@ra <- c(x@ra, as.numeric(tmp[,2]))
      x@ja <- c(x@ja, as.numeric(tmp[,1]))
    }
    i <- i + 1
    x@ia[i] <- x@ia[i-1] + length(s)
  }
  x@dimension <- c(i - 1, if (is.null(ncol)) max(x@ja) else max(ncol, x@ja))
  class(x) <- "matrix.csr"

  if (length(y))
    list (x = x, y = if (fac) as.factor(y) else y)
  else x
}

write.matrix.csr <- function (x, file="out.dat", y=NULL) {
  library(methods)
  if (!is.null(y) & (length(y) != nrow(x)))
    stop(paste("Length of y (=", length(y),
               ") does not match number of rows of x (=",
               nrow(x), ")!", sep=""))
  sink(file)
  for (i in 1:nrow(x)) {
    if (!is.null(y)) cat (y[i],"")
    for (j in x@ia[i]:(x@ia[i+1] - 1))
      cat(x@ja[j], ":", x@ra[j], " ", sep="")
    cat("\n")
  }
  sink()
}

na.fail.matrix.csr <- function(object, ...) {
  library(methods)
  if (any(is.na(object@ra)))
    stop("missing values in object") else return(object)
}






stft <- function(X, win=min(80,floor(length(X)/10)),
                 inc=min(24, floor(length(X)/30)), coef=64,
                 wtype="hanning.window")
{
  numcoef <- 2*coef
  if (win > numcoef)
  {
    win <- numcoef
    cat ("stft: window size adjusted to", win, ".\n")
  }
  numwin <- trunc ((length(X) - win) / inc)

  ## compute the windows coefficients
  wincoef <- eval(parse(text=wtype))(win)

  ## create a matrix Z whose columns contain the windowed time-slices
  z <- matrix (0, numwin + 1, numcoef)
  y <- z
  st <- 1
  for (i in 0:numwin)
  {
    z[i+1, 1:win] <- X[st:(st+win-1)] * wincoef
    y[i+1,] <- fft(z[i+1,])
    st <- st + inc
  }

  Y<- list (values = Mod(y[,1:coef]), windowsize=win, increment=inc,
            windowtype=wtype)
  class(Y) <- "stft"
  return(Y)
}



svm <- function (x, ...)
  UseMethod ("svm")

svm.formula <-
  function (formula, data=NULL, subset, na.action=na.fail, ...)
  {
    call <- match.call()
    if (!inherits(formula, "formula"))
      stop("method is only for formula objects")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 0
    x <- model.matrix(Terms, m)
    y <- model.extract(m, response)
    ret <- svm.default (x, y, ...)
    ret$call <- call
    ret$terms <- Terms
    if (!is.null(attr(m, "na.action")))
      ret$na.action <- attr(m, "na.action")
    class (ret) <- c("svm.formula", class(ret))
    return (ret)
  }

svm.default <-
  function (x,
            y         = NULL,
            scale     = TRUE,
            subset,
            na.action = na.fail,
            type      = NULL,
            kernel    = "radial",
            degree    = 3,
            gamma     = 1 / ncol(as.matrix(x)),
            coef0     = 0,
            cost      = 1,
            nu        = 0.5,
            class.weights = NULL,
            cachesize = 40,
            tolerance = 0.001,
            epsilon   = 0.1,
            shrinking = TRUE,
            cross     = 0,
            fitted    = TRUE,
            ...)
  {
    sparse  <- inherits(x, "matrix.csr")
    if (sparse) {
      if (!require(SparseM))
        stop("Need SparseM package for handling of sparse structures!")
    }

    if(is.null(degree)) stop("`degree' must not be NULL!")
    if(is.null(gamma)) stop("`gamma' must not be NULL!")
    if(is.null(coef0)) stop("`coef0' must not be NULL!")
    if(is.null(cost)) stop("`cost' must not be NULL!")
    if(is.null(nu)) stop("`nu' must not be NULL!")
    if(is.null(epsilon)) stop("`epsilon' must not be NULL!")
    if(is.null(tolerance)) stop("`tolerance' must not be NULL!")

    xhold   <- if (fitted) x else NA
    x.scale <- y.scale <- NULL

    if (sparse) {
      scale <- FALSE
      na.fail(x)
      if(!is.null(y)) na.fail(y)
      x <- t(t(x)) ## make shure that col-indices are sorted
    } else {
      formula <- inherits(x, "svm.formula")
      x <- if (is.vector(x)) t(t(x)) else as.matrix(x)
      if (!formula) {
        if (!missing(subset)) x <- x[subset,]
        x <- if (is.null(y))
          na.action(x)
        else {
          df <- data.frame(y, x)
          df <- na.action(df)
          y <- df[,1]
          as.matrix(df[,-1])
        }
      }

      if (scale) {
        co <- if (is.data.frame(x)) !sapply(x, var) else !apply(x, 2, var)
        if (any(co)) {
          scale <- FALSE
          warning(paste("Variable(s)",
                        paste("`",colnames(x)[co], "'", sep="", collapse=" and "),
                        "constant. Cannot scale data.")
          )
        } else {
          x <- scale(x)
          x.scale <- attributes(x)[c("scaled:center","scaled:scale")]
          if (is.numeric(y)) {
            y <- scale(y)
            y.scale <- attributes(y)[c("scaled:center","scaled:scale")]
            y <- as.vector(y)
          }
        }
      }
    }

    if (is.null(type)) type <-
      if (is.null(y)) "one-classification"
    else if (is.factor(y)) "C-classification"
    else "eps-regression"

    type <- pmatch (type, c("C-classification",
                            "nu-classification",
                            "one-classification",
                            "eps-regression",
                            "nu-regression"), 99) - 1

    if (type > 10) stop("wrong type specification!")

    kernel <- pmatch (kernel, c("linear",
                                "polynomial",
                                "radial",
                                "sigmoid"), 99) - 1

    if (kernel > 10) stop("wrong kernel specification!")

    if (!is.vector(y) && !is.factor (y) && type != 2) stop("y must be a vector or a factor.")
    if (type != 2 && length(y) != nrow(x)) stop("x and y don't match.")

    if (cachesize < 0.1) cachesize <- 0.1

    lev <- NULL
    weightlabels <- NULL
    # in case of classification: transform factors into integers
    if (type == 2) # one class classification --> set dummy
      y <- 1
    else
      if (is.factor(y)) {
        lev <- levels(y)
        y <- as.integer(y)
        if (!is.null(class.weights)) {
          if (is.null(names(class.weights)))
            stop("Weights have to be specified along with their according level names !")
          weightlabels <- match (names(class.weights),lev)
          if (any(is.na(weightlabels)))
            stop("At least one level name is missing or misspelled.")
        }
      } else {
        if (type < 3 && any(as.integer(y) != y))
          stop("dependent variable has to be of factor or integer type for classification mode.")
        lev <- unique(y)
      }

    nclass <- 2
    if (type < 2) nclass <- length(lev)

    cret <- .C ("svmtrain",
                # data
                as.double  (if (sparse) x$ra else t(x)),
                as.integer (nrow(x)), as.integer(ncol(x)),
                as.double  (y),

                # sparse index info
                as.integer (if (sparse) x$ia else 0),
                as.integer (if (sparse) x$ja else 0),

                # parameters
                as.integer (type),
                as.integer (kernel),
                as.double  (degree),
                as.double  (gamma),
                as.double  (coef0),
                as.double  (cost),
                as.double  (nu),
                as.integer (weightlabels),
                as.double  (class.weights),
                as.integer (length (class.weights)),
                as.double  (cachesize),
                as.double  (tolerance),
                as.double  (epsilon),
                as.integer (shrinking),
                as.integer (cross),
                as.integer (sparse),

                # results
                nclasses = integer  (1),
                nr       = integer  (1), # nr of support vectors
                index    = integer  (nrow(x)),
                labels   = integer  (nclass),
                nSV      = integer  (nrow(x)),
                rho      = double   (nclass * (nclass - 1) / 2),
                coefs    = double   (nrow(x) * (nclass - 1)),

                cresults = double   (cross),
                ctotal1  = double   (1),
                ctotal2  = double   (1),
                error    = character(1),

                PACKAGE = "e1071")

    if (nchar(cret$error))
      stop(paste(cret$error, "!", sep=""))

    ret <- list (
      call     = match.call(),
      type     = type,
      kernel   = kernel,
      cost     = cost,
      degree   = degree,
      gamma    = gamma,
      coef0    = coef0,
      nu       = nu,
      epsilon  = epsilon,
      sparse   = sparse,
      scaled   = scale,
      x.scale  = x.scale,
      y.scale  = y.scale,

      nclasses = cret$nclasses,            #number of classes
      levels   = lev,
      tot.nSV  = cret$nr,                  #total number of sv
      nSV      = cret$nSV[1:cret$nclasses],#number of SV in diff. classes
      labels   = cret$label[1:cret$nclasses],#labels of the SVs.
      SV       = t(t(x[cret$index[1:cret$nr],])), #copy of SV
      index    = cret$index[1:cret$nr],     #indexes of sv in x
      #constants in decision functions
      rho      = cret$rho[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
      #coefficiants of sv
      coefs    = if (cret$nr == 0) NULL else
        t(matrix(cret$coefs[1:((cret$nclasses - 1) * cret$nr)],
                 nrow = cret$nclasses - 1,
                 byrow = TRUE))
    )

    # cross-validation-results
    if (cross > 0)
      if (type > 2) {
        scale.factor     <- if (scale) crossprod(y.scale$"scaled:scale") else 1;
        ret$MSE          <- cret$cresults * scale.factor;
        ret$tot.MSE      <- cret$ctotal1  * scale.factor;
        ret$scorrcoeff   <- cret$ctotal2;
      } else {
        ret$accuracies   <- cret$cresults;
        ret$tot.accuracy <- cret$ctotal1;
      }

    class (ret) <- "svm"
    ret$fitted  <- if (fitted) predict(ret, xhold) else NA
    ret
  }

predict.svm <- function (object, newdata, ...) {
  if (missing(newdata))
    return(fitted(object))

  sparse <- inherits(newdata, "matrix.csr")
  if (object$sparse || sparse) {
    if (!require(SparseM))
      stop("Need SparseM package for handling of sparse structures!")
  }

  if (is.vector(newdata) || sparse) newdata <- t(t(newdata))
  if (!object$sparse) {
    if (inherits(object, "svm.formula")) {
      if(is.null(colnames(newdata)))
        colnames(newdata) <- colnames(object$SV)
      newdata <- model.matrix(delete.response(terms(object)),
                              as.data.frame(newdata), na.action = na.omit)
    } else if (!sparse) newdata <- as.matrix(newdata)
  }

  if (object$scaled)
    newdata <- scale(newdata,
                     center = object$x.scale$"scaled:center",
                     scale  = object$x.scale$"scaled:scale"
    )

  if (ncol(object$SV) != ncol(newdata)) stop ("test data does not match model !")

  ret <- .C ("svmpredict",
             # model
             as.double  (if (object$sparse) object$SV$ra else t(object$SV)),
             as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
             as.integer (if (object$sparse) object$SV$ia else 0),
             as.integer (if (object$sparse) object$SV$ja else 0),
             as.double  (as.vector(object$coefs)),
             as.double  (object$rho),
             as.integer (object$nclasses),
             as.integer (object$tot.nSV),
             as.integer (object$labels),
             as.integer (object$nSV),
             as.integer (object$sparse),

             # parameter
             as.integer (object$type),
             as.integer (object$kernel),
             as.double  (object$degree),
             as.double  (object$gamma),
             as.double  (object$coef0),

             # test matrix
             as.double  (if (sparse) newdata$ra else t(newdata)),
             as.integer (nrow(newdata)),
             as.integer (if (sparse) newdata$ia else 0),
             as.integer (if (sparse) newdata$ja else 0),
             as.integer (sparse),

             # decision-values
             ret = double(nrow(newdata)),

             PACKAGE = "e1071"
  )$ret

  if (is.character(object$levels))
    # classification: return factors
    factor (object$levels[ret], levels=object$levels)
  else if (object$type == 2)
    # one-class-classification: return TRUE/FALSE
    ret == 1
  else if (object$scaled)
    # return raw values, possibly scaled back
    ret * object$y.scale$"scaled:scale" + object$y.scale$"scaled:center"
  else
    ret
}

print.svm <- function (x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n")
  cat("Parameters:\n")
  cat("   SVM-Type: ", c("C-classification",
                         "nu-classification",
                         "one-classification",
                         "eps-regression",
                         "nu-regression")[x$type+1], "\n")
  cat(" SVM-Kernel: ", c("linear",
                         "polynomial",
                         "radial",
                         "sigmoid")[x$kernel+1], "\n")
  if (x$type==0 || x$type==3 || x$type==4)
    cat("       cost: ", x$cost, "\n")
  if (x$kernel==1)
    cat("     degree: ", x$degree, "\n")
  cat("      gamma: ", x$gamma, "\n")
  if (x$kernel==1 || x$kernel==3)
    cat("     coef.0: ", x$coef0, "\n")
  if (x$type==1 || x$type==2 || x$type==4)
    cat("         nu: ", x$nu, "\n")
  if (x$type==3)
    cat("    epsilon: ", x$epsilon, "\n\n")

  cat("\nNumber of Support Vectors: ", x$tot.nSV)
  cat("\n\n")

}

summary.svm <- function(object, ...)
  structure(object, class="summary.svm")

print.summary.svm <- function (x, ...) {
  print.svm(x)
  if (x$type<2) {
    cat(" (", x$nSV, ")\n\n")
    cat("\nNumber of Classes: ", x$nclasses, "\n\n")
    cat("Levels:", if(is.numeric(x$levels)) "(as integer)", "\n", x$levels)
  }
  cat("\n\n")
  if (x$type==2) cat("\nNumber of Classes: 1\n\n\n")

  if ("MSE" %in% names(x)) {
    cat(length (x$MSE), "-fold cross-validation on training data:\n\n", sep="")
    cat("Total Mean Squared Error:", x$tot.MSE, "\n")
    cat("Squared Correlation Coefficient:", x$scorrcoef, "\n")
    cat("Mean Squared Errors:\n", x$MSE, "\n\n")
  }
  if ("accuracies" %in% names(x)) {
    cat(length (x$accuracies), "-fold cross-validation on training data:\n\n", sep="")
    cat("Total Accuracy:", x$tot.accuracy, "\n")
    cat("Single Accuracies:\n", x$accuracies, "\n\n")
  }
  cat("\n\n")
}

scale.data.frame <- function(x, center = TRUE, scale = TRUE) {
  i <- sapply(x, is.numeric)
  if (ncol(x[, i, drop = FALSE])) {
    x[, i] <- tmp <- scale.default(x[, i, drop = FALSE], na.omit(center), na.omit(scale))
    if(center || !is.logical(center))
      attr(x, "scaled:center")[i] <- attr(tmp, "scaled:center")
    if(scale || !is.logical(scale))
      attr(x, "scaled:scale")[i]  <- attr(tmp, "scaled:scale")
  }
  x
}

plot.svm <- function(x, data, formula = NULL, fill = TRUE,
                     grid = 50, slice = list(), ...) {
  if (x$type < 3) {
    if (is.null(formula) && ncol(data) == 3) {
      formula <- formula(delete.response(terms(x)))
      formula[2:3] <- formula[[2]][2:3]
    }
    if(is.null(formula)) stop("missing formula.")
    if (fill) {
      sub <- model.frame(formula, data)
      xr <- seq(min(sub[,2]), max(sub[,2]), length = grid)
      yr <- seq(min(sub[,1]), max(sub[,1]), length = grid)
      l <- length(slice)
      if (l < ncol(data) - 3) {
        slnames <- names(slice)
        slice <- c(slice, rep(list(0), ncol(data) - 3 - l))
        names <- labels(delete.response(terms(x)))
        names(slice) <- c(slnames, names[! names %in% c(colnames(sub), slnames)])
      }
      lis <- c(list(yr), list(xr), slice)
      names(lis)[1:2] <- colnames(sub)
      new <- expand.grid(lis)[,labels(terms(x))]
      preds <- predict(x, new)
      filled.contour(xr, yr, matrix(codes(preds), nr = length(xr), byrow=TRUE),
                     plot.axes = {
                       axis(1)
                       axis(2)
                       colors <- codes(model.response(model.frame(x, data[-x$index,])))
                       points(formula, data = data[-x$index,], col = colors)
                       points(formula, data = data[x$index,], pch = "x", col = colors)
                     },
                     levels = 1:(length(unique(codes(preds)))+1),
                     key.axes = axis(4,
                                     1:length(unique(codes(preds)))+0.5,
                                     labels = levels(preds)[unique(preds)], las = 3
                     ),
                     plot.title = title(main = "SVM classification plot",
                                        xlab = names(lis)[1], ylab = names(lis)[2]),
                     ...
      )
    } else {
      plot(formula, data = data, type = "n", ...)
      colors <- codes(model.response(model.frame(m, data[-x$index,])))
      points(formula, data = data[-x$index,], col = colors)
      points(formula, data = data[x$index,], pch = "x", col = colors)
    }
  }
}
tune <- function(method, train.x, train.y = NULL, data = list(),
                 validation.x = NULL, validation.y = NULL,
                 ranges,
                 random = FALSE,
                 nrepeat = 1,
                 repeat.aggregate = min,
                 sampling = c("cross", "fix", "bootstrap"),
                 sampling.aggregate = mean,
                 cross = 10,
                 fix = 2 / 3,
                 nboot = 10,
                 boot.size = 9 / 10,
                 predict.func = predict,
                 best.model = TRUE,
                 performances = TRUE,
                 ...
) {
  ## internal helper functions
  resp <- function(formula, data)
    model.response(model.frame(formula, data))

  classAgreement <- function (tab) {
    n <- sum(tab)
    if (!is.null(dimnames(tab))) {
      lev <- intersect(colnames(tab), rownames(tab))
      p0 <- sum(diag(tab[lev, lev])) / n
    } else {
      m <- min(dim(tab))
      p0 <- sum(diag(tab[1:m, 1:m])) / n
    }
    p0
  }

  ## parameter handling
  sampling <- match.arg(sampling)
  method <- deparse(substitute(method))
  if (sampling == "cross") validation.x <- validation.y <- NULL
  useFormula <- is.null(train.y)
  if (useFormula && (is.null(data) || length(data) == 0))
    data <- model.frame(train.x)
  if (is.vector(train.x)) train.x <- t(t(train.x))

  ## prepare training indices
  if (!is.null(validation.x)) fix <- 1
  n <- nrow(if (useFormula) data else train.x)
  perm.ind <- sample(n)
  if (sampling == "cross") {
    if (cross > n)
      stop("`cross' must not exceed sampling size!")
    if (cross == 1)
      stop("`cross' must be greater than 1!")
  }
  train.ind <- if (sampling == "cross")
    tapply(1:n, cut(1:n, breaks = cross), function(x) perm.ind[-x])
  else if (sampling == "fix")
    list(perm.ind[1:trunc(n * fix)])
  else ## bootstrap
    lapply(1:nboot, function(x) sample(n, n * boot.size))

  ## find best model
  parameters <- expand.grid(ranges)
  p <- nrow(parameters)
  if (!is.logical(random)) {
    if (random < 1)
      stop ("random must be a strictly positive integer")
    if (random > p) random <- p
    parameters <- parameters[sample(1:p, random),]
  }
  model.errors <- c()

  ## - loop over all models
  for (para.set in 1:p) {
    sampling.errors <- c()

    ## - loop over all training samples
    for (sample in 1:length(train.ind)) {
      repeat.errors <- c()

      ## - repeat training `nrepeat' times
      for (reps in 1:nrepeat) {

        ## train one model
        model <- if (useFormula)
          do.call(method, c(list(train.x,
                                 data = data,
                                 subset = train.ind[[sample]]),
                            parameters[para.set,,drop = FALSE], ...
          )
          )
        else
          do.call(method, c(list(train.x[train.ind[[sample]],],
                                 y = train.y[train.ind[[sample]]]),
                            parameters[para.set,,drop = FALSE], ...
          )
          )

        ## predict validation set
        pred <- predict.func(model,
                             if (!is.null(validation.x))
                               validation.x
                             else if (useFormula)
                               data[-train.ind[[sample]],]
                             else
                               train.x[-train.ind[[sample]],]
        )

        ## compute performance measure
        true.y <- if (!is.null(validation.y))
          validation.y
        else if (useFormula)
          resp(train.x, data[-train.ind[[sample]],])
        else
          train.y[-train.ind[[sample]]]

        repeat.errors[reps] <- if (is.factor(true.y)) ## classification error
          1 - classAgreement(table(pred, true.y))
        else ## mean squared error
          crossprod(pred - true.y) / length(pred)
      }
      sampling.errors[sample] <- repeat.aggregate(repeat.errors)
    }
    model.errors[para.set] <- sampling.aggregate(sampling.errors)
  }

  ## return results
  best <- which.min(model.errors)
  structure(list(best.parameters  = parameters[best,,drop = FALSE],
                 best.performance = model.errors[best],
                 method           = method,
                 sampling         = switch(sampling,
                                           fix = "fixed training/validation set",
                                           bootstrap = "bootstrapping",
                                           cross = if (cross == n) "leave-one-out" else
                                             paste(cross,"-fold cross validation", sep="")
                 ),
                 performances     = if (performances) cbind(parameters, error = model.errors),
                 best.model       = if (best.model)
                   if (useFormula)
                     do.call(method, c(list(train.x, data = substitute(data)),
                                       parameters[best,,drop = FALSE], ...))
                 else
                   do.call(method, c(list(x = train.x,
                                          y = train.y),
                                     parameters[best,,drop = FALSE], ...))
  ),
  class = "tune"
  )
}

print.tune <- function(x, ...) {
  cat("\nParameter tuning of `", x$method, "':\n\n", sep="")
  cat("- sampling method:", x$sampling,"\n\n")
  cat("- best parameters:\n")
  print(x$best.parameters)
  cat("\n- best performance:", x$best.performance, "\n")
  cat("\n")
}

summary.tune <- function(object, ...)
  structure(object, class = "summary.tune")

print.summary.tune <- function(x, ...) {
  print.tune(x)
  if (!is.null(x$performances)) {
    cat("- Detailed performance results:\n")
    print(x$performances)
    cat("\n")
  }
}

plot.tune <- function(x,
                      type=c("contour","perspective"),
                      theta=60,
                      col="lightblue",
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      swapxy = FALSE,
                      transform.x = NULL,
                      transform.y = NULL,
                      transform.z = NULL,
                      color.palette = topo.colors,
                      nlevels = 20,
                      ...)
{
  if (is.null(x$performance))
    stop("Object does not contain detailed performance measures!")
  k <- ncol(x$performance)
  if (k > 3) stop("Cannot visualize more than 2 parameters")
  type = match.arg(type)

  if (is.null(main))
    main <- paste("Performance of `",x$method,"' on `",x$data,"'", sep="")

  if (k == 2)
    plot(x$performances, type = "b", main = main)
  else  {
    if (!is.null(transform.x))
      x$performances[,1] <- transform.x(x$performances[,1])
    if (!is.null(transform.y))
      x$performances[,2] <- transform.y(x$performances[,2])
    if (!is.null(transform.z))
      x$performances[,3] <- transform.z(x$performances[,3])
    if (swapxy)
      x$performances[,1:2] <- x$performances[,2:1]
    x <- xtabs(error~., data = x$performances)
    if (is.null(xlab)) xlab <- names(dimnames(x))[1 + swapxy]
    if (is.null(ylab)) ylab <- names(dimnames(x))[2 - swapxy]
    if (type == "perspective")
      persp(x=as.double(rownames(x)),
            y=as.double(colnames(x)),
            z=x,
            xlab=xlab,
            ylab=ylab,
            zlab="accuracy",
            theta=theta,
            col=col,
            ticktype="detailed",
            main = main,
            ...
      )
    else
      filled.contour(x=as.double(rownames(x)),
                     y=as.double(colnames(x)),
                     xlab=xlab,
                     ylab=ylab,
                     nlevels=nlevels,
                     color.palette = color.palette,
                     main = main,
                     x, ...)
  }
}

#############################################
## convenience functions for some methods
#############################################

tune.svm <- function(x, y = NULL, degree = NULL, gamma = NULL,
                     coef0 = NULL, cost = NULL, nu = NULL, ...) {
  ranges <- list(degree = degree, gamma = gamma,
                 coef0 = coef0, cost = cost, nu = nu)
  ranges[sapply(ranges, is.null)] <- NULL
  if (length(ranges) < 1)
    stop("No parameter range given.")
  tune(svm, train.x = x, train.y = y, ranges = ranges, ...)
}

best.svm <- function(x, ...)
  tune.svm(x, ..., best.model = TRUE)$best.model

tune.nnet <- function(x, y = NULL, size = NULL, decay = NULL, nrepeat = 5, trace = FALSE,
                      predict.func = function(...) predict(..., type="class"), ...) {
  require(nnet)
  ranges <- list(size = size, decay = decay)
  ranges[sapply(ranges, is.null)] <- NULL
  if (length(ranges) < 1)
    stop("No parameter range given.")
  tune(nnet, train.x = x, train.y = NULL, ranges = ranges, nrepeat = nrepeat, trace = trace,
       predict.func = predict.func, ...)
}

best.nnet <- function(x, ...)
  tune.nnet(x, ..., best.model = TRUE)$best.model

tune.randomForest <- function(x, y = NULL, nodesize = NULL, mtry = NULL, ntree = NULL, ...) {
  require(randomForest)
  ranges <- list(nodesize = nodesize, mtry = mtry, ntree = ntree)
  ranges[sapply(ranges, is.null)] <- NULL
  if (length(ranges) < 1)
    stop("No parameter range given.")
  tune(randomForest, train.x = x, ranges = ranges, ...)
}

best.randomForest <- function(x, ...)
  tune.randomForest(x, ..., best.model = TRUE)$best.model

knn.wrapper <- function(x, y, k = 1, l = 0, ...)
  list(train = x, cl = y, k = k, l = l, ...)

tune.knn <- function(x, y, k = NULL, l = NULL, ...) {
  require(class)
  ranges <- list(k = k, l = l)
  ranges[sapply(ranges, is.null)] <- NULL
  if (length(ranges) < 1)
    stop("No parameter range given.")
  tune(knn.wrapper,
       train.x = x, train.y = y, ranges = ranges,
       predict.func = function(x, ...) knn(train = x$train, cl = x$cl, k = x$k, l = x$l, ...),
       ...)
}

rpart.wrapper <- function(formula, minsplit=20, minbucket=round(minsplit/3), cp=0.01,
                          maxcompete=4, maxsurrogate=5, usesurrogate=2, xval=10,
                          surrogatestyle=0, maxdepth=30, ...)
  rpart(formula, control = rpart.control(minsplit=minsplit, minbucket=minbucket, cp=cp,
                                         maxcompete=maxcompete, maxsurrogate=maxsurrogate,
                                         usesurrogate=usesurrogate, xval=xval,
                                         surrogatestyle=surrogatestyle, maxdepth=maxdepth),
        ...
  )

tune.rpart <- function(formula, data, na.action = na.omit,
                       minsplit=NULL, minbucket=NULL, cp=NULL,
                       maxcompete=NULL, maxsurrogate=NULL, usesurrogate=NULL, xval=NULL,
                       surrogatestyle=NULL, maxdepth=NULL,
                       predict.func = NULL,
                       ...) {
  require(rpart)
  ranges <- list(minsplit=minsplit, minbucket=minbucket, cp=cp,
                 maxcompete=maxcompete, maxsurrogate=maxsurrogate,
                 usesurrogate=usesurrogate, xval=xval,
                 surrogatestyle=surrogatestyle, maxdepth=maxdepth)
  ranges[sapply(ranges, is.null)] <- NULL
  if (length(ranges) < 1)
    stop("No parameter range given.")

  predict.func <- if (is.factor(model.response(model.frame(formula, data))))
    function(...) predict(..., type="class")
  else
    predict
  tune(rpart.wrapper, train.x = formula, data = data, ranges = ranges,
       predict.func = predict.func, na.action = na.omit, ...)
}

best.rpart <- function(formula, ...)
  tune.rpart(formula, ..., best.model = TRUE)$best.model





.First.lib <- function(lib, pkg){

  library.dynam("e1071", pkg, lib)
}
