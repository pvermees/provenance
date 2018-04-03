#' Kolmogorov-Smirnov dissimilarity
#'
#' Returns the Kolmogorov-Smirnov dissimilarity between two samples
#'
#' @param x the first sample as a vector
#' @param y the second sample as a vector
#' @return a scalar value representing the maximum vertical distance
#' between the two cumulative distributions
#' @examples
#' data(Namib)
#' print(KS.diss(Namib$DZ$x[['N1']],Namib$DZ$x[['T8']]))
#' @export
KS.diss <- function(x,y) {
    xy <- c(x,y)
    cad1 <- stats::ecdf(x)
    cad2 <- stats::ecdf(y)
    max(abs(cad2(xy)-cad1(xy)))
}

#' Kuiper dissimilarity
#'
#' Returns the Kuiper dissimilarity between two samples
#'
#' @param x the first sample as a vector
#' @param y the second sample as a vector
#' @return a scalar value representing the sum of the maximum vertical
#' distances above and below the cumulative distributions of x and y
#' @examples
#' data(Namib)
#' print(Kuiper.diss(Namib$DZ$x[['N1']],Namib$DZ$x[['T8']]))
#' @export
Kuiper.diss <- function(x,y){
    xy <- c(x,y)
    cad1 <- stats::ecdf(x)
    cad2 <- stats::ecdf(y)
    d <- cad2(xy) - cad1(xy)
    M <- max(d)
    m <- min(d)
    if (M<0) M <- 0
    if (m>0) m <- 0
    M-m
}

#' Calculate the dissimilarity matrix between two \code{distributional} or
#' \code{compositional} datasets
#'
#' Calculate the dissimilarity matrix between two datasets of class
#' \code{distributional} or \code{compositional} using the Kolmogorov-Smirnov,
#' Sircombe-Hazelton, Aitchison or Bray Curtis distance
#' 
#' @param x an object of class \code{distributional} or
#'     \code{compositional}
#' @param method (optional) either "KS", "Kuiper", "SH", "aitchison"
#'     or "bray"
#' @examples
#' data(Namib)
#' print(round(100*diss(Namib$DZ)))
#' @return an object of class \code{diss}
#' @rdname diss
#' @export
diss <- function(x,method){ UseMethod("diss",x) }
#' @rdname diss
#' @export
diss.distributional <- function(x,method=NULL) {
    if (!is.null(method)) x$method <- method
    n <- length(x$x)
    d <- mat.or.vec(n,n)
    rownames(d) <- names(x$x)
    colnames(d) <- names(x$x)
    if (x$method=="SH") c2 <- getc2(x)
    for (i in 1:n){
        for (j in 1:n){
            if (x$method=="SH"){
                d[i,j] <- SH.diss(x,i,j,c.con=c2)
            } else if (x$method=="KS"){
                d[i,j] <- KS.diss(x$x[[i]],x$x[[j]])
            } else if (x$method=="Kuiper"){
                d[i,j] <- Kuiper.diss(x$x[[i]],x$x[[j]])
            }
        }
    }
    out <- stats::as.dist(d)
    class(out) <- append("diss",class(out))
    return(out)
}
#' @rdname diss
#' @export
diss.compositional <- function(x,method=NULL){
    if (!is.null(method)) x$method <- method
    if (x$method=="aitchison"){
        out <- stats::dist(CLR(x))
    } else {
        snames <- names(x)
        ns <- length(snames)
        d <- mat.or.vec(ns,ns)
        rownames(d) <- snames
        colnames(d) <- snames
        for (i in 1:ns){
            for (j in 1:ns){
                d[i,j] <- bray.diss(x$x[i,],x$x[j,])
            }
        }
        out <- stats::as.dist(d)
    }
    class(out) <- append("diss",class(out))
    return(out)
}
#' @rdname diss
#' @export
diss.counts <- function(x,method=NULL){
    if (!is.null(method)) x$method <- method
    snames <- names(x)
    ns <- length(snames)
    d <- mat.or.vec(ns,ns)
    X <- x$x / rowSums(x$x) # normalise data
    CC <- colSums(X)
    for (i in 1:ns){
        for (j in 1:ns){
            if (x$method=='bray'){
                d[i,j] <- bray.diss(x$x[i,],x$x[j,])
            } else { # chisq
                d[i,j] <- sqrt(sum( (1/CC)*(X[i,] - X[j,])^2 ))
            }
        }
    }
    rownames(d) <- snames
    colnames(d) <- snames
    out <- stats::as.dist(d)
    class(out) <- append("diss",class(out))
    return(out)
}

#' Bray-Curtis dissimilarity
#'
#' Calculates the Bray-Curtis dissimilarity between two samples
#' @param x a vector containing the first compositional sample
#' @param y a vector of length(x) containing the second compositional sample
#' @return a scalar value
#' @examples
#' data(Namib)
#' print(bray.diss(Namib$HM$x["N1",],Namib$HM$x["N2",]))
#' @export
bray.diss <- function(x,y){
    return(as.numeric(sum(abs(x-y))/sum(x+y)))
}

#' Multidimensional Scaling
#'
#' Performs classical or nonmetric Multidimensional Scaling analysis
#' of provenance data
#'
#' @param x an object of class \code{distributional},
#'     \code{compositional} or \code{diss}
#' @param classical boolean flag indicating whether classical (TRUE)
#'     or nonmetric (FALSE) MDS should be used
#' @param k the desired dimensionality of the solution
#' @param ... optional arguments to be passed onto \code{diss} (if
#'     \code{x} is of class \code{compositional} or
#'     \code{distributional}) or onto \code{cmdscale} or \code{isoMDS}
#'     (if \code{x} is of class \code{dist}).
#' @return an object of class \code{MDS}, i.e. a list containing the
#'     following items:
#'
#' \code{points}: a two column vector of the fitted configuration
#'
#' \code{classical}: a boolean flag indicating whether the MDS
#' configuration was obtained by classical (\code{TRUE}) or nonmetric
#' (\code{FALSE}) MDS.
#'
#' \code{diss}: the dissimilarity matrix used for the MDS analysis
#' 
#' \code{stress}: (only if \code{classical=TRUE}) the final stress
#' achieved (in percent)
#' @examples
#' data(Namib)
#' plot(MDS(Namib$Major,classical=TRUE))
#' @rdname MDS
#' @importFrom MASS isoMDS
#' @export
MDS <- function(x,...){ UseMethod("MDS",x) }
#' @rdname MDS
#' @export
MDS.compositional <- function(x,classical=FALSE,k=2,...){
    d <- diss.compositional(x,...)
    return(MDS.diss(d,classical=classical,k=k))
}
#' @rdname MDS
#' @export
MDS.counts <- function(x,classical=FALSE,k=2,...){
    d <- diss.counts(x,...)
    return(MDS.diss(d,classical=classical,k=k))
}
#' @rdname MDS
#' @export
MDS.distributional <- function(x,classical=FALSE,k=2,...){
    d <- diss.distributional(x,...)
    return(MDS.diss(d,classical=classical,k=k))
}
#' @rdname MDS
#' @export
MDS.diss <- function(x,classical=FALSE,k=2,...){
    out <- list() 
    if (classical){
        out$points <- stats::cmdscale(x,k=k)
    } else {
        out <- MASS::isoMDS(d=x,k=k,...)
    }
    out$classical <- classical
    out$diss <- x
    class(out) <- "MDS"
    return(out)
}

#' Centred logratio transformation
#'
#' Calculates Aitchison's centered logratio transformation for a
#' dataset of class \code{compositional} or a compositional data
#' matrix.
#' @param x an object of class \code{compositional} OR a matrix of
#'     numerical values
#' @param inverse perform the inverse inverse logratio transformation?
#' @param ... optional arguments
#' @return a matrix of CLR coordinates OR an object of class
#'     \code{compositional} (if \code{inverse=TRUE})
#' @examples
#' # The following code shows that applying provenance's PCA function
#' # to compositional data is equivalent to applying R's built-in
#' # princomp function to the CLR transformed data.
#' data(Namib)
#' plot(PCA(Namib$Major))
#' dev.new()
#' clrdat <- CLR(Namib$Major)
#' biplot(princomp(clrdat))
#' @rdname CLR
#' @export
CLR <- function(x,...){ UseMethod("CLR",x) }
#' @rdname CLR
#' @export
CLR.default <- function(x,inverse=FALSE,...){
    if (inverse){
        closure <-  rowSums(exp(x)) %*% matrix(1,nrow=1,ncol=length(x))
        out <- list()
        class(out) <- "compositional"
        out$x <- exp(x) / closure
    } else {
        g <- apply(log(x),1,mean)
        nc <- ncol(x)
        gg <- matrix(rep(g,nc),ncol=nc,byrow=FALSE)
        out <- log(x) - gg
    }
    return(out)
}
#' @rdname CLR
#' @export
CLR.compositional <- function(x,...){
    return(CLR(x$x,...))
}

#' Additive logratio transformation
#'
#' Calculates Aitchison's additive logratio transformation for a
#' dataset of class \code{compositional} or a compositional data
#' matrix.
#' @param x an object of class \code{compositional} OR a matrix of
#'     numerical values
#' @param inverse perform the inverse inverse logratio transformation?
#' @param ... optional arguments
#' @return a matrix of ALR coordinates OR an object of class
#'     \code{compositional} (if \code{inverse=TRUE}).
#' @examples
#' # logratio plot of trace element concentrations:
#' data(Namib)
#' alr <- ALR(Namib$Trace)
#' pairs(alr[,1:5])
#' title('log(X/Pb)')
#' @export
#' @rdname ALR
#' @export
ALR <- function(x,...){ UseMethod("ALR",x) }
#' @rdname ALR
#' @export
ALR.default <- function(x,inverse=FALSE,...){
    dat <- as.matrix(x)
    nr <- nrow(dat)
    nc <- ncol(dat)
    if (inverse){
        num <- matrix(1,nr,nc+1)
        num[,1:nc] <- exp(dat)
        den <- rowSums(num) %*% matrix(1,1,nc+1)
        out <- list()
        class(out) <- "compositional"
        out$x <- num / den
    } else {
        num <- log(dat[,1:(nc-1)])
        den <- log(dat[,nc]) %*% matrix(1,1,nc-1)
        out <- num - den
    }
    return(out)
}
#' @rdname ALR
#' @export
ALR.compositional <- function(x,...){
    return(ALR(x$x,...))
}

#' Principal Component Analysis
#'
#' Performs PCA of compositional data using a centred logratio
#' distance
#' @param x an object of class \code{compositional}
#' @param ... optional arguments to R's \code{princomp} function
#' @return an object of classes \code{PCA}, which is synonymous to the
#'     stats packages' \code{princomp} class.
#' @examples
#' data(Namib)
#' plot(MDS(Namib$Major,classical=TRUE))
#' dev.new()
#' plot(PCA(Namib$Major),asp=1)
#' print("This example demonstrates the equivalence of classical MDS and PCA")
#' @export
PCA <- function(x,...){
    if (methods::is(x,'compositional')){
        dat <- CLR(x)
    } else {
        dat <- x
    }
    pc <- stats::princomp(dat,...)
    class(pc) <- append("PCA",class(pc))
    return(pc)
}

#' Correspondence Analysis
#'
#' Performs Correspondence Analysis of point-counting data
#' @param x an object of class \code{counts}
#' @param nf number of correspondence factors (dimensions)
#' @param ... optional arguments to the \code{corresp} function of the
#'     \code{MASS} package
#' @return an object of classes \code{CA}, which is synonymous to the
#'     MASS packages' \code{correspondence} class.
#' @examples
#' data(Namib)
#' plot(CA(Namib$PT))
#' @export
CA <- function(x,nf=2,...){
    if (methods::is(x,'counts') |
        methods::is(x,'compositional')){
        X <- x$x
    } else {
        X <- x
    }
    out <- MASS::corresp(X,nf=nf,...)
    class(out) <- append("CA",class(out))
    return(out)
}

#' Get a subset of distributional data
#'
#' Return a subset of provenance data according to some specified
#' indices
#' @param x an object of class \code{distributional}
#' @param subset logical expression indicating elements or rows to
#'     keep: missing values are taken as false.
#' @param select a vector of sample names
#' @param ... optional arguments for the generic subset function
#' @return an object of class \code{distributional}
#' @seealso read.distributional
#' @examples
#' data(Namib)
#' coast <- c("N1","N2","T8","T13","N12","N13")
#' ZTRcoast <- subset(Namib$HM,select=coast,components=c('gt','cpx','ep'))
#' DZcoast <- subset(Namib$DZ,select=coast)
#' summaryplot(ZTRcoast,KDEs(DZcoast),ncol=2)
#' @name subset
#' @export
subset.distributional <- function(x,subset=NULL,select=NULL,...){
    out <- x
    if (!is.null(subset)){
        i <- which(subset,arr.ind=TRUE)
    } else if (!is.null(select)){
        i <- which(names(x) %in% select)
    } else {
        return(out)
    }
    if (length(x$err)==length(x$x)) out$err <- x$err[i]    
    out$x <- x$x[i]
    return(out)
}
#' @param components categories to keep
#' @rdname subset
#' @export
subset.compositional <- function(x,subset=NULL,components=NULL,select=NULL,...){
    out <- x
    if (!is.null(subset)){
        i <- which(subset,arr.ind=TRUE)
    } else if (!is.null(select)){
        i <- which(names(x) %in% select)
    } else {
        i <- 1:length(names(x))
    }
    if (!is.null(components)){
        j <- which(colnames(x$x) %in% components,arr.ind=TRUE)
    } else {
        j <- 1:length(colnames(x$x))
    }
    out$x <- x$x[i,j]
    if (methods::is(x,"SRDcorrected")){
        out$restoration <- x$restoration[i]
        for (sname in rownames(out$x)){
            out$restoration[[sname]] <- subset(x$restoration[[sname]],select=j)
        }
    }
    return(out)
}
#' @rdname subset
#' @export
subset.counts <- function(x,subset=NULL,components=NULL,select=NULL,...){
    out <- subset.compositional(x,subset=subset,select=select,components=components,...)
    if (methods::is(x,"ternary")){
        i <- which(rownames(x$raw) %in% rownames(out$x))
        j <- which(colnames(x$raw) %in% colnames(out$x))
        out$raw <- x$raw[i,j]
    }
    out
}

# returns list of dissimilarities between common items
getdisslist <- function(slist){
    dnames <- names(slist)
    lablist <- lapply(slist,function(x) names(x))
    commonlabels <- Reduce(intersect,lablist)
    for (name in dnames){
        slist[[name]] <- subset(slist[[name]],select=commonlabels)
    }
    disslist <- slist
    for (name in dnames){
        disslist[[name]] <- diss(slist[[name]])
    }
    return(disslist)
}

#' Generalised Procrustes Analysis of provenance data
#'
#' Given a number of input datasets, this function performs an MDS
#' analysis on each of these and the feeds the resulting
#' configurations into the \code{GPA()} function.
#'
#' @param ... a sequence of datasets of classes \code{distributional}
#' and \code{compositional}
#' @return an object of class \code{GPA}, i.e. a list containing the
#' following items:
#' 
#' \code{points}: a two column vector with the coordinates of the
#' group configuration
#'
#' \code{labels}: a list with the sample names
#' @author Pieter Vermeesch
#' @references Gower, J.C. (1975). Generalized Procrustes analysis,
#' Psychometrika, 40, 33-50.
#' @examples
#' data(Namib)
#' gpa <- procrustes(Namib$DZ,Namib$HM)
#' plot(gpa)
#' @importFrom MASS isoMDS
#' @seealso GPA
#' @export
procrustes <- function(...) {
    dnames <- sapply(match.call(expand.dots=TRUE)[-1], deparse)
    slist <- list(...)
    names(slist) <- dnames
    disslist <- getdisslist(slist)
    n <- length(labels(disslist[[1]]))
    m <- length(disslist)
    X <- array(dim=c(n,2,m))
    for (i in 1:m){
        md <- MDS(disslist[[i]],FALSE)
        if (md$stress < 0.05) md <- MDS(disslist[[i]],TRUE)
        X[,,i] <- md$points
    }
    result <- GPA(X)
    out <- list()
    out$points <- result
    out$labels <- labels(disslist[[1]])
    class(out) <- "GPA"
    return(out)
}

#  based on a Wikipedia algorithm
#' Generalised Procrustes Analysis of configurations
#'
#' Given a number of (2D) configurations, this function uses a
#' combination of transformations (reflections, rotations,
#' translations and scaling) to find a 'consensus' configuration which
#' best matches all the component configurations in a least-squares
#' sense.
#' 
#' @param X a list of dissimilarity matrices
#' @param scale boolean flag indicating if the transformation should include the scaling operation
#' @return a two column vector with the coordinates of the
#' group configuration
#' @seealso procrustes
#' @export
GPA <- function(X,scale=TRUE){
    if (length(dim(X))<3) {
        return(X)
    } else if (dim(X)[3]<3){
        return(procfit(X[,,1],X[,,2])$Yrot)
    } else {
        Y <- X # initialise fitted configurations
        refconf <- X[,,1] # reference configuration
        for (j in 1:100){
            for (i in 1:dim(X)[3]){
                Y[,,i] <- procfit(refconf,X[,,i])$Yrot
            }
            meanconf <- apply(Y,c(1,2),'mean')
            misfit <- sum((refconf-meanconf)^2)
            if (misfit < 1e-10){
                break
            } else {
                refconf <- meanconf
            }
        }
        return(refconf)
    }
}

# Procrustes analysis of two configurations
# based on the 'procrustes' function of the 'vegan' package
procfit <- function (X, Y, scale=TRUE, symmetric=FALSE, ...) {
    if (nrow(X) != nrow(Y)) 
        stop("Matrices have different number of rows: ", nrow(X), 
            " and ", nrow(Y))
    if (ncol(X) < ncol(Y)) {
        warning("X has fewer axes than Y: X adjusted to comform Y\n")
        addcols <- ncol(Y) - ncol(X)
        for (i in 1:addcols) X <- cbind(X, 0)
    }
    ctrace <- function(MAT) sum(MAT^2)
    c <- 1
    if (symmetric) {
        X <- scale(X, scale = FALSE)
        Y <- scale(Y, scale = FALSE)
        X <- X/sqrt(ctrace(X))
        Y <- Y/sqrt(ctrace(Y))
    }
    xmean <- apply(X, 2, mean)
    ymean <- apply(Y, 2, mean)
    if (!symmetric) {
        X <- scale(X, scale = FALSE)
        Y <- scale(Y, scale = FALSE)
    }
    XY <- crossprod(X, Y)
    sol <- svd(XY)
    A <- sol$v %*% t(sol$u)
    if (scale) {
        c <- sum(sol$d)/ctrace(Y)
    }
    Yrot <- c * Y %*% A
    b <- xmean - c * ymean %*% A
    R2 <- ctrace(X) + c * c * ctrace(Y) - 2 * c * sum(sol$d)
    reslt <- list(Yrot = Yrot, X = X, ss = R2, rotation = A, 
        translation = b, scale = c, xmean = xmean, symmetric = symmetric, 
        call = match.call())
    reslt$svd <- sol
    class(reslt) <- "procrustes"
    reslt
}

# calculate the trace of a matrix
tr <- function (m){
    if (!is.matrix(m) | (dim(m)[1] != dim(m)[2])) 
        stop("m must be a square matrix")
    return(sum(diag(m)))
}

get.data.names <- function(dlist){
    out <- c()
    for (d in dlist){
        out <- c(out,d$name)
    }
    out
}

# set minimum and maximum values of a dataset
setmM <- function(x,from=NA,to=NA,log=FALSE){
    if (is.na(from)) { from <- min(x); setm <- TRUE }
    else { setm <- FALSE }
    if (is.na(to)) { to <- max(x); setM <- TRUE }
    else { setM <- FALSE }
    if (setm) {
        if (log) { from <- from/2 }
        else {
            if (2*from-to<0) {from <- 0}
            else {from <- from-(to-from)/10}
        }
    }
    if (setM) {
        if (log) { to <- 2*to }
        else { to <- to+(to-from)/10 }
    }
    return(list(m=from,M=to))
}

#' @export
names.distributional <- function(x){
    return(names(x$x))
}
#' @export
names.compositional <- function(x){
    return(rownames(x$x))
}
#' @export
names.counts <- function(x){
    return(rownames(x$x))
}
#' @export
names.KDEs <- function(x){
    return(names(x$kdes))
}
#' @export
names.ternary <- function(x){
    return(rownames(x$x))
}

#' Calculate the number of grains required to achieve a desired level of sampling resolution
#'
#' Returns the number of grains that need to be analysed to decrease
#' the likelihood of missing any fraction greater than a given size
#' below a given level.
#' @param f the size of the smallest resolvable fraction (0<f<1)
#' @param p the probability that all n grains in the sample have missed
#' at least one fraction of size f
#' @param n, the number of grains in the sample
#' @return the number of grains needed to reduce the chance of missing
#' at least one fraction f of the total population to less than p
#' @references Vermeesch, Pieter. "How many grains are needed for a
#' provenance study?." Earth and Planetary Science Letters 224.3
#' (2004): 441-451.
#' @examples
#' # number of grains required to be 99% that no fraction greater than 5% was missed:
#' print(get.n(0.01))
#' # number of grains required to be 90% that no fraction greater than 10% was missed:
#' print(get.n(p=0.1,f=0.1))
#' @export
get.n <- function(p=0.05,f=0.05){
    n <- 1
    while(T){
        pp <- get.p(n,f)
        if (pp<p){ break }
        else {n <- n+1}
    }
    return(n)
}

#' Calculate the probability of missing a given population fraction
#'
#' For a given sample size, returns the likelihood of missing any
#' fraction greater than a given size
#' @param n the number of grains in the detrital sample
#' @param f the size of the smallest resolvable fraction (0<f<1)
#' @return the probability that all n grains in the sample have missed
#' at least one fraction of size f
#' @references Vermeesch, Pieter. "How many grains are needed for a
#' provenance study?." Earth and Planetary Science Letters 224.3
#' (2004): 441-451.
#' @examples
#' print(get.p(60))
#' print(get.p(117))
#' @export
get.p <- function(n,f=0.05){
    p <- 0
    M <- 1/f
    for (i in 1:M){
        p <- p + (-1)^(i-1) * choose(M,i)*(1-i*f)^n
    }
    return(p)
}

#' Calculate the largest fraction that is likely to be missed
#'
#' For a given sample size, returns the largest fraction which has
#' been sampled with p x 100 % likelihood.
#' @param n the number of grains in the detrital sample
#' @param p the required level of confidence
#' @return the largest fraction that is sampled with at least 100 x p%
#' certainty
#' @references Vermeesch, Pieter. "How many grains are needed for a
#' provenance study?." Earth and Planetary Science Letters 224.3
#' (2004): 441-451.
#' @examples
#' print(get.f(60))
#' print(get.f(117))
#' @export
get.f <- function(n,p=0.05){
    fmin <- 0
    fmax <- 1
    for (i in 1:100){
        f <- (fmax+fmin)/2
        if (get.p(n,f)<p) { fmax <- f }
        else { fmin <- f }
    }
    return((fmin+fmax)/2)
}

get.densities <- function(X,dtable){
    if (!(methods::is(X,"compositional") | methods::is(X,"counts")))
        stop("input is not of class compositional or counts")
    minerals <- colnames(X$x)
    i <- which(colnames(dtable) %in% colnames(X$x), arr.ind=TRUE)
    return(dtable[i])
}

ndim <- function(X){
    return(length(dim(X)))
}

sumcols <- function(X,x){
    if (ndim(X)==0){
        out <- X[x]
    } else {
        dat <- subset(X,select=x)
        out <- rowSums(dat)
    }    
    return(out)
}

#' Group components of a composition
#'
#' Adds several components of a composition together into a single
#' component
#' @param X a compositional dataset
#' @param ... a series of new labels assigned to strings or vectors of
#'     strings denoting the components that need amalgamating
#' @return an object of the same class as X with fewer components
#' @examples
#' data(Namib)
#' HMcomponents <- c("zr","tm","rt","TiOx","sph","ap","ep",
#'                   "gt","st","amp","cpx","opx")
#' am <- amalgamate(Namib$PTHM,feldspars=c("KF","P"),
#'                  lithics=c("Lm","Lv","Ls"),heavies=HMcomponents)
#' plot(ternary(am))
#' @rdname amalgamate
#' @export
amalgamate <- function(X,...){ UseMethod("amalgamate",X) }
#' @rdname amalgamate
#' @export
amalgamate.default <- function(X,...){
    groups <- list(...)
    ng <- length(groups)
    labels <- names(groups)
    out <- NULL
    for (i in 1:ng){
        colsum <- sumcols(X,groups[[i]])
        out <- cbind(out,colsum)
    }
    colnames(out) <- labels
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.compositional <- function(X,...){
    out <- X
    out$x <- amalgamate(X$x,...)
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.counts <- function(X,...){
    out <- X
    out$x <- amalgamate(X$x,...)
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.SRDcorrected <- function(X,...){
    out <- X
    out$x <- amalgamate.default(X$x,...)
    for (sname in names(X$restoration)){
        out$restoration[[sname]] <-
            amalgamate.default(X$restoration[[sname]],...)
    }
    return(out)
}

#' Combine samples of distributional data
#'
#' Lumps all single grain analyses of several samples together under a
#' new name
#' @param X a distributional dataset
#' @param ... a series of new labels assigned to strings or vectors of
#'     strings denoting the samples that need amalgamating
#' @return a distributional data object with fewer samples than X
#' @examples
#' data(Namib)
#' combined <- combine(Namib$DZ,
#'                     east=c('N3','N4','N5','N6','N7','N8','N9','N10'),
#'                     west=c('N1','N2','N11','N12','T8','T13'))
#' summaryplot(KDEs(combined))
#' @export
combine <- function(X,...){
    out <- X
    groups <- list(...)
    ng <- length(groups)
    labels <- names(groups)
    out$x <- list()
    out$err <- list()
    loadErr <- (length(X$err)>0)
    for (i in 1:ng){
        out$x[[labels[i]]] <- NULL
        for (g in groups[[i]]){
            out$x[[labels[i]]] <- c(out$x[[labels[i]]],X$x[[g]])
            if (loadErr) {
                out$err[[labels[i]]] <- c(out$err[[labels[i]]],X$err[[g]])
            }
        }
    }
    return(out)    
}

# from package mvtnorm:
dmvnorm <- function(x,mean=rep(0,p),sigma=diag(p),log=FALSE) {
    if (is.vector(x)) 
        x <- matrix(x, ncol = length(x))
    p <- ncol(x)
    if (!missing(mean)) {
        if (!is.null(dim(mean))) 
            dim(mean) <- NULL
        if (length(mean) != p) 
            stop("mean and sigma have non-conforming size")
    }
    if (!missing(sigma)) {
        if (p != ncol(sigma)) 
            stop("x and sigma have non-conforming size")
        if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
            check.attributes = FALSE)) 
            stop("sigma must be a symmetric matrix")
    }
    dec <- tryCatch(chol(sigma), error = function(e) e)
    if (inherits(dec, "error")) {
        x.is.mu <- colSums(t(x) != mean) == 0
        logretval <- rep.int(-Inf, nrow(x))
        logretval[x.is.mu] <- Inf
    }
    else {
        tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
        rss <- colSums(tmp^2)
        logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * 
            pi) - 0.5 * rss
    }
    names(logretval) <- rownames(x)
    if (log) 
        logretval
    else exp(logretval)
}
