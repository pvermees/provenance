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
        out <- as.compositional(exp(x) / closure)
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
#' @param x an object of class \code{compositional} OR a matrix of numerical values
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
        den <- as.matrix(rowSums(num)) %*% matrix(1,1,nc+1)
        out <- list()
        out$x <- num/den
        out$method <- 'KS'
        out$colmap <- 'rainbow'
        class(out) <- 'compositional'
    } else {
        num <- log(dat[,1:(nc-1),drop=FALSE])
        den <- log(subset(dat,select=nc)) %*% matrix(1,1,nc-1)
        out <- num - den
    }
    return(out)
}
#' @rdname ALR
#' @export
ALR.compositional <- function(x,...){
    return(ALR(x$x,...))
}

#' Get a subset of provenance data
#'
#' Return a subset of provenance data according to some specified
#' indices
#' @param x an object of class \code{distributional},
#'     \code{compositional}, \code{counts} or \code{varietal}.
#' @param subset logical expression indicating elements or rows to
#'     keep: missing values are taken as false.
#' @param select a vector of sample names
#' @param ... optional arguments for the generic subset function
#' @return an object of the same class as \code{x}
#' @seealso \link{amalgamate}, \link{combine}
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
        i <- match(select,names(x))
    } else {
        return(out)
    }
    if (length(x$err)==length(x$x)) out$err <- x$err[i]    
    out$x <- x$x[i]
    return(out)
}
#' @param components vector of categories (column names) to keep
#' @rdname subset
#' @export
subset.compositional <- function(x,subset=NULL,components=NULL,select=NULL,...){
    out <- x
    if (!is.null(subset)){
        i <- which(subset,arr.ind=TRUE)
    } else if (!is.null(select)){
        i <- match(select,names(x))
    } else {
        i <- 1:nrow(x$x)
    }
    if (!is.null(components)){
        j <- match(components,colnames(x$x))
    } else {
        j <- 1:ncol(x$x)
    }
    out$x <- x$x[i,j,drop=FALSE]
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
    out <- subset.compositional(x,subset=subset,select=select,
                                components=components,...)
    if (methods::is(x,"ternary")){
        i <- match(rownames(out$x),rownames(x$raw))
        j <- match(colnames(out$x),colnames(x$raw))
        out$raw <- x$raw[i,j,drop=FALSE]
    }
    out
}
#' @rdname subset
#' @export
subset.varietal <- function(x,subset=NULL,components=NULL,select=NULL,...){
    if (is.null(subset)){
        X <- x$x
    } else {
        X <- x$x[subset,]
    }
    if (!is.null(components)){
        X <- X[,components]
    }
    out <- x
    if (is.null(select)){
        out$x <- X
        out$snames <- NULL
        for (sname in x$snames){
            matches <- grepl(sname,rownames(X))
            if (any(matches)) out$snames <- c(out$snames,sname)
        }
    } else {
        out$x <- NULL
        out$snames <- select
        for (sname in select){
            matches <- grepl(sname,rownames(X))
            out$x <- rbind(out$x,X[matches,])
        }
    }
    out
}

# calculate the trace of a matrix
tr <- function (m){
    if (!is.matrix(m) | (dim(m)[1] != dim(m)[2])) 
        stop("m must be a square matrix")
    return(sum(diag(m)))
}

get.data.names <- function(dlist){
    out <- c()
    nd <- length(dlist)
    for (i in 1:nd){
        if (is.null(dlist[[i]]$name)) dname <- i
        else dname <- dlist[[i]]$name
        out <- c(out,dname)
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
    out <- names(x$x)
    if (is.null(out)) out <- 1:length(x$x)
    return(out)
}
#' @export
names.compositional <- function(x){
    out <- rownames(x$x)
    if (is.null(out)) out <- 1:nrow(x$x)
    return(out)
}
#' @export
names.counts <- function(x){
    out <- rownames(x$x)
    if (is.null(out)) out <- 1:nrow(x$x)
    return(out)
}
#' @export
names.varietal <- function(x){
    return(x$snames)
}

#' @export
names.KDEs <- function(x){
    out <- names(x$kdes)
    if (is.null(out)) out <- 1:length(x$kdes)
    return(out)
}
#' @export
names.ternary <- function(x){
    out <- rownames(x$x)
    if (is.null(out)) out <- 1:nrow(x$x)
    return(out)
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
#' @param x a compositional dataset
#' @param ... a series of new labels assigned to strings or vectors of
#'     strings denoting the components that need amalgamating
#' @return an object of the same class as \code{X} with fewer
#'     components
#' @examples
#' data(Namib)
#' HMcomponents <- c("zr","tm","rt","TiOx","sph","ap","ep",
#'                   "gt","st","amp","cpx","opx")
#' am <- amalgamate(Namib$PTHM,feldspars=c("KF","P"),
#'                  lithics=c("Lm","Lv","Ls"),heavies=HMcomponents)
#' plot(ternary(am))
#' @rdname amalgamate
#' @export
amalgamate <- function(x,...){ UseMethod("amalgamate",x) }
#' @rdname amalgamate
#' @export
amalgamate.default <- function(x,...){
    groups <- list(...)
    ng <- length(groups)
    labels <- names(groups)
    out <- NULL
    for (i in 1:ng){
        colsum <- sumcols(x,groups[[i]])
        out <- cbind(out,colsum)
    }
    colnames(out) <- labels
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.compositional <- function(x,...){
    out <- x
    out$x <- amalgamate(x$x,...)
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.counts <- function(x,...){
    out <- x
    out$x <- amalgamate(x$x,...)
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.SRDcorrected <- function(x,...){
    out <- x
    out$x <- amalgamate.default(x$x,...)
    for (sname in names(x$restoration)){
        out$restoration[[sname]] <-
            amalgamate.default(x$restoration[[sname]],...)
    }
    return(out)
}
#' @rdname amalgamate
#' @export
amalgamate.varietal <- function(x,...){
    out <- x
    out$x <- amalgamate(x$x,...)
    return(out)
}

#' Combine samples of distributional data
#'
#' Lumps all single grain analyses of several samples together under a
#' new name
#' @param x a distributional dataset
#' @param ... a series of new labels assigned to strings or vectors of
#'     strings denoting the samples that need amalgamating
#' @return a distributional data object with fewer samples than
#'     \code{x}
#' @examples
#' data(Namib)
#' combined <- combine(Namib$DZ,
#'                     east=c('N3','N4','N5','N6','N7','N8','N9','N10'),
#'                     west=c('N1','N2','N11','N12','T8','T13'))
#' summaryplot(KDEs(combined))
#' @export
combine <- function(x,...){
    out <- x
    groups <- list(...)
    ng <- length(groups)
    labels <- names(groups)
    out$x <- list()
    out$err <- list()
    loadErr <- (length(x$err)>0)
    for (i in 1:ng){
        out$x[[labels[i]]] <- NULL
        for (g in groups[[i]]){
            out$x[[labels[i]]] <- c(out$x[[labels[i]]],x$x[[g]])
            if (loadErr) {
                out$err[[labels[i]]] <- c(out$err[[labels[i]]],x$err[[g]])
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

removeNAcols <- function(x){
    bad <- apply(apply(x,2,is.na),2,all)
    subset(x,select=!bad)
}

resample <- function(x,nb=10){
    snames <- names(x$x)
    ns <- length(snames) # number of samples
    out <- x
    for (i in 1:nb){
        for (j in 1:ns){
            sname <- paste0(snames[j],'[',i,']')
            ng <- length(x$x[[j]]) # number of grains
            out$x[[sname]] <- sample(x$x[[j]],size=ng,replace=TRUE)
        }
    }
    out
}

#' Convert varietal to distributional data
#'
#' Convert an object of class \code{varietal} either to a list of
#' distributional objects by breaking it up into separate elements, or
#' to a single distributional object corresponding to the first
#' principal component.
#'
#' @param x an object of class \code{varietal}.
#' @param bycol logical. If \code{TRUE}, returns a list of
#'     distributional objects (one for each element). If \code{FALSE},
#'     returns a single distributional object (containing the PC1
#'     scores for each sample).
#' @param plot logical. If \code{TRUE}, shows the PCA biplot that is
#'     used when \code{bycol} is \code{FALSE}.
#' @examples
#' Apfile <- system.file("apatite.csv",package="provenance")
#' Ap <- read.varietal(fn=Apfile,snames=3)
#' varietal2distributional(Ap,bycol=FALSE,plot=TRUE)
#' @export
varietal2distributional <- function(x,bycol=FALSE,plot=FALSE){
    template <- list(x=list(),colmap='rainbow',method=x$method)
    class(template) <- 'distributional'
    if (bycol){
        out <- list()
        for (cn in colnames(x$x)){
            out[[cn]] <- template
            out[[cn]]$xlab <- 'concentration'
            for (sname in x$snames){
                matches <- grepl(sname,rownames(x$x))
                if (any(matches)){
                    out[[cn]]$name <- sname
                    out[[cn]]$x[[sname]] <- x$x[matches,cn]
                }
            }
            out[[cn]]$breaks <- getbreaks(out[[cn]]$x)
        }
    } else {
        pc <- PCA(x$x,scale.=TRUE)
        if (plot) plot(pc)
        out <- template
        out$name <- x$name
        out$rotation <- pc$rotation[,'PC1']
        out$breaks <- getbreaks(pc$x[,'PC1'])
        out$xlab <- 'PC1'
        for (sname in x$snames){
            matches <- grepl(sname,rownames(x$x))
            if (any(matches)){
                out$x[[sname]] <- pc$x[matches,'PC1']
            }
        }
    }
    out
}

# dat = vector of numbers
getbreaks <- function(dat){
    d <- unlist(dat)
    ns <- length(dat)
    ng <- length(d)
    nb <- log(ng/ns,base=2)+1
    seq(min(d),max(d),length.out=nb+1)
}

RStudio <- function(){
    Sys.getenv("RSTUDIO") == '1'
}
