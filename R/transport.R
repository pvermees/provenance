distmat <- function(x,y){
    m <- nrow(x)
    n <- nrow(y)
    xy <- x %*% t(y)
    xx <- matrix( rep(apply(x*x,1,sum), n), m, n, byrow=FALSE)
    yy <- matrix( rep(apply(y*y,1,sum), m), m, n, byrow=TRUE)
    out <- sqrt(pmax(xx+yy-2*xy, 0))
    return(out)
}

#' Wasserstein distance
#'
#' Returns the Wasserstein distance between two samples
#'
#' @param x the first sample as a vector
#' @param y the second sample as a vector
#' @param log logical. Take the lograthm of the data before
#'     calculating the distances?
#' @param package the name of the package that provides the 2D
#'     Wasserstein distance. Currently, this can be either
#'     \code{'transport'} or \code{T4transport}.
#' @param verbose logical. If \code{TRUE}, gives progress updates
#'     during the construction of the dissimilarity matrix.
#' @param ... optional arguments to the
#'     \code{transport::wasserstein()} or
#'     \code{T4transport::wasserstein()} functions. Warning: the
#'     latter function is very slow.
#' @author The default S3 method was written by Pieter Vermeesch,
#'     using modified code from Dominic Schuhmacher's \code{transport}
#'     package (\code{transport1d} function), as implemented in
#'     \code{IsoplotR}.
#' @return a scalar value
#' @examples
#' data(Namib)
#' print(Wasserstein.diss(Namib$DZ$x[['N1']],Namib$DZ$x[['T8']]))
#' @rdname Wasserstein.diss
#' @export
Wasserstein.diss <- function(x,...){ UseMethod("Wasserstein.diss",x) }
#' @rdname Wasserstein.diss
#' @export
Wasserstein.diss.default <- function(x,y,...){
    IsoplotR::diss(x,y,method="W2")
}
#' @rdname Wasserstein.diss
#' @export
Wasserstein.diss.distributional <- function(x,log=FALSE,...){
    diss.distributional(x,method="W2_1D",log=log,...)
}
#' @rdname Wasserstein.diss
#' @export
Wasserstein.diss.varietal <- function(x,package="transport",verbose=FALSE,...){
    snames <- names(x$x)
    ns <- length(snames)
    out <- matrix(0,ns,ns)
    rownames(out) <- colnames(out) <- snames
    for (snamei in snames){
        xi <- CLR(x$x[[snamei]])
        ni <- nrow(xi)
        for (snamej in snames){
            if (verbose){
                msg <- paste0('Comparing ',snamei,' with ',snamej)
            }
            xj <- CLR(x$x[[snamej]])
            if (!identical(snamei,snamej)){
                if (identical(package,"T4transport")){
                    W <- T4transport::wasserstein(X=xi,Y=xj,...)
                    out[snamei,snamej] <- W$distance
                } else if (identical(package,"transport")){
                    nj <- nrow(xj)
                    wi <- rep(1,ni)/ni
                    wj <- rep(1,nj)/nj
                    a <- transport::wpp(xi,mass=wi)
                    b <- transport::wpp(xj,mass=wj)
                    out[snamei,snamej] <- transport::wasserstein(a=a,b=b,...)
                } else {
                    stop("Unknown package")
                }
            }
        }
    }
    out <- stats::as.dist(out)
    class(out) <- append("diss",class(out))
    return(out)
}
