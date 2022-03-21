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
#' @param x an object of class \code{distributional},
#'     \code{compositional} or \code{counts}
#' @param method (optional) either "KS", "Kuiper", "SH", "aitchison",
#'     "bray" or "chisq"
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
    NN <- sum(x$x)
    RR <- rowSums(x$x)
    CC <- colSums(x$x)
    for (i in 1:ns){
        for (j in 1:ns){
            if (x$method=='bray'){
                d[i,j] <- bray.diss(x$x[i,],x$x[j,])
            } else { # chisq
                d[i,j] <- sqrt(sum( (NN/CC)*(x$x[i,]/RR[i] - x$x[j,]/RR[j])^2 ))
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
#' @param y a vector of \code{length(x)} containing the second
#'     compositional sample
#' @return a scalar value
#' @examples
#' data(Namib)
#' print(bray.diss(Namib$HM$x["N1",],Namib$HM$x["N2",]))
#' @export
bray.diss <- function(x,y){
    return(as.numeric(sum(abs(x-y))/sum(x+y)))
}

# returns list of normalised dissimilarities between common items
getdisslist <- function(slist){
    dnames <- names(slist)
    lablist <- lapply(slist,function(x) names(x))
    commonlabels <- Reduce(intersect,lablist)
    for (name in dnames){
        slist[[name]] <- subset(slist[[name]],select=commonlabels)
    }
    ns <- length(commonlabels)
    disslist <- slist
    for (name in dnames){
        dl <- diss(slist[[name]])
        # normalise according to pers. comm. by Jan de Leeuw
        disslist[[name]] <- dl*sqrt(ns*(ns-1)*0.5/sum(dl^2))
    }
    return(disslist)
}

