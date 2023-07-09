#' Kolmogorov-Smirnov dissimilarity
#'
#' Returns the Kolmogorov-Smirnov dissimilarity between two samples
#'
#' @param x the first sample as a vector
#' @param y the second sample as a vector
#' @param ... optional arguments
#' @return a scalar value representing the maximum vertical distance
#' between the two cumulative distributions
#' @examples
#' data(Namib)
#' print(KS.diss(Namib$DZ$x[['N1']],Namib$DZ$x[['T8']]))
#' @rdname KS.diss
#' @export
KS.diss <- function(x,...){ UseMethod("KS.diss",x) }
#' @rdname KS.diss
#' @export
KS.diss.default <- function(x,y,...){
    IsoplotR::diss(x,y,method="KS")
}
#' @rdname KS.diss
#' @export
KS.diss.distributional <- function(x,...){
    diss.distributional(x,method="KS",...)
}

#' Kuiper dissimilarity
#'
#' Returns the Kuiper dissimilarity between two samples
#'
#' @param x the first sample as a vector
#' @param y the second sample as a vector
#' @param ... optional arguments
#' @return a scalar value representing the sum of the maximum vertical
#' distances above and below the cumulative distributions of x and y
#' @examples
#' data(Namib)
#' print(Kuiper.diss(Namib$DZ$x[['N1']],Namib$DZ$x[['T8']]))
#' @rdname Kuiper.diss
#' @export
Kuiper.diss <- function(x,...){ UseMethod("Kuiper.diss",x) }
#' @rdname Kuiper.diss
#' @export
Kuiper.diss.default <- function(x,y,...){
    xy <- c(x,y)
    cad1 <- stats::ecdf(x)
    cad2 <- stats::ecdf(y)
    d <- cad2(xy) - cad1(xy)
    M <- max(d)
    m <- min(d)
    if (M<0) M <- 0
    if (m>0) m <- 0
    return(M-m)
}
#' @rdname Kuiper.diss
#' @export
Kuiper.diss.distributional <- function(x,...){
    diss.distributional(x,method="Kuiper",...)
}

#' Calculate the dissimilarity matrix between two datasets of class
#' \code{distributional}, \code{compositional}, \code{counts} or
#' \code{varietal}
#'
#' Calculate the dissimilarity matrix between two datasets of class
#' \code{distributional} or \code{compositional} using the Kolmogorov-Smirnov,
#' Sircombe-Hazelton, Aitchison or Bray-Curtis distance
#' 
#' @param x an object of class \code{distributional},
#'     \code{compositional} or \code{counts}
#' @param method if \code{x} has class \code{distributional}: either
#'     \code{"KS"}, \code{"Wasserstein"}, \code{"Kuiper"} or
#'     \code{"SH"};
#'
#' if \code{x} has class \code{compositional}: either
#' \code{"aitchison"} or \code{"bray"};
#'
#' if \code{x} has class \code{counts}: either \code{"chisq"} or
#' \code{"bray"};
#'
#' if \code{x} has class \code{varietal}: either \code{"KS"},
#' \code{"W2_1D"} or \code{"W2"}.
#' @param log logical. If \code{TRUE}, subjects the distributional
#'     data to a logarithmic transformation before calculating the
#'     Wasserstein distance.
#' @param verbose logical. If \code{TRUE}, gives progress updates
#'     during the construction of the dissimilarity matrix.
#' @param ... optional arguments
#' @details \code{"KS"} stands for the Kolmogorov-Smirnov statistic,
#'     \code{"W2_1D"} for the 1-dimensional Wasserstein-2 distance,
#'     \code{"Kuiper"} for the Kuiper statistic, \code{"SH"} for the
#'     Sircombe-Hazelton distance, \code{"aitchison"} for the
#'     Aitchison logratio distance, \code{"bray"} for the Bray-Curtis
#'     distance, \code{"chisq"} for the Chi-square distance, and "W2"
#'     for the 2-dimensional Wasserstein-2 distance.
#' @examples
#' data(Namib)
#' print(round(100*diss(Namib$DZ)))
#' @return an object of class \code{diss}
#' @seealso KS.diss bray.diss SH.diss Wasserstein.diss Kuiper.diss
#' @rdname diss
#' @export
diss <- function(x,method,...){ UseMethod("diss",x) }
#' @rdname diss
#' @export
diss.distributional <- function(x,method=NULL,log=FALSE,verbose=FALSE,...) {
    if (!is.null(method)) x$method <- method
    n <- length(x$x)
    d <- mat.or.vec(n,n)
    snames <- names(x$x)
    rownames(d) <- snames
    colnames(d) <- snames
    if (x$method=="SH") c2 <- getc2(x)
    for (i in 1:n){
        for (j in 1:n){
            if (verbose){
                msg <- paste0('Comparing ',snames[i],' with ',snames[j])
                print(msg)
            }
            if (x$method=="SH"){
                d[i,j] <- SH.diss(x,i,j,c.con=c2)
            } else if (x$method%in%c("W2","W2_1D","Wasserstein")){
                if (log){
                    a <- log(x$x[[i]])
                    b <- log(x$x[[j]])
                } else {
                    a <- x$x[[i]]
                    b <- x$x[[j]]
                }
                d[i,j] <- Wasserstein.diss(a,b)
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
diss.compositional <- function(x,method=NULL,...){
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
diss.counts <- function(x,method=NULL,...){
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
#' @rdname diss
#' @export
diss.varietal <- function(x,method=NULL,...){
    if (!is.null(method)) x$method <- method
    if (x$method=='W2'){
        out <- Wasserstein.diss(x,...)
    } else {
        xd <- varietal2distributional(x)
        out <- diss.distributional(xd,method=x$method)
    }
    return(out)
}

#' Bray-Curtis dissimilarity
#'
#' Calculates the Bray-Curtis dissimilarity between two samples
#' @param x a vector containing the first compositional sample
#' @param y a vector of \code{length(x)} containing the second
#'     compositional sample
#' @param ... optional arguments
#' @return a scalar value
#' @examples
#' data(Namib)
#' print(bray.diss(Namib$HM$x["N1",],Namib$HM$x["N2",]))
#' @rdname bray.diss
#' @export
bray.diss <- function(x,...){ UseMethod("bray.diss",x) }
#' @rdname bray.diss
#' @export
bray.diss.default <- function(x,y,...){
    return(as.numeric(sum(abs(x-y))/sum(x+y)))
}
#' @rdname bray.diss
#' @export
bray.diss.compositional <- function(x,...){
    diss.compositional(x,method="bray",...)
}

# returns list of normalised dissimilarities between common items
getdisslist <- function(slist){
    dnames <- names(slist)
    lablist <- lapply(slist,function(x) names(x))
    commonlabels <- Reduce(intersect,lablist)
    if (length(commonlabels)>0){
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
    } else {
        stop('This dataset contains no common sample names.')
    }
}
