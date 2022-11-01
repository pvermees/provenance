#' Multidimensional Scaling
#'
#' Performs classical or nonmetric Multidimensional Scaling analysis
#' of provenance data
#'
#' @param x an object of class \code{distributional},
#'     \code{compositional}, \code{counts}, \code{varietal} or
#'     \code{diss}
#' @param classical boolean flag indicating whether classical
#'     (\code{TRUE}) or nonmetric (\code{FALSE}) MDS should be used
#' @param k the desired dimensionality of the solution
#' @param ... optional arguments
#'
#' If \code{x} has class \code{distributional}, \code{...} is passed
#' on to \code{diss.distributional}.
#' 
#' If \code{x} has class \code{compositional}, \code{...} is passed on
#' to \code{diss.compositional}.
#'
#' If \code{x} has class \code{counts}, \code{...} is passed on to
#' \code{diss.counts}.
#'
#' If \code{x} has class \code{varietal}, \code{...} is passed on to
#' \code{diss.varietal}.
#'
#' Otherwise, \code{...} is passed on to \code{cmdscale} (if
#'     \code{classical=TRUE}) or \code{isoMDS} (if
#'     \code{classical=FALSE}).
#' 
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
MDS.default <- function(x,classical=FALSE,k=2,...){
    if (classical){
        out <- list()
        out$points <- stats::cmdscale(x,k=k,...)
    } else {
        out <- MASS::isoMDS(d=x,k=k,...)
    }
    out$classical <- classical
    out$diss <- x
    out$nb <- 0
    class(out) <- "MDS"
    return(out)

}
#' @rdname MDS
#' @export
MDS.compositional <- function(x,classical=FALSE,k=2,...){
    d <- diss.compositional(x,...)
    MDS.default(d,classical=classical,k=k)
}
#' @rdname MDS
#' @export
MDS.counts <- function(x,classical=FALSE,k=2,...){
    d <- diss.counts(x,...)
    MDS.default(d,classical=classical,k=k)
}
#' @param nb number of bootstrap resamples. If \code{nb>0}, then
#'     \code{plot.MDS(...)} will visualise the sampling uncertainty as
#'     polygons (inspired by Nordsvan et al. 2020). The bigger
#'     \code{nb}, the slower the calculations. \code{nb=10} seems a
#'     good compromise.
#' 
#' @references
#' Nordsvan, A.R., Kirscher, U., Kirkland, C.L., Barham, M. and
#' Brennan, D.T., 2020. Resampling (detrital) zircon age distributions
#' for accurate multidimensional scaling solutions. Earth-Science
#' Reviews, p.103149.
#' 
#' Vermeesch, P., 2013, Multi-sample comparison of detrital age
#' distributions. Chemical Geology v.341, 140-146,
#' doi:10.1016/j.chemgeo.2013.01.010
#' 
#' @rdname MDS
#' @export
MDS.distributional <- function(x,classical=FALSE,k=2,nb=0,...){
    if (nb>0) X <- resample(x,nb=nb)
    else X <- x
    d <- diss.distributional(X,...)
    out <- MDS.default(d,classical=classical,k=k)
    out$nb <- nb
    return(out)
}
#' @rdname MDS
#' @export
MDS.varietal <- function(x,classical=FALSE,k=2,nb=0,...){
    if (nb>0) X <- resample(x,nb=nb)
    else X <- x
    print('Constructing the dissimilarity matrix...')
    d <- diss.varietal(X,...)
    print('Computing the MDS configuration...')
    out <- MDS.default(d,classical=classical,k=k)
    out$nb <- nb
    return(out)
}

#' Principal Component Analysis
#'
#' Performs PCA of compositional data using a centred logratio
#' distance
#' @param x an object of class \code{compositional}
#' @param ... optional arguments to R's \code{princomp} function
#' @return an object of classes \code{PCA}, which is synonymous to the
#'     stats package's \code{prcomp} class.
#' @examples
#' data(Namib)
#' plot(MDS(Namib$Major,classical=TRUE))
#' dev.new()
#' plot(PCA(Namib$Major),asp=1)
#' print("This example demonstrates the equivalence of classical MDS and PCA")
#' @export
PCA <- function(x,...){
    if (methods::is(x,'compositional') |
        methods::is(x,'counts')){
        dat <- CLR.compositional(x)
    } else {
        dat <- x
    }
    pc <- stats::prcomp(dat,...)
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
#'     MASS package's \code{correspondence} class.
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

#' Generalised Procrustes Analysis of provenance data
#'
#' Given a number of input datasets, this function performs an MDS
#' analysis on each of these and the feeds the resulting
#' configurations into the \code{GPA()} function.
#'
#' @param ... a sequence of datasets of classes \code{distributional},
#'     \code{counts}, \code{compositional} and \code{varietal} OR a
#'     single object of class \code{varietal}.
#' @return an object of class \code{GPA}, i.e. a list containing the
#'     following items:
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
#' gpa1 <- procrustes(Namib$DZ,Namib$HM)
#' plot(gpa1)
#' 
#' data(SNSM)
#' gpa2 <- procrustes(SNSM$ap)
#' plot(gpa2)
#' @importFrom MASS isoMDS
#' @seealso GPA
#' @export
procrustes <- function(...) {
    slist <- list(...)
    names(slist) <- get.data.names(slist)
    if (length(slist)==1 & 'varietal' %in% class(slist[[1]])){
        slist <- varietal2distributional(slist[[1]],bycol=TRUE)
    }
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
#' translations and scaling) to find a `consensus' configuration which
#' best matches all the component configurations in a least-squares
#' sense.
#' 
#' @param X a list of dissimilarity matrices
#' @param scale boolean flag indicating if the transformation should
#'     include the scaling operation
#' @return a two column vector with the coordinates of the group
#'     configuration
#' @seealso procrustes
#' @export
GPA <- function(X,scale=TRUE){
    if (length(dim(X))<3) {
        return(X)
    } else if (dim(X)[3]<3){
        return(procfit(X[,,1],X[,,2])$Yhat)
    } else {
        Y <- X # initialise fitted configurations
        refconf <- X[,,1]
        misfit <- Inf
        for (j in 1:100){
            for (i in 1:dim(X)[3]){
                Y[,,i] <- procfit(X=refconf,Y=X[,,i])$Yhat
            }
            meanconf <- apply(Y,c(1,2),'mean')
            newmisfit <- sum((refconf-meanconf)^2)
            if (abs(newmisfit-misfit) < 1e-10){
                break
            } else {
                misfit <- newmisfit
                refconf <- procfit(X=refconf,Y=meanconf)$Yhat                
            }
        }
        return(refconf)
    }
}

# Procrustes analysis of two configurations
# based on the 'Procrustes' function of the 'smacof' package
procfit <- function(X, Y) {
    n <- dim(X)[1]
    E <- diag(1, nrow = n)
    eins <- rep(1, n)
    k <- 1/n
    Z <- E - k * eins %*% t(eins)
    C <- t(X) %*% Z %*% Y
    s <- svd(C)
    f <- diag(s$d)
    P <- s$u
    Q <- s$v
    T <- Q %*% t(P)
    streck <- sum(diag((t(X) %*% Z %*% Y %*% T)))/sum(diag((t(Y) %*% 
        Z %*% Y)))
    trans <- as.vector(k * t(X - streck * Y %*% T) %*% eins)
    Yhut <- streck * Y %*% T + eins %*% t(trans)
    colnames(Yhut) <- rownames(T) <- colnames(T) <- names(trans) <- colnames(Y)
    dX <- stats::dist(X)
    dY <- stats::dist(Y)
    dYhat <- stats::dist(Yhut)
    cong <- sum(dX * dY)/(sqrt(sum(dX^2)) * sqrt(sum(dY^2)))
    alien <- sqrt(1 - cong^2)
    pairdist <- sort(sqrt(rowSums((X - Yhut)^2)))
    res <- list(X = X, Y = Y, Yhat = Yhut, translation = trans, 
        dilation = streck, rotation = T, confdistX = dX, confdistY = dY, 
        confdistYhat = dYhat, congcoef = cong, aliencoef = alien, 
        pairdist = pairdist, call = match.call())
    class(res) <- "procrustes"
    return(res)
}
