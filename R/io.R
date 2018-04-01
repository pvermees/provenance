#' Read a .csv file with continuous (detrital zircon) data
#'
#' Reads a data table containing continuous data (e.g. detrital zircon
#' ages)
#' @param fname the path of a .csv file with the input data, arranged
#'     in columns.
#' @param errorfile the (optional) path of a .csv file with the
#'     standard errors of the input data, arranged by column in the
#'     same order as \code{fname}. Must be specified if the data are
#'     to be compared with the Sircombe-Hazelton dissimilarity.
#' @param method an optional string specifying the dissimilarity
#'     measure which should be used for comparing this with other
#'     datasets. Should be one of either \code{"KS"} (for
#'     Kolmogorov-Smirnov) or \code{"SH"} (for Sircombe and
#'     Hazelton). If \code{method = "SH"}, then \code{errorfile}
#'     should be specified. If \code{method = "SH"} and
#'     \code{errorfile} is unspecified, then the program will default
#'     back to the Kolmogorov-Smirnov dissimilarity.
#' @param xlab an optional string specifying the nature and units of
#'     the data.  This string is used to label kernel density
#'     estimates.
#' @param colmap an optional string with the name of one of R's
#'     built-in colour palettes (e.g., heat.colors, terrain.colors,
#'     topo.colors, cm.colors), which are to be used for plotting the
#'     data.
#' @return an object of class \code{distributional}, i.e. a list with
#'     the following items:
#' 
#' \code{x}: a named list of vectors containing the numerical data for
#' each sample
#' 
#' \code{err}: an (optional) named list of vectors containing the
#' standard errors of \code{x}
#'
#' \code{method}: either "KS" (for Kolmogorov-Smirnov), "Kuiper" (for
#' the Kuiper statistic) or "SH" (for Sircombe Hazelton)
#' 
#' \code{breaks}: a vector with the locations of the histogram bin edges
#' 
#' \code{xlab}: a string containing the label to be given to the
#' x-axis on all plots
#'
#' \code{colmap}: the colour map provided by the input argument
#'
#' @examples
#'     agefile <- system.file("DZ.csv",package="provenance")
#'     errfile <- system.file("DZerr.csv",package="provenance")
#'     DZ <- read.distributional(agefile,errfile)
#'     plot(KDE(DZ$x$N1))
#' @export
read.distributional <- function(fname,errorfile=NA,method="KS",
                                xlab="age [Ma]",colmap='rainbow') {
    out <- list()
    out$name <- basename(substr(fname,1,nchar(fname)-4))
    if (method=="SH" & is.na(errorfile)) method <- "KS"
    class(out) <- "distributional"
    out$method <- method
    out$x <- list()
    out$err <- list()
    out$colmap <- colmap
    dat <- utils::read.csv(fname,header=TRUE)
    ns = length(dat)
    for (i in 1:ns){
        out$x[[names(dat)[i]]] = dat[!is.na(dat[,i]),i]
    }
    if (!is.na(errorfile)){
        err <- utils::read.csv(errorfile,header=TRUE)
        for (i in 1:ns) {
            out$err[[names(dat)[i]]] = dat[!is.na(err[,i]),i]
        }
    }
    d <- unlist(out$x)
    ng <- length(d) # number of grains
    nb <- log(ng/ns,base=2)+1
    out$breaks <- seq(min(d),max(d),length.out=nb+1)
    out$xlab <- xlab
    return(out)
}

#' Read a .csv file with compositional data
#'
#' Reads a data table containing compositional data (e.g. chemical
#' concentrations)
#'
#' @param fname a string with the path to the .csv file
#' @param method either "bray" (for the Bray-Curtis distance) or
#'     "aitchison" (for Aitchison's central logratio distance). If
#'     omitted, the function defaults to 'aitchison', unless there are
#'     zeros present in the data.
#' @param colmap an optional string with the name of one of R's
#'     built-in colour palettes (e.g., heat.colors, terrain.colors,
#'     topo.colors, cm.colors), which are to be used for plotting the
#'     data.
#' @return an object of class \code{compositional}, i.e. a list with
#'     the following items:
#' 
#' \code{x}: a data frame with the samples as rows and the categories as columns
#'
#' \code{method}: either "aitchison" (for Aitchison's centred logratio
#' distance) or "bray" (for the Bray-Curtis distance)
#'
#' \code{colmap}: the colour map provided by the input argument
#' @examples
#'     fname <- system.file("Major.csv",package="provenance")
#'     Major <- read.compositional(fname)
#'     plot(PCA(Major))
#' @export
read.compositional <- function(fname,method=NULL,colmap='rainbow') {
    out <- list()
    out$name <- basename(substr(fname,1,nchar(fname)-4))
    class(out) <- "compositional"
    out$x <- utils::read.csv(fname,header=TRUE,row.names=1)
    if (is.null(method)){
        if (any(out$x==0)) { method <- "bray" }
        else { method <- "aitchison" }
    }
    out$method <- method
    out$colmap <- colmap
    if (any(out$x==0) & method=="aitchison"){
        stop(paste("This dataset contains zeros and is",
                   "incompatible with the 'aitchison' distance"))
    }
    return(out)
}

#' Read a .csv file with point-counting data
#'
#' Reads a data table containing point-counting data
#' (e.g. petrographic, heavy mineral, palaeontological or
#' palynological data)
#'
#' @param fname a string with the path to the .csv file
#' @param method either "chisq" (for the chi-square distance) or
#'     "bray" (for the Bray-Curtis distance)
#' @param colmap an optional string with the name of one of R's
#'     built-in colour palettes (e.g., heat.colors, terrain.colors,
#'     topo.colors, cm.colors), which are to be used for plotting the
#'     data.
#' @return an object of class \code{counts}, i.e. a list with the
#'     following items:
#' 
#' \code{x}: a data frame with the samples as rows and the categories
#' as columns
#'
#' \code{colmap}: the colour map provided by the input argument
#'
#' @examples
#'     fname <- system.file("HM.csv",package="provenance")
#'     Major <- read.counts(fname)
#'     #plot(PCA(HM))
#' @export
read.counts <- function(fname,method='chisq',colmap='rainbow'){
    out <- list()
    class(out) <- "counts"
    out$name <- basename(substr(fname,1,nchar(fname)-4))
    out$x <- utils::read.csv(fname,header=TRUE,row.names=1)
    out$method <- method
    out$colmap <- colmap
    return(out)
}

#' Read a .csv file with mineral and rock densities
#'
#' Reads a data table containing densities to be used for
#' hydraulic sorting corrections (minsorting and srd functions)
#'
#' @param fname a string with the path to the .csv file
#' @return a vector with mineral and rock densities
#' @examples
#' data(Namib,densities)
#' N8 <- subset(Namib$HM,select="N8")
#' distribution <- minsorting(N8,densities,phi=2,sigmaphi=1,medium="air",by=0.05)
#' plot(distribution)
#' @export
read.densities <- function(fname){
    return(utils::read.csv(fname,header=TRUE))
}

#' create a \code{data.frame} object
#'
#' Convert an object of class \code{compositional} to a
#' \code{data.frame} for use in the \code{robCompositions} package
#'
#' @param x an object of class \code{compositional}
#' @param ... optional arguments to be passed on to the generic function
#' @return a \code{data.frame}
#' @examples
#' data(Namib)
#' Major.frame <- as.data.frame(Namib$Major)
#' ## uncomment the next two lines to plot an error
#' ## ellipse using the robCompositions package:
#' # library(robCompositions)
#' # plot(pcaCoDa(Major.frame))
#' @export
as.data.frame.compositional <- function(x,...){
    nc <- ncol(as.matrix(x$x))
    if (nc==3) out <- data.frame(x$x[,c(2,3,1)],...)
    if (nc>3) out <- data.frame(x$x[,c(2,3,1,4:nc)],...)
    if (nc<3) out <- data.frame(x$x,...)
    return(out)
}

#' create a \code{data.frame} object
#'
#' Convert an object of class \code{counts} to a \code{data.frame} for
#' use in the \code{robCompositions} package
#'
#' @param x an object of class \code{counts}
#' @param ... optional arguments to be passed on to the generic function
#' @return a \code{data.frame}
#' @examples
#' data(Namib)
#' qfl <- ternary(Namib$PTHM,c('Q'),c('KF','P'),c('Lm','Lv','Ls'))
#' plot(qfl,type="QFL.dickinson")
#' qfl.frame <- as.data.frame(qfl)
#' ## uncomment the next two lines to plot an error
#' ## ellipse using the robCompositions package:
#' # library(robCompositions)
#' # pca <- pcaCoDa(qfl.frame)
#' # plot(pca,xlabs=rownames(qfl.frame))
#' @export
as.data.frame.counts <- function(x,...){
    nc <- ncol(as.matrix(x$x))
    if (nc==3) out <- data.frame(x$x[,c(2,3,1)],...)
    if (nc>3) out <- data.frame(x$x[,c(2,3,1,4:nc)],...)
    if (nc<3) out <- data.frame(x$x,...)
    return(out)
}

#' create an \code{acomp} object
#'
#' Convert an object of class \code{compositional} to an object of
#' class \code{acomp} for use in the \code{compositions} package
#'
#' @param x an object of class \code{compositional}
#' @return a \code{data.frame}
#' @examples
#' data(Namib)
#' qfl <- ternary(Namib$PT,c('Q'),c('KF','P'),c('Lm','Lv','Ls'))
#' plot(qfl,type="QFL.dickinson")
#' qfl.acomp <- as.acomp(qfl)
#' ## uncomment the next two lines to plot an error
#' ## ellipse using the compositions package: 
#' # library(compositions)
#' # ellipses(mean(qfl.acomp),var(qfl.acomp),r=2)
#' @export
as.acomp <- function(x){
    if (!(methods::is(x,"compositional")| methods::is(x,"counts"))){
        stop("not an object of class compositional or counts")
    }
    dat <- as.matrix(as.data.frame(x))
    out <- structure(dat)
    attributes(out) <- list(dim = dim(dat), dimnames=dimnames(dat), class="acomp")
    return(out)
}

as.compositional.matrix <- function(x,method=NULL,colmap='rainbow'){
    out <- list(x=NULL,method=method,colmap=colmap)
    class(out) <- "compositional"
    out$x <- as.matrix(x)
    nc <- ncol(out$x)
    if (nc==3) out$x <- out$x[,c(3,1,2)]
    if (nc>3) out$x <- out$x[,c(3,1,2,4:nc)]
    if (nc<3) out$x <- out$x
    return(out)
}

#' create a \code{compositional} object
#'
#' Convert an object of class \code{matrix}, \code{data.fram} or
#' \code{acomp} to an object of class \code{compositional}
#'
#' @param x an object of class \code{matrix}, \code{data.fram} or
#' \code{acomp}
#' @param method dissimilarity measure, either 'aitchison' for
#' Aitchison's CLR-distance or 'bray' for the Bray-Curtis distance.
#' @param colmap the colour map to be used in pie charts.
#' @return an object of class \code{compositional}
#' @examples
#' data(Namib)
#' PT.acomp <- as.acomp(Namib$PT)
#' PT.compositional <- as.compositional(PT.acomp)
#' print(Namib$PT$x - PT.compositional$x)
#' ## uncomment the following lines for an illustration of using this 
#' ## function to integrate the \code{provenance} package with \code{compositions}
#' # library(compositions)
#' # data(Glacial)
#' # a.glac <- acomp(Glacial)
#' # c.glac <- as.compositional(a.glac)
#' # summaryplot(c.glac,ncol=8)
#' @export
as.compositional <- function(x,method=NULL,colmap='rainbow'){
    if (methods::is(x,"acomp")){
        attr <- attributes(x)
        attributes(x) <- NULL
        x[!is.numeric(x)] <- NA
        y <- matrix(x,nrow=attr$dim[[1]],ncol=attr$dim[[2]],dimnames=attr$dimnames)
        print(print(attr$dim[[1]]))
        return(as.compositional.matrix(y))
    } else if (methods::is(x,"data.frame") | methods::is(x,"matrix")){
        y <- as.matrix(x)
        dimnames(y) <- dimnames(x)
        return(as.compositional.matrix(y,method,colmap))
    } else {
        stop(paste("cannot convert an object of class",class(x),
                   "into an object of class compositional"))
    }
}
