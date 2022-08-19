#' Plot a kernel density estimate
#'
#' Plots an object of class \code{KDE}
#' @param x an object of class \code{KDE}
#' @param pch the symbol used to show the samples. May be a vector.
#'     Set \code{pch = NA} to turn them off.
#' @param xlab the label of the x-axis
#' @param ylab the label of the y-axis
#' @param ... optional parameters to be passed on to the graphics
#'     object
#' @examples
#' data(Namib)
#' samp <- Namib$DZ$x[['N1']]
#' dens <- KDE(samp,from=0,to=3000)
#' plot(dens)
#' @seealso KDE
#' @method plot KDE
#' @export
plot.KDE <- function(x,pch='|',xlab="age [Ma]",ylab="",...){
    if (x$log) {
        graphics::plot(x$x,x$y,type='l',log="x",xlab=xlab,ylab=ylab,...)
    } else {
        graphics::plot(x$x,x$y,type='l',xlab=xlab,ylab=ylab,...)
    }
    graphics::points(x$ages,rep(graphics::par("usr")[3]/2,length(x$ages)),pch=pch)
    graphics::text(utils::tail(x$x,n=1),.9*max(x$y),paste0("n=",length(x$ages)),pos=2)
}

#' Plot one or more kernel density estimates
#'
#' Plots an object of class \code{KDEs}
#' @param x an object of class \code{KDEs}
#' @param sname optional sample name. If \code{sname=NA}, all samples
#'     are shown on a summary plot
#' @param annotate add a time axis?
#' @param pch symbol to be used to mark the sample points along the
#'     x-axis. Change to \code{NA} to omit.
#' @param ... optional parameters to be passed on to the
#'     \code{summaryplot} function
#' @examples
#' data(Namib)
#' kdes <- KDEs(Namib$DZ)
#' plot(kdes,ncol=2)
#' @seealso KDEs summaryplot
#' @method plot KDEs
#' @export
plot.KDEs <- function(x,sname=NA,annotate=TRUE,pch='|',...){
    if (is.na(sname)) {
        summaryplot(x,pch=pch,...)
    } else {
        if (x$themax>0){ # normalise
            M <- x$themax
        } else {
            M <- max(x$kdes[[sname]]$y)
        }
        if (annotate){
            plot.KDE(x$kdes[[sname]],pch=pch,ylim=c(0,M),...)
        } else {
            plot.KDE(x$kdes[[sname]],pch=pch,axes=FALSE,xlab="",ylab="",ylim=c(0,M),...)
        }
    }
}

#' Plot continuous data as histograms or cumulative age distributions
#'
#' Plot one or several samples from a \code{distributional} dataset as
#' a histogram or Cumulative Age Distributions (CAD).
#' @param x an object of class \code{distributional}
#' @param snames a string or a vector of string with the names of the
#'     samples that need plotting if \code{snames} is a vector, then
#'     the function will default to a CAD.
#' @param annotate boolean flag indicating whether the x- and y-axis
#'     should be labelled
#' @param CAD boolean flag indicating whether the data should be
#'     plotted as a cumulative age distribution or a histogram. For
#'     multi-sample plots, the function will override this value with
#'     \code{TRUE}.
#' @param pch an optional symbol to mark the sample points along the
#'     CAD
#' @param verticals boolean flag indicating if the horizontal lines of
#'     the CAD should be connected by vertical lines
#' @param colmap an optional string with the name of one of R's
#'     built-in colour palettes (e.g., heat.colors, terrain.colors,
#'     topo.colors, cm.colors), which are to be used for plotting the
#'     data.
#' @param ... optional arguments to the generic \code{plot} function
#' @examples
#' data(Namib)
#' plot(Namib$DZ,c('N1','N2'))
#' @method plot distributional
#' @export
plot.distributional <- function(x,snames=NULL,annotate=TRUE,CAD=FALSE,
                       pch=NA,verticals=TRUE,colmap=NULL,...){
    if (is.null(snames)) snames <- names(x)
    if (is.null(colmap)) { colmap <- x$colmap }
    n <- length(snames)
    col <- do.call(colmap,list(n))
    if (n>1){
        graphics::plot(stats::ecdf(x$x[[snames[1]]]),pch=pch,verticals=verticals,
             col=col[1],xlab=x$xlab,main="",...)
        for (i in 2:n){
            graphics::lines(stats::ecdf(x$x[[snames[i]]]),pch=pch,verticals=verticals,
                  col=col[i],xlab=x$xlab,main="",...)
        }
        graphics::legend("bottomright",legend=snames,lwd=1,col=col)
    } else {
        if (annotate){
            if (CAD) { graphics::plot(stats::ecdf(x$x[[snames]]),pch=pch,
                       verticals=verticals,col=col,main=snames,...) }
            else { graphics::hist(x$x[[snames]],x$breaks,col=col) }
        } else {
            if (CAD) { graphics::plot(stats::ecdf(x$x[[snames]]),pch=pch,verticals=verticals,
                       axes=FALSE,xlab="",ylab="",main="",col=col,...)
            } else { graphics::hist(x$x[[snames]],x$breaks, axes=FALSE,
                                    xlab="",ylab="",main="",col=col) }
        }
    }
}

#' Plot a pie chart
#'
#' Plots an object of class \code{compositional} as a pie chart
#' @param x an object of class \code{compositional}
#' @param sname the sample name
#' @param annotate a boolean flag controlling if the pies of the
#'     pie-chart should be labelled
#' @param colmap an optional string with the name of one of R's
#'     built-in colour palettes (e.g., heat.colors, terrain.colors,
#'     topo.colors, cm.colors), which are to be used for plotting the
#'     data.
#' @param ... optional parameters to be passed on to the graphics
#'     object
#' @examples
#' data(Namib)
#' plot(Namib$Major,'N1',colmap='heat.colors')
#' @method plot compositional
#' @export
plot.compositional <- function(x,sname,annotate=TRUE,colmap=NULL,...){
    i <- which(names(x) %in% sname)
    if (is.null(colmap)){ colmap <- x$colmap }
    col <- do.call(colmap,list(length(colnames(x$x))))
    if (annotate){
        graphics::pie(unlist(x$x[i,]),col=col,...)
    } else {
        graphics::pie(unlist(x$x[i,]),labels=NA,col=col,...)
    }
}
#' @method plot counts
#' @export
plot.counts <- function(x,sname,annotate=TRUE,colmap=NULL,...){
    plot.compositional(x,sname,annotate=annotate,colmap=colmap,...)
}

#' Plot a Procrustes configuration
#'
#' Plots the group configuration of a Generalised Procrustes Analysis
#'
#' @param x an object of class \code{GPA}
#' @param pch plot symbol
#' @param pos position of the sample labels relative to the plot
#'     symbols if pch != NA
#' @param col plot colour (may be a vector)
#' @param bg background colour (may be a vector)
#' @param cex relative size of plot symbols
#' @param ... optional arguments to the generic \code{plot} function
#' @examples
#' data(Namib)
#' GPA <- procrustes(Namib$DZ,Namib$HM)
#' coast <- c('N1','N2','N3','N10','N11','N12','T8','T13')
#' snames <- names(Namib$DZ)
#' bgcol <- rep('yellow',length(snames))
#' bgcol[which(snames %in% coast)] <- 'red'
#' plot(GPA,pch=21,bg=bgcol)
#' @seealso procrustes
#' @method plot GPA
#' @export
plot.GPA <- function(x,pch=NA,pos=NULL,col='black',bg='white',cex=1,...){
    graphics::plot(x$points[,1],x$points[,2],asp=1,pch=pch,col=col,bg=bg,cex=cex,...)
    if (!is.na(pch) & is.null(pos)) { pos <- 1 }
    graphics::text(x$points[,1],x$points[,2],x$labels,pos=pos,col=col,bg=bg,cex=cex)
}

#' Compositional biplot
#'
#' Plot the results of a principal components analysis as a biplot
#' @param x an object of class \code{PCA}
#' @param labelcol colour(s) of the sample labels (may be a vector).
#' @param vectorcol colour of the vector loadings for the variables
#' @param choices see the help pages of the generic \code{biplot}
#'     function.
#' @param scale see the help pages of the generic \code{biplot}
#'     function.
#' @param pc.biplot see the help pages of the generic \code{biplot}
#'     function.
#' @param ... optional arguments of the generic \code{biplot} function
#' @examples
#' data(Namib)
#' plot(PCA(Namib$Major))
#' @seealso PCA
#' @method plot PCA
#' @export
plot.PCA <- function(x, labelcol='black', vectorcol='red',
                     choices = 1L:2L, scale = 1, pc.biplot = FALSE, ...){
    if (length(choices) != 2L) 
        stop("length of choices must be 2")
    if (!length(scores <- x$x)) 
        stop(gettextf("object '%s' has no scores", deparse(substitute(x))), 
             domain = NA)
    if (is.complex(scores)) 
        stop("biplots are not defined for complex PCA")
    lam <- x$sdev[choices]
    n <- NROW(scores)
    lam <- lam * sqrt(n)
    if (scale < 0 || scale > 1) 
        warning("'scale' is outside [0, 1]")
    if (scale != 0) 
        lam <- lam^scale
    else lam <- 1
    if (pc.biplot) 
        lam <- lam/sqrt(n)
    biplotHelper(t(t(scores[,choices])/lam),
                 t(t(x$rotation[, choices])*lam),
                 labelcol=labelcol, vectorcol=vectorcol,...)
    invisible()
}

#' Point-counting biplot
#'
#' Plot the results of a correspondence analysis as a biplot
#' @param x an object of class \code{CA}
#' @param labelcol colour of the sample labels (may be a vector).
#' @param vectorcol colour of the vector loadings for the variables
#' @param ... optional arguments of the generic \code{biplot} function
#' @examples
#' data(Namib)
#' plot(CA(Namib$PT))
#' @seealso CA
#' @method plot CA
#' @export
plot.CA <- function(x,labelcol='black',vectorcol='red',...){
    X <- x$rscore[, 1L:2]
    X <- X %*% diag(x$cor[1L:2])
    Y <- x$cscore[, 1L:2]
    Y <- Y %*% diag(x$cor[1L:2])
    biplotHelper(X, Y, labelcol=labelcol, vectorcol=vectorcol,
                 xlab='Component 1', ylab='Component 2',...)
    invisible()
}

#' Plot an MDS configuration
#'
#' Plots the coordinates of a multidimensional scaling analysis as an
#' X-Y scatter plot or `map' and, if x$classical = FALSE, a Shepard
#' plot.
#' 
#' @param x an object of class \code{MDS}
#' @param nnlines if TRUE, draws nearest neighbour lines
#' @param pch plot character (see ?plot for details). May be a vector.
#' @param pos position of the sample labels relative to the plot
#'     symbols if \code{pch != NA}
#' @param cex relative size of plot symbols (see \code{?par} for
#'     details)
#' @param col plot colour (may be a vector)
#' @param bg background colour (may be a vector)
#' @param oma A vector of the form \code{c(bottom, left, top, right)}
#'     giving the size of the outer margins in lines of text.
#' @param mar A numerical vector of the form \code{c(bottom, left,
#'     top, right)} that gives the number of lines of margin to be
#'     specified on the four sides of the plot.
#' @param mgp The margin line (in \code{mex} units) for the axis
#'     title, axis labels and axis line. See \code{?par} for further
#'     details.
#' @param xpd A logical value or \code{NA}.  See \code{?par} for
#'     further details.
#' @param Shepard either:
#'
#' \code{0}: only plot the MDS configuration, do not show the Shepard plot
#'
#' \code{1}: only show the Shepard plot, do not plot the MDS configuration
#'
#' \code{2}: show both the MDS configuration and Shepard plot in separate
#' windows
#' 
#' @param ... optional arguments to the generic \code{plot} function
#' @seealso MDS
#' @method plot MDS
#' @examples
#' data(Namib)
#' mds <- MDS(Namib$DZ)
#' coast <- c('N1','N2','N3','N10','N11','N12','T8','T13')
#' snames <- names(Namib$DZ)
#' bgcol <- rep('yellow',length(snames))
#' bgcol[which(snames %in% coast)] <- 'red'
#' plot(mds,pch=21,bg=bgcol)
#' @export
plot.MDS <- function(x,nnlines=FALSE,pch=NA,pos=NULL,cex=1,
                     col='black',bg='white',oma=rep(1,4),
                     mar=rep(2,4),mgp=c(2,1,0),xpd=NA,Shepard=2,...){
    ns <- nrow(x$points)/(x$nb+1)
    k <- ncol(x$points)
    op <- graphics::par(mfrow=c(k-1,k-1), oma=oma, mar=mar, mgp=mgp, xpd=xpd)
    on.exit(graphics::par(op))
    if (!x$classical & Shepard>0){ # Shepard plot
        graphics::par(mfrow=c(k-1,k-1), oma=oma, mar=mar, mgp=mgp, xpd=xpd)
        for (i in 1:(k-1)){
            for (j in 2:k){
                if (i>=j){
                    graphics::plot.new() # empty plot
                } else {
                    ylab <- "Distance/Disparity"
                    if (k>2) ylab <- paste0(ylab,' (Dims ',i,' & ',j,')')
                    D <- as.matrix(x$diss)[1:ns,1:ns]
                    shep <- MASS::Shepard(stats::as.dist(D), x$points[1:ns,c(i,j)])
                    graphics::plot(shep,pch=20,xlab="Dissimilarity",ylab=ylab)
                    graphics::lines(shep$x, shep$yf, type="S")
                    if (i==1 & j==2)
                        graphics::title(paste0("Stress = ",x$stress))    
                }
            }
        }
        if (Shepard==1) return(invisible())
        if (!RStudio()) grDevices::dev.new()
    }
    for (i in 1:(k-1)){ # Configuration
        for (j in 2:k){
            if (i>=j){
                graphics::plot.new() # empty plot
            } else {
                xlab <- paste0('Dim ',i)
                ylab <- paste0('Dim ',j)
                graphics::plot(x$points[,c(i,j)], type='n',
                               asp=1, xlab=xlab, ylab=ylab, ...)
                if (nnlines) { # draw lines between closest neighbours
                    if (is.na(pch)) pch=21
                    plotlines(x$points[1:ns,c(i,j)],x$diss)
                }
                graphics::points(x$points[1:ns,c(i,j)],
                                 pch=pch, cex=cex, col=col, bg=bg)
                graphics::text(x$points[1:ns,c(i,j)], labels=labels(x$diss)[1:ns],
                               pos=pos, col=col, bg=bg)
                if (x$nb>0){ # bootstrap
                    for (l in 1:ns){
                        if (length(col)==ns) bcol <- col[l]
                        else bcol <- col
                        m <- ns + l + seq(from=1,to=ns*x$nb,by=ns) - 1
                        chi <- grDevices::chull(x$points[m,c(i,j)])
                        graphics::polygon(x$points[m[chi],c(i,j)],border=bcol)
                    }
                }
            }
        }
    }
}

#' Joint plot of several provenance datasets
#'
#' Arranges kernel density estimates and pie charts in a grid format
#' 
#' @param ... a sequence of datasets of class \code{compositional},
#'     \code{distributional}, \code{counts} or \code{KDEs}.
#' @param ncol the number of columns
#' @param pch (optional) symbol to be used to mark the sample points
#'     along the x-axis of the KDEs (if appropriate).
#' @return a summary plot of all the data comprised of KDEs for the
#'     datasets of class \code{KDEs}, pie charts for those of class
#'     \code{compositional} or \code{counts} and histograms for those
#'     of class \code{distributional}.
#' @examples
#' data(Namib)
#' KDEs <- KDEs(Namib$DZ,0,3000)
#' summaryplot(KDEs,Namib$HM,Namib$PT,ncol=2)
#' @seealso \link{KDEs}
#' @export
summaryplot <- function(...,ncol=1,pch=NA){
    oldpar <- graphics::par(no.readonly=T)
    dlist <- list(...)
    dnames <- get.data.names(dlist)
    names(dlist) <- dnames
    classes <- unlist(lapply(dlist,class))
    nd <- length(dlist) # number of datasets
    snames <- unique(unlist(lapply(dlist,names)))
    ns <- length(snames)
    w <- rep(1,nd) # column widths
    idist <- which( classes %in% c("KDEs","distributional") )
    w[idist] <- 2
    w <- rep(c(1,w),ncol)
    nppc <- ceiling(ns/ncol)
    np <- (nppc+1)*ncol*(nd+1) # number of subpanels
    graphics::layout(matrix(1:np,nppc+1,length(w)),w,rep(1,nppc+1))
    si <- ceiling(seq(from=0,to=ns,length.out=ncol+1)) # sample index
    graphics::par(xpd=TRUE,mar=c(0,0,0,0))#, mfcol=c(nppc,nd))
    for (i in 1:ncol){ # loop through columns
        for (j in (si[i]+1):(si[i+1])){
            sname <- snames[j]
            emptyplot()
            graphics::text(1,0.5,labels=sname,pos=2)
        }
        if (si[i+1]-si[i]<nppc) emptyplot()
        emptyplot()
        for (k in 1:nd){ # loop through datasets
            d <- dlist[[k]]
            for (j in (si[i]+1):si[i+1]){
                sname <- snames[j]
                if (sname %in% names(d)){
                    if (k %in% idist) graphics::plot(d,sname,annotate=FALSE,pch=pch)
                    else graphics::plot(d,sname,annotate=FALSE)
                } else {
                    emptyplot()
                }
            }
            if (si[i+1]-si[i]<nppc) emptyplot()
            if (i==ncol){
                ds <- grDevices::dev.size()[2]/(nppc+1)
                annotation(d,height=ds)
                if (nd>2) graphics::title(d$name,line=-1)
            } else {
                emptyplot()
            }
        }
    }
    graphics::par(oldpar)
}

#' Plot inferred grain size distributions
#'
#' Plot the grain size distributions of the different minerals under
#' consideration
#' @param x an object of class \code{minsorting}
#' @param cumulative boolean flag indicating whether the grain size
#'     distribution should be plotted as a density or cumulative
#'     probability curve.
#' @param components string or list of strings with the names of a
#'     subcomposition that needs plotting
#' @param ... optional parameters to be passed on to graphics::matplot
#'     (see \code{?par} for details)
#' @examples
#' data(endmembers,densities)
#' OPH <- subset(endmembers,select="ophiolite")
#' distribution <- minsorting(OPH,densities,phi=2,sigmaphi=1,
#'                            medium="air",by=0.05)
#' plot(distribution,components=c('F','px','opaques'))
#' @seealso minsorting
#' @method plot minsorting
#' @export
plot.minsorting <- function(x,cumulative=FALSE,components=NULL,...){
    if (is.null(components)){
        plotval <- x$mfract
    } else {
        plotval <- x$mfract[,components,drop=FALSE]
        colnames(plotval) <- components
    }
    gsize <- as.numeric(rownames(plotval)) # grain size fraction
    if (cumulative) plotval <- apply(plotval,2,'cumsum')
    ncat <- ncol(plotval)
    graphics::matplot(gsize,plotval, type="l", xlab="phi", ylab="%",
                      lty=1, col=grDevices::rainbow(ncat),...)
    graphics::legend("right", inset=.01, legend=colnames(plotval),
                     pch="-", col=grDevices::rainbow(ncat), horiz=FALSE, cex=0.7)
}

#' Plot an INDSCAL group configuration and source weights
#'
#' Given an object of class \code{INDSCAL}, generates two plots: the
#' group configuration and the subject weights. Together, these
#' describe a 3-way MDS model.
#' @param x an object of class \code{INDSCAL}
#' @param asp the aspect ratio of the plot
#' @param pch plot symbol (may be a vector)
#' @param pos position of the sample labels relative to the plot
#'     symbols if pch != NA
#' @param col plot colour (may be a vector)
#' @param bg background colour (may be a vector)
#' @param cex relative size of plot symbols
#' @param xlab a string with the label of the x axis
#' @param ylab a string with the label of the y axis
#' @param xaxt if = 'y', adds ticks to the x axis
#' @param yaxt if = 'y', adds ticks to the y axis
#' @param ... optional arguments to the generic plot function
#' @examples
#' data(Namib)
#' coast <- c('N1','N2','N3','N10','N11','N12','T8','T13')
#' snames <- names(Namib$DZ)
#' pch <- rep(21,length(snames))
#' pch[which(snames %in% coast)] <- 22
#' plot(indscal(Namib$DZ,Namib$HM),pch=pch)
#' @seealso indscal
#' @method plot INDSCAL
#' @export
plot.INDSCAL <- function(x,asp=1,pch=NA,pos=NULL,col='black',
                         bg='white',cex=1,xlab="X",ylab="Y",
                         xaxt='n',yaxt='n',...){
    if (!any(is.na(pch)) && !any(is.null(pos))) { pos <- 1 }
    X <- unlist(lapply(x$cweights,function(foo) foo[1,1]))
    Y <- unlist(lapply(x$cweights,function(foo) foo[2,2]))
    graphics::plot(X,Y,asp=1,pch=pch[1],cex=cex,...)
    graphics::text(X,Y,names(x$cweights),pos=pos,cex=cex)
    graphics::title('Source Weights')
    if (!RStudio()) grDevices::dev.new()    
    graphics::plot(x$gspace,asp=asp,pch=pch,col=col,bg=bg,cex=cex,
                   xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,...)
    graphics::text(x$gspace,labels=rownames(x$gspace),pos=pos,col=col,bg=bg,cex=cex)
    graphics::title('Group Configuration')
}

# a function to plot the nearest neighbour lines
plotlines <- function(conf,diss) {
    # rank the samples according to their pairwise proximity
    i = t(apply(as.matrix(diss),1,function(x) order(x))[2:3,])
    # coordinates for the lines
    x1 = as.vector(conf[i[,1],1]) # calculate (x,y)-coordinates ...
    y1 = as.vector(conf[i[,1],2]) # ... of nearest neighbours
    x2 = as.vector(conf[i[,2],1]) # calculate (x,y)-coordinates ...
    y2 = as.vector(conf[i[,2],2]) # ... of second nearest neighbours
    for (j in 1:nrow(conf)) {
        graphics::lines(c(conf[j,1],x1[j]),c(conf[j,2],y1[j]),lty=1) # solid line
        graphics::lines(c(conf[j,1],x2[j]),c(conf[j,2],y2[j]),lty=2) # dashed line
    }
}

# annotation of various plots used in summaryplot function
annotation <- function(x,...){ UseMethod("annotation",x) }
annotation.default <- function(x,...){
    stop('x not of class KDEs, compositional or distributional in annotation(x)')
}
annotation.KDEs <- function(x,height=NULL,...){
    oldpar <- graphics::par()
    if (is.null(height)){ graphics::par(mar=c(2,0,0,0)) }
    else { graphics::par(mai=c(height/2,0,0,0)) }
    if (x$log) {
        graphics::plot(c(x$from,x$to),c(0,1),type='n',log='x',
                       axes=FALSE,xlab="",ylab="",...)
    } else {
        graphics::plot(c(x$from,x$to),c(0,1),type='n',
                       axes=FALSE,xlab="",ylab="",...)
    }
    graphics::Axis(side=1)
    if (x$log){
        middle <- sqrt(x$from*x$to)
    } else {
        middle <- .5*(x$from+x$to)
    }
    graphics::text(x=middle,.5,label=x$xlab)
    graphics::par(mai=oldpar$mai)
    graphics::par(mar=oldpar$mar)
}
annotation.compositional <- function(x,height=NULL,...){
    labels <- colnames(x$x)
    comp <- rep(1,length(labels))
    col <- do.call(x$colmap,list(length(labels)))
    graphics::pie(comp,labels=labels,col=col,xpd=NA,...)
}
annotation.counts <- function(x,height=NULL,...){
    annotation.compositional(x,height=NULL,...)
}
annotation.distributional <- function(x,height=NULL,...){
    oldpar <- graphics::par()
    if (is.null(height)){ graphics::par(mar=c(2,0,0,0)) }
    else { graphics::par(mai=c(height/2,0,0,0)) }
    m <- min(x$breaks)
    M <- max(x$breaks)
    graphics::plot(c(m,M),c(0,1),type='n',axes=FALSE,xlab="",ylab="",...)
    graphics::Axis(side=1)
    graphics::text(x=.5*(m+M),.5,label=x$xlab)
    graphics::par(mai=oldpar$mai)
    graphics::par(mar=oldpar$mar)
}

ternary.ticks <- function(ticks=seq(0,1,0.25),ticklength=0.02){
    for (tick in ticks){
        xtick <- xyz2xy(c(tick,1-tick-ticklength,ticklength),
                        c(tick,1-tick,0),
                        c(tick-ticklength,1-tick,ticklength))
        ytick <- xyz2xy(c(ticklength,tick,1-tick-ticklength),
                        c(0,tick,1-tick),
                        c(ticklength,tick-ticklength,1-tick))
        ztick <- xyz2xy(c(1-tick-ticklength,ticklength,tick),
                        c(1-tick,0,tick),
                        c(1-tick,ticklength,tick-ticklength))
        if (tick>0 & tick<1){
            graphics::lines(xtick)
            graphics::lines(ytick)
            graphics::lines(ztick)
        }
    }
}
ternary.grid <- function(ticks=seq(0,1,0.25),
                         col='cornflowerblue',lty=2,lwd=1){
    oldpar <- graphics::par('col','lty','lwd')
    graphics::par(col=col,lty=lty,lwd=lwd)
    for (tick in ticks){
        xline <- xyz2xy(c(tick,1-tick,0),c(tick,0,1-tick))
        yline <- xyz2xy(c(0,tick,1-tick),c(1-tick,tick,0))
        zline <- xyz2xy(c(1-tick,0,tick),c(0,1-tick,tick))
        if (tick>0 & tick<1){
            graphics::lines(xline)
            graphics::lines(yline)
            graphics::lines(zline)
        }
    }
    graphics::par(oldpar)
}
ternary.lines <- function(type='empty',col='cornflowerblue',
                          ticks=seq(0,1,0.25),lty=2,lwd=1){
    oldpar <- graphics::par('col','lty','lwd')
    graphics::par(col=col,lty=lty,lwd=lwd)
    thelabels <- c('x','y','z')
    if (type=='QFL.descriptive'){
        xy1 <- xyz2xy(c(90,10,0),c(0,10,90))
        xy2 <- xyz2xy(c(90,0,10),c(0,90,10))
        xy3 <- xyz2xy(c(10,90,0),c(10,0,90))
        xy4 <- xyz2xy(c(50,50,0),c(10,10,80))
        xy5 <- xyz2xy(c(80,10,10),c(0,50,50))
        xy6 <- xyz2xy(c(10,80,10),c(50,0,50))
        graphics::lines(xy1); graphics::lines(xy2);
        graphics::lines(xy3); graphics::lines(xy4);
        graphics::lines(xy5); graphics::lines(xy6);
        graphics::text(xyz2xy(c(50,32,18)),labels=expression(paste('lF',bold('Q'))))
        graphics::text(xyz2xy(c(50,18,32)),labels=expression(paste('fL',bold('Q'))))
        graphics::text(xyz2xy(c(32,50,18)),labels=expression(paste('lQ',bold('F'))))
        graphics::text(xyz2xy(c(32,18,50)),labels=expression(paste('fQ',bold('L'))))
        graphics::text(xyz2xy(c(18,50,32)),labels=expression(paste('qL',bold('F'))))
        graphics::text(xyz2xy(c(18,32,50)),labels=expression(paste('qF',bold('L'))))
        graphics::text(xyz2xy(c(65,30,5)),
                       labels=expression(paste('feldspatho-',bold('quartzose'))),srt=60)
        graphics::text(xyz2xy(c(30,65,5)),
                       labels=expression(paste('quartzo-',bold('feldspathic'))),srt=60)
        graphics::text(xyz2xy(c(65,5,30)),
                       labels=expression(paste('litho-',bold('quartzose'))),srt=-60)
        graphics::text(xyz2xy(c(30,5,65)),
                       labels=expression(paste('quartzo-',bold('lithic'))),srt=-60)
        graphics::text(xyz2xy(c(5,30,65)),
                       labels=expression(paste('feldspatho-',bold('lithic'))))
        graphics::text(xyz2xy(c(5,65,30)),
                       labels=expression(paste('litho-',bold('feldspathic'))))
        thelabels <- c('Q','F','L')
    } else if (type=='QFL.folk'){
        xy1 <- xyz2xy(c(90,10,0),c(90,0,10))
        xy2 <- xyz2xy(c(75,25,0),c(75,0,25))
        xy3 <- xyz2xy(c(0,25,75),c(75,25/4,25*3/4))
        xy4 <- xyz2xy(c(0,75,25),c(75,25*3/4,25/4))
        xy5 <- xyz2xy(c(0,50,50),c(90,5,5))
        graphics::lines(xy1); graphics::lines(xy2);
        graphics::lines(xy3); graphics::lines(xy4); graphics::lines(xy5)
        graphics::text(xyz2xy(c(20,-1,2)),labels='quartzarenite',adj=0)
        graphics::text(xyz2xy(c(40,-2,10)),labels='sublitharenite',adj=0)
        graphics::text(xyz2xy(c(40,10,-2)),labels='subarkose',adj=1)
        graphics::text(xyz2xy(c(30,70,10)),labels='arkose',srt=68)
        graphics::text(xyz2xy(c(30,50,30)),labels='lithic arkose',srt=85)
        graphics::text(xyz2xy(c(30,30,50)),labels='feldspathic litharenite',srt=-85)
        graphics::text(xyz2xy(c(30,10,70)),labels='litharenite',srt=-68)
        thelabels <- c('Q','F','L')
    } else if (type=='QFL.dickinson') {
        xy1 <- xyz2xy(c(97,0,3),c(0,85,15))
        xy2 <- xyz2xy(c(25,0,75),c(51.6,40,8.4))
        graphics::lines(xy1); graphics::lines(xy2)
        graphics::text(xyz2xy(c(20,40,40)),labels='magmatic arc')
        graphics::text(xyz2xy(c(55,15,30)),labels='recycled orogen')
        graphics::text(xyz2xy(c(50,45,5)),labels='continental block',srt=65)
        thelabels <- c('Q','F','L')
    } else if (type=='QmFLt.dickinson') {
        xy1 <- xyz2xy(c(89,0,11),c(0,77,23))
        xyi1 <- ternary_intersection(c(89,0,11),c(0,77,23),
                                     c(80,20,0),c(0,13,87))
        xy2 <- rbind(xyi1,xyz2xy(c(0,13,87),test=TRUE))
        graphics::lines(xy1); graphics::lines(xy2)
        graphics::text(xyz2xy(c(20,40,40)),labels='magmatic arc')
        graphics::text(xyz2xy(c(42,11,47)),labels='recycled orogen',srt=-58)
        graphics::text(xyz2xy(c(45,47,8)),labels='continental block',srt=64)
        thelabels <- c('Qm','F','Lt')
    }
    graphics::par(oldpar)
    return(thelabels)
}

# X is an object of class 'ternary'
plotpath <- function(X){
    for (sname in names(X$restoration)){
        Y <- X$restoration[[sname]]
        xy <- xyz2xy(Y)
        graphics::lines(xy[,1],xy[,2])
        graphics::points(utils::tail(xy[,1],n=1),utils::tail(xy[,2],n=1))
    }
}

xyz2xy <- function(...,test=FALSE){
    xyz <- rbind(...)
    x <- 0.5*(xyz[,1]+2*xyz[,3])/rowSums(xyz)
    y <- sin(pi/3)*xyz[,1]/rowSums(xyz)
    cbind(x,y)
}

ternary_intersection <- function(xyz11,xyz12,xyz21,xyz22){
    xyz1 <- rbind(xyz11,xyz12)
    xyz2 <- rbind(xyz21,xyz22)
    xy1 <- xyz2xy(xyz1)
    xy2 <- xyz2xy(xyz2)
    b1 <- (xy1[2,'y']-xy1[1,'y'])/(xy1[2,'x']-xy1[1,'x'])
    b2 <- (xy2[2,'y']-xy2[1,'y'])/(xy2[2,'x']-xy2[1,'x'])
    xi <- (xy2[1,'y']-xy1[1,'y']+b1*xy1[1,'x']-b2*xy2[1,'x'])/(b1-b2)
    yi <- xy1[1,'y']+b1*(xi-xy1[1,'x'])
    c(xi,yi)
}

emptyplot <- function(){
    graphics::plot(c(0,1),c(0,1),type='n',axes=FALSE,xlab="",ylab="")
}

# save plot d as a pdf file with name f
saveplot <- function(f, d){
    grDevices::dev.set(d)
    grDevices::dev.copy2pdf(file=f)
}

# modified from stats::biplot.default
biplotHelper <- function(x, y, var.axes = TRUE, labelcol='black', vectorcol='red',
                         cex = rep(graphics::par("cex"), 2), pch=NULL,
                         xlabs = NULL, ylabs = NULL, expand = 1,
                         xlim = NULL, ylim = NULL,
                         arrow.len = 0.1, main = NULL, sub = NULL,
                         xlab = NULL, ylab = NULL, ...){
    n <- nrow(x)
    p <- nrow(y)
    if (missing(xlabs)) {
        xlabs <- dimnames(x)[[1L]]
        if (is.null(xlabs)) 
            xlabs <- 1L:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
    if (missing(ylabs)) {
        ylabs <- dimnames(y)[[1L]]
        if (is.null(ylabs)) 
            ylabs <- paste("Var", 1L:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2L]])
    if (length(cex) == 1L) 
        cex <- c(cex, cex)
    unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)),
                                    abs(max(x, na.rm = TRUE)))
    rangx1 <- unsigned.range(x[, 1L])
    rangx2 <- unsigned.range(x[, 2L])
    rangy1 <- unsigned.range(y[, 1L])
    rangy2 <- unsigned.range(y[, 2L])
    if (missing(xlim) && missing(ylim)) 
        xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if (missing(xlim)) 
        xlim <- rangx1
    else if (missing(ylim)) 
        ylim <- rangx2
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(graphics::par(op))
    op <- graphics::par(pty = "s")
    if (!is.null(main)) 
        op <- c(op, graphics::par(mar = graphics::par("mar") + c(0, 0, 1, 0)))
    graphics::plot(y, axes = FALSE, type = "n", xlim = xlim * ratio,
                   ylim = ylim * ratio, xlab = "", ylab = "", ...)
    graphics::axis(3, col = vectorcol, ...)
    graphics::axis(4, col = vectorcol, ...)
    graphics::box(col = 'black')
    graphics::text(y, labels = ylabs, cex = cex[2L], col = vectorcol, ...)
    if (var.axes) 
        graphics::arrows(0, 0, y[, 1L] * 0.8, y[, 2L] * 0.8,
                         col = vectorcol, length = arrow.len)
    graphics::par(new = TRUE)
    grDevices::dev.hold()
    on.exit(grDevices::dev.flush(), add = TRUE)
    graphics::plot(x, type = "n", xlim = xlim, ylim = ylim, col = labelcol, 
                   xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    if (is.null(pch))
        graphics::text(x, labels=xlabs, cex = cex[1L], col = labelcol, ...)
    else
        graphics::points(x, pch=pch, cex = cex[1L], bg = labelcol, ...)
    invisible()
}
