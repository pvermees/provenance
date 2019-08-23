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
#'     should be labeled
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
#'     pie-chart should be labeled
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
#' @param ... optional arguments of the \code{biplot} function
#' @examples
#' data(Namib)
#' plot(PCA(Namib$Major))
#' @seealso PCA
#' @method plot PCA
#' @export
plot.PCA <- function(x,...){
    stats::biplot(x,...)
}

#' Point-counting biplot
#'
#' Plot the results of a correspondence analysis as a biplot
#' @param x an object of class \code{CA}
#' @param ... optional arguments of the \code{biplot} function
#' @examples
#' data(Namib)
#' plot(CA(Namib$PT))
#' @seealso CA
#' @method plot CA
#' @export
plot.CA <- function(x,...){
    X <- x$rscore[, 1L:2]
    X <- X %*% diag(x$cor[1L:2])
    Y <- x$cscore[, 1L:2]
    Y <- Y %*% diag(x$cor[1L:2])
    stats::biplot(X,Y,xlab='Component 1',ylab='Component 2',...)
}

#' Plot an MDS configuration
#'
#' Plots the coordinates of a multidimensional scaling analysis as an
#' X-Y scatter plot or 'map' and, if x$classical = FALSE, a Shepard
#' plot.
#' 
#' @param x an object of class \code{MDS}
#' @param nnlines if TRUE, draws nearest neighbour lines
#' @param pch plot character (see ?plot for details). May be a vector.
#' @param pos position of the sample labels relative to the plot
#'     symbols if pch != NA
#' @param cex relative size of plot symbols (see ?par for details)
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
                     mar=rep(2,4),mgp=c(2,1,0),xpd=NA,...){
    k <- ncol(x$points)
    graphics::par(mfrow=c(k-1,k-1), oma=oma, mar=mar, mgp=mgp, xpd=xpd)
    for (i in 1:(k-1)){
        for (j in 2:k){
            if (i>=j){
                graphics::plot.new() # empty plot
            } else {
                xlab <- paste0('Dim ',i)
                ylab <- paste0('Dim ',j)
                graphics::plot(x$points[,c(i,j)], type='n', asp=1, xlab=xlab, ylab=ylab, ...)
                if (nnlines) { # draw lines between closest neighbours
                    if (is.na(pch)) pch=21
                    plotlines(x$points[,c(i,j)],x$diss)
                }
                graphics::points(x$points[,c(i,j)], pch=pch, cex=cex, col=col, bg=bg)
                graphics::text(x$points[,c(i,j)], labels=labels(x$diss), pos=pos, col=col, bg=bg)    
            }
        }
    }
    if (!x$classical){
        grDevices::dev.new()
        graphics::par(mfrow=c(k-1,k-1), oma=oma, mar=mar, mgp=mgp, xpd=xpd)
        for (i in 1:(k-1)){
            for (j in 2:k){
                if (i>=j){
                    graphics::plot.new() # empty plot
                } else {
                    ylab <- "Distance/Disparity"
                    if (k>2) ylab <- paste0(ylab,' (Dims ',i,' & ',j,')')
                    shep <- MASS::Shepard(x$diss, x$points[,c(i,j)])
                    graphics::plot(shep,pch=20,xlab="Dissimilarity",ylab=ylab)
                    graphics::lines(shep$x, shep$yf, type="S")
                    if (i==1 & j==2)
                        graphics::title(paste0("Stress = ",x$stress))    
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
#'     \code{KDEs}, or \code{distributional}
#' @param ncol the number of columns
#' @param pch (optional) symbol to be used to mark the sample points
#'     along the x-axis of the KDEs (if appropriate).
#' @return a summary plot of all the data comprised of KDEs for the
#'     datasets of class \code{KDEs}, pie charts for those of class
#'     \code{compositional} and histograms for those of class
#'     \code{distributional}.
#' @examples
#' data(Namib)
#' KDEs <- KDEs(Namib$DZ,0,3000)
#' summaryplot(KDEs,Namib$HM,Namib$PT,ncol=2)
#' @seealso KDEs
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
#'     (see ?par for details)
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
    graphics::plot(x$gspace,asp=asp,pch=pch,col=col,bg=bg,cex=cex,
                   xlab=xlab,ylab=ylab,xaxt=xaxt,yaxt=yaxt,...)
    graphics::text(x$gspace,labels=rownames(x$gspace),pos=pos,col=col,bg=bg,cex=cex)
    graphics::title('Group Configuration')
    X <- unlist(lapply(x$cweights,function(foo) foo[1,1]))
    Y <- unlist(lapply(x$cweights,function(foo) foo[2,2]))
    grDevices::dev.new()
    graphics::plot(X,Y,asp=1,pch=pch[1],cex=cex,...)
    graphics::text(X,Y,names(x$cweights),pos=pos,cex=cex)
    graphics::title('Source Weights')
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
    graphics::pie(comp,labels=labels,col=col,...)
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
        xtick <- xyz2xy(matrix(c(tick,1-tick-ticklength,ticklength,
                                 tick,1-tick,0,
                                 tick-ticklength,1-tick,ticklength),
                               ncol=3,byrow=TRUE))
        ytick <- xyz2xy(matrix(c(ticklength,tick,1-tick-ticklength,
                                 0,tick,1-tick,
                                 ticklength,tick-ticklength,1-tick),
                               ncol=3,byrow=TRUE))
        ztick <- xyz2xy(matrix(c(1-tick-ticklength,ticklength,tick,
                                 1-tick,0,tick,
                                 1-tick,ticklength,tick-ticklength),
                               ncol=3,byrow=TRUE))
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
        xline <- xyz2xy(matrix(c(tick,1-tick,0,
                                 tick,0,1-tick),
                               ncol=3,byrow=TRUE))
        yline <- xyz2xy(matrix(c(0,tick,1-tick,
                                 1-tick,tick,0),
                               ncol=3,byrow=TRUE))
        zline <- xyz2xy(matrix(c(1-tick,0,tick,
                                 0,1-tick,tick),
                               ncol=3,byrow=TRUE))
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
        xy1 <- xyz2xy(matrix(c(90,0,10,10,0,90),ncol=3))
        xy2 <- xyz2xy(matrix(c(90,0,0,90,10,10),ncol=3))
        xy3 <- xyz2xy(matrix(c(10,10,90,0,0,90),ncol=3))
        xy4 <- xyz2xy(matrix(c(50,10,50,10,0,80),ncol=3))
        xy5 <- xyz2xy(matrix(c(80,0,10,50,10,50),ncol=3))
        xy6 <- xyz2xy(matrix(c(10,50,80,0,10,50),ncol=3))
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
        xy1 <- xyz2xy(matrix(c(90,90,10,0,0,10),ncol=3))
        xy2 <- xyz2xy(matrix(c(75,75,25,0,0,25),ncol=3))
        xy3 <- xyz2xy(matrix(c(0,75,25,25/4,75,25*3/4),ncol=3))
        xy4 <- xyz2xy(matrix(c(0,75,75,25*3/4,25,25/4),ncol=3))
        xy5 <- xyz2xy(matrix(c(0,90,50,5,50,5),ncol=3))
        graphics::lines(xy1); graphics::lines(xy2);
        graphics::lines(xy3); graphics::lines(xy4); graphics::lines(xy5)
        graphics::text(xyz2xy(c(20,-1,2)),labels='quartarenite',adj=0)
        graphics::text(xyz2xy(c(40,-2,10)),labels='sublitharenite',adj=0)
        graphics::text(xyz2xy(c(40,10,-2)),labels='subarkose',adj=1)
        graphics::text(xyz2xy(c(30,70,10)),labels='arkose',srt=68)
        graphics::text(xyz2xy(c(30,50,30)),labels='lithic arkose',srt=85)
        graphics::text(xyz2xy(c(30,30,50)),labels='feldspathic litharenite',srt=-85)
        graphics::text(xyz2xy(c(30,10,70)),labels='litharenite',srt=-68)
        thelabels <- c('Q','F','L')
    } else if (type=='QFL.dickinson') {
        xy1 <- xyz2xy(matrix(c(97,0,0,85,3,15),ncol=3))
        xy2 <- xyz2xy(matrix(c(25,51.6,0,40,75,8.4),ncol=3))
        graphics::lines(xy1); graphics::lines(xy2)
        graphics::text(xyz2xy(c(20,40,40)),labels='magmatic arc')
        graphics::text(xyz2xy(c(55,15,30)),labels='recycled orogen')
        graphics::text(xyz2xy(c(50,45,5)),labels='continental block',srt=65)
        thelabels <- c('Q','F','L')
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

xyz2xy <- function(xyz){
    if (methods::is(xyz,"matrix")){
        n <- nrow(xyz)
        x <- xyz[,1]
        y <- xyz[,2]
        z <- xyz[,3]
    } else {
        n <- 1
        x <- xyz[1]
        y <- xyz[2]
        z <- xyz[3]
    }
    xy <- matrix(0,nrow=n,ncol=2)
    xy[,1] <- 0.5*(x+2*z)/(x+y+z)
    xy[,2] <- sin(pi/3)*x/(x+y+z)
    return(xy)
}

emptyplot <- function(){
    graphics::plot(c(0,1),c(0,1),type='n',axes=FALSE,xlab="",ylab="")
}

# save plot d as a pdf file with name f
saveplot <- function(f, d){
    grDevices::dev.set(d)
    grDevices::dev.copy2pdf(file=f)
}
