#' Define a ternary composition
#'
#' Create an object of class \code{ternary}
#' @param X an object of class \code{compositional} OR a matrix or
#'     data frame with numerical data
#' @param x string or a vector of strings indicating the variables
#'     making up the first subcomposition of the ternary system. If
#'     omitted, the first component of X is used instead.
#' @param y second (set of) variables
#' @param z third (set of) variables
#' @return an object of class \code{ternary}, i.e. a list containing:
#'
#' x: a three column matrix (or vector) of ternary compositions.
#'
#' and (if X is of class \code{SRDcorrected})
#'
#' restoration: a list of intermediate ternary compositions inherited
#' from the SRD correction
#'
#' @seealso restore
#' @examples
#' data(Namib)
#' tern <- ternary(Namib$PT,c('Q'),c('KF','P'),c('Lm','Lv','Ls'))
#' plot(tern,type="QFL")
#' @export
ternary <- function(X,x=NA,y=NA,z=NA){
    out <- list()
    if (any(class(X) %in% c("compositional","counts"))) dat <- X$x
    else dat <- X
    if (ndim(dat)>1) cnames <- colnames(dat)
    else cnames <- names(dat)
    hasnames <-  (!is.null(cnames))
    if (all(is.na(x)) & hasnames) x <- cnames[1] else x <- 1
    if (all(is.na(y)) & hasnames) y <- cnames[2] else y <- 2
    if (all(is.na(z)) & hasnames) z <- cnames[3] else z <- 3
    out$raw <- cbind(sumcols(dat,x),sumcols(dat,y),sumcols(dat,z))
    colnames(out$raw) <- c(x,y,z)
    class(out) <- append("ternary",class(X))
    arg <- deparse(substitute(x))
    out$x <- ternaryclosure(out$raw,x,y,z)
    if (methods::is(X,"SRDcorrected")){
        out$restoration <- list()
        snames <- names(X$restoration)
        for (sname in snames){
            out$restoration[[sname]] <-
                ternaryclosure(X$restoration[[sname]],x,y,z)
        }
        class(out) <- append("SRDcorrected",class(out))
    }
    return(out)
}

# X is a matrix or vector
# x, y, z is an index or (vector) of string(s)
ternaryclosure <- function(X,x=1,y=2,z=3){
    xlab <- sumlabels(X,x)
    ylab <- sumlabels(X,y)
    zlab <- sumlabels(X,z)
    out <- cbind(sumcols(X,x),sumcols(X,y),sumcols(X,z))
    den <- rowSums(out)
    out <- apply(out,2,'/',den)
    if (methods::is(out,"matrix")) {
        colnames(out) <- c(xlab,ylab,zlab)
    } else {
        names(out) <- c(xlab,ylab,zlab)
    }
    return(out)
}

#' Plot a ternary diagram
#'
#' Plots triplets of compositional data on a ternary diagram
#' @param x an object of class \code{ternary}, or a three-column data
#'     frame or matrix
#' @param type adds annotations to the ternary diagram, one of either
#'     \code{empty}, \code{grid}, \code{QFL.descriptive},
#'     \code{QFL.folk} or \code{QFL.dickinson}
#' @param pch plot character, see \code{?par} for details (may be a
#'     vector)
#' @param pos position of the sample labels relative to the plot
#'     symbols if pch != NA
#' @param labels vector of strings to be added to the plot symbols
#' @param showpath if \code{x} has class \code{SRDcorrected}, and
#'     \code{showpath}==TRUE, the intermediate values of the SRD
#'     correction will be plotted on the ternary diagram as well as
#'     the final composition
#' @param col colour to be used for the background lines (if
#'     applicable)
#' @param ticks vector of tick values between 0 and 1
#' @param ticklength number between 0 and 1 to mark the length of the
#'     ticks
#' @param lty line type for the annotations (see \code{type})
#' @param lwd line thickness for the annotations
#' @param bg background colour for the plot symbols (may be a vector)
#' @param ... optional arguments to the generic \code{points} function
#' @examples
#' data(Namib)
#' tern <- ternary(Namib$PT,'Q',c('KF','P'),c('Lm','Lv','Ls'))
#' plot(tern,type='QFL.descriptive',pch=21,bg='red',labels=NULL)
#' @seealso ternary
#' @method plot ternary
#' @export
plot.ternary <- function(x,type='grid',pch=NA,pos=NULL,
                         labels=names(x),showpath=FALSE,bg=NA,
                         col='cornflowerblue',ticks=seq(0,1,0.25),
                         ticklength=0.02,lty=2,lwd=1,...){
    graphics::plot(c(0,1),c(0,1),type='n',xaxt='n',yaxt='n',
                   xlab='',ylab='',asp=1,bty='n',...)
    ternary.ticks(ticks=ticks,ticklength=ticklength)
    if (type=='empty'){
        cornerlabels <- colnames(x$x)
    } else if (type=='grid'){
        cornerlabels <- colnames(x$x)
        ternary.grid(ticks=ticks,col=col,lty=lty,lwd=lwd)
    } else {
        cornerlabels <- ternary.lines(type=type,col=col,lty=lty,lwd=lwd)
    }
    corners <- xyz2xy(matrix(c(1,0,0,1,0,1,0,0,0,0,1,0),ncol=3))
    graphics::lines(corners)
    graphics::text(corners[1:3,],labels=cornerlabels,pos=c(3,1,1))
    xy <- xyz2xy(x$x)
    if (is.null(pch)) return()
    if (is.na(pch) && is.null(labels)){ pch <- 1 }
    if (!is.na(pch) && is.null(pos)){ pos <- 1 }
    if (!is.na(pch)) graphics::points(xy,pch=pch,bg=bg,...)
    if (!is.null(labels)){ graphics::text(xy,labels=labels,pos=pos) }
    if (showpath & methods::is(x,'SRDcorrected')) plotpath(x)
}
#' Ternary point plotting
#'
#' Add points to an existing ternary diagram
#' @param x an object of class \code{ternary}, or a three-column data
#'     frame or matrix
#' @param ... optional arguments to the generic \code{points} function
#' @aliases lines text
#' @examples
#' tern <- ternary(Namib$PT,'Q',c('KF','P'),c('Lm','Lv','Ls'))
#' plot(tern,pch=21,bg='red',labels=NULL)
#' # add the geometric mean composition as a yellow square:
#' gmean <- ternary(exp(colMeans(log(tern$x))))
#' points(gmean,pch=22,bg='yellow')
#' @method points ternary
#' @export
points.ternary <- function(x,...){
    xy <- xyz2xy(x$x)
    graphics::points(xy,...)
}
#' Ternary line plotting
#' 
#' Add lines to an existing ternary diagram
#' @param x an object of class \code{ternary}, or a three-column data
#'     frame or matrix
#' @param ... optional arguments to the generic \code{lines} function
#' @examples
#' tern <- ternary(Namib$PT,'Q',c('KF','P'),c('Lm','Lv','Ls'))
#' plot(tern,pch=21,bg='red',labels=NULL)
#' middle <- matrix(c(0.01,0.49,0.01,0.49,0.98,0.02),2,3)
#' lines(ternary(middle))
#' @method lines ternary
#' @export
lines.ternary <- function(x,...){
    xy <- xyz2xy(x$x)
    graphics::lines(xy,...)
}
#' Ternary text plotting
#'
#' Add text an existing ternary diagram
#' @param x an object of class \code{ternary}, or a three-column data
#'     frame or matrix
#' @param labels a character vector or expression specifying the text
#'     to be written
#' @param ... optional arguments to the generic \code{text} function
#' @examples
#' tern <- ternary(Namib$PT,'Q',c('KF','P'),c('Lm','Lv','Ls'))
#' plot(tern,pch=21,bg='red',labels=NULL)
#' # add the geometric mean composition as a text label:
#' gmean <- ternary(exp(colMeans(log(tern$x))))
#' text(gmean,labels='geometric mean')
#' @method text ternary
#' @export
text.ternary <- function(x,labels=1:nrow(x$x),...){
    xy <- xyz2xy(x$x)
    graphics::text(xy,labels=labels,...)
}
#' Ternary confidence ellipse
#' 
#' plot a logistic \eqn{100(1-\alpha)\%} confidence region around the
#'     data or around its mean.
#' @param x an object of class \code{ternary}
#' @rdname ternary.ellipse
#' @export
ternary.ellipse <- function(x,...){ UseMethod("ternary.ellipse",x) }
#' @param alpha cutoff level for the confidence ellipse
#' @param population show the standard deviation of the entire
#'     population or the standard error of the mean?
#' @param ... optional formatting arguments
#' @rdname ternary.ellipse
#' @export
ternary.ellipse.default <- function(x,alpha=0.05,population=TRUE,...){
    uv <- ALR(x)
    u <- subset(uv,select=1)
    v <- subset(uv,select=2)
    n <- length(u)
    m <- 2
    df1 <- m
    df2 <- n-m
    if (population) k <- n+1
    else k <- 1
    hk <- k*(n-1)*stats::qf(1-alpha,df1,df2)/(n*(n-m))
    S <- stats::cov(uv)
    VW2VT <- svd(S)
    V <- VW2VT$u
    W2 <- diag(VW2VT$d)
    w1 <- sqrt(VW2VT$d[1])
    w2 <- sqrt(VW2VT$d[2])
    YbarY <- rbind(mean(u)-u,mean(v)-v)
    res <- 100
    g1 <- w1*sqrt(hk)*cos(seq(from=0,to=2*pi,length.out=res))
    g2 <- w2*sqrt(hk)*sin(seq(from=0,to=2*pi,length.out=res))
    G <- rbind(g1,g2)
    ell <- list()
    class(ell) <- 'compositional'
    ell$x <- t(V%*%G + rbind(rep(mean(u),res),rep(mean(v),res)))
    XYZ <- ALR(ell,inverse=TRUE)
    xy <- xyz2xy(XYZ$x)
    graphics::lines(xy,...)    
}
#' @rdname ternary.ellipse
#' @export
ternary.ellipse.compositional <- function(x,alpha=0.05,population=TRUE,...){
    ternary.ellipse.default(x,alpha=alpha,population=population,...)
}
ternary.ellipse.counts <- function(x,alpha=0.05,population=TRUE,...){
    fit13 <- central(x,components=colnames(x$raw)[c(1,3)])
    fit23 <- central(x,components=colnames(x$raw)[c(2,3)])
    pars <- c(log(fit13$ratio),log(fit23$ratio),
              fit13$sigma^2,fit23$sigma^2,0)
    if (fit13$sigma>0.01 & fit23$sigma>0.01){
        sol <- optimise(LL.ternary.random.effects,interval=c(-0.99,0.99),
                        dat=x$raw,pars=pars,tol=0.05)
        b1 <- pars[1]
        b2 <- pars[2]
        E <- matrix(0,2,2)
        E[1,1] <- pars[3]
        E[2,2] <- pars[4]
        E[1,2] <- sol$minimum*sqrt(E[1,1]*E[2,2])
        E[2,1] <- E[1,2]
        ell <- IsoplotR::ellipse(b1,b2,E,alpha=alpha)
        XYZ <- ALR(ell,inverse=TRUE)
        xy <- xyz2xy(XYZ$x)
        graphics::lines(xy,...)
        pars[5] <- E[1,2]
    } else if (fit13$sigma>0.01 & fit23$sigma<0.01){
        pars[4] <- 0
        m13 <- qnorm(alpha/2,mean=log(fit13$ratio),sd=fit13$sigma)
        M13 <- qnorm(1-alpha/2,mean=log(fit13$ratio),sd=fit13$sigma)
        np <- 20 # number of points
        u <- seq(m13,M13,length.out=np)
        v <- rep(log(fit23$ratio),np)
        uv <- cbind(u,v)
        XYZ <- ALR(uv,inverse=TRUE)
        xy <- xyz2xy(XYZ$x)
        graphics::lines(xy,...)
    } else if (fit13$sigma<0.01 & fit23$sigma>0.01){
        pars[3] <- 0
        m23 <- qnorm(alpha/2,mean=log(fit23$ratio),sd=fit23$sigma)
        M23 <- qnorm(1-alpha/2,mean=log(fit23$ratio),sd=fit23$sigma)
        np <- 20 # number of points
        u <- rep(log(fit13$ratio),np)
        v <- seq(m23,M23,length.out=np)
        uv <- cbind(u,v)
        XYZ <- ALR(uv,inverse=TRUE)
        xy <- xyz2xy(XYZ$x)
        graphics::lines(xy,...)
    } else {
        pars[3:4] <- 0
        XYZ <- ALR(log(c(fit13$ratio,fit23$ratio)),inverse=TRUE)
        xy <- xyz2xy(XYZ$x)
        graphics::points(xy,...)
    }
    pars
}

# pars = mu1, mu2, var1, var2, cov12
# dat = [n x 3] matrix of counts
LL.ternary.random.effects <- function(rho,dat,pars){
    mu <- pars[1:2]
    E <- matrix(0,2,2)
    E[1,1] <- pars[3]
    E[2,2] <- pars[4]
    E[1,2] <- rho*sqrt(E[1,1]*E[2,2])
    E[2,1] <- E[1,2]
    LL <- 0
    for (i in 1:nrow(dat)){
        LL <- LL + get.LL.sample(dat[i,],mu,E)
    }
    LL
}

get.LL.sample <- function(nn,mu,E){
    lnfact <- sum(1:sum(nn)) - sum(1:nn[1]) - sum(1:nn[2]) - sum(1:nn[3])
    p.int <- 
        integrate(function(b2) { 
            sapply(b2, function(b2) {
                integrate(function(b1) get.p.sample(b1,b2,nn,mu,E),
                          -Inf, Inf)$value
            })
        }, -Inf, Inf)$value
    - lnfact - log(p.int)
}

get.p.sample <- function(b1,b2,nn,mu,E){
    logbfact <- b1*nn[1] + b2*nn[2] - sum(nn) *
                log( exp(b1) + exp(b2) + 1 )
    mfact <- dmvnorm(cbind(b1,b2),mean=mu,sigma=E)
    exp(logbfact) * mfact
}
