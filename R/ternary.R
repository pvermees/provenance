#' Define a ternary composition
#'
#' Create an object of class \code{ternary}
#' @param X an object of class \code{compositional} OR a matrix or
#'     data frame with numerical data
#' @param x string/number or a vector of strings/numbers indicating the
#'     variables/indices making up the first subcomposition of the ternary system.
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
ternary <- function(X,x=1,y=2,z=3){
    out <- list()
    if (any(class(X) %in% c("compositional","counts"))){
        dat <- X$x
        nms <- names(X)
        class(out) <- class(X)
    } else {
        dat <- X
        nms <- rownames(X)
    }
    if (ndim(dat)>1){
        cnames <- colnames(dat)
        if (is.null(cnames)) cnames <- 1:ncol(dat)
    } else {
        cnames <- names(dat)
        if (is.null(cnames)) cnames <- 1:length(dat)
    }
    if (is.numeric(x)) x <- cnames[x]
    if (is.numeric(y)) y <- cnames[y]
    if (is.numeric(z)) z <- cnames[z]
    out$raw <- ternary.amalgamate(dat,x,y,z)
    if (nrow(out$raw)>1) rownames(out$raw) <- nms
    class(out) <- append("ternary",class(out))
    out$x <- ternaryclosure(out$raw)
    if (methods::is(X,"SRDcorrected")){
        out$restoration <- list()
        snames <- names(X$restoration)
        for (sname in snames){
            ternary_restoration <-
                ternary.amalgamate(X$restoration[[sname]],x,y,z)
            out$restoration[[sname]] <-
                ternaryclosure(ternary_restoration)
        }
        class(out) <- append("SRDcorrected",class(out))
    }
    return(out)
}

ternary.amalgamate <- function(dat,x,y,z){
    out <- cbind(sumcols(dat,x),sumcols(dat,y),sumcols(dat,z))
    xlab <- paste(x,collapse='+')
    ylab <- paste(y,collapse='+')
    zlab <- paste(z,collapse='+')
    colnames(out) <- c(xlab,ylab,zlab)
    out
}

# X is a 3-column matrix or vector
ternaryclosure <- function(X){
    den <- rowSums(X)
    out <- apply(X,2,'/',den)
    if (methods::is(out,"matrix")) {
        colnames(out) <- colnames(X)
    } else {
        names(out) <- names(X)
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
    if (type!='empty')
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
    if (is.na(pch) && is.null(pos) && showpath && methods::is(x,'SRDcorrected')){ pos <- 1 }
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
#' data(Namib)
#' tern <- ternary(Namib$Major,'CaO','Na2O','K2O')
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
#' plot a \eqn{100(1-\alpha)\%} confidence region around the data or
#'     around its mean.
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
    E <- stats::cov(uv)
    ell <- errellipse(n=nrow(x$x),mu=c(mean(u),mean(v)),Sigma=E,
                      alpha=alpha,population=population)
    XYZ <- ALR(ell,inverse=TRUE)
    xy <- xyz2xy(XYZ$x)
    graphics::lines(xy,...)
}
#' @examples
#' data(Namib)
#' tern <- ternary(Namib$Major,'CaO','Na2O','K2O')
#' plot(tern)
#' ternary.ellipse(tern)
#' @rdname ternary.ellipse
#' @export
ternary.ellipse.compositional <- function(x,alpha=0.05,population=TRUE,...){
    ternary.ellipse.default(x,alpha=alpha,population=population,...)
}
#' @rdname ternary.ellipse
#' @export
ternary.ellipse.counts <- function(x,alpha=0.05,population=TRUE,...){
    pars <- rep(0,5)
    fit <- central.multivariate(x)
    pars[1] <- log(fit['theta',1]) - log(1-fit['theta',1]) # mu[1]
    pars[2] <- log(fit['theta',2]) - log(1-fit['theta',2]) # mu[2]
    pars[3] <- fit['sigma',1]
    pars[4] <- fit['sigma',2]
    if (abs(pars[3])>0.01 & abs(pars[4])>0.01){ # draw ellipse
        if (TRUE){ # approximate but fast
            rho <- stats::cor(log(x$raw[,1:2]+0.5)-
                              log(x$raw[,3]+0.5))[1,2]
        } else { # accurate but slow and unstable
            init <- cor(log(x$raw[,1:2]+0.5)-log(x$raw[,3]+0.5))[1,2]
            message(paste0('Warning: calculating an error ellipse for \n',
                           'point-counting data may take a minute ...'))
            rho <- stats::optim(init,LL.ternary.random.effects.cor,
                                method="L-BFGS-B",lower=-0.99,upper=0.99,
                                pars=pars,dat=x$raw,
                                control=list(factr=1e14))$par
        }
        covariance <- rho*pars[3]*pars[4]
        b1 <- pars[1]
        b2 <- pars[2]
        E <- matrix(0,2,2)
        E[1,1] <- pars[3]^2
        E[2,2] <- pars[4]^2
        E[1,2] <- covariance
        E[2,1] <- E[1,2]
        pars[5] <- E[1,2]
        ell <- errellipse(n=nrow(x$x),mu=c(b1,b2),Sigma=E,
                          alpha=alpha,population=population)
        XYZ <- ALR(ell,inverse=TRUE)
    } else if (pars[3]<0.01 & pars[4]<0.01){ # draw point
        pars[3:4] <- 0
        XYZ <- ALR(pars[1:2],inverse=TRUE)
    } else { # draw line
        np <- 20 # number of points
        if (pars[3]>0.01){
            pars[4] <- 0
            if (population) sigma <- pars[3]
            else sigma <- fit['err',1]/(pars[1]*(1-pars[1]))
            m1 <- stats::qnorm(alpha/2,mean=pars[1],sd=sigma)
            M1 <- stats::qnorm(1-alpha/2,mean=pars[1],sd=sigma)
            u <- seq(m1,M1,length.out=np)
            v <- rep(pars[2],np)
        } else if (pars[4]>0.01){
            pars[3] <- 0
            if (population) sigma <- pars[4]
            else sigma <- fit['err',2]/(pars[2]*(1-pars[2]))
            m2 <- stats::qnorm(alpha/2,mean=pars[2],sd=sigma)
            M2 <- stats::qnorm(1-alpha/2,mean=pars[2],sd=sigma)
            u <- rep(pars[1],np)
            v <- seq(m2,M2,length.out=np)
        }
        uv <- cbind(u,v)
        XYZ <- ALR(uv,inverse=TRUE)
    }
    xy <- xyz2xy(XYZ$x)
    graphics::lines(xy,...)
    invisible(pars)
}

errellipse <- function(n,mu,Sigma,alpha=0.05,population=TRUE){
    m <- 2
    df1 <- m
    df2 <- n-m
    if (population) k <- n+1
    else k <- 1
    hk <- k*m*(n-1)*stats::qf(1-alpha,df1,df2)/(n*(n-m))
    VW2VT <- svd(Sigma)
    V <- VW2VT$u
    W2 <- diag(VW2VT$d)
    w1 <- sqrt(VW2VT$d[1])
    w2 <- sqrt(VW2VT$d[2])
    res <- 100
    g1 <- w1*sqrt(hk)*cos(seq(from=0,to=2*pi,length.out=res))
    g2 <- w2*sqrt(hk)*sin(seq(from=0,to=2*pi,length.out=res))
    G <- rbind(g1,g2)
    ell <- list()
    class(ell) <- 'compositional'
    ell$x <- t(V%*%G + rbind(rep(mu[1],res),rep(mu[2],res)))
    ell
}

# pars = mu1, mu2, var1, var2, cov12
# dat = [n x 3] matrix of counts
LL.ternary.random.effects.cor <- function(rho,pars,dat){
    pars[5] <- rho*pars[3]*pars[4]
    LL.ternary.random.effects(pars,dat)
}
LL.ternary.random.effects <- function(pars,dat){
    mu <- pars[1:2]
    E <- matrix(0,2,2)
    E[1,1] <- pars[3]^2
    E[2,2] <- pars[4]^2
    E[1,2] <- pars[5]
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
        stats::integrate(function(b2) { 
            sapply(b2, function(b2) {
                stats::integrate(function(b1) get.p.sample(b1,b2,nn,mu,E),
                                 -Inf, Inf)$value
            })
        }, -Inf, Inf)$value
    if (p.int==0) log.p.int <- -1000
    else log.p.int <- log(p.int)
    - lnfact - log.p.int
}
get.p.sample <- function(b1,b2,nn,mu,E){
    logbfact <- b1*nn[1] + b2*nn[2] - sum(nn) *
                log( exp(b1) + exp(b2) + 1 )
    mfact <- dmvnorm(cbind(b1,b2),mean=mu,sigma=E,log=TRUE)
    exp(logbfact + mfact)
}
