#' Visualise point-counting data on a radial plot
#'
#' Implementation of a graphical device developed by Rex Galbraith to
#' display several estimates of the same quantity that have different
#' standard errors.
#'
#' @details
#'
#' The radial plot (Galbraith, 1988, 1990) is a graphical device that
#' was specifically designed to display heteroscedastic data, and is
#' constructed as follows.  Consider a set of dates
#' \eqn{\{t_1,...,t_i,...,t_n\}} and uncertainties
#' \eqn{\{s[t_1],...,s[t_i],...,s[t_n]\}}. Define \eqn{z_i = z[t_i]}
#' to be a transformation of \eqn{t_i} (e.g., \eqn{z_i = log[t_i]}),
#' and let \eqn{s[z_i]} be its propagated analytical uncertainty
#' (i.e., \eqn{s[z_i] = s[t_i]/t_i} in the case of a logarithmic
#' transformation). Create a scatterplot of \eqn{(x_i,y_i)} values,
#' where \eqn{x_i = 1/s[z_i]} and \eqn{y_i = (z_i-z_\circ)/s[z_i]},
#' where \eqn{z_\circ} is some reference value such as the mean. The
#' slope of a line connecting the origin of this scatterplot with any
#' of the \eqn{(x_i,y_i)}s is proportional to \eqn{z_i} and, hence,
#' the date \eqn{t_i}.  These dates can be more easily visualised by
#' drawing a radial scale at some convenient distance from the origin
#' and annotating it with labelled ticks at the appropriate
#' angles. While the angular position of each data point represents
#' the date, its horizontal distance from the origin is proportional
#' to the precision. Imprecise measurements plot on the left hand side
#' of the radial plot, whereas precise age determinations are found
#' further towards the right. Thus, radial plots allow the observer to
#' assess both the magnitude and the precision of quantitative data in
#' one glance.
#'
#' @param x an object of class \code{counts}
#' @param num index or name of the numerator variable
#' @param den index or name of the denominator variable
#' @param from minimum limit of the radial scale
#' @param to maximum limit of the radial scale
#' @param t0 central value
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output.
#' @param show.numbers boolean flag (\code{TRUE} to show sample
#'     numbers)
#' @param pch plot character (default is a filled circle)
#' @param levels a vector with additional values to be displayed as
#'     different background colours of the plot symbols.
#' @param clabel label of the colour legend
#' @param bg a vector of two background colours for the plot symbols.
#'     If \code{levels=NA}, then only the first colour is used. If
#'     \code{levels} is a vector of numbers, then \code{bg} is used to
#'     construct a colour ramp.
#' @param title add a title to the plot?
#' @param ... additional arguments to the generic \code{points}
#'     function
#' @references
#' Galbraith, R.F., 1988. Graphical display of estimates
#' having differing standard errors. Technometrics, 30(3),
#' pp.271-281.
#'
#' Galbraith, R.F., 1990. The radial plot: graphical assessment of
#' spread in ages. International Journal of Radiation Applications and
#' Instrumentation. Part D. Nuclear Tracks and Radiation Measurements,
#' 17(3), pp.207-214.
#'
#' Galbraith, R.F. and Laslett, G.M., 1993. Statistical models for
#' mixed fission track ages. Nuclear Tracks and Radiation
#' Measurements, 21(4), pp.459-470.
#' @examples
#' data(Namib)
#' radialplot(Namib$PT,components=c('Q','P'))
#' @export
radialplot <- function(x,num=1,den=2,from=NA,to=NA,t0=NA,
                       sigdig=2,show.numbers=FALSE,pch=21,
                       levels=NA,clabel="",
                       bg=c("white","red"),title=TRUE,...){
    ncol <- ncol(x$x)
    if (all(is.numeric(num))) num <- colnames(x$x)[num]
    if (all(is.numeric(den))) den <- colnames(x$x)[den]
    label <- paste0('central ',num,'/',den,'-ratio')
    dat <- subset(x,components=c(num,den))
    X <- x2zs(dat$x)
    X$transformation <- 'arctan'
    pcol <- IsoplotR:::set.ellipse.colours(ns=nrow(x$x),levels=levels,col=bg)
    IsoplotR:::radial.plot(X,show.numbers=show.numbers,pch=pch,
                           levels=levels,clabel=clabel,bg=pcol,...)
    fit <- central(dat)
    ratio <- fit['theta',1]/fit['theta',2]
    err <- fit['err',1]*ratio
    rounded.ratio <- IsoplotR:::roundit(ratio,err,sigdig=sigdig)
    line1 <- substitute(a~'='~b%+-%c~(1*sigma),
                        list(a=label,
                             b=rounded.ratio[1],
                             c=rounded.ratio[2]))
    line2 <- substitute('MSWD ='~a~', p('*chi^2*')='~b,
                        list(a=signif(fit['mswd',1],sigdig),
                             b=signif(fit['p.value',1],sigdig)))
    line3 <- substitute('dispersion ='~a*'%',
                        list(a=signif(100*fit['sigma',1],sigdig)))
    graphics::mtext(line1,line=2)
    graphics::mtext(line2,line=1)
    graphics::mtext(line3,line=0)
}

x2zs <- function(x){
    out <- list()
    n <- x[,1]
    m <- x[,2]
    out$z <- atan(sqrt((n+3/8)/(m+3/8)))
    out$s <- sqrt(1/(n+m+1/2))/2
    out$z0 <- atan(sqrt(sum(n,na.rm=TRUE)/
                        sum(m,na.rm=TRUE)))
    out$from <- min(tan(out$z)^2)
    out$to <- max(tan(out$z)^2)
    out$xlab <- paste(colnames(x),collapse='+')
    out
}

#' Calculate central compositions
#'
#' Computes the geometric mean composition of a continuous mixture of
#' point-counting data.
#'
#' @details The central composition assumes that the observed
#'     point-counting distribution is the combination of two sources
#'     of scatter: counting uncertainty and true geological
#'     dispersion.
#'
#' @param x an object of class \code{counts}
#' @param ... optional arguments
#' @return an \code{[5 x n]} matrix with \code{n} being the number
#' of categories and the rows containing:
#'
#' \describe{
#' \item{theta}{ the `central' composition. }
#' \item{err}{ the standard error for the central composition. }
#' \item{sigma}{ the overdispersion parameter, i.e. the coefficient of
#'               variation of the underlying logistic normal
#'               distribution. \code{central} computes a continuous
#'               mixture model for each component (column)
#'               separately. Covariance terms are not reported.}
#' \item{mswd}{ the mean square of the weighted deviates, a.k.a.
#'              reduced chi-square statistic.}
#' \item{p.value}{ the p-value for age homogeneity }
#' }
#' @export
central <- function(x,...){
    if ("ternary" %in% class(x)) dat <- x$raw
    else dat <- x$x
    ns <- nrow(dat)
    nc <- ncol(dat)
    out <- matrix(0,5,nc)
    colnames(out) <- colnames(dat)
    rownames(out) <- c('theta','err','sigma','mswd','p.value')
    for (i in 1:nc){
        Nsj <- subset(dat,select=i)
        Nij <- rowSums(subset(dat,select=-i))
        out[,i] <- central_helper(Nsj,Nij)
    }
    out
}
# calculates proportions relative to the last component
central.multivariate <- function(x,...){
    if ("ternary" %in% class(x)) dat <- x$raw
    else dat <- x$x
    ns <- nrow(dat)
    nc <- ncol(dat)
    out <- matrix(0,5,nc-1)
    rownames(out) <- c('theta','err','sigma','mswd','p.value')
    Nij <- subset(dat,select=nc) # use last column as denominator
    for (i in 1:(nc-1)){ # loop over all but the last column
        Nsj <- subset(dat,select=i)
        out[,i] <- central_helper(Nsj,Nij)
    }
    colnames(out) <- colnames(dat)[-nc]
    out
}

central_helper <- function(Nsj,Nij){
    sigma <- 0.15 # convenient starting value
    Ns <- sum(Nsj)
    Ni <- sum(Nij)
    mj <- Nsj+Nij
    ispos <- (mj>0)
    pj <- 0*Nsj
    pj[ispos] <- Nsj[ispos]/mj[ispos]
    pj[!ispos] <- Ns/(Ns+Ni)
    theta <- Ns/sum(mj)
    #nn <- length(Nsj)
    for (i in 1:30){ # from page 49 of Galbraith (2005)
        wj <- mj/(theta*(1-theta)+(mj-1)*(theta*(1-theta)*sigma)^2)
        sigma <- sigma * sqrt(sum((wj*(pj-theta))^2)/sum(wj))
        #sigma <- sigma * sqrt(sum(wj*(pj-theta)^2)/(nn-1))
        theta <- sum(wj*pj)/sum(wj)
    }
    if (sigma>2){ # the point iteration method may be fast but it doesn't
        ts <- central_helper_ML(Nsj,Nij) # work well for very dispersed 
        theta <- ts[1]                   # datasets
        sigma <- ts[2]
        wj <- mj/(theta*(1-theta)+(mj-1)*(theta*(1-theta)*sigma)^2)
    }
    num <- (Nsj*Ni-Nij*Ns)^2
    Chi2 <- sum(num[ispos]/mj[ispos])/(Ns*Ni)
    # remove one d.o.f. for theta (for homogeneity test)
    df <- length(Nsj)-1
    mswd <- Chi2/df    
    p.value <- 1-stats::pchisq(Chi2,df)
    err <- 1/(sqrt(sum(wj))*(1-theta)^2)
    c(theta,err,sigma,mswd,p.value)
}

# maximum likelihood alternative to central_helper(Nsj,Nij)
central_helper_ML <- function(Nsj,Nij){
    init <- pilot(Nsj,Nij)
    mulogsigma <- stats::optim(init,LL.random.effects,Nsj=Nsj,Nij=Nij)$par
    mu <- mulogsigma[1]
    sigma <- exp(mulogsigma[2])
    theta <- exp(mu)/(1+exp(mu))
    c(theta,sigma)
}
pilot <- function(Nsj,Nij){
    logits <- log(Nsj+0.5)-log(Nij+0.5)
    mu <- mean(logits)
    logsigma <- log(stats::sd(logits))
    c(mu,logsigma)
}
LL.random.effects <- function(mls,Nsj,Nij){
    mu <- mls[1]
    sigma <- exp(mls[2])
    LL <- 0
    for (u in 1:length(Nsj)){
        LL <- LL + log(pyumu(Nsj[u],Nij[u],mu,sigma))
    }
    -LL
}
pyumu <- function(Ns,Ni,mu,sigma){
    stats::integrate(function(b) integrand(b,Ns,Ni,mu,sigma),lower=-Inf,upper=Inf)$value
}
integrand <- function(b,Ns,Ni,mu,sigma){
    lognum <- b*Ns
    logden <- -(Ns+Ni)*log(1+exp(b))
    logout <- lognum + logden - 0.5*((b-mu)/sigma)^2 - log(sigma) - pi
    exp(logout)
}
