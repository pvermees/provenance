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
#' @param components a vector specifying a subcomposition
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
#' @param alpha cutoff value for confidence intervals
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
#'
#' @rdname radialplot
#' @export
radialplot.counts <- function(x,components=NA,from=NA,to=NA,t0=NA,
                              sigdig=2,show.numbers=FALSE,pch=21,
                              levels=NA,clabel="",
                              bg=c("white","red"),title=TRUE,
                              alpha=0.05,...){
    if (is.na(components)) components <- colnames(x$x)[1:2]
    label <- paste0('central ',components[1],'/',components[2],'-ratio')
    dat <- x$x[,components]
    X <- x2zs(dat)
    X$transformation <- 'arctan'
    IsoplotR:::radial.plot(X,show.numbers=show.numbers,pch=pch,
                           levels=levels,clabel=clabel,bg=bg,...)
    fit <- central(x,components=components)
    rounded.ratio <- IsoplotR:::roundit(fit$ratio,fit$err,sigdig=sigdig)
    line1 <- substitute(a~'='~b%+-%c~'%',
                        list(a=label,
                             b=rounded.ratio[1],
                             c=rounded.ratio[2]))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,sigdig),
                             b=signif(fit$p.value,sigdig)))
    line3 <- substitute('dispersion ='~a~'%',
                        list(a=signif(100*fit$sigma,sigdig)))
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

#' @rdname central
#' @export
central.counts <- function(x,components=NA,...){
    if (all(is.na(components))) components <- colnames(x$x)[1:2]
    out <- list()
    sigma <- 0.15 # convenient starting value
    valid <- which(rowSums(x$x[components])>0)
    Nsj <- x$x[valid,components[1]]
    Nij <- x$x[valid,components[2]]
    Ns <- sum(Nsj)
    Ni <- sum(Nij)
    num <- (Nsj*Ni-Nij*Ns)^2
    den <- Nsj+Nij
    Chi2 <- sum(num/den)/(Ns*Ni)
    mj <- Nsj+Nij
    pj <- Nsj/mj
    theta <- Ns/sum(mj)
    for (i in 1:30){ # page 49 of Galbraith (2005)
        wj <- mj/(theta*(1-theta)+(mj-1)*(theta*(1-theta)*sigma)^2)
        sigma <- sigma * sqrt(sum((wj*(pj-theta))^2)/sum(wj))
        theta <- sum(wj*pj)/sum(wj)
    }
    # remove two d.o.f. for mu and sigma
    out$df <- length(Nsj)-2
    # add back one d.o.f. for homogeneity test
    out$mswd <- Chi2/(out$df+1)
    out$p.value <- 1-stats::pchisq(Chi2,out$df+1)
    out$ratio <- theta/(1-theta)
    out$err <- sqrt(1/(sum(wj)*(theta*(1-theta))^2))
    out$sigma <- sigma
    out
}
