#' Calculate the number of grains required to achieve a desired level of sampling resolution
#'
#' Returns the number of grains that need to be analysed to decrease
#' the likelihood of missing any fraction greater than a given size
#' below a given level.
#' @param f the size of the smallest resolvable fraction (0<f<1)
#' @param p the probability that all n grains in the sample have missed
#' at least one fraction of size \code{f}
#' @param n, the number of grains in the sample
#' @return the number of grains needed to reduce the chance of missing
#' at least one fraction f of the total population to less than \code{p}
#' @references Vermeesch, Pieter. "How many grains are needed for a
#' provenance study?." Earth and Planetary Science Letters 224.3
#' (2004): 441-451.
#' @examples
#' # number of grains required to be 99% that no fraction greater than 5% was missed:
#' print(get.n(0.01))
#' # number of grains required to be 90% that no fraction greater than 10% was missed:
#' print(get.n(p=0.1,f=0.1))
#' @export
get.n <- function(p=0.05,f=0.05){
    n <- 1
    while(TRUE){
        pp <- get.p(n,f)
        if (pp<p){ break }
        else {n <- n+1}
    }
    return(n)
}

#' Calculate the probability of missing a given population fraction
#'
#' For a given sample size, returns the likelihood of missing any
#' fraction greater than a given size
#' @param n the number of grains in the detrital sample
#' @param f the size of the smallest resolvable fraction
#'     (0<\code{f}<1)
#' @return the probability that all \code{n} grains in the sample have
#'     missed at least one fraction of size \code{f}
#' @references Vermeesch,
#'     Pieter. "How many grains are needed for a provenance study?."
#'     Earth and Planetary Science Letters 224.3 (2004): 441-451.
#' @examples
#' print(get.p(60))
#' print(get.p(117))
#' @export
get.p <- function(n,f=0.05){
    if (f<0.03){ # bins are approximately independent
        p <- 1-(1-(1-f)^n)^(1/f)
    } else { # bins are not independent
        p <- 0
        M <- 1/f
        if (M*f == 1) A <- 0
        else A <- 1
        if (M*f < 1) B <- 1
        else B <- 0
        for (i in 1:(M-A)){
            p1 <- choose(M-A,i)*(1-i*f)^n
            p2 <- B*choose(M-1,i-1)*((M-i)*f)^n
            p <- p + (-1)^(i-1) * (p1 + p2)
        }
    }
    return(p)
}

#' Calculate the largest fraction that is likely to be missed
#'
#' For a given sample size, returns the largest fraction which has
#' been sampled with (1-p) x 100 \% likelihood.
#' @param n the number of grains in the detrital sample
#' @param p the required level of confidence
#' @return the largest fraction that is sampled with at least (1-p) x
#'     100\% certainty
#' @references
#' Vermeesch, Pieter. "How many grains are needed for a provenance study?"
#' Earth and Planetary Science Letters 224.3 (2004): 441-451.
#' @examples
#' print(get.f(60))
#' print(get.f(117))
#' @export
get.f <- function(n,p=0.05){
    misfit <- function(f,n,p){ get.p(n,f) - p }
    stats::uniroot(misfit,interval=c(0,1),n=n,p=p)$root
}
