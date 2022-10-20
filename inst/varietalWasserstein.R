distmat <- function(x,y){
    m <- nrow(x)
    n <- nrow(y)
    xy <- x %*% t(y)
    xx <- matrix( rep(apply(x*x,1,sum), n), m, n, byrow=FALSE)
    yy <- matrix( rep(apply(y*y,1,sum), m), m, n, byrow=TRUE)
    out <- sqrt(pmax(xx+yy-2*xy, 0))
    return(out)
}

Wasserstein.diss.varietal <- function(x,package="approxOT"){
    snames <- names(x)
    ns <- length(snames)
    out <- matrix(0,ns,ns)
    rownames(out) <- colnames(out) <- snames
    for (snamei in snames){
        matches <- grepl(snamei,rownames(x$x))
        xi <- CLR(as.matrix(x$x[matches,]))
        ni <- nrow(xi)
        for (snamej in snames){
            matches <- grepl(snamej,rownames(x$x))
            xj <- CLR(as.matrix(x$x[matches,]))
            if (!identical(snamei,snamej)){
                nj <- nrow(xj)
                wi <- rep(1,ni)/ni
                wj <- rep(1,nj)/nj
                if (identical(package,"approxOT")){
                    d <- distmat(x=xi,y=xj)
                    out[snamei,snamej] <-
                        approxOT::wasserstein(a=wi,b=wj,cost=d,method="exact")
                } else if (identical(package,"transport")){
                    a <- transport::wpp(xi,mass=wi)
                    b <- transport::wpp(xj,mass=wj)
                    out[snamei,snamej] <- transport::wasserstein(a=a,b=b)
                } else {
                    stop("Unknown package")
                }
            }
        }
    }
    out <- as.dist(out)
    class(out) <- append("diss",class(out))
    return(out)
}
