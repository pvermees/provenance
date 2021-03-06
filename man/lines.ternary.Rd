% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ternary.R
\name{lines.ternary}
\alias{lines.ternary}
\title{Ternary line plotting}
\usage{
\method{lines}{ternary}(x, ...)
}
\arguments{
\item{x}{an object of class \code{ternary}, or a three-column data
frame or matrix}

\item{...}{optional arguments to the generic \code{lines} function}
}
\description{
Add lines to an existing ternary diagram
}
\examples{
tern <- ternary(Namib$PT,'Q',c('KF','P'),c('Lm','Lv','Ls'))
plot(tern,pch=21,bg='red',labels=NULL)
middle <- matrix(c(0.01,0.49,0.01,0.49,0.98,0.02),2,3)
lines(ternary(middle))
}
