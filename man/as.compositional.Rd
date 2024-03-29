% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{as.compositional}
\alias{as.compositional}
\title{create a \code{compositional} object}
\usage{
as.compositional(x, method = NULL, colmap = "rainbow")
}
\arguments{
\item{x}{an object of class \code{matrix}, \code{data.frame} or
\code{acomp}}

\item{method}{dissimilarity measure, either \code{"aitchison"} for
Aitchison's CLR-distance or \code{"bray"} for the Bray-Curtis
distance.}

\item{colmap}{the colour map to be used in pie charts.}
}
\value{
an object of class \code{compositional}
}
\description{
Convert an object of class \code{matrix}, \code{data.frame} or
\code{acomp} to an object of class \code{compositional}
}
\examples{
data(Namib)
PT.acomp <- as.acomp(Namib$PT)
PT.compositional <- as.compositional(PT.acomp)
print(Namib$PT$x - PT.compositional$x)
## uncomment the following lines for an illustration of using this 
## function to integrate 'provenance' with 'compositions'
# library(compositions)
# data(Glacial)
# a.glac <- acomp(Glacial)
# c.glac <- as.compositional(a.glac)
# summaryplot(c.glac,ncol=8)
}
