# provenance

**provenance** bundles a number of established statistical methods to
  facilitate the visual interpretation of large datasets in
  sedimentary geology. Includes functionality for adaptive kernel
  density estimation, principal component analysis, correspondence
  analysis, multidimensional scaling, generalised procrustes analysis
  and individual differences scaling using a variety of dissimilarity
  measures. Univariate provenance proxies, such as single-grain ages
  or (isotopic) compositions are compared with the Kolmogorov-Smirnov,
  Kuiper or Sircombe-Hazelton L2 distances. Categorical provenance
  proxies such as chemical compositions are compared with the
  Aitchison and Bray-Curtis distances, and point-counting data with
  the chi-square distance. Also included are tools to plot
  compositional and point-counting data on ternary diagrams, to
  calculate the sample size required for specified levels of
  statistical precision, and to assess the effects of hydraulic
  sorting on detrital compositions. Includes an intuitive query-based
  user interface for users who are not proficient in R..

## Prerequisites

You must have **R** installed on your system (see
[http://r-project.org](http://r-project.org)).  Additionally, to
install provenance from Github, you also need the **devtools**
package.  This can be installed by typing the following code at the R
command line prompt:

```
install.packages('devtools')
```

## Installation

The most recent stable version of provenance is available from **CRAN** at
[https://cran.r-project.org/package=provenance](https://cran.r-project.org/package=provenance)
and can be installed on your system as follows:

```
install.packages('provenance')
```

Alternatively, to install the current development version of
provenance from Github, type:

```
library(devtools)
install_github('pvermees/provenance')
```

## Further information

See [http://provenance.london-geochron.com](http://provenance.london-geochron.com)

[Vermeesch, P., Resentini, A. and Garzanti, E., 2016, An R package for
statistical provenance analysis, Sedimentary Geology, 336,
14-25](http://www.ucl.ac.uk/~ucfbpve/papers/VermeeschSedGeol2016/)

[Vermeesch, P., 2018, Statistical models for point-counting
data. Earth and Planetary Science Letters (in
press)](http://www.ucl.ac.uk/~ucfbpve/papers/VermeeschEPSL2018/)

## Author

[Pieter Vermeesch](http://pieter.london-geochron.com)

## License

This project is licensed under the GPL-3 License