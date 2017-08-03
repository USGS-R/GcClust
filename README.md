GcClust
=======

[![status](https://img.shields.io/badge/USGS-Research-blue.svg)](https://owi.usgs.gov/R/packages.html#research)

GcClust is a software package for statistical clustering of regional geochemical data, and similar data such as regional mineralogical data. The clustering procedure partitions the field samples for a data set into two clusters. Each cluster is partitioned again to create two sub-clusters, and so on, generating a hierarchy of clusters. The clustering method is based on a Bayesian finite mixture model. The model parameters are estimated with Hamiltonian Monte Carlo sampling of the posterior probability density function.

The functions that implement the clustering, the documentation of the functions, and a user’s guide that describes each step of the clustering are bundled together in an “R package,” which is the unit of sharable code for the R statistical programming language.

If you are interested in obtaining the source code for the package, then it can be downloaded from this GIT repository. However, if you only want to use the package, then you can obtain the installation files (as well as the user’s guide formatted as a USGS publication) from <https://pubs.er.usgs.gov/publication/tm7C13>. When the digital object identifier has been registered, then the url will be <http://dx.doi.org/10.3133/tm7c13>.

<!-- README.md is generated from README.Rmd. Please edit that file -->
