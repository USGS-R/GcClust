% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GcClusterFunctions.R
\name{plotClusters}
\alias{plotClusters}
\title{Plot the field samples as clusters}
\usage{
plotClusters(gcData, condProbs1, probIntervals = c(0, 0.1, 0.5, 0.9, 1),
  symbolIndices = c(16, 16, 16, 16), symbolSizes = c(1/3, 1/3, 1/3, 1/3),
  symbolColors = c("red", "yellow", "green", "blue"))
}
\arguments{
\item{gcData}{List containing the geochemical and related data. This container is
described in the package documentation.}

\item{condProbs1}{A matrix containing the Monte Carlo samples of the
conditional probabilities. This matrix is returned by function
\code{\link{calcCondProbs1}}, for which the documentation includes a
complete description of container \code{condProbs1}.}

\item{probIntervals}{Vector containing intervals of conditional
probability. All field samples within an given interval are plotted the same
way.}

\item{symbolIndices}{Vector containing the indices of the plotting symbols
for the conditional probability intervals.}

\item{symbolSizes}{Vector containing the relative sizes of the plotting symbols
for the conditional probability intervals.}

\item{symbolColors}{Vector containing the colors of the plotting symbols
for the conditional probability intervals.}
}
\description{
Plot the locations of the field samples on a previously-plotted
map. The attributes of each location symbol (for example, color) indicate the
cluster to which the field sample belongs. That is, the attributes indicate
the conditional probability that the field sample is associated with
a particular pdf in the finite mixture model.
}
\details{
The conditional probabilities indicate the extend
to which the field samples are associated with the first pdf in the finite
mixture model. The conditional probabilities
in container \code{condProbs1} are Monte Carlo samples,
and their medians are used to assign plotting attributes for the field
samples.

The plotting attributes are specified by arguments \code{probIntervals},
\code{symbolIndices}, \code{symbolSizes}, and \code{symbolColors}. To
understanding the specification of these attributes, consider their default
values, which pertain to four probability intervals.

Vector \code{probIntervals} has elements 0, 0.1, 0.5, 0.9, and 1. These five
elements specify four probability intervals: [0,0.1], [0.1,0.5],
[0.5,0.9] and
[0.9,1]. Notice that the first and last elements of \code{probIntervals}
are 0 and 1 respectively. The probability intervals are used to
classify the field samples based upon their associated conditional
probabilities:
\itemize{
 \item If the conditional probability of a field sample is
 within the interval [0,0.1], then the field sample is classified as
 "strongly associated with pdf 2" and is assigned the color red. This
 color is consistent with the colors used in functions
 \code{\link{plotStdCompMeans}} and \code{\link{plotCompMeans}}.
 \item If the conditional probability of a field sample is
 within the interval [0.1,0.5], then the field sample is classified as
 "moderately associated with pdf 2" and is assigned the color yellow.
 \item If the conditional probability of a field sample is
 within the interval [0.5,0.9], then the field sample is classified as
 "moderately associated with pdf 1" and is assigned the color green.
 \item If the conditional probability of a field sample is
 within the interval [0.9,1], then the field sample is classified as
 "strongly associated with pdf 1" and is assigned the color blue. This
 color is consistent with the colors used in functions
 \code{\link{plotStdCompMeans}} and \code{\link{plotCompMeans}}.
}

Because, in this explanation, vector \code{probIntervals} specifies four
probability intervals, arguments for vectors \code{symbolIndices},
\code{symbolSizes}, and \code{symbolColors} must have four elements.
The symbol indices, symbol sizes, and symbol colors are described in
Murrell (2006, p. 55-56, 68, 69).

This function adds symbols to a map that has already been plotted.
}
\examples{
\dontrun{
map('state', fill = TRUE, col = "gray60", border = "white")
plotClusters(concentrationData, condProbs1)
}

}
\references{
Murrell, P., 2006, R graphics: Chapman & Hall / CRC.
}

