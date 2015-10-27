#' @title Colorado geochemical data
#'
#' @description
#' In 2006, soil samples were collected at 960 sites (1 site per 280 square
#' kilometers) throughout the state of Colorado. These samples were collected
#' from a depth of 0-15 centimeters and, following a near-total multi-acid
#' digestion, were analyzed for a suite of 44 major and trace
#' elements.
#'
#' @format
#' The data is organized as a list with three elements:
#' \itemize{
#' \item The first list element is the geochemical data and is
#' called "concData". The data are from
#' Smith and others (2010); selected summary statistics of these
#' data are presented in Ellefsen and others (2014).
#'
#' These data have been editted to make them suitable for
#' analysis. First, field sample "06co437" has been deleted because
#' it has an anomalously high copper concentration.
#' Second, elements silver (Ag), cesium (Cs), mercury (Hg),
#' tellurium (Te), and selenium (Se) have been deleted because 98%, 77%,
#' 69%, 94%, and 45%, respectively, of the measured concentrations are below
#' the lower limit of determination. (The lower limit of determination is
#' the threshold for left-censoring.) The lower limits of determination are
#' listed in Table 1 of Smith et al., 2010.
#'
#' Third, elements antimony (Sb), arsenic (As), bismuth (Bi), cadmium (Cd),
#' indium (In), phosphorous (P), and sulfur (S) have 0.208%, 0.417%, 1.15%,
#' 5.62%, 14.1%, 0.521%, and 3.02%, respectively, of the measured
#' concentrations below the lower limit of determination.
#' These left-censored concentrations were assigned concentrations equal to
#' 0.65 times their respective lower limit of determination.
#' Finally, element concentrations were
#' scaled, as necessary, so that the units for all concentrations are "mg/kg"
#' or equivalently "ppm".
#'
#' After this editting, there are 959 field samples for which
#' 39 element concentrations are reported.
#' The locations of the field samples are specified by lattitude and longitude
#' using the WGS84 datum.
#'
#' Both the element concentrations for the field samples
#' and the locations of the field samples are
#' stored in SpatialPointsDataFrame from the "sp" library. In this
#' container, the element concentrations are stored in an R data frame. It has
#' 959 rows, corresponding to the 959 field samples, and the row names are
#' the field sample names. It has 39 columns, corresponding to the 39
#' elements, and the column names are the abbreviated element names
#' (for example, "Al"). If explanations of these common abbreviations are
#' needed, they may be found in Smith and others (2010).
#'
#' \item The second list element is an array that identifies the concentrations
#' that are censored. Its name is "censorIndicators."
#'
#' The structure of the array is identical to the structure of the data frame
#' containing the concentrations, so that identifying the status of a
#' reported concentration is easy.
#' The array elements have two possible values: "no" meaning that the
#' concentration is not censored and "left" meaning that the concentration
#' is left-censored and now has an imputed value.
#'
#' \item The third list element is a scalar and is called "kappa".
#'
#' For each field sample, the reported geochemical concentrations in container
#' "concentrationData" plus all unreported geochemical concentrations
#' (for example, Si) must equal 1000000 mg/kg or equivalently 1000000 ppm.
#' Thus, kappa equals 1000000.
#' }
#'
#' @source \url{http://pubs.usgs.gov/ds/520/}
#'
#' @references
#' Ellefsen, K.J., Smith, D.B., Horton, J.D., 2014, A modified procedure
#' for mixture-model clustering of regional geochemical data: Applied
#' Geochemistry, vol. 51, p. 315-326,
#' doi: http://dx.doi.org/10.1016/j.apgeochem.2014.10.011.
#'
#' Smith, D.B., Ellefsen, K.J., and Kilburn, J.E., 2010,
#' Geochemical data for Colorado soils-Results from the 2006
#' state-scale geochemical survey: U.S. Geological Survey,
#' Data Series 520, 9 p.
#'
"CoGeochemData"

