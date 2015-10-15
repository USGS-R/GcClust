
#' @title Calculate sample statistics pertinent to the simplex
#'
#' @description Calculate the sample center, total variation matrix, and
#' the metric variance.
#'
#' @param gcData   SpatialPointsDataFrame storing the geochemical
#'                 concentrations and locations of each field sample
#' @param kappa    Constant-sum value for the geochemical concentrations.
#'                 Typically \code{kappa} equals 1,000,000.
#'
#' @details
#' The sample statistics are described in
#' section 5.2 of Pawlowsky-Glahn et al. (2015).
#'
#' @return A list with three elements is returned. Variable D is the number
#' of geochemical concentrations reported for each sample.
#' @return \item{sampleCenter}{Vector of length D.}
#' @return \item{variationMatrix}{Matrix of dimension D x D.}
#' @return \item{metricVariance}{Scalar. The metric variance is also called
#' the total sample variance.}
#'
#' @references
#' Pawlowsky-Glahn, V., Egozcue, J.J., and Tolosana-Delgado, R., 2015, Modeling
#' and analysis of compositional data: John Wiley and Sons, Ltd.
#'
#' @examples
#' \dontrun{
#' simplexStats <- calcSimplexStats( gcData, kappa )
#' }
#'
#' @export
calcSimplexStats <- function( gcData, kappa ) {

  concData <- as.matrix(gcData@data)

  sampleCenter <- CoDA::CalcCompCenter( concData, kappa=kappa )
  variationMatrix <- CoDA::CalcVariationMatrix( concData )
  metricVariance <- CoDA::CalcTotalVariance( variationMatrix )

  return( list(sampleCenter = sampleCenter,
               variationMatrix = variationMatrix,
               metricVariance = metricVariance))
}


#' @title Transform the geochemical data
#'
#' @description Transform the geochemical data first with the
#' isometric log-ratio transform and then with the robust principal
#' component transform.
#'
#' @param gcData   SpatialPointsDataFrame storing the geochemical
#'                 concentrations and location for each field sample.
#'
#' @param alpha         Fraction of the data used for the robust principal
#'                      component transform (See Details).
#'
#' @details
#' The geochemical data are stored in SpatialPointsDataFrame \code{gcData},
#' which is part of the sp package.
#'
#' The chemical concentrations are transformed twice: First, the concentrations
#' are transformed to isometric log-ratio (ilr) coordinates. This
#' transformation is described in Pawlowsky-Glahn et al. (2015, p. 36-38).
#' Second, the ilr coordinates are transformed to robust principal coordinates.
#' This tranformation is decribed in Filzmoser et al. (2009).
#' The transformation requires mean vector and covariance matrix; robust
#' values for these two statistics are calculated
#' are calculated using function covMcd from
#' package robustbase. An argument to function covMcd is \code{alpha}, which is
#' the fraction of the ilr-transformed coordinates
#' that are used to calculate the two statistics.
#'
#' @return A list with seven elements is returned. The elements are vectors
#' and matrices for which the dimensions depend on N, the number of field
#' samples, and D, the number of geochemical concentrations reported for
#' each sample.
#' @return \item{Psi}{Contrast matrix that is used for the ilr transformation.
#' The matrix dimensions are (D-1) x D.}
#' @return \item{ilrCoefs}{Matrix of ilr coefficients (coordinates) resulting
#' from ilr transformation of the geochemical concentrations. The
#' matrix dimensions are N x (D-1).}
#' @return \item{robustIlrCenter}{Vector containing the robust mean of the
#' ilr coordinates. The vector dimension is D-1.}
#' @return \item{robustEigenvectors}{Matrix containing the eigenvectors of the
#' robust covariance matrix for the ilr coordinates. The matrix dimension
#' is (D-1) x (D-1).}
#' @return \item{robustEigenvalues}{Vector containing the eigenvalues of
#' the robust covariance matrix for the ilr coordinates. The vector
#' dimension is D-1.}
#' @return \item{robustPCs}{Matrix containing the robust principal components.
#' The matrix dimension is N x (D-1).}
#' @return \item{alpha}{Scalar containing the input argument \code{alpha}.}
#'
#' @references
#' Filzmoser, P., Hron, K., and Reimann, C., 2009, Principal component
#' analysis for compositional data with outliers: Environmetrics, v. 20,
#' p. 621-632.
#'
#' Pawlowsky-Glahn, V., Egozcue, J.J., and Tolosana-Delgado, R., 2015, Modeling
#' and analysis of compositional data: John Wiley and Sons, Ltd.
#'
#' @examples
#' \dontrun{
#' transData <- gcTransform(X)
#' }
#'
#' @export
transformGcData <- function(gcData, alpha = 0.98) {

  X <- as.matrix(gcData@data)

  Psi <- CoDA::CalcPsiMatrix2( ncol(X) )
  ilrCoefs <- CoDA::CalcIlrCoefs( X, t(Psi) )

  mcdResult <- robustbase::covMcd( ilrCoefs, alpha=alpha )

  robustIlrCenter <- mcdResult$center
  robustIlrCov <- mcdResult$cov

  # Implements part of equation 10 in Filzmoser et al. (2009).
  # ilrCoefs corresponds to matrix Z
  # the second term (right side) corresponds to matrix 1T(Z)'
  # centeredIlrCoefs corresponds to matrix Z - 1T(Z)'
  centeredIlrCoefs <- ilrCoefs -  rep.int( 1, nrow( ilrCoefs ) ) %o% robustIlrCenter

  # Implements equation 11 in Filzmoser et al. (2009).
  # svdResult$u corresponds to matrix Gz
  # svdResult$d corresponds to the diagonal of matrix Lz (which is a diagonal matrix)
  svdResult <- svd( robustIlrCov )

  robustEigenvectors <- svdResult$u
  robustEigenvalues <- svdResult$d   # namely, robust variances

  # Implements the rest of equation 10 in Filzmoser et al. (2009).
  # centeredIlrCoefs %*% svdResult$u  corresponds to [Z - 1T(Z)']Gz
  # robustPCs corresponds to Z*
  robustPCs <- centeredIlrCoefs %*% robustEigenvectors

  return( list( Psi = Psi,
                ilrCoefs = ilrCoefs,
                robustIlrCenter = robustIlrCenter,
                robustEigenvectors = robustEigenvectors,
                robustEigenvalues = robustEigenvalues,
                robustPCs = robustPCs,
                alpha = alpha ))
}

#' @title Plot the scree plot
#'
#' @description Plot the variances associated with the (robust) principal
#' components.
#'
#' @param transData     List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}; the
#' documentation for function transformGcData includes a complete
#' description of container \code{transData}.
#' @param relOffset  Scalar specifying the relative distance that each
#' percentage of cummulative variance (see Details) is offset from the top
#' of its respective bar.
#' @param size  Scalar specifying the text size for the percentage of
#' cummulative variance (see Details).
#'
#' @details
#' In principal component analysis, this plot is called a "scree plot."
#'
#' The variance for each component is represented by a bar, and above each
#' bar is the cummulative percentage of total variance. The total variance
#' is the sum of the variances for each component. The percentage of variance
#' associated a component is its variance divided by the total variance
#' times 100. The cummulative percentage of total variance for component k
#' is the sum of the percentages from component 1 to component k.
#'
#' @references
#' Johnson, R.A., and Wichern, D.W., 2007, Applied multivariate statistical
#' analysis (6th ed.): Pearson Prentice Hall.
#'
#' @examples
#' \dontrun{
#' plotVariances(transData)
#' }
#'

#' @export
plotVariances <- function(transData, relOffset = 0.04, size = 3) {

  variances <- transData$robustEigenvalues
  tmp <- 100 * cumsum( variances ) / sum( variances )
  # cPTV: cummulative percentage of total variance
  varInfo <- data.frame(component = 1:length(variances),
                        variances = variances,
                        cPTV = as.character( round( tmp, digits=2 )),
                        stringsAsFactors = FALSE )

  offset <- relOffset * variances[1]

  p <- ggplot2::ggplot( varInfo,
                        ggplot2::aes(x = component, y = variances),
                        environment = environment() ) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::geom_text(ggplot2::aes(x = component,
                                    y = variances+offset,
                                    label = cPTV),
                       angle = 90, size = size ) +
    ggplot2::xlab("Component") +
    ggplot2::ylab("Variance")

  print(p)
}



#' @title Optimize the posterior pdf using package mclust
#'
#' @description Optimize the posterior probability density function of the
#' finite mixture model. The model implementation and the optimization are
#' performed with package mclust.
#'
#' @param transData     List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}; the
#' documentation for function transformGcData includes a complete
#' description of container \code{transData}.
#' @param nPCS          Number of principal components that are used in the
#'                      finite mixture model (See Details).
#' @param nCpuCores      Number of central processing unit (cpu) cores that
#'                       are used
#'                      for parallel computations (See Details).
#' @param nSoln         Number of solutions (See Details).
#'
#' @details
#' The parameters in the finite mixture model are estimated
#' from the robust principal components, which are
#' calculated with function gcTransform and stored as a matrix within
#' \code{transData}. All matrix elements must contain values;
#' that is, missing values are not allowed.
#'
#' The number of principal components, \code{nPCS}, means that
#' principal components 1, 2, ..., nPCs are used in the finite mixture model.
#' That is, all higher-order principal components are not used.
#'
#' The maximum of the posterior probability density function (pdf) is calculated
#' repetitively with different starting points.
#' The number of repetitive solutions is specified by argument \code{nSoln}.
#' In practice, it is usually best to use the default value for
#' argument \code{nSoln}.
#'
#' The repetitive solutions are calculated in parallel on multiple central
#' processing units (cpu's). Of course, the number of requested cpu's must be
#' less than or equal to the actual number.
#'
#' @return A list with the following components is returned.
#' @return \item{nSoln}{Number of solutions.}
#' @return \item{theLogLikelihoods}{The logarithms of the likelihoods for
#' the multiple solutions.}
#' @return \item{bestClusterResult}{The clustering result with the highest
#' log-likelihood. The structure is identical to the structure returned by
#' function Mclust in package mclust.}
#' @return \item{worstClusterResult}{The clustering result with the lowest
#' log-likelihood. The structure is identical to the structure returned by
#' function Mclust in package mclust.}
#' @return \item{nErrors}{Sometimes the function Mclust fails, and the
#' total number of such failures is nErrors.}

#' @references
#' Ellefsen, K.J., Smith, D.B., Horton, J.D., 2014, A modified procedure for
#' mixture-model clustering of regional geochemical data: Applied Geochemistry,
#' vol. 51, p. 315-326, doi: http://dx.doi.org/10.1016/j.apgeochem.2014.10.011.
#'
#' Fraley, C., Raftery, A.E., 2002. Model-based clustering, discriminant
#' analysis, and density estimation: Journal of the American Statistical
#' Association, vol. 97, no. 458, p. 611-631.
#'
#' @examples
#' \dontrun{
#' fmmMode_mclust <- optFmm_mclust(transData, nPCS)
#' }
#'
#' @export
optFmm_mclust <- function(transData, nPCS, nCpuCores = 4, nSoln = NULL) {

  CalcClusters_Em <- function( theData, nPDFs, sampleSize, sampleSpace, nSoln ) {
    # Make mclust available to the processors. This function won't work otherwise.
    # So, I'm violating the principles in "R packages", p. 34, 82-84
    require(mclust, quietly = TRUE)
    theLogLikelihoods <- vector( mode="numeric" )
    nErrors <- 0
    for( i in 1:nSoln ) {
      S <- sample( sampleSpace, size = sampleSize )

      clusterResult <- tryCatch( mclust::Mclust( theData, nPDFs, modelNames=c("VVV"), initialization=list(subset=S) ), error=function(e){ e } )

      if( inherits( clusterResult, "error" ) ) {
        nErrors <- nErrors + 1
      } else {
        theLogLikelihoods <- append( theLogLikelihoods, clusterResult$loglik )
        if( max( theLogLikelihoods ) == clusterResult$loglik ) {
          bestClusterResult <- clusterResult
        }
        if( min( theLogLikelihoods ) == clusterResult$loglik ) {
          worstClusterResult <- clusterResult
        }
      }
    }

    return( list( theLogLikelihoods = theLogLikelihoods,
                  bestClusterResult = bestClusterResult,
                  worstClusterResult = worstClusterResult,
                  nErrors = nErrors ) )
  }

  if(nCpuCores > parallel::detectCores())
    stop("The number of requested cpu's must be <= the number of actual cpu's.")

  cl <- parallel::makeCluster( nCpuCores )
  parallel::clusterSetRNGStream( cl, 123 )

  nRows <- nrow( transData$robustPCs )
  sampleSize <- as.integer( 0.75 * nRows )
  sampleSpace <- 1:nRows
  nPDFs <- 2

  if(is.null(nSoln)) {
    nSoln <- 400 + ( nPDFs - 2 ) * 150 + ( nPCs - 4 ) * 30
  }

  nSolnPerCore <- as.integer( nSoln / nCpuCores )

  tmpResult <- parallel::clusterCall( cl, CalcClusters_Em,
                                      transData$robustPCs[,1:nPCs], nPDFs,
                                      sampleSize, sampleSpace, nSolnPerCore)

  parallel::stopCluster( cl )

  theLogLikelihoods <- tmpResult[[1]]$theLogLikelihoods
  bestClusterResult <- tmpResult[[1]]$bestClusterResult
  worstClusterResult <- tmpResult[[1]]$worstClusterResult
  nErrors <- tmpResult[[1]]$nErrors

  for( i in 2:nCpuCores ) {

    theLogLikelihoods <- append( theLogLikelihoods, tmpResult[[i]]$theLogLikelihoods )

    if( bestClusterResult$loglik < tmpResult[[i]]$bestClusterResult$loglik )
      bestClusterResult <- tmpResult[[i]]$bestClusterResult

    if( worstClusterResult$loglik > tmpResult[[i]]$worstClusterResult$loglik )
      worstClusterResult <- tmpResult[[i]]$worstClusterResult

    nErrors <- nErrors + tmpResult[[i]]$nErrors
  }

  return( list(
    nSoln = nSoln,
    theLogLikelihoods = theLogLikelihoods,
    bestClusterResult = bestClusterResult,
    worstClusterResult = worstClusterResult,
    nErrors = nErrors
  ))

}

#' @title Optimize the posterior pdf using package rstan.
#'
#' @description Optimize the posterior probability density function of the
#' finite mixture model. The model implementation and the optimization are
#' performed with package rstan.
#'
#' @param transData     List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}; the
#' documentation for function transformGcData includes a complete
#' description of container \code{transData}.
#' @param nPCS          Number of principal components that are used in the
#'                      finite mixture model (See Details).
#' @param tauBounds     Vector of length 2, containing the lower and upper
#'                      bounds for tau (See Details).
#'
#' @details
#' The parameters in the finite mixture model are estimated
#' from the robust principal components, which are
#' calculated with function gcTransform and stored as a matrix within
#' \code{transData}. All matrix elements must contain values;
#' that is, missing values are not allowed.
#'
#' The number of principal components, \code{nPCS}, means that
#' principal components 1, 2, ..., nPCs are used in the finite mixture model.
#' That is, all higher-order principal components are not used.
#'
#' In the finite mixture model, the standard deviations for
#' a pdf are stored in a vector tau. For each element in this vector, the
#' prior pdf is uniform; its lower and upper bounds are stored in
#' tauBounds[1] and tauBounds[2], respectively.
#'
#' Unlike function optFmm_mclust, only one solution is calculated with
#' function optFmm_rstan.
#'
#' The finite mixture model, which has two pdfs, involves five parameters:
#' \describe{
#'  \item{theta}{Proportion of population associated with pdf 1.}
#'  \item{mu1}{Mean vector for pdf 1.}
#'  \item{mu2}{Mean vector for pdf 2.}
#'  \item{Sigma1}{Covariance matrix for pdf 1.}
#'  \item{Sigma2}{Covariance matrix for pdf 2.}
#' }
#' There are also two additional variables associated with the model:
#' \describe{
#'  \item{g}{Conditional probability that a field sample is associated
#'  with pdf 1.}
#'  \item{log_lik}{The logarithm of the likelihood function.}
#' }
#'
#' @return A list that is returned from function optimizing within the
#' rstan package. If the optimization successfully locates a mode in
#' the posterior pdf, then the list has two elements:
#' @return \item{par}{List comprising the point estimates for all
#' model parameters and the two additional variables associated with the model.
#' Although there are many model parameters, only five are important and
#' are described in Details; the others may be safely ignored.}
#' @return \item{value}{The logarithm of the posterior pdf at the mode. Note
#' that this variable is not the logarithm of the likelihood function.}
#'
#' @references
#' Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A.,
#' and Rubin, D.B., 2014, Bayesian data analysis (3rd ed.),
#' CRC Press.
#'
#' Stan Development Team, 2015, Stan Modeling Language - User’s Guide
#' and Reference Manual, Version 2.6.0, available on line at
#' http://mc-stan.org/ (last accessed October 2015).
#'
#' @examples
#' \dontrun{
#' fmmMode_rstan <- optFmm_rstan(transData$robustPCs, nPCS,
#'                               tauBounds = c(0.001,20))
#' }
#'
#' @export
optFmm_rstan <- function(transData, nPCS, tauBounds = c(0.001, 50)) {

  # Load rstan because it calls internal functions.
  # So, I'm violating the principles in "R packages", p. 34, 82-84
  require(rstan, quietly = TRUE)

  stanData <- list( M = nPCs,
                    N = nrow(transData$robustPCs),
                    Z = transData$robustPCs[,1:nPCs],
                    tauBounds = tauBounds )

  pdfMode <- rstan::optimizing( GcClust:::sm, data = stanData, init="0",
                                refresh = 1, algorithm = "BFGS",
                                hessian = FALSE, as_vector = FALSE )
  return(pdfMode)
}

#' @title Sample the posterior pdf
#'
#' @description Sample the posterior probability density function of the
#' finite mixture model. The model implementation and the sampling are
#' performed with package rstan.
#'
#' @param transData     List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}; the
#' documentation for function transformGcData includes a complete
#' description of container \code{transData}.
#' @param nPCS          Number of principal components that are used in the
#'                      finite mixture model (See Details).
#' @param tauBounds     Vector of length 2, containing the lower and upper
#'                      bounds for tau (See Details).
#' @param nCpuCores      Number of central processing units (cpu's) that are
#'                      used
#'                      for parallel computations (See Details).
#' @param nChainsPerCore Number of chains that each cpu core calculates
#'                      (See Details).
#'
#' @details
#' The parameters in the finite mixture model are estimated
#' from the robust principal components, which are
#' calculated with function gcTransform and stored as a matrix within
#' \code{transData}. All matrix elements must contain values;
#' that is, missing values are not allowed.
#'
#' The number of principal components, \code{nPCS}, means that
#' principal components 1, 2, ..., nPCs are used in the finite mixture model.
#' That is, all higher-order principal components are not used.
#'
#' In the finite mixture model, the standard deviations for
#' a pdf are stored in a vector tau. For each element in this vector, the
#' prior pdf is uniform; its lower and upper bounds are stored in
#' tauBounds[1] and tauBounds[2], respectively.
#'
#' The posterior probability density function (pdf) is sampled to yield a chain.
#' The sampling for each chain is performed in parallel on multiple central
#' processing units (cpu's). Of course, the number of requested cpu's
#' \code{nCpuCores} must be less than or equal to the actual number.
#' The total number of
#' chains is \code{nCpuCores} times \code{nChainsPerCore}.
#'
#' Each chain has 1000 samples, of which 500 are warm-up and 500 are post
#' warm-up. These values are fixed.
#'
#' The finite mixture model, which has two pdfs, involves five parameters:
#' \describe{
#'  \item{theta}{Proportion of population associated with pdf 1.}
#'  \item{mu1}{Mean vector for pdf 1.}
#'  \item{mu2}{Mean vector for pdf 2.}
#'  \item{Sigma1}{Covariance matrix for pdf 1.}
#'  \item{Sigma2}{Covariance matrix for pdf 2.}
#' }
#' There are also two additional variables associated with the model:
#' \describe{
#'  \item{g}{Conditional probability that a field sample is associated
#'  with pdf 1.}
#'  \item{log_lik}{The logarithm of the likelihood function.}
#' }
#'
#' The optimzation may fail to converge to a mode of the posterior pdf because
#' the likelihood function is unbounded (Marin et al., 2005).
#' Consequently, checking the results of the optimization is important.
#'
#' @return A list with the following components is returned.
#' @return \item{stanSamples}{List for which the size equals the total number
#' of chains. Each element is itself a list containing the samples for
#' one chain. The samples are of the five parameters in the finite mixture
#' model and of the two additional variables, which are described in Details.
#' Only the 500 samples after warm-up are included, and their order is permuted.}
#' @return \item{stanSelectedTraces}{List for which the size equals the total number
#' of chains. Each element is itself a list containing the samples for
#' one chain. The samples are of \code{theta}, \code{mu1[1]} (the element 1
#' of vector \code{mu1}), \code{Sigma1[1,1]} (the element 1,1
#' of matrix \code{Sigma1}), \code{mu2[1]} (the element 1
#' of vector \code{mu2}), and \code{Sigma2[1,1]} (the element 1,1
#' of matrix \code{Sigma2}). Only the 500 samples after warm-up are included,
#' and their order is not permuted.}
#' @return \item{samplingParams}{List with 3 elements: nChains (the number
#' of chains), nTotalSamplesPerChain (the total number of samples per
#' chain), and nSamplesPerChain (the number of samples per chain after
#' warm-up.)}
#'
#' @references
#' Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A.,
#' and Rubin, D.B., 2014, Bayesian data analysis (3rd ed.),
#' CRC Press.
#'
#' Marin, J.M., Mengersen, K. and Robert, C.P., 2005,
#' Bayesian modelling and inference on mixtures of distributions,
#' in Handbook of Statistics 25, D. Dey and C.R. Rao (eds.),
#' Elsevier-Sciences.
#'
#' Stan Development Team, 2015, Stan Modeling Language - User’s Guide
#' and Reference Manual, Version 2.6.0, available on line at
#' http://mc-stan.org/ (last accessed October 2015).
#'
#' @examples
#' \dontrun{
#' fmmSamples <- sampleFmm(transData, nPCs, tauBounds = c(0.001,20),
#'                         nChainsPerCore = 5, nCpuCores = 7)
#' }
#'
#' @export
sampleFmm <- function(transData, nPCs, tauBounds = c(0.001, 50),
                      nSamplesPerChain = 1000,
                      nRepetitions = 2, nCpuCores = 4, procDir = NULL) {

  rstanParallelSampler <- function(stanData, nSamplesPerChain,
                                   rng_seed, nCpuCores, iRep, procDir ) {

    CL <- parallel::makeCluster(nCpuCores)

    parallel::clusterExport(cl = CL,
                            c("stanData", "nSamplesPerChain", "rng_seed",
                              "iRep", "procDir"),
                            envir=environment())

    sflist <- parallel::parLapply(CL, 1:nCpuCores, fun = function(cid) {

      # Make rstan available to the processors. This function won't work
      # otherwise. So, I'm violating the principles in "R packages",
      # p. 34, 82-84
      require(rstan, quietly = TRUE)

      rawSamples <- rstan::sampling(GcClust:::sm, data=stanData, init = "0",
                                    control = list(stepsize = 0.0001),
                                    chains = 1, iter = nSamplesPerChain,
                                    seed = rng_seed, chain_id = cid,
                                    pars=c("theta", "mu1", "mu2",
                                           "Sigma1", "Sigma2", "log_lik", "g"),
                                    save_dso = FALSE)

      fileName <- paste("RawSamples", iRep, "-", cid, ".dat", sep = "")
      save( rawSamples, file = paste(procDir, "\\", fileName, sep = ""))
      return(fileName)

    } )

    parallel::stopCluster(CL)
    return(sflist)
  }

  if(nCpuCores > parallel::detectCores())
    stop("The number of requested cpu's must be <= the number of actual cpu's.")

  stanData <- list( M = nPCs,
                    N = nrow(transData$robustPCs),
                    Z = transData$robustPCs[,1:nPCs],
                    tauBounds = tauBounds)

  fileNames <- NULL
  for(i in 1:nRepetitions) {
    rng_seed <- sample.int(.Machine$integer.max,1)
    tmp <- rstanParallelSampler(stanData, nSamplesPerChain, rng_seed,
                         nCpuCores, i, procDir )
    fileNames <- c(fileNames, unlist(tmp))
  }

  return(list(nChains = nRepetitions * nCpuCores,
              nSamplesPerChain = nSamplesPerChain,
              nSamplesPWU = nSamplesPerChain/2,
              fileNames = fileNames))
}

#' @title Plot the likelihoods of the modes
#'
#' @description Plot selected traces for each chain to assess whether
#' within-chain label switching has occurred.
#'
#' @param fmmSamples
#' List containing samples of the posterior pdf and
#' related information. This list is return by function
#' \code{\link{sampleFmm}}; the documentation for function sampleFmm includes
#' a complete description of container \code{fmmSamples}.
#'
#' @details
#' Three plots are generated:
#' \itemize{
#'  \item The first plot comprises two traces: A trace for element [1]
#'  of the mean vector for pdf 1 (mu1[1]), and another trace for
#'  element [1] of the mean vector for pdf 2 (mu2[1]).
#'  \item The second plot also comprises two traces: A trace for element [1,1]
#'  of the covariance matrix for pdf 1 (Sigma1[1,1]), and another trace for
#'  element [1,1] of the covariance matrix for pdf 2 (Sigma2[1,1]).
#'  \item The third plot has one trace of model proportion associated with
#'  pdf 1 (theta)
#' }
#'
#' @examples
#' \dontrun{
#' plotTraces(fmmSamples)
#' }
#'
#' @export
plotModeLikelihoods <- function(fmmMode_mclust, fmmMode_rstan,
                                binwidthScale = 30,
                                rugColour = "green",
                                rugSize = 1) {

  mclustData <- data.frame(logLikelihood = fmmMode_mclust$theLogLikelihoods)
  theRange <- range(mclustData$logLikelihood)
  p <- ggplot2::ggplot(mclustData,
                       ggplot2::aes(x = logLikelihood),
                       environment = environment()) +
    ggplot2::geom_histogram(binwidth = (theRange[2] - theRange[1])/binwidthScale) +
    ggplot2::xlab("Log-likelihood") +
    ggplot2::ylab("Count")

  if(!is.na(fmmMode_rstan$par$log_lik)) {
    rstanData <- data.frame(logLikelihood = fmmMode_rstan$par$log_lik)
    p <- p + ggplot2::geom_rug(data = rstanData, colour=rugColour, size = rugSize )
  }

  print(p)

}



#' @title Plot selected traces
#'
#' @description Plot selected traces for each chain to assess whether
#' within-chain label switching has occurred.
#'
#' @param fmmSamples
#' List containing samples of the posterior pdf and
#' related information. This list is return by function
#' \code{\link{sampleFmm}}; the documentation for function sampleFmm includes
#' a complete description of container \code{fmmSamples}.
#'
#' @details
#' Three plots are generated:
#' \itemize{
#'  \item The first plot comprises two traces: A trace for element [1]
#'  of the mean vector for pdf 1 (mu1[1]), and another trace for
#'  element [1] of the mean vector for pdf 2 (mu2[1]).
#'  \item The second plot also comprises two traces: A trace for element [1,1]
#'  of the covariance matrix for pdf 1 (Sigma1[1,1]), and another trace for
#'  element [1,1] of the covariance matrix for pdf 2 (Sigma2[1,1]).
#'  \item The third plot has one trace of model proportion associated with
#'  pdf 1 (theta)
#' }
#'
#' @examples
#' \dontrun{
#' plotTraces(fmmSamples)
#' }
#'
#' @export
plotSelectedTraces <- function(samplePars, procDir) {

  N <- length(samplePars$fileNames)
  for(i in 1:N){

    devAskNewPage( ask = TRUE )

    load( paste(procDir, "\\", samplePars$fileNames[i], sep = ""))
    theSamples <- rstan::extract(rawSamples, permuted = FALSE)

    df1 <- data.frame(indices = 1:samplePars$nSamplesPWU,
                     y = theSamples[, 1, "theta"] )
    p1 <- ggplot2::ggplot( df1,
                           ggplot2::aes(x = df1$indices, y = df1$y),
                           environment = environment() ) +
      ggplot2::geom_line() +
      ggplot2::xlab("Sample index") +
      ggplot2::ylab("theta")

    df2 <- data.frame(indices = 1:samplePars$nSamplesPWU,
                     y1 = theSamples[, 1, "mu1[1]"],
                     y2 = theSamples[, 1, "mu2[1]"])
    p2 <- ggplot2::ggplot( df2,
                           ggplot2::aes(x = df2$indices),
                           environment = environment() ) +
      ggplot2::geom_line(ggplot2::aes(y = df2$y1, colour = "1")) +
      ggplot2::geom_line(ggplot2::aes(y = df2$y2, colour = "2")) +
      ggplot2::scale_color_manual("Pdf", values = c("1" = "blue", "2" = "red")) +
      ggplot2::xlab("Sample index") +
      ggplot2::ylab("Element [1] of mean vectors")

    df3 <- data.frame(indices = 1:samplePars$nSamplesPWU,
                      y1 = theSamples[, 1, "Sigma1[1,1]"],
                      y2 = theSamples[, 1, "Sigma2[1,1]"])
    p3 <- ggplot2::ggplot( df3,
                           ggplot2::aes(x = df3$indices),
                           environment = environment() ) +
      ggplot2::geom_line(ggplot2::aes(y = df3$y1, colour = "1")) +
      ggplot2::geom_line(ggplot2::aes(y = df3$y2, colour = "2")) +
      ggplot2::scale_color_manual("Pdf", values = c("1" = "blue", "2" = "red")) +
      ggplot2::xlab("Sample index") +
      ggplot2::ylab("Element [1, 1] of variance matrices")

    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,1)))
    print(p1, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
    print(p2, vp=grid::viewport(layout.pos.row=2, layout.pos.col=1))
    print(p3, vp=grid::viewport(layout.pos.row=3, layout.pos.col=1))

    # devAskNewPage( ask = options()$device.ask.default )
  }
}

#' @title Plot point statistics
#'
#' @description Plot point statistics for each chain, the mode calculated with
#' function mclust, and the mode calculated with rstan.
#'
#' @param fmmSamples
#' List containing samples of the posterior pdf and
#' related information. This list is return by function
#' \code{\link{sampleFmm}}; the documentation for function sampleFmm includes
#' a complete description of container \code{fmmSamples}.
#' @param fmmMode_mclust  List containing information about the optimization
#' of the posterior pdf for which the model and the optimization are
#' performed by package mclust. This list is returned by
#' function \code{\link{optFmm_mclust}}; the documentation for
#' function optFmm_mclust includes a complete description of container
#' \code{fmmMode_mclust}.
#'
#' @param fmmMode_rstan List containing information about the optimization
#' of the posterior pdf for which the model and the optimization are
#' performed by package rstan. This list is returned by
#' function \code{\link{optFmm_rstan}}; the documentation for
#' function optFmm_rstan includes a complete description of container
#' \code{fmmMode_rstan}.
#' @param excludedChains Vector with the indices of the chains for which
#' the point statistics are not used to calculate plot ranges. See Details.
#'
#' @details
#' The point statistics are calculated and plotted for five model
#' parameters and two addtional variables:
#' \itemize{
#'  \item The first element of the mean vector for pdf 1, which is
#' designated "mu1[1]."
#'  \item The first element of the mean vector for pdf 2, which is
#' designated "mu2[1]."
#'  \item The first element of the covariance matrix for pdf 1, which is
#' designated "Sigma1[1,1]."
#'  \item The first element of the covariance matrix for pdf 2, which is
#' designated "Sigma2[1,1]."
#'  \item The proportion in the finite mixture model associated with pdf 1,
#' which is designated "theta."
#'  \item The logarithm of the likelihood, which is
#' designated "Log-likelihood."
#'  \item The logarithm of the posterior pdf, which is
#' designated "Log-Posterior Pdf."
#' }
#'
#' To prevent a point statistic from being included in the calculation of
#' the plot range, the index of the associated chain is specified in
#' argument \code{excludedChains}.
#'
#' @examples
#' \dontrun{
#' plotPointStats(fmmSamples, fmmMode_mclust, fmmMode_rstan,
#'                excludedChains = NULL)
#' }
#'
#' @export
plotPointStats <- function(samplePars, fmmMode_mclust, fmmMode_rstan,
                           excludedChains = NULL ) {

  N <- length(samplePars$fileNames)

  statNames <-
    c("mu1[1]","mu2[1]","Sigma1[1,1]","Sigma2[1,1]","theta","log_lik",
      "lp__")

  probs <- c(0.025, 0.50, 0.975)
  theQuantiles <- array(NA_real_,
                        dim = c(length(probs), length(statNames), N),
                        dimnames = list(probs, statNames, paste("Chain",1:N)) )

  for(i in 1:N) {
    load( paste(procDir, "\\", samplePars$fileNames[i], sep = ""))
    chainStats <- rstan::summary(rawSamples, pars = statNames, probs = probs)$c_summary

    # excludes the first two columns, which contain the mean and sd
    theQuantiles[, , i] <- t(chainStats[statNames, -(1:2), 1])
  }

  Internal1 <- function( X, excludedChains,
                         parMode_mclust, parMode_rstan, yLabel ) {
    if(is.null(excludedChains)) {
      yRange <- range( X, parMode_mclust, parMode_rstan, na.rm = TRUE )
    } else {
      yRange <- range( X[, -excludedChains], parMode_mclust,
                       parMode_rstan, na.rm = TRUE )
    }

    solutions <- c( as.character(1:ncol(X)), "mclust", "rstan")

    df <- data.frame(x = factor(solutions, levels = solutions),
                      y = c( X["0.5", ], parMode_mclust, parMode_rstan ),
                      ymin = c( X["0.025",], parMode_mclust, parMode_rstan ),
                      ymax = c( X["0.975",], parMode_mclust, parMode_rstan ) )
    p <- ggplot2::ggplot(df,
                          ggplot2::aes(x = x, y = y, ymin = ymin, ymax = ymax),
                          environment = environment() ) +
      ggplot2::geom_pointrange() +
      ggplot2::ylim(yRange[1], yRange[2]) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, vjust=0.125)) +
      ggplot2::xlab("Solution") +
      ggplot2::ylab(yLabel)

    return(p)
  }

  Internal2 <- function( X, Y, excludedChains,
                         parMode_mclust, parMode_rstan, yLabel ) {
    if(is.null(excludedChains)) {
      yRange <- range( X, Y, parMode_mclust, parMode_rstan, na.rm = TRUE )
    } else {
      yRange <- range( X[, -excludedChains], Y[, -excludedChains],
                       parMode_mclust, parMode_rstan, na.rm = TRUE )
    }

    solutions <- c( as.character(1:ncol(X)), "mclust", "rstan")

    df <- data.frame(x = factor(solutions, levels = solutions),
                     y1 = c( X["0.5", ], parMode_mclust[1], parMode_rstan[1] ),
                     ymin1 = c( X["0.025",], parMode_mclust[1], parMode_rstan[1] ),
                     ymax1 = c( X["0.975",], parMode_mclust[1], parMode_rstan[1] ),
                     y2 = c( Y["0.5", ], parMode_mclust[2], parMode_rstan[2] ),
                     ymin2 = c( Y["0.025",], parMode_mclust[2], parMode_rstan[2] ),
                     ymax2 = c( Y["0.975",], parMode_mclust[2], parMode_rstan[2] ) )
    p <- ggplot2::ggplot(df, environment = environment() ) +
      ggplot2::geom_pointrange(ggplot2::aes(x = x, y = y1, ymin = ymin1, ymax = ymax1),
                               colour = "blue") +
      ggplot2::geom_pointrange(ggplot2::aes(x = x, y = y2, ymin = ymin2, ymax = ymax2),
                               colour = "red") +
      ggplot2::ylim(yRange[1], yRange[2]) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, vjust=0.125)) +
      ggplot2::xlab("Solution") +
      ggplot2::ylab(yLabel)

    return(p)
  }

  fmmMode_mclust <- fmmMode_mclust$bestClusterResult

  p1 <- Internal2(theQuantiles[, "mu1[1]", ],
                  theQuantiles[, "mu2[1]", ],
                  excludedChains,
                  fmmMode_mclust$parameters$mean[1,1:2],
                  c(fmmMode_rstan$par$mu1[1],fmmMode_rstan$par$mu2[1]),
                 "Element [1] of the mean vectors")
  p2 <- Internal2(theQuantiles[, "Sigma1[1,1]", ],
                  theQuantiles[, "Sigma2[1,1]", ],
                  excludedChains,
                  fmmMode_mclust$parameters$variance$sigma[1,1,1:2],
                  c(fmmMode_rstan$par$Sigma1[1,1],fmmMode_rstan$par$Sigma2[1,1]),
                 "Element [1,1] of the variance matrices")
  p3 <- Internal1(theQuantiles[, "theta", ], excludedChains,
               fmmMode_mclust$parameters$pro[1],
               fmmMode_rstan$par$theta,
               "theta")
  p4 <- Internal1(theQuantiles[, "log_lik", ], excludedChains,
               fmmMode_mclust$loglik,
               fmmMode_rstan$par$log_lik,
               "Log-likelihood")

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2,2)))
  print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p3, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p4, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))

}

#' @title Combine selected chains
#'
#' @description Combine selected chains, creating one set of samples for each
#' parameter in the finite mixture model.
#'
#' @param fmmSamples      List, which includes the Monte Carlo chains and
#'                        other data, returned by function
#'                        \code{\link{sampleFmm}}; the
#'                        documentation for fmmSample includes a complete
#'                        description of fmmSamples.
#' @param selectedChains  Dataframe listing the indices of the chains that
#'  are combined. (See Details).
#'
#' @details
#' Argument \code{selectedChains} is a dataframe with two columns. The first
#' column, for which the heading is "Chain," comprises the indices of the chains
#' that will be combined. The second column, for which the heading is
#' "isSwitched", comprises logical values (that is, TRUE or FALSE), which
#' indicate whether the variables in the associated chain are switched.
#'
#' To understand \code{selectedChains} consider two examples. (1) Assume that
#' row 1 of \code{selectedChains} comprises the values 2 and FALSE. These
#' values
#' indicate that chain 2 will be combined and that the variables in the
#' finite mixture model are not switched. That is, variable mu1 in the
#' chain is assigned to mu1, variable mu2 in the chain is assiged to mu2, and
#' so on. (2) Assume that
#' row 2 of \code{selectedChains} comprises the values 4 and TRUE. These
#' values
#' indicate that chain 4 will be combined and that the variables in the
#' finite mixture model are switched. That is, variable mu1 in the
#' chain is assigned to mu2, variable mu2 in the chain is assiged to mu1, and
#' so on.
#'
#' @return A list with the following components is returned.
#' @return \item{theta}{Vector containing samples of the proportion
#' associated with pdf 1 that is within the finite mixture model.}
#' @return \item{mu1}{Matrix containing Monte Carlo
#' samples of the mean vector for pdf 1. The first dimension of the matrix
#' pertains to indices of the samples; the second dimension pertains to
#' indices of the mean vector.}
#' @return \item{mu2}{Matrix containing Monte Carlo
#' samples of the mean vector for pdf 2. Its structure is identical to
#' that for \code{mu1}.}
#' @return \item{Sigma1}{Array containing Monte Carlo samples of the
#' covariance matrix
#' for pdf 1. The first dimension of the array pertains to indices of
#' the samples; the second and the third dimensions pertain to indices of
#' the covariance matrix.}
#' @return \item{Sigma2}{Array containing samples of the covariance matrix
#' for pdf 2. Its structure is identical to that for \code{Sigma1}.}
#' @return \item{g}{Matrix containing Monte Carlo samples of the
#' conditional probabilites that the field samples are associated with pdf 1.
#' The first dimension of the matrix
#' pertains to indices of the Monte Carlo samples;
#' the second dimension pertains to
#' indices of the field samples. }
#'
#' @examples
#' \dontrun{
#' combinedChains <- combineChains(fmmSamples$stanSamples, selectedChains)
#' }
#'
#' @export
combineChains <- function(samplePars, selectedChains) {

  sfList <- vector(mode = "list")
  for(k in 1:nrow(selectedChains)) {
    iChain <- selectedChains[k, "Chain"]
    load( paste(procDir, "\\", samplePars$fileNames[iChain], sep = ""))

    if(selectedChains[k, "isSwitched"] == TRUE) {

      rawSamples@sim$samples[[1]]$theta <- 1 - rawSamples@sim$samples[[1]]$theta

      N <- rawSamples@par_dims$mu1
      for(i in 1:N) {
        var1 <- paste("rawSamples@sim$samples[[1]]$\"mu1[", i, "]\"", sep="")
        tmp1 <- eval(parse(text = var1))

        var2 <- paste("rawSamples@sim$samples[[1]]$\"mu2[", i, "]\"", sep="")
        tmp2 <- eval(parse(text = var2))

        eval(parse(text = paste(var1, " <- tmp2")))

        eval(parse(text = paste(var2, " <- tmp1")))

        for(j in 1:N) {
          var1 <- paste("rawSamples@sim$samples[[1]]$\"Sigma1[", i, ",", j, "]\"", sep="")
          tmp1 <- eval(parse(text = var1))

          var2 <- paste("rawSamples@sim$samples[[1]]$\"Sigma2[", i, ",", j, "]\"", sep="")
          tmp2 <- eval(parse(text = var2))

          eval(parse(text = paste(var1, " <- tmp2")))

          eval(parse(text = paste(var2, " <- tmp1")))

        }
      }

      M <- rawSamples@par_dims$g
      for(i in 1:M) {
        var <- paste("rawSamples@sim$samples[[1]]$\"g[", i, "]\"", sep="")
        tmp <- eval(parse(text = var))
        eval(parse(text = paste(var, " <- 1 - tmp")))
      }


    }
    sfList[[k]] <- rawSamples

  }

  return(rstan::sflist2stanfit(sfList))
}

#' @title Calculate test quantities for model checking
#'
#' @description Calculate test quantities for posterior predictive checking
#' of the parameters in the finite mixture model. The test quantities pertain
#' to the observed data that have undergone both the isometric log-ratio
#' transform and the robust, principal component transform.
#'
#' @param transData     List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}; the
#' documentation for function transformGcData includes a complete
#' description of container \code{transData}.
#' @param nPCS          Number of principal components that are used in the
#'                      finite mixture model (See Details).
#' @param combinedChains     List in which each element contains the
#' Monte Carlo samples of a parameter from the finite mixture model.
#' This list is return by function \code{\link{combineChains}}; the
#'                        documentation for function combineChains includes a
#'                        complete description of container
#'                        \code{combinedChains}.
#'
#' @details
#' The number of principal components is specified by \code{nPCS}. For example,
#' if \code{nPCS} equals 7, then components 1, 2, ..., 7 are used in the
#' finite mixture model. These components are in columns 1 through 7 of
#' matrix \code{robustPCs}.
#'
#' @return A list with the following test quantities is returned.
#' @return \item{mu1}{Vector containing the test quantity for the mean for
#' pdf 1.}
#' @return \item{mu2}{Vector containing the test quantity for the mean for
#' pdf 1.}
#' @return \item{Sigma1}{Matrix containing the test quantity for the covariance
#' matrix for pdf 1.}
#' @return \item{Sigma2}{Matrix containing the test quantity for the covariance
#' matrix for pdf 1.}
#'
#' @examples
#' \dontrun{
#' obsTestQuantities <- calcObsTestQuantities(transData, nPCS, combinedChains)
#' }
#'
#'
#' @export
calcObsTestQuantities <- function(transData, nPCS, combinedChains) {

  # associated with pdf 1
  condProb1 <- rstan::extract(combinedChains, pars="g")$g
  S1 <- cov.wt(transData$robustPCs[,1:nPCs], wt = colMeans(condProb1) )

  # associated with pdf 2
  condProb2 <- 1 - condProb1
  S2 <- cov.wt(transData$robustPCs[,1:nPCs], wt = colMeans(condProb2) )

  return(list(mu1 = S1$center,
              mu2 = S2$center,
              Sigma1 = S1$cov,
              Sigma2 = S2$cov))

}


#' @title Plot test quantities for model checking
#'
#' @description The test quantities are the elements of the mean vector,
#' the standard deviations from the diagonal of the covariance matrix, or
#' the correlation matrix derived from the diagonal matrix.
#' These statistics pertain to one of the two pdfs in the finite mixture model.
#' One set of test quantities
#' are calcuated from the data, which have undergone both the isometric
#' log-ratio transform and the robust, principal component transform. The other
#' set of test quantities are calculated from the Monte Carlo samples.
#' This function plots the test quantities so that they can be compared;
#' this comparison is a graphical posterior predictive check of the model.
#'
#' @param combinedChains    List in which each element contains the
#' Monte Carlo samples of a parameter from the finite mixture model.
#' This list is return by function \code{\link{combineChains}}; the
#'                        documentation for function combineChains includes a
#'                        complete description of container
#'                        \code{combinedChains}.
#'
#' @param obsTestQuantities     List containing the test quantities for
#' posterior predicted checks of the finit mixture model.
#' This list is return by function \code{\link{calcObsTestQuantities}}; the
#'                        documentation for function calcObsTestQuantities
#'                        includes a
#'                        complete description of container
#'                        \code{obsTestQuantities}.
#'
#' @param pdf   Index for the pdf in the finite mixture model. It must be
#' either 1 or 2.
#' @param testQuantity   The test quantity that will be plotted. It must be
#' "mean", "sd" (for standard deviations), or "cor" (for the correlation
#' matrix).
#' @param intervalPercentage Interval for the distributions of the test
#' quantity. Typical values might be 50, 90, or 95. Not used when
#' argument testQuantity is "cor".
#' @param arePValuesPlotted   Logical variable indicating whether the p-values
#' are plotted.
#'
#' @details
#' The p-value for the test quantity is defined as the probability that the
#' replicated could be more extreme that the observed data, as measured
#' by the test quantity (Gelman et al., 2014, p. 146). It mathematical
#' terms,
#'     pvalue = Pr( T.rep >= T)
#' where T.rep is the test quantity calculated for the replicated data and
#' T is the test quantity calculated for the observed data.
#'
#' Frequently, the p-value is close to 1 when T is in the left tail of the
#' distribution for T.rep. This can confuse the interpretation of the p-value.
#' In such situations, the mathematical definition is modified slightly:
#'    pvalue = Pr( T.rep < T)
#' (Gelman et al., 2014, p. 148). Consequently, the calculated p-value is
#' always less than 0.5, and it may be interpreted in the standard way.
#'
#' If the testQuantity is either "mean" or "sd", then the p-values are
#' scaled. The reason is that, if the a p-value is expressed as "0.12",
#' then four characters are needed to print it. The consequence is that
#' the p-values, which are printed along the top of the plot, would overlap
#' one another, making them difficult to read. Instead, the
#' p-values are multiplied by 100 and then rounded to the nearest whole number.
#' Consequently the scaled p-values require only two characters to print them.
#' (There is one exception: If the scaled p-value is 100, then three
#' characters are needed to print it. But this exception should occur rarely.)
#'
#' @references
#' Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B.,
#' Vehtari, A., and Rubin, D.B., 2014, Bayesian data analysis (3rd ed.):
#' CRC Press.
#'
#' @examples
#' \dontrun{
#' plotTestComparison( combinedChains, obsTestQuantities,
#'                      pdf = 1, testQuantity = "mean",
#'                      arePValuesPlotted = TRUE )
#'
#' plotTestComparison( combinedChains, obsTestQuantities,
#'                      pdf = 1, testQuantity = "sd",
#'                      arePValuesPlotted = TRUE )
#'
#' plotTestComparison( combinedChains, obsTestQuantities,
#'                      pdf = 1, testQuantity = "cor",
#'                      arePValuesPlotted = TRUE )
#' }
#'
#' @export
plotTestComparison <- function( combinedChains, obsTestQuantities,
                                pdf = 1, testQuantity = "mean",
                                intervalPercentage = 95,
                                arePValuesPlotted = FALSE ) {

  ExtractSd <- function(Sigma) {
    tmp <- dim(Sigma)
    nMcSamples <- tmp[1]
    M <- tmp[2]
    sampSd <- matrix(NA_real_, nrow = nMcSamples, ncol = M)
    for(i in 1:nMcSamples) {
      sampSd[i, ] <- sqrt(diag(Sigma[i, ,]))
    }
    return(sampSd)
  }

  PlotModelCheck <- function( testStat, sampModelStat, intervalPercentage,
                              arePValuesPlotted, plotTitle ){

    origPar <- par( lend="butt" )
    on.exit(par(origPar), add = TRUE )

    M <- ncol(sampModelStat)
    nMcSamples <- nrow(sampModelStat)

    tailPercentage <- 0.5*(100.0-intervalPercentage)
    interval <- c(tailPercentage,100.0-tailPercentage)/100.0

    yLimits <- range(sampModelStat, testStat )

    plot( c(1,M), yLimits, type="n",
          xlim=c(0.5,M+0.5), ylim=yLimits,
          xlab="Component", ylab="ilr and pc-transformed concentration (no units)",
          main=plotTitle )
    tmp <- par()$usr
    rect(tmp[1], tmp[3], tmp[2], tmp[4], col="gray90")

    for(j in 1:M) {

      lines( c(j-0.5,j+0.5), rep.int(testStat[j],2), lwd=3, col="orange" )

      theQuantiles <- quantile(sampModelStat[, j], probs=c(interval,0.50),
                               names = FALSE)
      lines( c(j,j), theQuantiles[1:2], lwd=1 )
      lines( c(j-0.5,j+0.5), rep.int(theQuantiles[3],2), lwd=1 )
    }

    if(arePValuesPlotted) {
      pValues <- vector(mode="numeric", length=M )
      for(j in 1:M) {
        # defined in Gelman et al., p. 146
        tmp <- sum(sampModelStat[, j] > testStat[j]) / nMcSamples
        pValues[j] <- ifelse(tmp <= 0.5, tmp, 1-tmp)
      }
      scaledPValues <- round(pValues*100, digits=0)

      mtext( scaledPValues,
             side=3, line=0, at=1:M, cex=0.6,
             col=ifelse(scaledPValues < 10, "red", "black") )

    }
  }

  PlotCorCheck <- function( testSigma, sampModelSigma,
                              arePValuesPlotted, plotTitle ){

    Internal1 <- function(testSigma, sampModelSigma) {
      X <- cov2cor(testSigma)
      Y <- cov2cor(apply(sampModelSigma, c(2,3), median))
      Z <- X
      Z[lower.tri(Z)] <- Y[lower.tri(Y)]
      Z <- reshape2::melt(Z)
      w <- ggplot2::ggplot(Z, ggplot2::aes(Var2, Var1))+
        ggplot2::geom_tile(data=Z, ggplot2::aes(fill=value), color="white")+
        ggplot2::scale_fill_gradient2(low="blue", high="red", mid="white",
                                      midpoint=0, limit=c(-1,1),
                                      name="Correlation\n(Pearson)")+
        ggplot2::xlab("Component") +
        ggplot2::ylab("Component") +
        ggplot2::coord_equal()

      return(w)
    }

    Internal2 <- function(testSigma, sampModelSigma) {
      testCor <- cov2cor(testSigma)

      theDimensions <- dim(sampModelSigma)
      nMcSamples <- theDimensions[1]
      M <- theDimensions[2]

      sampModelCor <- array(NA_real_, dim=theDimensions)
      for(i in 1:nMcSamples) {
        sampModelCor[i, , ] <- cov2cor(sampModelSigma[i, , ])
      }


      pValues <- matrix(NA_real_, nrow = M, ncol = M )
      for(i in 1:(M-1)) {
        for(j in (i+1):M) {
          tmp <- sum(sampModelCor[ , i, j ] > testCor[i,j]) / nMcSamples
          pValues[i,j] <- ifelse(tmp <= 0.5, tmp, 1-tmp)
        }
      }

      Z <- reshape2::melt(pValues)
      Z <- na.omit(Z)

      w <- ggplot2::ggplot(Z, ggplot2::aes(Var2, Var1))+
        ggplot2::geom_tile(data=Z, ggplot2::aes(fill=value), color="white")+
        ggplot2::scale_fill_gradient(limit=c(min(Z$value),0.5),
                                     high = "#132B43", low = "#56B1F7",
                                      name="P-value")+
        ggplot2::xlab("Component") +
        ggplot2::ylab("Component") +
        ggplot2::coord_equal()

      return(w)


    }

    x <- Internal1(testSigma, sampModelSigma)

    w <- Internal2(testSigma, sampModelSigma)

    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,2)))
    print(x, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
    print(w, vp=grid::viewport(layout.pos.row=1, layout.pos.col=2))

  }

  if(pdf == 1 && testQuantity == "mean") {
    PlotModelCheck( obsTestQuantities$mu1,
                    rstan::extract(combinedChains, pars="mu1")$mu1,
                    intervalPercentage, arePValuesPlotted,
                    plotTitle = "Mean for pdf 1")
  } else if(pdf == 1 && testQuantity == "sd") {
    PlotModelCheck( sqrt(diag(obsTestQuantities$Sigma1)),
                    ExtractSd(rstan::extract(combinedChains, pars="Sigma1")$Sigma1),
                    intervalPercentage, arePValuesPlotted,
                    plotTitle = "Sd for pdf 1")
  } else if(pdf == 1 && testQuantity == "cor") {
    PlotCorCheck( obsTestQuantities$Sigma1,
                  rstan::extract(combinedChains, pars="Sigma1")$Sigma1,
                  arePValuesPlotted,
                  plotTitle = "Correlation for pdf 1")
  } else if(pdf == 2 && testQuantity == "mean") {
    PlotModelCheck( obsTestQuantities$mu2,
                    rstan::extract(combinedChains, pars="mu2")$mu2,
                    intervalPercentage, arePValuesPlotted,
                    plotTitle = "Mean for pdf 2")
  } else if(pdf == 2 && testQuantity == "sd") {
    PlotModelCheck( sqrt(diag(obsTestQuantities$Sigma2)),
                    ExtractSd(rstan::extract(combinedChains, pars="Sigma2")$Sigma2),
                    intervalPercentage, arePValuesPlotted,
                    plotTitle = "Sd for pdf 2")
  } else if(pdf == 2 && testQuantity == "cor") {
    PlotCorCheck( obsTestQuantities$Sigma2,
                  rstan::extract(combinedChains, pars="Sigma2")$Sigma2,
                  arePValuesPlotted,
                  plotTitle = "Correlation for pdf 1")
  } else {
    msg <- paste( "Argument pdf must be 1 or 2, ",
                  "and argument testQuantity must be ",
                  "mean, sd, or cor.", sep="")
    stop(msg)
  }

}


#' @title Back transform to the simplex
#'
#' @description Back transform the mean vectors and covariance matrices
#' of the finite mixture model
#' to the correpsonding quantities for simplex
#' (that is, compostional means and variation matrices).
#'
#' @param gcData   SpatialPointsDataFrame storing the geochemical
#'                 concentrations and locations of each field sample
#' @param kappa    Constant-sum value for the geochemical concentrations.
#'                 Typically \code{kappa} equals 1,000,000.
#' @param nPCS          Number of principal components that were used in the
#'                      finite mixture model.
#'
#' @param transData     List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}; the
#' documentation for function transformGcData includes a complete
#' description of container \code{transData}.
#' @param combinedChains List in which each element contains the
#' Monte Carlo samples of a parameter from the finite mixture model.
#' This list is return by function \code{\link{combineChains}}; the
#'                        documentation for function combineChains includes a
#'                        complete description of container
#'                        \code{combinedChains}.
#'
#'
#'
#' @return A list with the following components is returned.
#' @return \item{compMean1}{Matrix containing Monte Carlo
#' samples of the compositional mean vector for pdf 1.
#' The first dimension of the matrix
#' pertains to indices of the samples; the second dimension pertains to
#' indices of the compositional mean vector.}
#' @return \item{compMean2}{Matrix containing Monte Carlo
#' samples of the mean vector for pdf 2. Its structure is identical to
#' that for \code{compMean1}.}
#' @return \item{varMatrix1}{Array containing Monte Carlo samples of the
#' variation matrix
#' for pdf 1. The first dimension of the array pertains to indices of
#' the samples; the second and the third dimensions pertain to indices of
#' the variation matrix.}
#' @return \item{varMatrix2}{Array containing samples of the variation matrix
#' for pdf 2. Its structure is identical to that for \code{varMatrix1}.}
#' @return \item{theta}{Vector containing samples of the proportion
#' associated with pdf 1 that is within the finite mixture model.}
#' @return \item{g}{Matrix containing Monte Carlo samples of the
#' conditional probabilites that the field samples are associated with pdf 1.
#' The first dimension of the matrix
#' pertains to indices of the Monte Carlo samples;
#' the second dimension pertains to
#' indices of the field samples. }
#'
#' @references
#' Hron, K., Filzmoser, P. (2015) Exploring Compositional Data with the Robust
#' Compositional Biplot. In: Carpita, M., Brentari, E., Qannari, M. (eds.)
#' Advances in Latent Variables. Methods, Models and Applications.
#' Springer, Heidelberg, pp 219-226.
#'
#' Pawlowsky-Glahn, V., Egozcue, J.J., and Tolosana-Delgado, R., 2015, Modeling
#' and analysis of compositional data: John Wiley and Sons, Ltd.
#'
#' @examples
#' \dontrun{
#' simplexModPar <- backTransform( gcData, kappa, nPCs,
#'                                 transData, combinedChains)
#' }
#'
#' @export

backTransform <- function(gcData, kappa, nPCs,
                          transData, combinedChains) {

  elementNames <- names(gcData)

  nMcSamples <- nrow(combinedChains$mu1)

  # dimension within the simplex
  D <- length(elementNames)

  # subset now, so that it won't be done repeatedly
  rEV <- transData$robustEigenvectors[,1:nPCs]
  rEVt <- t(rEV)

  # undo the rotation and translation, which are associated
  # with the principal component transformation
  tmp.center <- rep.int( 1, nMcSamples ) %o% transData$robustIlrCenter

  mu1.tmp <- combinedChains$mu1 %*%  rEVt + tmp.center
  mu2.tmp <- combinedChains$mu2 %*%  rEVt + tmp.center

  # Back-transformation to concentrations
  compMean1 <- CoDA::CalcInvIlr( mu1.tmp, t(transData$Psi), kappa=kappa )
  compMean2 <- CoDA::CalcInvIlr( mu2.tmp, t(transData$Psi), kappa=kappa )

  colnames(compMean1) <- elementNames
  colnames(compMean2) <- elementNames

  varMatrix1 <- array( NA_real_, dim=c(nMcSamples,D,D),
                       dimnames = list(NULL,elementNames,elementNames))
  varMatrix2 <- varMatrix1

  for(i in 1:nMcSamples) {

    # undo the rotation, which is associated
    # with the principal component transformation
    Sigma1.tmp <- rEV %*% combinedChains$Sigma1[i, , ] %*% rEVt
    Sigma2.tmp <- rEV %*% combinedChains$Sigma2[i, , ] %*% rEVt

    # Back-transformation to concentrations/simplex
    varMatrix1[i, , ] <- CoDA::InvCovTransform( Sigma1.tmp, transData$Psi )
    varMatrix2[i, , ] <- CoDA::InvCovTransform( Sigma2.tmp, transData$Psi )
  }

  return( list( compMean1 = compMean1,
                compMean2 = compMean2,
                varMatrix1 = varMatrix1,
                varMatrix2 = varMatrix2,
                theta = combinedChains$theta,
                g = combinedChains$g ))

}

#' @title Plot standardized compositional means
#'
#' @description Plot standardized compositional means for the pdfs in the
#' finite mixture model. The standardization is based on the geochemical
#' concentrations that were used in the finite mixture model.
#'
#' @param simplexStats   SpatialPointsDataFrame storing the geochemical
#'                 concentrations that were used in the finite mixture model.
#' @param kappa    Constant-sum value for the geochemical concentrations.
#'                 Typically \code{kappa} equals 1,000,000.
#' @param simplexModPar
#' List containing Monte Carlo samples of the parameters in the finite mixture
#' model.
#' These parameters are expressed in terms of their equivalent values in the
#' simplex.
#' This list is return by function \code{\link{backTranform}}; the
#'                        documentation for function backTranform includes a
#'                        complete description of container
#'                        \code{simplexModPar}.
#' @param elementOrder Vector specifying the order in which the elements are
#' plotted.
#' @param intervalPercentage Interval for the distributions of the standardized
#' compositional means. Typical values are 50, 90, or 95.
#' @param symbolSize The size of the plotting symbol.
#'
#' @details
#' ???
#'
#' @references
#' Pawlowsky-Glahn, V., Egozcue, J.J., and Tolosana-Delgado, R., 2015, Modeling
#' and analysis of compositional data: John Wiley and Sons, Ltd.
#'
#' @examples
#' \dontrun{
#' plotStdCompMeans( simplexStats, kappa, simplexModPar, elementOrder)
#' }
#'
#' @export
plotStdCompMeans <- function(simplexStats, kappa, simplexModPar, elementOrder,
                          intervalPercentage = 95, symbolSize = 0.75) {

  Internal1 <- function(compMeans, kappa, elementOrder,
                        center, metricVariance,
                        interval, pdf) {

    # reorder the compositional means and then standardize them
    tmp <- compMeans[, elementOrder]
    tmp <- CoDA::Perturb(center[elementOrder]^(-1), tmp, kappa = kappa)
    stdCompMeans <- CoDA::Power(tmp, 1/sqrt(metricVariance), kappa = kappa)

    a <- t( apply(stdCompMeans, 2, quantile, probs=c(interval, 0.50 ),
                  names = FALSE ) )

    # The code for "Elements" ensures that the chemical elements are
    # plotted in the order specified by elementOrder.
    b <- data.frame( qlower = a[, 1],
                     qupper = a[, 2],
                     qmid = a[, 3],
                     Elements = factor(elementOrder, levels = elementOrder),
                     Pdf = rep.int( as.character(pdf), nrow(a) ))
    return(b)

  }

  # barycenter for the standardized concentrations
  # Strictly, the barycenter is a vector in which all elements are equal. Here
  # variable barycenter is a scalar, which records the value of one element.
  barycenter <- kappa / length(simplexStats$sampleCenter)

  tailPercentage <- 0.5*(100.0-intervalPercentage)
  interval <- c(tailPercentage,100.0-tailPercentage)/100.0

  df1 <- Internal1(simplexModPar$compMean1, kappa, elementOrder,
                   simplexStats$sampleCenter, simplexStats$metricVariance, interval, 1 )
  df2 <- Internal1(simplexModPar$compMean2, kappa, elementOrder,
                   simplexStats$sampleCenter, simplexStats$metricVariance, interval, 2 )

  compData <- rbind(df1, df2)

  # The environment is set to the local environment, so that local variables
  # (namely, barycenter and symbolsize) can be found. Otherwise, ggplot
  # will search the global environment.
  w <- ggplot2::ggplot(compData,
                       ggplot2::aes(x = Elements, y = qmid,
                                    ymin = qlower, ymax = qupper,
                                    colour = Pdf),
                       environment = environment() ) +
    ggplot2::scale_colour_manual(values=c("blue","red")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=barycenter)) +
    ggplot2::geom_pointrange(size=symbolSize) +
    ggplot2::ylab("Standardized concentrations (no units)")

  print(w)

}

#' @title Plot compositional means
#'
#' @description Plot the compositional means for the pdfs in the
#' finite mixture model.
#'
#' @param simplexModPar
#' List containing Monte Carlo samples of the parameters in the finite mixture
#' model.
#' These parameters are expressed in terms of their equivalent values in the
#' simplex.
#' This list is return by function \code{\link{backTranform}}; the
#'                        documentation for function backTranform includes a
#'                        complete description of container
#'                        \code{simplexModPar}.
#' @param elementOrder Vector specifying the order in which the elements are
#' plotted.
#' @param symbolSize The size of the plotting symbol.
#'
#' @details
#' The plot does not include intervals from the distributions of the
#' mean components, because the intervals are smaller than the symbol size.
#'
#' @examples
#' \dontrun{
#' plotCompMeans( simplexModPar, elementOrder)
#' }
#'
#' @export
plotCompMeans <- function(simplexModPar, elementOrder, symbolSize = 2) {

  Internal1 <- function(compMeans, elementOrder, pdf) {

    a <- apply(compMeans[, elementOrder], 2, median )

    # The code for "Elements" ensures that the chemical elements are
    # plotted in the order specified by elementOrder.
    b <- data.frame( qmid = a,
                     Elements = factor(elementOrder, levels = elementOrder),
                     Pdf = rep.int( as.character(pdf), length(a) ))
    return(b)
  }

  df1 <- Internal1(simplexModPar$compMean1, elementOrder, 1 )
  df2 <- Internal1(simplexModPar$compMean2, elementOrder, 2 )

  compData <- rbind(df1, df2)

  # The environment is set to the local environment, so that local variables
  # (namely, symbolsize) can be found. Otherwise, ggplot
  # will search the global environment.
  w <- ggplot2::ggplot(compData,
                       ggplot2::aes(x = Elements, y = qmid, colour = Pdf),
                       environment = environment() ) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values=c("blue","red")) +
    ggplot2::geom_point(size=symbolSize) +
    ggplot2::ylab("Concentration (mg/kg)")

  print(w)

}


#' @title Plot both variation matrices, scaled by the square root
#'
#' @description Plot both variation matrices, scaled by the square root.
#' Because a variation matrix is symmetric and its diagonal is zero, The
#' two variation matrices from the two pdfs are combined and displayed
#' in a single plot.
#'
#' @param simplexModPar
#' List containing Monte Carlo samples of the parameters in the finite mixture
#' model.
#' These parameters are expressed in terms of their equivalent values in the
#' simplex.
#' This list is return by function \code{\link{backTranform}}; the
#'                        documentation for function backTranform includes a
#'                        complete description of container
#'                        \code{simplexModPar}.
#' @param elementOrder Vector specifying the order in which the elements are
#' plotted.
#' @param colorScale    Character string specifying the color scale for
#' the plot. The choices are either "spectrum" (default) and "rainbow."
#'
#' @details
#' In the plot, the upper left side represents one triangle from the
#' variation matrix for pdf 1, and the lower right side represents one triangle
#' from the variation matrix for pdf 2.
#'
#' The pixels in the plot represent scaled variances of the log-ratios
#' between the respective chemical elements. The scaling is desirable because
#' it reduces the large range of the variances, making it easier to
#' visualize all of the variances together. The scaling function is the
#' square root; so, the pixels in the plot strickly represent standard devations
#' of the log-ratios between the respective chemical elements.
#'
#' @references
#' Pawlowsky-Glahn, V., Egozcue, J.J., and Tolosana-Delgado, R., 2015, Modeling
#' and analysis of compositional data: John Wiley and Sons, Ltd.
#'
#' @examples
#' \dontrun{
#' plotSqrtVarMatrix( simplexModPar, elementOrder, colorScale = "rainbow" )
#' }
#'
#' @export
plotSqrtVarMatrices <- function( simplexModPar, elementOrder,
                                colorScale = "spectrum" ) {

  D <- dim(simplexModPar$varMatrix1)[3]  # D is the standard notation, and is concise.


  medVarMatrix1 <- apply(simplexModPar$varMatrix1, c(2,3), median)
  medVarMatrix2 <- apply(simplexModPar$varMatrix2, c(2,3), median)

  # reorder
  medVarMatrix1 <- medVarMatrix1[elementOrder, elementOrder]
  medVarMatrix2 <- medVarMatrix2[elementOrder, elementOrder]

  Z <- matrix( NA_real_, nrow=D, ncol=D, dimnames=list(elementOrder,elementOrder) )
  Z[upper.tri(Z)] <- medVarMatrix1[upper.tri(medVarMatrix1)]
  Z[lower.tri(Z)] <- medVarMatrix2[lower.tri(medVarMatrix2)]

  Z <- sqrt(Z)

  Z <- reshape2::melt(Z)
  Z <- na.omit(Z)

  w <- ggplot2::ggplot(Z, ggplot2::aes(Var2, Var1)) +
    ggplot2::geom_tile(data=Z, ggplot2::aes(fill=value), color="white") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,colour = "black")) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(colour = "black")) +
    ggplot2::theme(panel.background = ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank()) +
    ggplot2::xlab("Element") +
    ggplot2::ylab("Element") +
    ggplot2::coord_equal()

  if(colorScale == "spectrum") {
    myPalette <- colorRampPalette( c("blue", "green", "yellow", "orange", "red"),
                                   space="rgb", interpolate="linear" )

    w <- w + ggplot2::scale_fill_gradientn(limit=range(Z$value),
                                           colours = myPalette(10),
                                           name="Std dev")
  } else if(colorScale == "rainbow") {
    w <- w + ggplot2::scale_fill_gradientn(limit=range(Z$value),
                                           colours = rev(colorspace::rainbow_hcl(7)),
                                           name="Std dev")
  } else {
    stop("The color scale is not recognized.")
  }

  print(w)

}


#' @title Plot the field samples as clusters
#'
#' @description Plot the locations of the field samples on a previously-plotted
#' map. The attributes of each location symbol (for example, color) indicate
#' the conditional probability that the field sample is associated with
#' the pdfs in the finite mixture model. That is, the attributes indicate the
#' cluster to which the field sample belongs.
#'
#' @param concData   SpatialPointsDataFrame storing the geochemical
#'                 concentrations and locations of each field sample
#' @param simplexModPar
#' List containing Monte Carlo samples of the parameters in the finite mixture
#' model.
#' These parameters are expressed in terms of their equivalent values in the
#' simplex.
#' This list is return by function \code{\link{backTranform}}; the
#'                        documentation for function backTranform includes a
#'                        complete description of container
#'                        \code{simplexModPar}.
#' @param probIntervals Vector containing intervals of conditional
#' probability. All field samples within an given interval are plotted the same
#' way.
#' @param symbolIndices Vector containing the indices of the plotting symbols
#' for the conditional probability intervals.
#' @param symbolSizes Vector containing the relative sizes of the plotting symbols
#' for the conditional probability intervals.
#' @param symbolColors Vector containing the colors of the plotting symbols
#' for the conditional probability intervals.
#'
#' @details
#' A field sample is assigned plotting attributes based upon the
#' median of the conditional probabilities, which are
#' calculated during the Monte Carlo sampling of the posterior pdf. These
#' conditional probabilities, and hence their median, indicate the extend
#' to which the field sample is associated with the first pdf in the finite
#' mixture model.
#'
#' The plotting attributes are specified by arguments \code{probIntervals},
#' \code{symbolIndices}, \code{symbolSizes}, and \code{symbolColors}. To
#' understanding the specification of these attributes, consider their default
#' values, which pertain to three probability intervals.
#'
#' Vector \code{probIntervals} has elements 0, 0.1, 0.9, and 1. These four
#' elements specify three probability intervals: [0,0.1], [0.1,0.9], and
#' [0.9,1]. Notice that the first and last elements of \code{probIntervals}
#' are 0 and 1 respectively. The first interval, [0,0.1] is interpreted to
#' indicate the
#' extent to which the field sample is associated with the second pdf
#' in the finite mixture model.
#' The third interval, [0.9,1] is interpreted to indicate the
#' extent to which the field sample is associated with the first pdf
#' in the finite mixture model. The second (middle) interval, [0.1,0.9]
#' is interpreted to indicate the extent to which the field sample is
#' associated with both pdfs in the finite mixture model.
#'
#' Because, in this explanation, vector \code{probIntervals} specifies three
#' probability intervals, arguments for vectors \code{symbolIndices},
#' \code{symbolSizes}, and \code{symbolColors} must have three elements.
#' The symbol indices, symbol sizes, and symbol colors are described in
#' Murrell (2006, p. 55-56, 68, 69).
#'
#' This function adds symbols to a map that has already been plotted. Obviously,
#' the user must plot that map before calling this funcion.
#'
#' @references
#' Murrell, P., 2006, R graphics: Chapman & Hall / CRC.
#'
#' @examples
#' \dontrun{
#' map('state', fill = TRUE, col = "gray60", border = "white")
#' plotSqrtVarMatrix( simplexModPar, elementOrder, colorScale = "rainbow" )
#' }
#'
#' @export
plotClusters <- function(concData, simplexModPar,
                        probIntervals = c( 0, 0.1, 0.9, 1.0 ),
                        symbolIndices = c( 16, 16, 16 ),
                        symbolSizes = c( 1/3, 1/2, 1/3 ),
                        symbolColors = c( "red", "green", "blue" ) ) {

  locations <- coordinates( concData )

  # median of conditional probabilites
  g_median <- apply(simplexModPar$g, 2, median)

  for (i in 1:(length(probIntervals)-1)) {

    areInInterval <- probIntervals[i] <= g_median & g_median <= probIntervals[i+1]

    if( sum(areInInterval) == 0 ) next

    locations.sp <- sp::SpatialPoints( locations[areInInterval,],
                                       proj4string=CRS( "+proj=longlat +ellps=WGS84" ) )
    sp::plot( locations.sp, add=TRUE, pch=symbolIndices[i],
              col=symbolColors[i], cex=symbolSizes[i] )

  }

}

#' @title Split the geochemical data
#'
#' @description The geochemical data, which have been clustered, are split
#' into two groups according to the clusters.
#'
#' @param concentrationData
#' SpatialPointsDataFrame storing the geochemical
#' concentrations and locations of each field sample. These data have been
#' clustered.
#'
#' @param censorIndicators
#' Matrix specifying which samples in \code{concentrationData} are
#' censored and, if so, have imputed values. Each column pertains to an
#' element,
#' and each row pertains to a field sample. Matrix elements may be "no",
#' "left", or "right", indicating that the concentration is "not censored",
#'  is "left-censored", or is "right-censored".
#'
#' @param simplexModPar
#' List containing Monte Carlo samples of the parameters in the finite mixture
#' model. These parameters are expressed in terms of their equivalent values in the
#' simplex. This list is return by function \code{\link{backTranform}}; the
#' documentation for function backTranform includes a
#' complete description of container \code{simplexModPar}.
#'
#' @param threshold
#' The threshold used to split the data into two groups. (See details.)
#'
#' @details
#' For each
#' field sample, the median of the Monte Carlo samples of conditional probability
#' is calculated. If this median is between 1-threshold and 1, then the
#' field sample is associated with pdf 1 in the finite mixture model. However,
#' if this median is between 0 and threshold, then the field samples is
#' associated with pdf 2 in the finite mixture model. This criterion is
#' used to split the field samples into two groups.
#'
#' The variable threshold must be less than 0.5 and greater than 0.
#'
#' @return A list with the following components is returned.
#' @return \item{concentrationData1}{SpatialPointsDataFrame storing the
#' geochemical concentrations and locations of each field sample
#' associated with pdf 1.}
#' @return \item{censorIndicators1}{Matrix specifying which
#' samples in \code{concentrationData1} are censored and, if so, have
#' imputed values. See the description for argument \code{censorIndicators}.}
#' @return \item{concentrationData2}{Data structure that is equivalent
#' to concentrationsData1, except that the field samples are associated with
#' pdf 2.}
#' @return \item{concentrationData2}{Data structure that is equivalent
#' to censorIndicators1, except that the field samples are associated with
#' pdf 2.}
#'
#' @examples
#' \dontrun{
#' tmp <- splitData(concentrationData, censorIndicators,
#'                  simplexModPar, threshold = 0.10 )
#' }
#'
#' @export
splitData <- function(concentrationData, censorIndicators,
                      simplexModPar, threshold = 0.10 ) {

  if(threshold <= 0.0 || 0.5 <= threshold)
    stop("Argument threshold must be > 0 and < 0.50.")

  g_median <- apply(simplexModPar$g, 2, median)

  areInPdf1 <- 1.0-threshold <= g_median & g_median <= 1.0
  areInPdf2 <- 0.0 <= g_median & g_median <= threshold

  return(list(
    concentrationData1 = concentrationData[areInPdf1, ],
    censorIndicators1 = censorIndicators[areInPdf1, ],
    concentrationData2 = concentrationData[areInPdf2, ],
    censorIndicators2 = censorIndicators[areInPdf2, ]))
}


# #' @title Plot autocorrelation functions
# #'
# #' @description Plot autocorrelation functions for selected traces
# #' from each chain.
# #'
# #' @param fmmSamples
# #' List containing samples of the posterior pdf and
# #' related information. This list is return by function
# #' \code{\link{sampleFmm}}; the documentation for function sampleFmm includes
# #' a complete description of container \code{fmmSamples}.
# #'
# #' @details
# #' Autocorrelation functions are plotted for the following model parameters:
# #' \itemize{
# #'  \item Element [1] of the mean vector for pdf 1 (mu1[1]).
# #'  \item Element [1,1] of the covariance matrix for pdf 1 (Sigma1[1,1]).
# #'  \item Model proportion associated with pdf 1 (theta)
# #' }
# #'
# #' @examples
# #' \dontrun{
# #' plotAcfs(fmmSamples)
# #' }
# #'
# #' @export
# plotAcfs <- function(fmmSamples) {
#
#   selectedTraces <- fmmSamples$stanSelectedTraces
#
#   N <- length(selectedTraces)
#
#   origPar <- par( mfrow=c(3,1))
#   on.exit(par(origPar), add=TRUE)
#
#   for(i in 1:N){
#     devAskNewPage( ask = TRUE )
#
#     acf(selectedTraces[[i]][, "mu1[1]"], main="mu1[1]")
#
#     acf(selectedTraces[[i]][, "Sigma1[1,1]"], main="Sigma1[1,1]")
#
#     acf(selectedTraces[[i]][, "theta"], main="theta")
#
#     title(main=paste("Chain", i), outer=TRUE, line=-1)
#
#     devAskNewPage( ask = options()$device.ask.default )
#
#   }
# }

# #' @title Check combined chains
# #'
# #' @description Check that the chains have been properly combined.
# #'
# #' @param combinedChains List in which each element contains the
# #' Monte Carlo samples of a parameter from the finite mixture model.
# #' This list is return by function \code{\link{combineChains}}; the
# #'                        documentation for function combineChains includes a
# #'                        complete description of container
# #'                        \code{combinedChains}.
# #'
# #' @details
# #' Histograms are shown for five model parameters:
# #' (1) The first element of the mean vector for pdf 1, which is
# #' designated "mu1[1]."
# #' (2) The first element of the mean vector for pdf 2, which is
# #' designated "mu2[1]."
# #' (3) The first element of the covariance matrix for pdf 1, which is
# #' designated "Sigma1[1,1]."
# #' (4) The first element of the covariance matrix for pdf 2, which is
# #' designated "Sigma2[1,1]."
# #' (5) The proportion in the finite mixture model associated with pdf 1,
# #' which is designated "theta."
# #'
# #' If the histograms show a bimodal distribution, then it is likely that
# #' the chains were improperly combined.
# #'
# #' @examples
# #' \dontrun{
# #' checkCombinedChains(combinedChains)
# #' }
# #'
# #' @export
#
# checkCombinedChains <- function(combinedChains) {
#
#   devAskNewPage( ask = TRUE )
#
#   origPar <- par(mfrow=c(1,2))
#   hist(combinedChains$mu1[,1], freq=FALSE,
#        col = "gray", border = "white",
#        xlab = "mu1[1]", main = "Pdf 1")
#   hist(combinedChains$mu2[,1], freq=FALSE,
#        col = "gray", border = "white",
#        xlab = "mu2[1]", main = "Pdf 2")
#   par(origPar)
#
#   origPar <- par(mfrow=c(1,2))
#   hist(combinedChains$Sigma1[,1,1], freq=FALSE,
#        col = "gray", border = "white",
#        xlab = "Sigma1[1,1]", main = "Pdf1")
#   hist(combinedChains$Sigma2[,1,1], freq=FALSE,
#        col = "gray", border = "white",
#        xlab = "Sigma2[1,1]", main = "Pdf 2")
#   par(origPar)
#
#   hist(combinedChains$theta, freq=FALSE,
#        col = "gray", border = "white",
#        xlab = "theta", main = "")
#
#   devAskNewPage( ask = options()$device.ask.default )
#
# }
#
# # chains a matrix of chains
# # each column is a separate chain.
# #
# # the returned value is a list with two elements,
# # the potential scale reduction (rHat) and the number of
# # effective samples (nEff)
# #
# # At first glance it may seem that the calculation of these
# # two measures should be in separate functions. However,
# # the calculation of nEff requires intermediate results from
# # the calucation of rHat.
# CalcMeasures <- function ( chains ) {
#
#   ##
#   ## Estimate the potential scale reduction
#   ##
#
#   # The variables are documented on p. 284-285 of Gelman et al. (3rd ed.)
#   m <- ncol( chains )
#   n <- nrow( chains )
#   theColMeans <- colMeans( chains )         # psiBar,j
#   theGlobalMean <- mean( theColMeans )      # psiBar..
#
#   # between sequence variance
#   B <- n / ( m - 1 ) * sum( (theColMeans-theGlobalMean)^2 )
#
#   s2 <- vector( mode="numeric", length=m )
#   for( j in 1:m ) {
#     s2[j] <- var(chains[,j])
#   }
#
#   # within sequence variance
#   W <- mean( s2 )
#
#   # marginal posterior variance, eqn. 11.3
#   V <- ( n - 1 ) / n * W + 1 / n * B
#
#   # potential scale reduction, eqn. 11.4
#   rHat <- sqrt( V / W )
#
#   # edition 2 of Gelman et al., p. 298
#   # effective number of independent draws, eqn. 11.4
#   # nEff <- min( floor( m * n * V / B ), m * n )
#
#   ##
#   ## Estimate the effective sample size
#   ##
#
#   # See edition 3 of Gelman et al. p. 286-287
#
#   # unnumbered eqn just above eqn 11.7
#   Vt <- vector(mode = "numeric", length = n)
#   Vt[1] <- 0.0
#   for(t in 2:n) {
#     Vt[t] <- sum( colSums( diff(chains, lag = (t-1))^2 ) )  / (m*(n-(t-1)))
#   }
#
#   # eqn 11.7
#   rho_t <- 1 - Vt / ( 2 * V)
#
#   #   # This following method of calculating the mean correlation coefficient
#   #   # differs from that on page 286, but the results are almost identical.
#   #   a <- acf(chains, plot=FALSE)
#   #   nLags <- dim(a$acf)[1]
#   #   rho_t <- vector(mode="numeric", length=nLags)
#   #   for(i in 1:nLags) {
#   #     partialSum <- 0.0
#   #     for(j in 1:m) {
#   #       partialSum <- partialSum + a$acf[i,j,j]
#   #     }
#   #     rho_t[i] <- partialSum / m
#   #   }
#
#
#   # One stopping criterion is in Gelman et al., p. 286-287. It is defined by
#   # the condition:
#   #           if( rho_t[t] + rho_t[t+1] < 0.0 ) break
#   #
#   # Another, slightly different, stopping criterion is presented in the
#   # rstan manual in the subsection "Estimation of effective sample size",
#   # section number 54.4, in Stan Modeling Language - User's Guide and
#   # Reference Manual, Stan version 2.6.0. It is defined by
#   # the condition:
#   #           if( rho_t[t] < 0.0 ) break
#   #
#   # The index starts at 2, which corresponds to a lag of 1.
#   sum_rho_t <- 0.0
#   for(t in 2:n) {
#     # if( rho_t[t] + rho_t[t+1] < 0.0 ) break
#     if( rho_t[t] < 0.0 ) break
#     sum_rho_t <- sum_rho_t + rho_t[t]
#   }
#
#   # eqn. 11.8, p. 287
#   nEff <- as.integer( m * n / ( 1 + 2 * sum_rho_t ) )
#
#   return( c(rHat=rHat, nEff=nEff ) )
# }
#
# # split the chains are described by Gelman, p. 284
# #
# # combinedChains is a vector: the first "nSamplesPerChain" samples are
# # from the first chain, the second "nSamplesPerChain" samples are from the
# # second chain, and so on
# SplitChains <- function( combinedChains, nSamplesPerChain ) {
#   X <- matrix(combinedChains, nrow = nSamplesPerChain)
#
#   nNewRows <- as.integer( nSamplesPerChain / 2 )
#   a <- X[1:nNewRows, ]
#   b <- X[(1:nNewRows)+nNewRows, ]
#   return(cbind(a, b))
# }
#
#
# #' @title Calculate convergence measures of the Monte Carlo chains
# #'
# #' @description Calculate the potential scale reduction and the effective
# #' number of simulations draws for each model parameter in the finite
# #' mixture model.
# #'
# #' @param combinedChains     List in which each element contains the
# #' Monte Carlo samples of a parameter from the finite mixture model.
# #' This list is return by function \code{\link{combineChains}}; the
# #'                        documentation for function combineChains includes a
# #'                        complete description of container
# #'                        \code{combinedChains}.
# #' @param fmmSamples       List, which includes the Monte Carlo chains and
# #'                        other data, returned by function
# #'                        \code{\link{sampleFmm}}; the
# #'                        documentation for fmmSample includes a complete
# #'                        description of fmmSamples.
# #'
# #' @details
# #' ???
# #'
# #' @return A list containing these 7 elements.
# #' \tabular{lrr}{
# #' Element name \tab Role in finite mixture model \tab Structure of rHat and nEff \cr
# #' theta \tab Proportion associated with \tab scalar \cr
# #' mu1 \tab Mean vector for pdf 1 \tab vector of length N \cr
# #' mu2 \tab Mean vector for pdf 2 \tab vector of length N \cr
# #' tau1 \tab Vector of standard deviations for pdf 1 \tab vector of length N \cr
# #' tau2 \tab Vector of standard deviations for pdf 2 \tab vector of length N \cr
# #' Corr1 \tab Matrix of correlations for pdf 1 \tab Matrix of size NxN. \cr
# #' Corr2 \tab Matrix of correlations for pdf 2 \tab Matrix of size NxN.
# #' }
# #' Each of the seven elements comprise another
# #' list with two elements. The first element pertains to
# #' the potential scale reduction (rHat), and the second element pertains to the
# #' number of effective samples (nEff). The structures of these two elements
# #' are identical and depend
# #' on the variable to which they pertain. The dimension of the two pdfs
# #' in the finite mixture model is designated as N.
# #'
# #' Potential scale reduction and the number of effective samples are
# #' described by Gelman et al. (2014, p. 284-287).
# #'
# #' @references
# #' Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A.,
# #' and Rubin, D.B., 2014, Bayesian data analysis (3rd ed.),
# #' CRC Press.
# #'
# #' @examples
# #' \dontrun{
# #' convergenceMeasures <- calcConvMeasures(combinedChains, fmmSamples)
# #' }
# #'
# #' @export
# calcConvMeasures <- function(combinedChains, fmmSamples) {
#
#   Internal1 <- function(x, N, nSamplesPerChain) {
#     rHat <- vector(mode = "numeric", length = N)
#     nEff <- vector(mode = "numeric", length = N)
#     for(i in 1:N) {
#       tmp <- CalcMeasures( SplitChains(x[, i], nSamplesPerChain) )
#       rHat[i] <- tmp["rHat"]
#       nEff[i] <- tmp["nEff"]
#     }
#     return(list(rHat = rHat, nEff = nEff))
#   }
#
#   Internal2 <- function(X, N, nSamplesPerChain) {
#     rHat <- matrix(NA_real_, nrow = N, ncol = N)
#     nEff <- matrix(NA_real_, nrow = N, ncol = N)
#     for(i in 1:(N-1)) {
#       for(j in (i+1):N) {
#         tmp <- CalcMeasures( SplitChains(X[, i, j], nSamplesPerChain) )
#         rHat[i,j] <- tmp["rHat"]
#         nEff[i,j] <- tmp["nEff"]
#         rHat[j,i] <- rHat[i,j]
#         nEff[j,i] <- nEff[i,j]
#       }
#     }
#     return(list(rHat = rHat, nEff = nEff))
#   }
#
#   Internal3 <- function(X) {
#     Y <- array( NA_real_, dim = dim(X) )
#     nSamples <- dim(X)[1]
#     for(i in 1:nSamples) {
#       Y[i, , ] <- cov2cor(X[i, ,])
#     }
#     return(Y)
#   }
#
#   Internal4 <- function(X) {
#     tmp <- dim(X)
#     Y <- matrix(NA_real_, nrow = tmp[1], ncol = tmp[2])
#     for(i in 1:tmp[2]) {
#       Y[, i] <- sqrt(X[, i, i])
#     }
#     return(Y)
#   }
#
#   nSamplesPerChain <- fmmSamples$samplingParams$nSamplesPerChain
#   N <- ncol(combinedChains$mu1)  # N is also the number of principal components
#
#   measures <- vector(mode = "list")
#
#   tmp <- CalcMeasures( SplitChains(combinedChains$theta, nSamplesPerChain) )
#   measures$theta <- list(rHat = unname(tmp["rHat"]), nEff = unname(tmp["nEff"]))
#
#   measures$mu1 <- Internal1(combinedChains$mu1, N, nSamplesPerChain)
#   measures$mu2 <- Internal1(combinedChains$mu2, N, nSamplesPerChain)
#
#   measures$tau1 <- Internal1(Internal4(combinedChains$Sigma1), N, nSamplesPerChain)
#   measures$tau2 <- Internal1(Internal4(combinedChains$Sigma2), N, nSamplesPerChain)
#   measures$Corr1 <- Internal2(Internal3(combinedChains$Sigma1), N, nSamplesPerChain)
#   measures$Corr2 <- Internal2(Internal3(combinedChains$Sigma2), N, nSamplesPerChain)
#
#   return( measures )
#
# }
#
# #' @title Plot convergence measures
# #'
# #' @description Plot the potential scale reduction and the number of
# #' effective samples for each parameter in the finite mixture model.
# #'
# #' @param convMeasures List containing the potential scale reductions and
# #' numbers of effective samples for every model parameter. This list is
# #' returned by function \code{\link{calcConvMeasures}}; the documentation
# #' for function calcConvMeasures includes a complete description of
# #' container \code{convMeasures}.
# #'
# #' @param combinedChains List in which each element contains the
# #' Monte Carlo samples of a parameter from the finite mixture model.
# #' This list is return by function \code{\link{combineChains}}; the
# #'                        documentation for function combineChains includes a
# #'                        complete description of container
# #'                        \code{combinedChains}.
# #'
# #' @details
# #' The potential scale reductions and
# #' numbers of effective samples are displayed as bar charts and matrices,
# #' depending upon the structure of the associated model parameters.
# #'
# #' When a bar chart shows potential scale reductions, there are two, horizontal
# #' green lines, which approximately delimit an optimal range. When a bar chart
# #' shows numbers of effective samples, there is one, horizontal line, which
# #' is the maximum (optimal) number of effective samples.
# #'
# #' @examples
# #' \dontrun{
# #' plotConvMeasures(convMeasures, combinedChains)
# #' }
# #'
# #' @export
# plotConvMeasures <- function(convMeasures, combinedChains) {
#
#   # For vectors
#   InternalPlot1 <- function(rHat, nEff, nSamples, paramName) {
#
#     df <- data.frame( component = 1:length(rHat),
#                       rHat = rHat,
#                       nEff = nEff )
#     a <- ggplot2::ggplot( df,
#                           ggplot2::aes(x = component, y = rHat),
#                           environment = environment()) +
#       ggplot2::geom_bar(stat = "identity") +
#       ggplot2::geom_hline(ggplot2::aes(yintercept = 0.99), col="green") +
#       ggplot2::geom_hline(ggplot2::aes(yintercept = 1.05), col="green") +
#       ggplot2::xlab("Component") +
#       ggplot2::ylab("Potential scale reduction (no units)") +
#       ggplot2::ggtitle(paramName)
#
#     b <- ggplot2::ggplot( df,
#                           ggplot2::aes(x = component, y = nEff),
#                           environment = environment()) +
#       ggplot2::geom_bar(stat = "identity") +
#       ggplot2::geom_hline(ggplot2::aes(yintercept = nSamples), col="green") +
#       ggplot2::xlab("Component") +
#       ggplot2::ylab("Number of effective samples (no units)")+
#       ggplot2::ggtitle(paramName)
#
#     grid::grid.newpage()
#     grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,1)))
#     print(a, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
#     print(b, vp=grid::viewport(layout.pos.row=2, layout.pos.col=1))
#
#   }
#
#   # for matrices
#   InternalPlot2 <- function(rHat, nEff, nSamples, paramName){
#
#     Z1 <- reshape2::melt(rHat)
#     Z1 <- na.omit(Z1)
#     p1 <- ggplot2::ggplot(Z1,
#                           ggplot2::aes(Var2, Var1),
#                           environment = environment())+
#       ggplot2::geom_tile(data=Z1, ggplot2::aes(fill=value), color="white")+
#       ggplot2::scale_fill_gradient2(low="blue", high="red", mid="white",
#                                     midpoint=1,
#                                     name="Potential\nscale\nreduction\n(no units)")+
#       ggplot2::xlab("Component") +
#       ggplot2::ylab("Component") +
#       ggplot2::coord_equal() +
#       ggplot2::ggtitle(paramName)
#
#     Z2 <- reshape2::melt(nEff)
#     Z2 <- na.omit(Z2)
#
#     p2 <- ggplot2::ggplot(Z2,
#                           ggplot2::aes(Var2, Var1),
#                           environment = environment())+
#       ggplot2::geom_tile(data=Z2, ggplot2::aes(fill=value), color="white")+
#       ggplot2::scale_fill_gradient(limit=c(min(Z2$value),nSamples),
#                                    high = "#132B43", low = "#56B1F7",
#                                    name="Number\nof effective\nsamples\n(no units)")+
#       ggplot2::xlab("Component") +
#       ggplot2::ylab("Component") +
#       ggplot2::coord_equal()
#
#     grid::grid.newpage()
#     grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,2)))
#     print(p1, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
#     print(p2, vp=grid::viewport(layout.pos.row=1, layout.pos.col=2))
#
#   }
#
#
#   nSamples <- length(combinedChains$theta)
#
#   devAskNewPage( ask = TRUE )
#
#   InternalPlot1(convMeasures$theta$rHat, convMeasures$theta$nEff, nSamples, "theta")
#   InternalPlot1(convMeasures$mu1$rHat, convMeasures$mu1$nEff, nSamples, "mu1")
#   InternalPlot1(convMeasures$mu2$rHat, convMeasures$mu2$nEff, nSamples, "mu2")
#   InternalPlot1(convMeasures$tau1$rHat, convMeasures$tau1$nEff, nSamples, "tau1")
#   InternalPlot1(convMeasures$tau2$rHat, convMeasures$tau2$nEff, nSamples, "tau2")
#
#   InternalPlot2(convMeasures$Corr1$rHat, convMeasures$Corr1$nEff, nSamples,
#                 "Correlation matrix 1")
#   InternalPlot2(convMeasures$Corr2$rHat, convMeasures$Corr2$nEff, nSamples,
#                 "Correlation matrix 2")
#
#   devAskNewPage( ask = options()$device.ask.default )
#
# }
