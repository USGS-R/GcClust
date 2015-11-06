
#' @title Calculate sample statistics pertinent to the simplex
#'
#' @description Calculate the sample center, total variation matrix, and
#' the metric variance.
#'
#' @param gcData
#' List containing the geochemical and related data. This container is
#' described in the package documentation.
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
#' simplexStats <- calcSimplexStats(gcData)
#' }
#'
#' @export
calcSimplexStats <- function(gcData) {

  concData <- as.matrix(gcData$concData@data)

  sampleCenter <- CalcCompCenter( concData, kappa = gcData$constSumValue )
  variationMatrix <- CalcVariationMatrix( concData )
  metricVariance <- CalcTotalVariance( variationMatrix )

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
#' @param gcData
#' List containing the geochemical and related data. This container is
#' described in the package documentation.
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

  X <- as.matrix(gcData$concData@data)

  Psi <- CalcPsiMatrix2( ncol(X) )
  ilrCoefs <- CalcIlrCoefs( X, t(Psi) )

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
#' This list is return by function \code{\link{transformGcData}}, for which the
#' documentation includes a complete description of container \code{transData}.
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
#' plotEdaVar(transData)
#' }
#'
#' @export
plotEdaVar <- function(transData, relOffset = 0.04, size = 3) {

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

#' @title Plot the principal component distributions
#'
#' @description Plot the distributions of the principal components as
#' boxplots and violin plots. At the top of the violin plots, print
#' the calculated standard deviations of principal components.
#'
#' @param transData     List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}, for which the
#' documentation includes a complete description of container \code{transData}.
#' @param relOffset  Scalar specifying the relative distance that the
#' standard deviations are are offset from the top of the highest violin figure.
#' @param size  Scalar specifying the text size for the standard deviations.
#'
#' @examples
#' \dontrun{
#' plotEdaDist(transData)
#' }
#'
#' @export
plotEdaDist <- function(transData, relOffset = 0.04, size = 3) {

  robustPCs <- transData$robustPCs
  colnames(robustPCs) <- 1:ncol(robustPCs)

  df <- reshape2::melt(robustPCs)[,-1]   # the first col contains the row numbers
  colnames(df) <- c("Component", "Value")

  df.sd <- data.frame(Component = colnames(robustPCs),
                      sd = round( apply(robustPCs, 2, sd), digits = 2) )

  p1 <- ggplot2::ggplot( df,
                         ggplot2::aes(x = factor(Component), y = Value),
                         environment = environment() ) +
    ggplot2::geom_boxplot() +
    ggplot2::xlab("Principal component")

  p2 <- ggplot2::ggplot( df,
                         ggplot2::aes(x = factor(Component), y = Value),
                         environment = environment() ) +
    ggplot2::geom_violin(scale = "width", fill = "grey50") +
    ggplot2::xlab("Principal component")
#     ggplot2::geom_text(ggplot2::aes(x = factor(Component),
#                                     y = max(df$Value) * (1 + relOffset),
#                                     label = sd), data = df.sd,
#                        size = size, vjust = 0, hjust = 1, angle = 90 )



  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,1)))
  print(p1, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
  print(p2, vp=grid::viewport(layout.pos.row=2, layout.pos.col=1))

}

#' @title Plot the correlations among the principal components
#'
#' @description Plot the correlation matrix of the principal components and
#' the plot a histogram of the correlations.
#'
#' @param transData     List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}, for which the
#' documentation includes a complete description of container \code{transData}.
#'
#' @examples
#' \dontrun{
#' plotEdaCorr(transData)
#' }
#'
#' @export
plotEdaCorr <- function(transData) {

  corMatrix <- cor(transData$robustPCs)

  p1 <- ggplot2::ggplot(reshape2::melt(corMatrix),
                        ggplot2::aes(Var2, Var1),
                        environment = environment())+
    ggplot2::geom_tile(ggplot2::aes(fill=value), color="white")+
    ggplot2::scale_fill_gradient2(low="blue", high="red", mid="white",
                                  midpoint=0, limit=c(-1,1),
                                  name="Correlation\n(Pearson)")+
    ggplot2::xlab("Component") +
    ggplot2::ylab("Component") +
    ggplot2::coord_equal()

  df <- data.frame(x = corMatrix[lower.tri(corMatrix)])
  p2 <- ggplot2::ggplot( df,
                         ggplot2::aes(x = x),
                         environment = environment()) +
    ggplot2::geom_histogram() +
    ggplot2::xlab("Correlation (Pearson)") +
    ggplot2::ylab("Count")

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,1)))
  print(p1, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
  print(p2, vp=grid::viewport(layout.pos.row=2, layout.pos.col=1))

}




#' @title Sample the posterior pdf
#'
#' @description Sample the posterior probability density function for the
#' finite mixture model. The model implementation and the sampling are
#' performed with package rstan.
#'
#' @param transData
#' List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}, for which the
#' documentation includes a complete description of container \code{transData}.
#' @param nPCS          Number of principal components that are used in the
#'                      finite mixture model (See Details).
#' @param tauBounds     Vector of length 2, containing the lower and upper
#'                      bounds for tau (See Details).
#' @param nWuSamples    Number of warm-up samples in each chain (See details).
#' @param nPwuSamples   Number of post warm-up samples in each chain (See details).
#' @param nChainsPerCore Number of chains that each core computes (See Details).
#' @param nCpuCores     Number of central processing units (cpu's) that are
#'                      used (See Details).
#' @param procDir      Directory into which the results are written
#'                      (See Details).
#'
#' @details
#' The parameters in the finite mixture model are estimated
#' from the robust principal components, which are stored as a matrix within
#' \code{transData}.
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
#' The number of samples per chain \code{nSamplesPerChain} equals the
#' number of warm-up samples plus the number of post warm-up samples. Both
#' groups have equal size.
#'
#' The sampling is performed in parallel on multiple central
#' processing units (cpu's). Of course, the number of requested cpu's
#' \code{nCpuCores} must be less than or equal to the actual number. The
#' parallel computations are organized in the following manner:
#' \code{nCpuCores} cpu's
#' are requested from the operating system; each cpu computes
#' \code{nChainsPerCore} individual chains, and each chain is written
#' to its own file in directory \code{procDir}.
#'
#' @return The returned values may be conveniently divided into two groups.
#' First, the samples in a single chain are stored in \code{stanfit} object,
#' which is described in the rstan documentation. The samples are of
#' the following model parameters:
#' \describe{
#'  \item{theta}{Proportion of population associated with pdf 1.}
#'  \item{mu1}{Mean vector for pdf 1.}
#'  \item{mu2}{Mean vector for pdf 2.}
#'  \item{tau1}{Standard deviation vector for pdf 1.}
#'  \item{tau2}{Standard deviation vector for pdf 2.}
#'  \item{L_Omega1}{Lower triangle of the correlation matrix for pdf 1.}
#'  \item{L_Omega2}{Lower triangle of the correlation matrix for pdf 2.}
#'  \item{log_lik}{The logarithm of the likelihood function.}
#' }
#' The \code{stanfit} object is written to a file in directory \code{procDir}.
#' The file names have the form "RawSamples?-?.dat" for which the first
#' ? is replaced by the repetition number and the second ? is replaced by the
#' cpu number. The important point is that the file name is unique.
#'
#' Second, a list with 4 elements is returned by the function:
#' @return \item{nChains}{Number of chains.}
#' @return \item{nWuSamples}{Number of warm-up samples.}
#' @return \item{nPwuSamples}{Number of post warm-up samples.}
#' @return \item{nSamplesPerChain}{See argument \code{nSamplesPerChain}.}
#' @return \item{fileNames}{Vector containing the names of the files
#' with the \code{stanfit} objects.}
#'
#' The total number of samples per chain equals \code{nWuSamples} plus
#' \code{nPwuSamples}. Only the post warm-up samples are used for interference.
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
#' fmmSamples <- sampleFmm(transData, nPCs,
#'                         nCpuCores = 7, nChainsPerCore = 5,
#'                         procDir = "Process121")
#' }
#'
#' @export
sampleFmm <- function(transData, nPCs,
                      nWuSamples = 500,
                      nPwuSamples = 500,
                      nChainsPerCore = 2,
                      nCpuCores = parallel::detectCores(),
                      procDir = ".") {


  rstanParallelSampler <- function(stanData, nWuSamples, nPwuSamples,
                                   nChainsPerCore, nCpuCores, procDir ) {

    CL <- parallel::makeCluster(nCpuCores)

    parallel::clusterExport(cl = CL,
                            c("stanData", "nWuSamples", "nPwuSamples",
                              "nChainsPerCore", "procDir"),
                            envir=environment())

    fnlist <- parallel::parLapply(CL, 1:nCpuCores, fun = function(cid) {

      # Make rstan available to the processors. This function won't work
      # otherwise. So, I'm violating the principles in "R packages",
      # p. 34, 82-84
      require(rstan, quietly = TRUE)

      fileNames <- vector(mode = "character", length = nChainsPerCore)

      for(i in 1:nChainsPerCore) {

        rng_seed <- sample.int(.Machine$integer.max,1)

        gen_inits <- function() {
          areInGrp1 <- sample(c(TRUE,FALSE), size = stanData$N,
                              prob = c(0.3, 0.7), replace = TRUE)
          return(list(
            theta = runif(1, min = 0.35, max = 0.65),
            mu1 = apply(stanData$Z[areInGrp1,], 2, mean ),
            mu2 = apply(stanData$Z[!areInGrp1,], 2, mean ),
            tau1 = apply(stanData$Z[areInGrp1,], 2, sd ),
            tau2 = apply(stanData$Z[!areInGrp1,], 2, sd ),
            L_Omega1 = diag(stanData$M),
            L_Omega2 = diag(stanData$M)
          ))
        }


        # Sampling from finite mixture model can be difficult --- the chains
        # tend to get stuck. I have read two different suggestions in the
        # emails for the rstan users group.
        #
        # 1. Suggestions from Michael Betacourt "This is ... a failure of
        # convergence. The problem is the extreme curvature near zero which
        # causes the gradients, and hence the HMC transitions, to fail."
        # "An ... option is to decrease the step size (or increase the
        # target average acceptance probability in the adaptation).
        # Alternatively you can try to find a parameterization that reduces
        # the change in curvature.  Lastly  you can reduce the curvature by
        # smoothing out your model ..."
        # "... just set adapt delta=0.9 (or 0.95, 0.99, increasing until
        # n_divergence is always zero).  The other defaults will change
        # (or not) as necessary automatically.  In general increasing the
        # target average acceptance probability may lead to smaller
        # efficiency but it will always lead to a more robust sampler
        # so it’s a very safe knob to turn."
        #
        # The user found that this suggestion did not work well.
        #
        # I too tried
        # this suggestion and also switched to a smooth prior pdf for tau1 and
        # tau2. I found that the number of stuck chains was signficantly
        # reduced.I don't know which of the two actions had this effect.
        # However, one chain did get stuck. In some
        # test of short chains, I found that the initial values of the chains
        # (for a particular parameter) varied a lot; some where quite large or
        # small. Is init = "0" not working? As adapt_delta increases,
        # the sampling time increases. If adapt_delta is 0.95 or larger, the
        # sampling time can be relatively long.
        #
        # https://groups.google.com/forum/?fromgroups#!searchin/stan-users/stuck$20chains/stan-users/xEG18UzoCGo/EeXeqPox3sAJ
        #
        # 2. Suggestion from another consultant: Reduce parameter stepsize
        # to a small value, say 0.0001. Betacourt did not comment on this
        # suggestion, so I suspect that he doesn't favor it. In some tests,
        # I found that the initial values of the chains (for a particular
        # parameter) were identical. This behavior is expected because
        # init = "0". Nevertheless, this suggestion solves the problem.

        rawSamples <- rstan::sampling(GcClust:::sm, data = stanData,
                                      init = gen_inits,
                                      # control = list(stepsize = 0.00001),
                                      control = list(stepsize = 0.0001),
                                      # control = list(adapt_delta = 0.95),
                                      chains = 1,
                                      iter = nWuSamples + nPwuSamples,
                                      warmup = nWuSamples,
                                      seed = rng_seed, chain_id = cid,
                                      pars=c("theta", "mu1", "mu2",
                                             "tau1", "tau2",
                                             "L_Omega1", "L_Omega2", "log_lik"),
                                      save_dso = FALSE)

        fileNames[i] <- paste("RawSamples", cid, "-", i, ".dat", sep = "")
        save( rawSamples, file = paste(procDir, "\\", fileNames[i], sep = "") )

      }
      return(fileNames)
    } )

    parallel::stopCluster(CL)
    return(unlist(fnlist))
  }


  if(nCpuCores > parallel::detectCores())
    stop("The number of requested cpu's must be <= the number of actual cpu's.")

  stanData <- list( M = nPCs,
                    N = nrow(transData$robustPCs),
                    Z = transData$robustPCs[,1:nPCs])

  fileNames <- rstanParallelSampler(stanData, nWuSamples, nPwuSamples,
                                    nChainsPerCore, nCpuCores, procDir )

  return(list(nChains = nChainsPerCore * nCpuCores,
              nWuSamples = nWuSamples,
              nPwuSamples = nPwuSamples,
              fileNames = fileNames))
}


#' @title Plot selected traces
#'
#' @description Plot selected traces for each chain to assess whether
#' within-chain label switching has occurred.
#'
#' List containing, among other things, the names of the files in which the
#' \code{stanfit} objects are stored. These objects contain the samples
#' of the posterior pdf. This list is return by function
#' \code{\link{sampleFmm}}, for which the documentation includes
#' a complete description of container \code{samplePars}.
#'
#' @param procDir
#' Directory containing the files with the \code{stanfit} objects.
#'
#' @details
#' A set of three plots are generated for each chain. The set comprises
#' \itemize{
#'  \item A plot of one trace: the model proportion associated with
#'  pdf 1 (theta)
#'  \item A plot of two traces: A trace for element [1]
#'  of the mean vector for pdf 1 (mu1[1]), and another trace for
#'  element [1] of the mean vector for pdf 2 (mu2[1]).
#'  \item A plot of two traces: A trace for element [1]
#'  of the standard deviation vector for pdf 1 (tau1[1]), and another trace for
#'  element [1] of the standard deviation vector for pdf 2 (tau2[1]).
#' }
#'
#' @examples
#' \dontrun{
#' plotTraces(samplePars, procDir)
#' }
#'
#' @export
plotSelectedTraces <- function(samplePars, procDir = ".") {

  oldSetting <- devAskNewPage( ask = TRUE )
  on.exit(oldSetting, add = TRUE)

  sampleIndices <- samplePars$nWuSamples + 1:samplePars$nPwuSamples
  N <- length(samplePars$fileNames)
  for(i in 1:N){

    load( paste(procDir, "\\", samplePars$fileNames[i], sep = ""))
    theSamples <- rstan::extract(rawSamples, permuted = FALSE)

    df <- data.frame(indices = sampleIndices,
                     theta = theSamples[, 1, "theta"],
                     mu1 = theSamples[, 1, "mu1[1]"],
                     mu2 = theSamples[, 1, "mu2[1]"],
                     tau1 = theSamples[, 1, "tau1[1]"],
                     tau2 = theSamples[, 1, "tau2[1]"] )

    p1 <- ggplot2::ggplot( df,
                           ggplot2::aes(x = indices, y = theta),
                           environment = environment() ) +
      ggplot2::geom_line() +
      ggplot2::xlab("Sample index") +
      ggplot2::ylab("theta") +
      ggplot2::ggtitle(paste("Chain ", i, sep = ""))

    p2 <- ggplot2::ggplot( df,
                           ggplot2::aes(x = indices),
                           environment = environment() ) +
      ggplot2::geom_line(ggplot2::aes(y = mu1, colour = "1")) +
      ggplot2::geom_line(ggplot2::aes(y = mu2, colour = "2")) +
      ggplot2::scale_color_manual("Pdf", values = c("1" = "blue", "2" = "red")) +
      ggplot2::xlab("Sample index") +
      ggplot2::ylab("Element 1 of mean vectors")

    p3 <- ggplot2::ggplot( df,
                           ggplot2::aes(x = indices),
                           environment = environment() ) +
      ggplot2::geom_line(ggplot2::aes(y = tau1, colour = "1")) +
      ggplot2::geom_line(ggplot2::aes(y = tau2, colour = "2")) +
      ggplot2::scale_color_manual("Pdf", values = c("1" = "blue", "2" = "red")) +
      ggplot2::xlab("Sample index") +
      ggplot2::ylab("Element 1 of std dev vectors")

    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,1)))
    print(p1, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
    print(p2, vp=grid::viewport(layout.pos.row=2, layout.pos.col=1))
    print(p3, vp=grid::viewport(layout.pos.row=3, layout.pos.col=1))

  }
}

#' @title Plot point statistics
#'
#' @description Plot point statistics for each chain, the mode calculated with
#' function mclust, and the mode calculated with rstan.
#'
#' @param samplePars
#' List containing, among other things, the names of the files in which the
#' \code{stanfit} objects are stored. These objects contain the samples
#' of the posterior pdf. This list is return by function
#' \code{\link{sampleFmm}}, for which the documentation includes
#' a complete description of container \code{samplePars}.
#'
#' @param procDir
#' Directory containing the files with the \code{stanfit} objects.
#'
#' @param excludedChains Vector with the indices of the chains for which
#' the point statistics are not used to calculate plot ranges. See Details.
#'
#' @details
#' The point statistics are calculated and plotted for element 1 of the
#' two mean vectors, element 1 of the two standard deviation vectors,
#' theta (the proportion associated with pdf1), and the log-likelihood.
#'
#' To prevent a point statistic associated with a sampling chain
#' from being included in the calculation of
#' the plot range, the index of the associated chain is specified in
#' argument \code{excludedChains}.
#'
#' The calculation of the point statistics for every chain requires a
#' minute or two.
#'
#' @examples
#' \dontrun{
#' plotPointStats(samplePars)
#' }
#'
#' @export
plotPointStats <- function(samplePars, procDir = ".",
                           excludedChains = NULL ) {

  N <- length(samplePars$fileNames)

  statNames <-
    c("mu1[1]","mu2[1]","tau1[1]","tau2[1]","theta","log_lik", "lp__")

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

  Internal1 <- function( X, excludedChains, yLabel ) {
    if(is.null(excludedChains)) {
      yRange <- range( X, na.rm = TRUE )
    } else {
      yRange <- range( X[, -excludedChains], na.rm = TRUE )
    }

    solutions <- as.character(1:ncol(X))

    df <- data.frame(x = factor(solutions, levels = solutions),
                      y = X["0.5", ],
                      ymin = X["0.025",],
                      ymax = X["0.975",] )
    p <- ggplot2::ggplot(df,
                          ggplot2::aes(x = x, y = y, ymin = ymin, ymax = ymax),
                          environment = environment() ) +
      ggplot2::geom_pointrange() +
      ggplot2::ylim(yRange[1], yRange[2]) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, vjust=0.125)) +
      ggplot2::xlab("Chain") +
      ggplot2::ylab(yLabel)

    return(p)
  }

  Internal2 <- function( X, Y, excludedChains, yLabel ) {
    if(is.null(excludedChains)) {
      yRange <- range( X, Y, na.rm = TRUE )
    } else {
      yRange <- range( X[, -excludedChains], Y[, -excludedChains],
                       na.rm = TRUE )
    }

    solutions <- as.character(1:ncol(X))

    df <- data.frame(x = factor(solutions, levels = solutions),
                     y1 = X["0.5", ],
                     ymin1 = X["0.025",],
                     ymax1 = X["0.975",],
                     y2 = Y["0.5", ],
                     ymin2 = Y["0.025",],
                     ymax2 = Y["0.975",] )
    p <- ggplot2::ggplot(df, environment = environment() ) +
      ggplot2::geom_pointrange(ggplot2::aes(x = x, y = y1, ymin = ymin1, ymax = ymax1),
                               colour = "blue") +
      ggplot2::geom_pointrange(ggplot2::aes(x = x, y = y2, ymin = ymin2, ymax = ymax2),
                               colour = "red") +
      ggplot2::ylim(yRange[1], yRange[2]) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90, vjust=0.125)) +
      ggplot2::xlab("Chain") +
      ggplot2::ylab(yLabel)

    return(p)
  }

  p1 <- Internal2(theQuantiles[, "mu1[1]", ], theQuantiles[, "mu2[1]", ],
                  excludedChains, "Element 1 of the mean vectors")
  p2 <- Internal2(theQuantiles[, "tau1[1]", ], theQuantiles[, "tau2[1]", ],
                  excludedChains, "Element 1 of the standard deviation vectors")
  p3 <- Internal1(theQuantiles[, "theta", ], excludedChains, "theta")
  p4 <- Internal1(theQuantiles[, "log_lik", ], excludedChains, "Log-likelihood")

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
#' @param samplePars
#' List containing, among other things, the names of the files in which the
#' \code{stanfit} objects are stored. These objects contain the samples
#' of the posterior pdf. This list is return by function
#' \code{\link{sampleFmm}}, for which the documentation includes
#' a complete description of container \code{samplePars}.
#'
#' @param selectedChains  Dataframe listing the indices of the chains that
#'  are combined. (See Details).
#'
#' @param procDir
#' Directory containing the files with the \code{stanfit} objects.
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
#' indicate that chain 2 will be included in the combination and that
#' the variables in the
#' finite mixture model are not switched. That is, variable mu1 in the
#' chain is assigned to mu1, variable mu2 in the chain is assiged to mu2, and
#' so on. (2) Assume that
#' row 2 of \code{selectedChains} comprises the values 4 and TRUE. These
#' values
#' indicate that chain 4 will be included in the combination and
#' that the variables in the
#' finite mixture model are switched. That is, variable mu1 in the
#' chain is assigned to mu2, variable mu2 in the chain is assiged to mu1, and
#' so on.
#'
#' @return The returned value is a \code{stanfit} object,
#' which is described in the rstan documentation. This object comprises
#' the multiple chains, which have samples of
#' the following model parameters:
#' \describe{
#'  \item{theta}{Proportion of population associated with pdf 1.}
#'  \item{mu1}{Mean vector for pdf 1.}
#'  \item{mu2}{Mean vector for pdf 2.}
#'  \item{tau1}{Standard deviation vector for pdf 1.}
#'  \item{tau2}{Standard deviation vector for pdf 2.}
#'  \item{L_Omega1}{Lower triangle of the correlation matrix for pdf 1.}
#'  \item{L_Omega2}{Lower triangle of the correlation matrix for pdf 2.}
#'  \item{log_lik}{The logarithm of the likelihood function.}
#' }
#'
#' @references
#  Stan Development Team, 2015, Stan Modeling Language - User’s Guide
#' and Reference Manual, Version 2.6.0, available on line at
#' http://mc-stan.org/ (last accessed October 2015).
#'
#' @examples
#' \dontrun{
#' combinedChains <- combineChains(samplePars, selectedChains)
#' }
#'
#' @export
combineChains <- function(samplePars, selectedChains, procDir = ".") {

  sfList <- vector(mode = "list")
  for(k in 1:nrow(selectedChains)) {
    iChain <- selectedChains[k, "Chain"]
    load( paste(procDir, "\\", samplePars$fileNames[iChain], sep = ""))

    if(selectedChains[k, "isSwitched"] == TRUE) {

      rawSamples@sim$samples[[1]]$theta <- 1 - rawSamples@sim$samples[[1]]$theta

      N <- rawSamples@par_dims$mu1
      for(i in 1:N) {

        # mean vectors
        var1 <- paste("rawSamples@sim$samples[[1]]$\"mu1[", i, "]\"", sep="")
        tmp1 <- eval(parse(text = var1))

        var2 <- paste("rawSamples@sim$samples[[1]]$\"mu2[", i, "]\"", sep="")
        tmp2 <- eval(parse(text = var2))

        eval(parse(text = paste(var1, " <- tmp2")))

        eval(parse(text = paste(var2, " <- tmp1")))

        # standard deviation vectors
        var1 <- paste("rawSamples@sim$samples[[1]]$\"tau1[", i, "]\"", sep="")
        tmp1 <- eval(parse(text = var1))

        var2 <- paste("rawSamples@sim$samples[[1]]$\"tau2[", i, "]\"", sep="")
        tmp2 <- eval(parse(text = var2))

        eval(parse(text = paste(var1, " <- tmp2")))

        eval(parse(text = paste(var2, " <- tmp1")))

        # lower triangles of the cholesky decompositions of the correlation matrices
        for(j in 1:N) {
          var1 <- paste("rawSamples@sim$samples[[1]]$\"L_Omega1[", i, ",", j, "]\"", sep="")
          tmp1 <- eval(parse(text = var1))

          var2 <- paste("rawSamples@sim$samples[[1]]$\"L_Omega2[", i, ",", j, "]\"", sep="")
          tmp2 <- eval(parse(text = var2))

          eval(parse(text = paste(var1, " <- tmp2")))

          eval(parse(text = paste(var2, " <- tmp1")))

        }
      }

    }
    sfList[[k]] <- rawSamples

  }

  return(rstan::sflist2stanfit(sfList))
}

#' @title Calculate conditional probabilities
#'
#' @description Calculate, for each field sample, Monte Carlo samples of
#' the conditional probability that the field sample is
#' associated with the first probability density function in the finite
#' mixture model.
#'
#' @param transData     List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}, for which the
#' documentation includes a complete description of container \code{transData}.
#' @param nPCS          Number of principal components that were used in the
#'                      finite mixture model.
#' @param combinedChains
#' A \code{stanfit} object containg multiple Monte Carlo chains. This
#' object is return by function \code{\link{combineChains}}, for which the
#' documentation includes a complete description of container
#' \code{combinedChains}.
#'
#' @details
#' The conditional probabilities are calculated using the formula in Gelman
#' et al. (2014, p. 540 bottom).
#'
#' @return A matrix containing the Monte Carlo samples of the
#' conditional probabilites. The number of matrix columns equals the
#' number of field samples, and the number of matrix rows equals the
#' number of Monte Carlo samples in \code{combinedChains}. That is,
#' column j of the matrix contains Monte Carlo samples of the conditional
#' probability that field sample j is associated with the first probability
#' density function in the finite mixture model.
#'
#' @references
#' Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A.,
#' and Rubin, D.B., 2014, Bayesian data analysis (3rd ed.),
#' CRC Press.
#'
#' @examples
#' \dontrun{
#' condProbs1 <- calcCondProbs1(transData, nPCs, combinedChains)
#' }
#'
#'
#' @export
calcCondProbs1 <- function(transData, nPCs, combinedChains) {

  theSamples <- rstan::extract(combinedChains)

  nMCSamples <- length(theSamples$theta)
  nFldSamples <- nrow(transData$robustPCs)

  condProb <- matrix(NA_real_, ncol = nFldSamples, nrow = nMCSamples )

  for(i in 1:nMCSamples){

    L_Sigma1 <- theSamples$tau1[i, ] * theSamples$L_Omega1[i, , ]
    Sigma1 <- L_Sigma1 %*% t(L_Sigma1)
    ps1 <- log(theSamples$theta[i]) +
      mvtnorm::dmvnorm( transData$robustPCs[,1:nPCs],
                        mean = theSamples$mu1[i, ],
                        sigma = Sigma1, log = TRUE )

    L_Sigma2 <- theSamples$tau2[i, ] * theSamples$L_Omega2[i, , ]
    Sigma2 <- L_Sigma2 %*% t(L_Sigma2)
    ps2 <- log(1-theSamples$theta[i]) +
      mvtnorm::dmvnorm( transData$robustPCs[,1:nPCs],
                        mean = theSamples$mu2[i, ],
                        sigma = Sigma2, log = TRUE )

    log_denom <- log(exp(ps1) + exp(ps2))
    condProb[i,] <- exp(ps1 - log_denom)

  }
  return(condProb)
}

#' @title Calculate observed test statistics for model checking
#'
#' @description Calculate test statistics for posterior
#' predictive checking
#' of the parameters in the finite mixture model. The test statistics pertain
#' to the observed data that have undergone both the isometric log-ratio
#' transform and the robust, principal component transform.
#'
#' @param transData     List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}, for which the
#' documentation includes a complete description of container \code{transData}.
#' @param nPCS          Number of principal components that were used in the
#'                      finite mixture model.
#' @param condProbs1
#' A matrix containing the Monte Carlo samples of the
#' conditional probabilites. This matrix is returned by function
#' \code{\link{calcCondProbs}}, for which the documentation includes a
#' complete descriptoin of container \code{condProbs1}.
#'
#' @details
#' Test statistics are defined in Gelman et al. (2014, p. 145) and
#' are used for posterior predictive checking of the finite mixture model.
#'
#' @return A list with the following test quantities is returned.
#' @return \item{mu1}{Mean vector for pdf 1.}
#' @return \item{mu2}{Mean vector for pdf 2.}
#' @return \item{tau1}{Standard deviation vector for pdf 1.}
#' @return \item{tau2}{Standard deviation vector for pdf 2.}
#' @return \item{Corr1}{Correlation matrix for pdf 1.}
#' @return \item{Corr2}{Correlation matrix for pdf 2.}
#'
#' @references
#' Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B.,
#' Vehtari, A., and Rubin, D.B., 2014, Bayesian data analysis (3rd ed.):
#' CRC Press.
#'
#' @examples
#' \dontrun{
#' obsTestStats <- calcObsTestStats(transData, nPCS, condProbs1)
#' }
#'
#'
#' @export
calcObsTestStats <- function(transData, nPCS, condProbs1) {

  # associated with pdf 1
  S1 <- cov.wt(transData$robustPCs[,1:nPCs], wt = colMeans(condProbs1) )

  # associated with pdf 2
  condProbs2 <- 1 - condProbs1
  S2 <- cov.wt(transData$robustPCs[,1:nPCs], wt = colMeans(condProbs2) )

  return(list(mu1 = S1$center,
              mu2 = S2$center,
              tau1 = sqrt(diag(S1$cov)),
              tau2 = sqrt(diag(S2$cov)),
              Corr1 = cov2cor(S1$cov),
              Corr2 = cov2cor(S2$cov)))

}


#' @title Plot model check --- means and standard deviations
#'
#' @description Perform posterior predictive checking of the
#' finite mixture model. The checks are performed for
#' the mean vector and the standard deviation vector of both pdfs.
#' The test quantities for the observed data were calculated by function
#' \code{\link{obsTestQuantities}}. The test quantities for the
#' replicated data are the Monte Carlo samples for the parameters in the
#' finite mixture model.
#' This function plots the test quantities so that they can be compared.
#'
#' @param combinedChains
#' A \code{stanfit} object containg multiple Monte Carlo chains. This
#' object is return by function \code{\link{combineChains}}, for which the
#' documentation includes a complete description of container
#' \code{combinedChains}.
#'
#' @param obsTestQuantities
#' List containing the test quantities for the observed data.
#' This list is return by function \code{\link{calcObsTestQuantities}},
#' for which the documentation includes a
#' complete description of container \code{obsTestQuantities}.
#'
#' @param intervalPercentage Interval for the distributions of the test
#' quantity. Typical values might be 50, 90, or 95.
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
#' @references
#' Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B.,
#' Vehtari, A., and Rubin, D.B., 2014, Bayesian data analysis (3rd ed.):
#' CRC Press.
#'
#' @examples
#' \dontrun{
#' plotModelCheck_MS( combinedChains, obsTestQuantities)
#' }
#'
#' @export
plotModelCheck_MS <- function( combinedChains, obsTestQuantities,
                                intervalPercentage = 95) {

  Internal1 <- function(Tobs, Trep, intervalPercentage, plotTitle){

    nComponents <- ncol(Trep)
    nMcSamples <- nrow(Trep)

    tailPercentage <- 0.5*(100.0-intervalPercentage)
    probs <- c(tailPercentage, 50, 100.0-tailPercentage)/100.0
    theQuantiles <- t(apply(Trep, 2, quantile, probs = probs, names = FALSE))

    pValues <- vector(mode="numeric", length=nComponents )
    for(j in 1:nComponents) {
      # defined in Gelman et al., p. 146
      tmp <- sum(Trep[, j] > Tobs[j]) / nMcSamples
      pValues[j] <- ifelse(tmp <= 0.5, tmp, 1-tmp)
    }
    pValues <- round(pValues, digits=2)

    df <- data.frame(x = 1:nComponents,
                     Trep_50 = theQuantiles[, 2],
                     Trep_min = theQuantiles[, 1],
                     Trep_max = theQuantiles[, 3],
                     T_obs = Tobs,
                     pValues = pValues)

    p <- ggplot2::ggplot(df,
                         ggplot2::aes(x = x),
                         environment = environment()) +
      ggplot2::geom_point(ggplot2::aes(y = Tobs),
                          shape = 16, size = 3, colour="#FF6C91") +
      ggplot2::geom_pointrange(ggplot2::aes(y = Trep_50, ymin = Trep_min,
                                            ymax = Trep_max), shape = 3) +
      ggplot2::xlab("Component") +
      ggplot2::ylab("Transformed concentration (no units)") +
      ggplot2::ggtitle(plotTitle)

    p <- p + ggplot2::geom_text(ggplot2::aes(y = max(Trep_max),
                                             label = pValues),
                                size = 3, vjust = 0, angle = 90,
                                colour = ifelse(pValues < 0.10, "red", "black"))
    return(p)
  }


  p_mu1 <- Internal1( obsTestQuantities$mu1,
                      rstan::extract(combinedChains, pars="mu1")$mu1,
                      intervalPercentage,
                      plotTitle = "Mean for pdf 1")
  p_tau1 <- Internal1( obsTestQuantities$tau1,
                       rstan::extract(combinedChains, pars="tau1")$tau1,
                       intervalPercentage,
                       plotTitle = "Sd for pdf 1")
  p_mu2 <- Internal1( obsTestQuantities$mu2,
                      rstan::extract(combinedChains, pars="mu2")$mu2,
                      intervalPercentage,
                      plotTitle = "Mean for pdf 2")
  p_tau2 <- Internal1( obsTestQuantities$tau2,
                       rstan::extract(combinedChains, pars="tau2")$tau2,
                       intervalPercentage,
                       plotTitle = "Sd for pdf 2")

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,2)))
  print(p_mu1, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
  print(p_mu2, vp=grid::viewport(layout.pos.row=1, layout.pos.col=2))
  print(p_tau1, vp=grid::viewport(layout.pos.row=2, layout.pos.col=1))
  print(p_tau2, vp=grid::viewport(layout.pos.row=2, layout.pos.col=2))

}

#' @title Plot model check --- correlation matrices
#'
#' @description Perform posterior predictive checking of the
#' finite mixture model. The checks are performed for
#' the correlation matrices of both pdfs.
#' The test quantities for the observed data were calculated by function
#' \code{\link{obsTestQuantities}}. The test quantities for the
#' replicated data are the Monte Carlo samples for the parameters in the
#' finite mixture model.
#' This function plots the test quantities so that they can be compared.
#'
#' @param combinedChains
#' A \code{stanfit} object containg multiple Monte Carlo chains. This
#' object is return by function \code{\link{combineChains}}, for which the
#' documentation includes a complete description of container
#' \code{combinedChains}.
#'
#' @param obsTestQuantities
#' List containing the test quantities for the observed data.
#' This list is return by function \code{\link{calcObsTestQuantities}},
#' for which the documentation includes a
#' complete description of container \code{obsTestQuantities}.
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
#' The plot for the correlation matrix has two parts, which are appear on the
#' left and the right sides of the graphic. The left side presents a comparison
#' of the correlation matrices:  Below the diagonal, which is solid red, is
#' one triangle of the correlation matrix that is calculated from the field data.
#' Above the diagonal is one triangle of the correlation matrix that is the
#' median correlation matrix calculated from the Monte Carlo samples.
#' The right side presents the associated p-values.
#' The color scale ranges from the smallest calculated p-value to 0.5,
#' which is the largest possible p-value (based on the previous definition).
#'
#' @references
#' Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B.,
#' Vehtari, A., and Rubin, D.B., 2014, Bayesian data analysis (3rd ed.):
#' CRC Press.
#'
#' @examples
#' \dontrun{
#' plotModelCheck_C(combinedChains, obsTestQuantities)
#' }
#'
#' @export
plotModelCheck_C <- function( combinedChains, obsTestQuantities) {

  Internal1 <- function(obsCorr, repCorr, plotTitle) {
    Y <- apply(repCorr, c(2,3), median)

    M <- nrow(obsCorr)
    Z <- matrix(1.0, nrow = M, ncol = M, dimnames = list( 1:M, 1:M))
    Z[upper.tri(Z)] <- obsCorr[upper.tri(obsCorr)]
    Z[lower.tri(Z)] <- Y[lower.tri(Y)]

    X <- reshape2::melt(Z)
    X <- na.omit(X)

    w <- ggplot2::ggplot(X, ggplot2::aes(Var2, Var1))+
      ggplot2::geom_tile(data=X, ggplot2::aes(fill=value), color="white")+
      ggplot2::scale_fill_gradient2(low="blue", high="red", mid="white",
                                    midpoint=0, limit=c(-1,1),
                                    name="Correlation\n(Pearson)")+
      ggplot2::xlab("Component") +
      ggplot2::ylab("Component") +
      ggplot2::coord_equal() +
      ggplot2::ggtitle(plotTitle)

    return(w)
  }

  Internal2 <- function(obsCorr, repCorr, plotTitle) {

    theDimensions <- dim(repCorr)
    nMcSamples <- theDimensions[1]
    M <- theDimensions[2]

    pValues <- matrix(NA_real_, nrow = M, ncol = M, dimnames = list( 1:M, 1:M) )
    for(i in 1:(M-1)) {
      for(j in (i+1):M) {
        tmp <- sum(repCorr[ , i, j ] > obsCorr[i,j]) / nMcSamples
        pValues[i,j] <- ifelse(tmp <= 0.5, tmp, 1-tmp)
      }
    }

    # The rational for the following code is presented in function Internal1.
    pValues <- pValues[rev(rownames(pValues)),]

    Z <- reshape2::melt(pValues)
    Z <- na.omit(Z)

    w <- ggplot2::ggplot(Z, ggplot2::aes(Var2, Var1))+
      ggplot2::geom_tile(data=Z, ggplot2::aes(fill=value), color="white")+
      ggplot2::scale_fill_gradient(limit=c(min(Z$value),0.5),
                                   high = "#132B43", low = "#56B1F7",
                                   name="P-value")+
      ggplot2::xlab("Component") +
      ggplot2::ylab("Component") +
      ggplot2::coord_equal() +
      ggplot2::ggtitle(plotTitle)

    return(w)
  }

  calcCorrMatrixSamples <- function(L_Omega) {

    corrMatrixSamples <- array(NA_real_, dim=dim(L_Omega) )
    nMcSamples <- dim(L_Omega)[1]
    for(i in 1:nMcSamples) {
      corrMatrixSamples[i, , ] <- L_Omega[i, , ] %*% t(L_Omega[i, , ])
    }
    return(corrMatrixSamples)
  }

  L_Omega1 <- rstan::extract(combinedChains, pars="L_Omega1")$L_Omega1
  corrMatrixSamples <- calcCorrMatrixSamples(L_Omega1)

  p1 <- Internal1(obsTestQuantities$Corr1, corrMatrixSamples,
                  "Correlation matrices for pdf 1")
  p2 <- Internal2(obsTestQuantities$Corr1, corrMatrixSamples,
                  "P-values, pdf 1" )

  L_Omega2 <- rstan::extract(combinedChains, pars="L_Omega2")$L_Omega2
  corrMatrixSamples <- calcCorrMatrixSamples(L_Omega2)

  p3 <- Internal1(obsTestQuantities$Corr2, corrMatrixSamples,
                  "Correlation matrices for pdf 2")
  p4 <- Internal2(obsTestQuantities$Corr2, corrMatrixSamples,
                  "P-values, pdf 2" )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,2)))
  print(p1, vp=grid::viewport(layout.pos.row=1, layout.pos.col=1))
  print(p2, vp=grid::viewport(layout.pos.row=1, layout.pos.col=2))
  print(p3, vp=grid::viewport(layout.pos.row=2, layout.pos.col=1))
  print(p4, vp=grid::viewport(layout.pos.row=2, layout.pos.col=2))

}


#' @title Back transform to the simplex
#'
#' @description Back transform the mean vectors and covariance matrices
#' of the finite mixture model to the correpsonding quantities for the simplex
#' (that is, compostional means and variation matrices).
#'
#' @param gcData
#' List containing the geochemical and related data. This container is
#' described in the package documentation.
#'
#' @param nPCS          Number of principal components that were used in the
#'                      finite mixture model.
#'
#' @param transData
#' List containing the transformed geochemical
#' concentrations and related information.
#' This list is return by function \code{\link{transformGcData}}, for which the
#' documentation includes a complete
#' description of container \code{transData}.
#'
#' @param combinedChains
#' A \code{stanfit} object containg multiple Monte Carlo chains. This
#' object is return by function \code{\link{combineChains}}, for which the
#' documentation includes a complete description of container
#' \code{combinedChains}.
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
#' simplexModPar <- backTransform( gcData, nPCs, transData, combinedChains)
#' }
#'
#' @export

backTransform <- function(gcData, nPCs, transData, combinedChains) {

  elementNames <- names(gcData$concData)

  # dimension within the simplex
  D <- length(elementNames)

  theSamples <- rstan::extract(combinedChains)

  nMcSamples <- length(theSamples$theta)

  # subset now, so that it won't be done repeatedly
  rEV <- transData$robustEigenvectors[,1:nPCs]
  rEVt <- t(rEV)

  # undo the rotation and translation, which are associated
  # with the principal component transformation
  tmp.center <- rep.int( 1, nMcSamples ) %o% transData$robustIlrCenter

  mu1.tmp <- theSamples$mu1 %*%  rEVt + tmp.center
  mu2.tmp <- theSamples$mu2 %*%  rEVt + tmp.center

  # Back-transformation to concentrations
  compMean1 <- CalcInvIlr( mu1.tmp, t(transData$Psi), kappa = gcData$constSumValue )
  compMean2 <- CalcInvIlr( mu2.tmp, t(transData$Psi), kappa = gcData$constSumValue )

  colnames(compMean1) <- elementNames
  colnames(compMean2) <- elementNames

  varMatrix1 <- array( NA_real_, dim=c(nMcSamples, D, D),
                       dimnames = list(NULL, elementNames, elementNames))
  varMatrix2 <- varMatrix1

  for(i in 1:nMcSamples) {

    # undo the rotation, which is associated
    # with the principal component transformation

    L_Sigma1 <- theSamples$tau1[i, ] * theSamples$L_Omega1[i, , ]
    Sigma1 <- L_Sigma1 %*% t(L_Sigma1)
    Sigma1.tmp <- rEV %*% Sigma1 %*% rEVt

    L_Sigma2 <- theSamples$tau2[i, ] * theSamples$L_Omega2[i, , ]
    Sigma2 <- L_Sigma2 %*% t(L_Sigma2)
    Sigma2.tmp <- rEV %*% Sigma2 %*% rEVt

    # Back-transformation to concentrations within the simplex
    varMatrix1[i, , ] <- InvCovTransform( Sigma1.tmp, transData$Psi )
    varMatrix2[i, , ] <- InvCovTransform( Sigma2.tmp, transData$Psi )
  }

  return( list( compMean1 = compMean1,
                compMean2 = compMean2,
                varMatrix1 = varMatrix1,
                varMatrix2 = varMatrix2))

}

#' @title Plot standardized compositional means
#'
#' @description Plot standardized compositional means for the pdfs in the
#' finite mixture model. The standardization is based on the geochemical
#' concentrations that were used in the finite mixture model.
#'
#' @param simplexModPar
#' List containing Monte Carlo samples of the selected parameters
#' in the finite mixture model.
#' These parameters (namely, the mean vector and covariance matrix for
#' each pdf) are expressed in terms of their equivalent values in the
#' simplex (namely, the compositional mean vector and the variation matrix
#' for each pdf).
#' This list is return by function \code{\link{backTranform}}, for which the
#' documentation includes a complete description of this container.
#'
#' @param simplexStats
#' List containing statistics for the simplex, which will be used for the
#' standardization. This list is return by function
#' \code{\link{calcSimplexStats}}, for which the
#' documentation includes a complete description of this container.
#'
#' @param gcData
#' List containing the geochemical and related data. This container is
#' described in the package documentation.
#'
#' @param elementOrder Vector specifying the order in which the elements are
#' plotted.
#'
#' @param intervalPercentage Interval for the distributions of the standardized
#' compositional means. Typical values are 50, 90, or 95.
#'
#' @param symbolSize The size of the plotting symbol.
#'
#' @details
#' The standardized compositional mean is a vector. The Monte Carlo samples
#' of this vector comprise Monte Carlo samples of each vector element.
#' The later must be summarized so that they can be visualized.
#' To this end, the 0.025, 0.5, and 0.975 quantiles of the samples for
#' each vector element is computed. (The 0.025 and 0.975 quantiles correspond
#' to the default value of argument intervalPercentage.) The compositional
#' operation of closure is not applied to these quantiles.
#' Each plot symbol represents the distribution for a vector element.
#' The vertical line within a symbol represents the range of the 0.025 and
#' 0.975 quantiles; the dot within a symbol represents the 0.5 quantile.
#'
#' @references
#' Pawlowsky-Glahn, V., Egozcue, J.J., and Tolosana-Delgado, R., 2015, Modeling
#' and analysis of compositional data: John Wiley and Sons, Ltd.
#'
#' @examples
#' \dontrun{
#' plotStdCompMeans(simplexModPar, simplexStats, gcData, elementOrder)
#' }
#'
#' @export
plotStdCompMeans <- function(simplexModPar, simplexStats, gcData, elementOrder,
                          intervalPercentage = 95, symbolSize = 0.75) {

  Internal1 <- function(compMeans, kappa, elementOrder,
                        center, metricVariance,
                        interval, pdf) {

    # reorder the compositional means and then standardize them
    tmp <- compMeans[, elementOrder]
    tmp <- Perturb(center[elementOrder]^(-1), tmp, kappa = kappa)
    stdCompMeans <- Power(tmp, 1/sqrt(metricVariance), kappa = kappa)

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
  barycenter <- gcData$constSumValue / length(simplexStats$sampleCenter)

  tailPercentage <- 0.5*(100.0-intervalPercentage)
  interval <- c(tailPercentage,100.0-tailPercentage)/100.0

  df1 <- Internal1(simplexModPar$compMean1, gcData$constSumValue, elementOrder,
                   simplexStats$sampleCenter, simplexStats$metricVariance,
                   interval, 1 )
  df2 <- Internal1(simplexModPar$compMean2, gcData$constSumValue, elementOrder,
                   simplexStats$sampleCenter, simplexStats$metricVariance,
                   interval, 2 )

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
#' List containing Monte Carlo samples of the selected parameters
#' in the finite mixture model.
#' These parameters (namely, the mean vector, the standard deviation vector,
#'  and the correlation matrix for
#' each pdf) are expressed in terms of their equivalent values in the
#' simplex (namely, the compositional mean vector and the variation matrix
#' for each pdf).
#' This list is return by function \code{\link{backTransform}}, for which the
#' documentation includes a complete description of this container.
#'
#' @param elementOrder
#' Vector specifying the order in which the elements are plotted.
#'
#' @param symbolSize
#' The size of the plotting symbol.
#'
#' @details
#' The compositional mean is a vector. The Monte Carlo samples of this
#' vector comprise Monte Carlo samples of each vector element. The later
#' must be summarized so that they can be visualized. To this end, the
#' median of the samples for each vector element is computed. The compositional
#' operation of closure is not applied to the medians.
#'
#' Although it is possible to compute, say, the 95% credible interval for
#' each vector element, the interval is smaller than the plot symbol for the
#' median.
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
#' List containing Monte Carlo samples of the selected parameters
#' in the finite mixture model.
#' These parameters (namely, the mean vector and covariance matrix for
#' each pdf) are expressed in terms of their equivalent values in the
#' simplex (namely, the compositional mean vector and the variation matrix
#' for each pdf).
#' This list is return by function \code{\link{backTranform}}, for which the
#' documentation includes a complete description of this container.
#' @param elementOrder
#' Vector specifying the order in which the elements are plotted.
#' @param colorScale
#' Character string specifying the color scale for
#' the plot. The choices are either "spectrum" (default) and "rainbow."
#'
#' @details
#' In the plot, the upper triangle is the upper triangle from the
#' variation matrix for pdf 1, and the lower triange is the lower triangle
#' from the variation matrix for pdf 2.
#'
#' The pixels represent scaled variances of the log-ratios
#' between the respective chemical elements. The scaling is desirable because
#' it reduces the large range of the variances, making it easier to
#' visualize all of the variances together. The scaling function is the
#' square root; so, the pixels in the plot strictly represent standard devations
#' of the log-ratios between the respective chemical elements.
#'
#' @references
#' Pawlowsky-Glahn, V., Egozcue, J.J., and Tolosana-Delgado, R., 2015, Modeling
#' and analysis of compositional data: John Wiley and Sons, Ltd.
#'
#' @examples
#' \dontrun{
#' plotSqrtVarMatrix(simplexModPar, elementOrder, colorScale = "rainbow")
#' }
#'
#' @export
plotSqrtVarMatrices <- function(simplexModPar, elementOrder,
                                colorScale = "spectrum" ) {

  # D is the standard notation, and is concise.
  D <- dim(simplexModPar$varMatrix1)[3]


  medVarMatrix1 <- apply(simplexModPar$varMatrix1, c(2,3), median)
  medVarMatrix2 <- apply(simplexModPar$varMatrix2, c(2,3), median)

  # reorder
  medVarMatrix1 <- medVarMatrix1[elementOrder, elementOrder]
  medVarMatrix2 <- medVarMatrix2[elementOrder, elementOrder]

  Z <- matrix( NA_real_, nrow=D, ncol=D, dimnames=list(elementOrder,elementOrder) )
  Z[upper.tri(Z)] <- medVarMatrix1[upper.tri(medVarMatrix1)]
  Z[lower.tri(Z)] <- medVarMatrix2[lower.tri(medVarMatrix2)]

  # If matrix Z were plotted in its current configuration,
  # the first row would be at the bottom, the second row
  # would be the second from the bottom, and so on. Consequently,
  # the diagonal would extend from the lower-left corner to the
  # upper-right corner. So, the matrix is plotted in an unfamilar way.
  #
  # The problem is corrected by reversing the rows.
  Z <- Z[rev(rownames(Z)),]

  Z <- sqrt(Z)

  Z <- reshape2::melt(Z)
  Z <- na.omit(Z)

  w <- ggplot2::ggplot(Z, ggplot2::aes(Var2, Var1)) +
    ggplot2::geom_tile(data=Z, ggplot2::aes(fill=value), color="white") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90,vjust=0.25,colour = "black")) +
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
#' @param gcData
#' List containing the geochemical and related data. This container is
#' described in the package documentation.
#'
#' @param condProbs1
#' A matrix containing the Monte Carlo samples of the
#' conditional probabilites. This matrix is returned by function
#' \code{\link{calcCondProbs}}, for which the documentation includes a
#' complete descriptoin of container \code{condProbs1}.
#'
#' @param probIntervals Vector containing intervals of conditional
#' probability. All field samples within an given interval are plotted the same
#' way.
#'
#' @param symbolIndices Vector containing the indices of the plotting symbols
#' for the conditional probability intervals.
#'
#' @param symbolSizes Vector containing the relative sizes of the plotting symbols
#' for the conditional probability intervals.
#'
#' @param symbolColors Vector containing the colors of the plotting symbols
#' for the conditional probability intervals.
#'
#' @details
#' The conditional probabilities indicate the extend
#' to which the field samples are associated with the first pdf in the finite
#' mixture model. The conditional probabilities
#' in container \code{condProbs1} are Monte Carlo samples,
#' and their medians are used to assign plotting attributes for the field
#' samples.
#'
#' The plotting attributes are specified by arguments \code{probIntervals},
#' \code{symbolIndices}, \code{symbolSizes}, and \code{symbolColors}. To
#' understanding the specification of these attributes, consider their default
#' values, which pertain to four probability intervals.
#'
#' Vector \code{probIntervals} has elements 0, 0.1, 0.5, 0.9, and 1. These five
#' elements specify four probability intervals: [0,0.1], [0.1,0.5],
#' [0.5,0.9] and
#' [0.9,1]. Notice that the first and last elements of \code{probIntervals}
#' are 0 and 1 respectively. The probability intervals are used to
#' classify the field samples based upon their associated conditional
#' probabilities:
#' \itemize{
#'  \item If the conditional probability of a field sample is
#'  within the interval [0,0.1], then the field sample is classified as
#'  "strongly associated with pdf 2" and is assigned the color red. This
#'  color is consistent with the colors used in functions
#'  \code{\link{plotStdCompMeans}} and \code{\link{plotCompMeans}}.
#'  \item If the conditional probability of a field sample is
#'  within the interval [0.1,0.5], then the field sample is classified as
#'  "moderately associated with pdf 2" and is assigned the color yellow.
#'  \item If the conditional probability of a field sample is
#'  within the interval [0.5,0.9], then the field sample is classified as
#'  "moderately associated with pdf 1" and is assigned the color green.
#'  \item If the conditional probability of a field sample is
#'  within the interval [0.9,1], then the field sample is classified as
#'  "strongly associated with pdf 1" and is assigned the color blue. This
#'  color is consistent with the colors used in functions
#'  \code{\link{plotStdCompMeans}} and \code{\link{plotCompMeans}}.
#' }
#'
#' Because, in this explanation, vector \code{probIntervals} specifies four
#' probability intervals, arguments for vectors \code{symbolIndices},
#' \code{symbolSizes}, and \code{symbolColors} must have four elements.
#' The symbol indices, symbol sizes, and symbol colors are described in
#' Murrell (2006, p. 55-56, 68, 69).
#'
#' This function adds symbols to a map that has already been plotted.
#'
#' @references
#' Murrell, P., 2006, R graphics: Chapman & Hall / CRC.
#'
#' @examples
#' \dontrun{
#' map('state', fill = TRUE, col = "gray60", border = "white")
#' plotClusters(concentrationData, condProbs1)
#' }
#'
#' @export
plotClusters <- function(gcData, condProbs1,
                        probIntervals = c( 0, 0.1, 0.5, 0.9, 1.0 ),
                        symbolIndices = c( 16, 16, 16, 16 ),
                        symbolSizes = c( 1/3, 1/3, 1/3, 1/3 ),
                        symbolColors = c( "red", "yellow", "green", "blue" ) ) {

  locations <- coordinates( gcData$concData )

  # median of conditional probabilites
  g <- apply(condProbs1, 2, median)

  for (i in 1:(length(probIntervals)-1)) {

    areInInterval <- probIntervals[i] <= g & g <= probIntervals[i+1]

    if( sum(areInInterval) == 0 ) next

    locations.sp <- sp::SpatialPoints( locations[areInInterval, ],
                                       proj4string = CRS( "+proj=longlat +ellps=WGS84" ) )
    sp::plot( locations.sp, add=TRUE, pch=symbolIndices[i],
              col=symbolColors[i], cex=symbolSizes[i] )

  }

}

#' @title Split the geochemical data
#'
#' @description The geochemical data, which have been clustered, are split
#' into two groups based on their conditional probabilites.
#'
#' @param gcData
#' List containing the geochemical and related data. This container is
#' described in the package documentation.
#'
#' @param condProbs1
#' A matrix containing the Monte Carlo samples of the
#' conditional probabilites. This matrix is returned by function
#' \code{\link{calcCondProbs}}, for which the documentation includes a
#' complete descriptoin of container \code{condProbs1}.
#'
#' @param threshold
#' The threshold used to split the data into two groups. (See details.)
#'
#' @details
#' For each
#' field sample, the median of the Monte Carlo samples of conditional
#' probability
#' is calculated. If this median is between 1-threshold and 1, then the
#' field sample is associated with pdf 1 in the finite mixture model. However,
#' if this median is between 0 and threshold, then the field samples is
#' associated with pdf 2 in the finite mixture model. This criterion is
#' used to split the field samples into two groups.
#'
#' The variable threshold must be greater than 0 and less than 0.5.
#'
#' @return A list with two components is returned.
#' @return \item{gcData1}{List containing the geochemical and related data that
#' are associated with pdf 1.
#' This container is described in the package documentation.}
#' @return \item{gcData2}{List containing the geochemical and related data that
#' are associated with pdf 2.
#' This container is described in the package documentation.}
#'
#' @examples
#' \dontrun{
#' theSplits <- splitGcData(gcData, condProb1, threshold = 0.10 )
#' }
#'
#' @export
splitGcData <- function(gcData, condProbs1, threshold = 0.10 ) {

  if(threshold <= 0.0 || 0.5 <= threshold)
    stop("Argument threshold must be > 0 and < 0.50.")

  g <- apply(condProbs1, 2, median)

  areInPdf1 <- 1.0-threshold <= g & g <= 1.0
  areInPdf2 <- 0.0 <= g & g <= threshold

  gcData1 <- list(concData = gcData$concData[areInPdf1, ],
                 censorIndicators = gcData$censorIndicators[areInPdf1, ],
                 constSumValue = gcData$constSumValue )

  gcData2 <- list(concData = gcData$concData[areInPdf2, ],
                 censorIndicators = gcData$censorIndicators[areInPdf2, ],
                 constSumValue = gcData$constSumValue )

  return( list(gcData1 = gcData1, gcData2 = gcData2))

}


