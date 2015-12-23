/*
Description: Finite mixture modeling of geochemical data. The model has only two
probability density functions (pdfs)). The geochemical data have
been transformed with the isometric, log-ratio transform and the
robust, principal components transform.

Arguments
M             number of principal components
N             number of field samples
Z             principal components
priorParams   parameters for the prior pdfs. See Details.

Details
Argument priorParams is a vector with 4 elements.
1. The first element pertains to prior pdf for theta, which is the proportion
in the finite mixture model. The prior pdf is a beta pdf, which has two shape
parameters. The shape parameters must be equal so that the beta pdf is
symmetric with respect to 0.5. Also, the shape parameters must be greater than
1, so that the pdf has a peak at 0.5 and equals 0 at both 0 and 1. Both shape
parameters are specified by the single value, priorParams[1].
2. The second element pertains to the prior pdf for the elements of mu1 and
mu2, which are the mean vectors in the finite mixture model. The prior pdf
is a normal pdf, for which the mean is 0 and the standard devation is
specified by priorParams[2]. Of course, priorParams[2] must be greater than 0.
3. The third element pertains to the prior pdf for the elements of tau1 and
tau2, which are the standard deviation vectors in the finite mixture model.
The prior pdf is a truncated Cauchy pdf: The Cauchy pdf before truncation has
a center of 0 and a scale specified by priorParams[3]. The truncation point is
0, which removes negative values of the random variable. Of course,
priorParams[3] must be greater than 0.
4. The fourth element pertains to the prior pdf for Omega1 and Omega2,
which are the correlation matrices in the finite mixture model.
The prior pdf is the LKJ correlation distribution, which has one shape
parameter. When the shape parameter is greater than 1, the LKJ correlation
distribution has a mode corresponding to the identity matrix; as the shape
parameter increases, the LKJ correlation distribution becomes increasingly
concentrated about this mode. The chosen shape parameter always should be
greater than 1. The shape parameter is specified by priorParams[4].

The previously described constraints on the parameters for the prior pdfs
are not checked.
*/
data {
  int<lower=1> M ;
  int<lower=1> N ;
  matrix[N,M] Z ;
  vector[4] priorParams ;
}
parameters {
  /*
  theta       proportion associated with the first pdf
              (The proportion associated with the second pdf is 1-theta.)
  mu1         mean vector for pdf 1
  mu2         mean vector for pdf 2
  tau1        standard deviation vector for pdf 1
  tau2        standard deviation vector for pdf 2
  L_Omega1    Cholesky factorization of the correlation matrix for pdf 1
  L_Omega2    Cholesky factorization of the correlation matrix for pdf 2
  */

  real<lower=0,upper=1> theta ;
  vector[M] mu1 ;
  vector[M] mu2 ;
  vector<lower=0>[M] tau1 ;
  vector<lower=0>[M] tau2 ;
  cholesky_factor_corr[M] L_Omega1 ;
  cholesky_factor_corr[M] L_Omega2 ;
}
transformed parameters {
  /*
  L_Sigma1    Cholesky factorization of the covariance matrix for pdf 1
  L_Sigma2    Cholesky factorization of the covariance matrix for pdf 2
  */
  cholesky_factor_cov[M] L_Sigma1 ;
  cholesky_factor_cov[M] L_Sigma2 ;

  L_Sigma1 <- diag_pre_multiply( tau1, L_Omega1 ) ;
  L_Sigma2 <- diag_pre_multiply( tau2, L_Omega2 ) ;

}
model {
  // prior pdfs ----------

  theta ~ beta( priorParams[1], priorParams[1] ) ;

  mu1 ~ normal( 0, priorParams[2] ) ;
  mu2 ~ normal( 0, priorParams[2] ) ;

  tau1 ~ cauchy( 0, priorParams[3] ) ;
  tau2 ~ cauchy( 0, priorParams[3] ) ;

  // implies L_Omega1 * L_Omega1' ~ lkj_corr(*)
  L_Omega1 ~ lkj_corr_cholesky( priorParams[4] ) ;
  // implies L_Omega2 * L_Omega2' ~ lkj_corr(*)
  L_Omega2 ~ lkj_corr_cholesky( priorParams[4] ) ;

  // likelihood
  for (n in 1:N) {
    increment_log_prob(log_mix( theta,
      multi_normal_cholesky_log( Z[n]', mu1, L_Sigma1 ),
      multi_normal_cholesky_log( Z[n]', mu2, L_Sigma2 )));
  }

}
generated quantities {
  real log_lik ;               // logarithm of the likelihood

  log_lik <- 0. ;
  for (n in 1:N) {
    log_lik <- log_lik + log_mix(theta,
      multi_normal_cholesky_log( Z[n]', mu1, L_Sigma1 ),
      multi_normal_cholesky_log( Z[n]', mu2, L_Sigma2 )) ;
  }
}
