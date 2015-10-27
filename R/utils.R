
# These functions apply various operations to compositional data.
# They are not exported.

Close <- function( X, kappa=1.0 )
{
  if( is.matrix( X ) ) {
    Y <- X
    for( i in 1:nrow(Y) ) {
      Y[i,] <- kappa * X[i,] / sum( X[i,] )
    }
    return( Y )
  } else if( is.vector( X ) ){
    Y <- kappa * X / sum( X )
    return( Y )
  } else {
    cat( sprintf( "ERROR: Function Close\n" ) )
    cat( sprintf( "Argument X is neither a matrix nor a vector.\n" ) )
    return( 0 )
  }
}

Perturb <- function ( a, X, kappa=1.0 )
{
  if( is.matrix( X ) ) {
    if( length(a) != ncol(X) ) {
      cat( sprintf( "ERROR: Function Perturb\n" ) )
      cat( sprintf( "The length of the perturbing vector = %d\n", length(a) ) )
      cat( sprintf( "The number of columns in the matrix = %d\n", ncol(X) ) )
    }
    Y <- matrix( a, nrow=nrow(X), ncol=ncol(X), byrow=TRUE ) * X
    return( Close( Y, kappa ) )
  } else if( is.vector( X ) ) {
    if( length(a) != length(X) ) {
      cat( sprintf( "ERROR: Function Perturb\n" ) )
      cat( sprintf( "The length of the perturbing vector = %d\n", length(a) ) )
      cat( sprintf( "The length of the vector = %d\n", length(X) ) )
    }
    Y <- a * X
    return( Close( Y, kappa ) )
  } else {
    cat( sprintf( "ERROR: Function Perturb\n" ) )
    cat( sprintf( "Argument X is neither a matrix nor a vector.\n" ) )
    return( 0 )
  }
}


Power <- function( X, alpha, kappa=1.0 )
{
  return( Close( X^alpha, kappa ) )
}


ContainsZeros <- function( X )
{
  return( any( X == 0, na.rm=TRUE ) )
}

# X must be a matrix
CalcIlrCoefs <- function( X, V )
{
  if( ContainsZeros( X ) )
    return( 0 )

  if( !is.matrix( X ) ) {
    # error message goes here
    return(0)
  }

  Z <- log( X ) %*% V

  colnames( Z ) <- colnames( V )

  return( Z )
}

# Y is a matrix containing the ilr coefficients
CalcInvIlr <- function( Y, V, kappa=1.0 )
{
  return( Close( exp( Y %*% t ( V ) ), kappa=kappa ) )
}


# Calculate the psi matrix, using a Gram-Schmidt technique.
# See example 4.4.1, p. 26
# D is the dimension.

CalcPsiMatrix2 <- function( D )
{
  psi <- matrix( 0, ncol=D, nrow=(D-1) )

  for( i in 1:(D-1) )
  {
    for( j in 1:(D-i) )
    {
      tmp <- ( D - i ) * ( D - i + 1 )
      psi[i,j] <- sqrt( 1.0 / tmp )
    }
    psi[i,D-i+1] <- - sqrt( ( D-i ) / (D-i+1) )
  }

  return( psi )

}

CalcCompCenter <- function( X, kappa=1.0 )
{
  if( ContainsZeros( X ) )
    return( 0 )

  h <- exp( colMeans( log( X ), na.rm=TRUE ) )
  names( h ) <- colnames( X )
  return( Close( h, kappa ) )
}


# Definition 5.2, p. 38

CalcVariationMatrix <- function( X )
{
  D <- ncol( X )
  T <- matrix( 0, nrow=D, ncol=D )
  colnames( T ) <- colnames( X )
  rownames( T ) <- colnames( X )
  #    n <- nrow( X )
  for( i in 1:(D-1) )
  {
    for( j in (i+1):D )
    {
      T[i,j] <- var( log( X[,i]/X[,j] ), na.rm=TRUE ) # Function var (in R) computes the variance with 1/(n-1)
      #            T[i,j] <- (n-1)/n * var( log( X[,i]/X[,j] ) ) # An alternative
      T[j,i] <- T[i,j]
    }
  }

  return( T )
}

# argument T is the variation matrix

CalcTotalVariance <- function( T )
{
  return( sum( T ) / ( 2.0 * ncol(T) ) )
}


# Transform the covariance matrix (in ilr coordinates) to the variation matrix
# (in simplex/concentration coordinates)
# from a talk by Karel Hron

InvCovTransform <- function( Sy, Psi ) {

  X <- t(Psi) %*% Sy %*% Psi
  D <- ncol(X)
  diagX <- diag( diag(X), D, D )
  J <- matrix(1, nrow=D, ncol=D)
  A <- J %*% diagX
  variationMatrix <- 0.5 * ( A + t(A) - 2 * X )
  return(variationMatrix)
}


