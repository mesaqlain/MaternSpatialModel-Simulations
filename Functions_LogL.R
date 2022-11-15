#' Log Likelihood Function
#'
#' This function calculates the value of log likelihood at given parameter values.
#' @param parms_to_estimate s0 (mar variance) and phi0 (dispersion parameter).
#' @param re random effects
#' @param Q Precision Matrix for the mesh
#' @param T.mat Block diagonal matrix with X, Z, 0, and I.
#' @param mu mean of the distribution.
#' @param eta vector of fixed effects and random effects.
#' @param family family of response variable
#' @return Log likelihood.
#' @examples
#' create.coords_mat(c(1,1,1,2,2,2),c(1,2,3,1,2,3))
logL <- function(parms_to_estimate, re, det.Q, Q, T.mat, mu, eta, family, y.response){
  n <- length(eta) #No. of observations
  ys<-y.response
  s0 <- parms_to_estimate[1]
  phi0 <- parms_to_estimate[2]
  wt <- (1/phi0)*as.numeric((family$mu.eta(eta))^2/family$variance(mu))
  Sigma_a <- bdiag(diag(wt),Q/s0)  ## DO this outside. Multiply with wt vector diag(1/sqrt(s0))
  dim(Sigma_a)
  
  CSI <- chol(Sigma_a)

  T1<- CSI%*%T.mat

  dim.Q <- dim(Q)[1]
  
  det.Q_new <-  dim.Q*log(1/s0) ## Decrease computation time?? 
  det.TT <- (determinant(crossprod(T1),logarithm=TRUE)) ##Try crossprod()function here.
  l <- sum(dnorm(ys,mean=eta,sd=sqrt(phi0),log=TRUE)) -
    (0.5)*t(as.vector(re))%*%(Q/s0)%*%as.vector(re) +
    0.5*(det.Q_new) -
    0.5*(det.TT$modulus[1])
  return(as.numeric(-1*l))
}
