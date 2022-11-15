# Version 21st of September 2020 from MSQ
#########################################
# This script contains the function to estimate fixed effects, random effects, marginal variance, and
# dispersion parameter. It also has the function to optimize the log likelihood.
# Version 3 Changes:
#     Create mesh is outside the function. This is because
#     for our three stages of testing, we want the same mesh, not a varying one.
# Version 5 Changes:
#     estimate.fun()  has a new argument, lattice_spacing.
#     Line 53 changed accordingly to incorporate the new argument.
#     Commented out the verbose messages - turn back on to debug if needed.
# Version 6 Changes:
#     Added N as an argumente to Q. Added range and alpha as arguments to Z.
#     Added uid as an argument to the estimation process which is used in Q.
#     Added a line to calculate det(Q) that is used as an argument in the logL function.
#     
#########################################


#' Create Coordinates Matrix Function
#'
#' This function returns a coordinates matrix from x and y vectors.
#' @param y.response Response variable.
#' @param x.coord Vector of x-coordinates of the observation locations.
#' @param y.coord Vector of y-coordinates of the observation locations.
#' @param range Range parameter.
#' @param alpha Alpha parameter, Alpha = nu + 1.
#' @param family Family of response variable.
#' @param init Vector of initial values of fixed effects, dispersion parameter, and marginal variance).
#' @param mesh Coordinates matrix of the mesh, created by create.regular.mesh() function.
#' @param uid Index of unique location id's for constructing the Z matrix.
#' @return Fixed effects, random effects, marginal variance, and dispersion parameter, upon convergence.
#' @examples
#' create.coords_mat(c(1,1,1,2,2,2),c(1,2,3,1,2,3))
estimate.fun <- function(y.response, x.coord, y.coord, range, alpha, family, init, mesh,
                         lattice_spacing, uid){
  # Create Model Frame
  model_frame <- model.frame(y.response ~ 1)
  
  # Design Matrices X
  X <- model.matrix(model_frame)

  # Create vector for response, y 
  y.resp <- as.numeric(model.response(model_frame))
  N <- length(y.resp) #No. of observations

  # Create Z matrix
  Z <- create.Z(x.coord, y.coord, mesh, range, alpha, uid)
  
  #Create the T matrix. (Lee, Neldar, Pawitan 2006, p154)
  T.mat <-rbind(cbind(X,Z),cbind(Matrix(0,nrow=ncol(Z),ncol=ncol(X)),diag(ncol(Z))))
  
  ######## Initial Values ######## 
  b0 <- init[1]
  phi0 <- 1
  s0 <- init[3]
  r0 <- c(rep(0,ncol(Z)))
  my.fam <- family
  
  # Create Q matrix for mesh
  #message(sprintf("Construct Q matrix."))
  #start_time <- proc.time()
  Q <- create.Q_matrix_augmented(create.Q_matrix(alpha, range, mesh, lattice_spacing), uid)
  #print(proc.time()-start_time)
  
  #message(sprintf("Calculaing determinant of Q"))
  #start_time <- proc.time()
  det.Q <- (determinant(Q,logarithm=TRUE)) ## Decrease computation time?? 
  #print(proc.time()-start_time)
  
  #message(sprintf("Calculating eta."))
  eta <-(X%*%b0) + (Z%*%r0)
  eta <- eta[,1]
  #dim(eta) #2237x1 as expected
  
  HL.correction <- 0
  mu <- family$linkinv(eta) #g(mu) = eta, so mu = g.inv(eta)
  
  # Column header for the parameters as they are being estimated
  message("Iter\tBeta\tTau\tPhi")
  
  for(i in 1:200){
    #message(sprintf("Iteration number: %i", i))
    
    #message(sprintf("Calculating working response."))
    #Working Reponse, pg 44 in text Lee, Neldar, Pawitan 2006
    z <- eta + (y.resp - mu) / family$mu.eta(eta) 
    z <- z - HL.correction
    #dim(z) #2237x1
    # Sigma_a
    
    #message(sprintf("Calculating weights."))
    wt <- (1/phi0)*as.numeric((family$mu.eta(eta))^2/family$variance(mu))
    
    #message(sprintf("Constructing Sigma_a matrix."))
    Si <- bdiag(diag(wt),Q/s0)
    #dim(Si) # (5273 x 5273) 
    #Makes sense, 2237x2237 for the sigma on top left block (covariance matrix for the observations)
    #and 3036x3036 for D in the bottom right block (covariance matrix for the random effects)
    
    #Construct T (pg 154 in text Lee, Neldar, Pawitan 2006)
    #dim(T) 5273 x 3046 , as we expected. 
    #Because on the left side we have a vector of y and psi(m) that is 5273x1.
    #y is 2237x1 (observations) and psi(m) is 3036x1 (0s for random effects)
    #On the right side, we have T%*%delta. delta is a vector of beta and random effects.
    #beta is 10x1 and random effects are 3036x1, so delta is 3046x1.
    #T should therefore be #5273 x 3046 , as we expected.
    
    #message(sprintf("Taking Cholesky of Sigma_a."))
    #start_time <- proc.time()
    CSI <- chol(Si) #Take Cholesky of Sigma Inverse 
    #print(proc.time()-start_time)
    #dim(CSI) #5273x5273
    
    #message(sprintf("Taking Cholesky of T."))
    T1 <- CSI%*%T.mat #T* 
    #dim(T) # 5273x3046
    
    #message(sprintf("Constructing vector for responses, y_a."))
    ys <- CSI%*%c(as.numeric(z),rep(0,ncol(Z)))
    
    #message(sprintf("Performing QR decomposition of the cholesky of T."))
    QR1 <- qr(T1) #QR decomposition
    
    #message(sprintf("Finding solutions from the QR decomposition."))
    b1 <- qr.coef(QR1,ys)
    
    #message(sprintf("Creating vector for fixed effects."))
    fe <- b1[1:ncol(X)] #Fixed effects
    
    #message(sprintf("Creating vector for random effects."))
    re <- b1[-(1:ncol(X))] #Random effects
    
    #message(sprintf("Updating eta."))
    eta1 <- (X%*%fe) + (Z%*%re)
    
    #message(sprintf("Checking for convergence."))
    if(sum((eta-eta1)^2)<1e-6*sum(eta1^2)) break #Convergence criterion
    
    #s0_opt <-optimize(logL, interval=c(0.00001,100),maximum=TRUE)
    #message(sprintf("Optimizing log likelihood."))
    opt_values <-nlminb(start = c(s0), #SHOULD I CHANGE THIS TO BE AN ARGUMENT??
                    logL_binom, 
                    lower=0.001,
                    upper=10,
                    control = list(rel.tol=1e-6),
                    phi = phi0,
                    re = re,
                    det.Q = det.Q,
                    Q = Q,
                    T.mat = T.mat,
                    mu = mu,
                    eta = eta,
                    family = family,
                    y.response = y.response)
    
    #message(sprintf("Optimized."))
    #s0 <- s0_opt$maximum
    #s0_obj <- s0_opt$objective
    #message(sprintf("Updating marginal variance, fixed effects, and random effects."))
    s0 <- opt_values$par[1]
    #phi0 <- opt_values$par[2]
    b0 <- fe
    r0 <- re
    
    eta <- eta1
    eta <- eta[,1]
    
    mu <- family$linkinv(eta) 
    HL.correction <- HL11(fv = mu, Z = Z, family = family, phi = phi0, s = s0)
    
    
    message(sprintf("%i\t%.4f\t%.4f", i, b0, s0))
  }
  
  #print(i)
  message(sprintf("Converged in %i iterations", i))
  return(list(beta = b0, phi = phi0, mar.var = s0, ran.eff = re))
}
