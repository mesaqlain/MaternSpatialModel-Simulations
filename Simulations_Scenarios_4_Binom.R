# Simulation Scenario 4
#########################################
# This script contains the codes to run simulation scenario 4 in the paper 
# Spatial Matérn Models on a Lattice using the Hierarchical Likelihood Framework
# 100 observation locations are simulated before estimations are performed.
# The response variable is binomial in this case.
# Later part of the code contains same estimations performed using spaMM instead
#########################################

rm(list=ls())

#Where the functions are located
# setwd("~/")

# Source of functions - CHECK VERSIONS
source("Functions_CreateQ.R") #Source file where the functions to create Q are
source("Functions_Estimation_Binomial.R") #Source file where the functions to estimate the parameters are
source("Functions_Q_and_Matern.R") #Source file where the functions to create regular lattice is
source("Functions_LogL_Binomial.R") #Source file for the Log likelihood function
source("Functions_CreateZ.R") #Source file where the functions to create Z with 3 corners are
#source("Functions_CreateZ_4corners.R") #Source file where the functions to create Z with 4 corners
#For 4 nearest corners in the mesh comment out line 19 and uncomment line 20
source("HL_Correction.R") #HL Correction file

# Load necessary packages
library(Matrix)
library(lgcp) #For circulant() function
library(spaMM)

################################
##### SIMULATE DATA ############
################################

## First, create random points
set.seed(12345) # seed used in paper
N.obs <- 100 #No. of observations

x.coords <- runif(N.obs, 0, 10)
y.coords <- runif(N.obs, 0, 5)
plot(x.coords,y.coords)

# Parameters for creating the mesh
lattice_spacing <- 1 # in paper varied between 0.5 and 1
extension_points <- 8

## We use the same mesh as in Regular Lattice case.
mesh <- create.regular.mesh(x.coords, y.coords, lattice_spacing, extension_points)
dim(mesh) ##594 x 2
## Check how the mesh and observations look like
plot_mesh_and_obs(x.coords, y.coords, mesh, -10, 20, -10, 15)


## Matern Parameters for simulating the gaussian responses
range <- 2 # in paper varied between 0.5, 1, and 2
alpha <- 2

# Calculate Matern Covariance
M <- get.matern_cov(d=as.matrix(dist(cbind(x.coords,y.coords))), range = range, alpha = alpha)
#View(M)
class(M)

################################
##### ESTIMATION STEPS #########
################################

repeated_obs <- 3 # repeated observation in each location
x.coordinates <- rep(x.coords,rep(repeated_obs,length(x.coords))) #x-coordinate location of response
y.coordinates <- rep(y.coords,rep(repeated_obs,length(y.coords))) #y-coordinate location of response
range.values <- range #Range of Matern
alpha.values <- alpha #Alpha, nu+1
family <- binomial(link=logit) #Response family
initial.values <- c(b0 = -0.8, phi0 = 1, s0 = 1) #Inital values of fixed effects, dispersion parameter, and marginal variance)
uid <-rep(1:length(x.coords), each = repeated_obs) ## unique location IDs


# Number of simulations.
number_of_iterations <- 2
# how many iterations to perform, for paper we used up to 200
iter <- seq(1,number_of_iterations) 

# Estimate
alpha_vector<- c()
range_vector <- c()
beta_vector <- c()
s0_vector <- c()
phi_vector <- c()
iter_vector <- c()
for (l in (1:length(iter))){
  message(sprintf("Simulation number: %i", l))
  
  v <- rnorm(N.obs) #Some continuous variable
  #v2 <- rep(v,rep(3,length(v)))
  sim.re <- rep((t(chol(M))%*%v), each = repeated_obs) #Random Effects
  lin.pred <- initial.values[1] + sim.re #Linear Predictor, with fixed effect = 2
  pr = 1/(1+exp(-lin.pred))         # pass through an inv-logit function
  Y = rbinom(N.obs*repeated_obs,1,pr)      # bernoulli response variable
  
  
  for (i in (1:length(alpha.values))){
    for (j in (1:length(range.values))){
      alpha_vector <- c(alpha_vector, alpha.values[i])
      range_vector <- c(range_vector, range.values[j])
      iter_vector <- c(iter_vector, iter[l])
      est <- estimate.fun(Y, x.coordinates, y.coordinates, 
                          range.values[j], alpha.values[i], 
                          family, 
                          initial.values, 
                          mesh,
                          lattice_spacing,
                          uid)
      # Fixed Effect
      est.beta <- est$beta
      beta_vector <- c(beta_vector, est.beta)
      # Mar. Variance
      est.s0 <- est$mar.var
      s0_vector <- c(s0_vector, est.s0)
      # Dispersion Parameter
      #est.phi <- est$phi
      #phi_vector <- c(phi_vector, est.phi)
      # Create Z matrix
      Z.mat <- create.Z(x.coordinates, y.coordinates, mesh, range, alpha, uid)
      # Random Effect
      est.re.ran <- Z.mat%*%(est$ran.eff)

      
    }
  }
}

#Uncomment the following lines to save the results in a csv

results_data <- data.frame(iter_vector,
                           beta_vector,
                           s0_vector)


# This gives the estimate of fixed effect, random effect and dispersion parameter 
# for range in line 54. For different range, change line 54 and rerun code from that line.
colMeans(results_data)
# Standard deviation of results
sapply(results_data, sd, na.rm = TRUE)

# uncomment to export to csv
#write.csv(results_data, 
#         "results_data.csv", row.names = F, quote = F)



#################################
##### ESTIMATION STEPS spaMM ####
#################################


#Fix parameters for CorrHLfit
nu <- 1
rho <- sqrt(8*nu)/range.values
fixed_params <- list(nu = nu, rho = rho)


# Estimate
alpha_vector<- c()
range_vector <- c()
beta_vector <- c()
s0_vector <- c()
phi_vector <- c()
iter_vector <- c()
for (l in (1:length(iter))){
  message(sprintf("Simulation number: %i", l))
  
  v <- rnorm(N.obs) #Some continuous variable
  #v2 <- rep(v,rep(3,length(v)))
  sim.re <- rep((t(chol(M))%*%v), each = repeated_obs) #Random Effects
  lin.pred <- initial.values[1] + sim.re #Linear Predictor, with fixed effect = 2
  pr = 1/(1+exp(-lin.pred))         # pass through an inv-logit function
  Y = rbinom(N.obs*repeated_obs,1,pr)      # bernoulli response variable
  
  
  for (i in (1:length(alpha.values))){
    for (j in (1:length(range.values))){
      alpha_vector <- c(alpha_vector, alpha.values[i])
      range_vector <- c(range_vector, range.values[j])
      iter_vector <- c(iter_vector, iter[l])
      
      dataset <- data.frame(Y, x.coordinates, y.coordinates)
      
      model.spamm <- corrHLfit(Y ~ 1 + Matern(1|x.coordinates + y.coordinates),
                               data = dataset,
                               family = binomial(),
                               rand.family = gaussian(),
                               ranFix = fixed_params,
                               #init.corrHLfit = list(rho = 0.5, nu = 1),
                               #lower = list(rho = 0.001, nu = 0.0001),
                               #upper = list(rho = 10, nu = 3),
                               HLmethod="ML")
      
      beta <- model.spamm$fixef #intercept
      beta_vector <- c(beta_vector, beta)
      # Mar. Variance
      tau <- model.spamm$lambda[1] #lambda (tau)
      s0_vector <- c(s0_vector, tau)
      # Dispersion Parameter
      phi <- model.spamm$phi #phi
      phi_vector <- c(phi_vector, phi)
      # Create Z matrix
      # Random Effect
      #est.re.ran <- Z.mat%*%(est$ran.eff)
      
      #file_name <- paste("alpha-", alpha.values[i], 
      #                   "_range-", range.values[j], 
      #                   "_iteration-", iter[l],
      #                   ".pdf", sep="")
      #Plot simulated random  effects against estimated random effects
      #pdf(file_name)
      #plot(v2,est.re.ran) 
      #dev.off()
      ran.eff <- model.spamm$eta
      
    }
  }
}


results_binomial_spaMM <- data.frame(iter_vector,
                                            beta_vector,
                                            s0_vector)
# Results
colMeans(results_binomial_spaMM)
# Standard deviation of results
sapply(results_binomial_spaMM, sd, na.rm = TRUE)
