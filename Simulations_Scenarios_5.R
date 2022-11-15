# Simulation Scenario 5 
#########################################
# This script contains the codes to run simulation scenario 5 in the paper 
# Spatial Matérn Models on a Lattice using the Hierarchical Likelihood Framework
# In Scenario 5(a) observation locations is fixed, but number of reps/observations in each location is varied
# In Scenario 5(b) observation locations are varied, but number of reps/observations in each location is fixed.
# HL Method is compared to spaMM estimation
#########################################

rm(list=ls()) #Clear environment

#Where the functions are located, uncomment and specify
# setwd("~/")

# Source of functions
source("Functions_CreateQ.R") #Source file where the functions to create Q are
source("Functions_Estimation.R") #Source file where the functions to estimate the parameters are
source("Functions_Q_and_Matern.R") #Source file where the functions to create regular lattice is
source("Functions_LogL.R") #Source file for the Log likelihood function
source("Functions_CreateZ.R") #Source file where the functions to create Z with 3 corners are
#source("Functions_CreateZ_4corners.R") #Source file where the functions to create Z with 4 corners
#For 4 nearest corners in the mesh comment out line 19 and uncomment line 20

# Load necessary packages
library(Matrix)
library(lgcp) #For circulant() function
library(spaMM)

reps <- 3 # for 5(a), varies between 3, 4, 5, and 6. For 5(b) keep fixed

# Load necessary packages
library(Matrix)
library(lgcp) #For circulant() function
library(spaMM)


################################
##### SIMULATE DATA ############
################################

## First, create a regular lattice that forms our observation locations

## First, create random points
set.seed(12345)
#A small simulation example
N.obs <- 500 #No. of observations # For 5(a) keep fixed 1000. For 5(b) vary between 500, 1000, 1500, and 2000

x.coords <- runif(N.obs, 0, 10)
y.coords <- runif(N.obs, 0, 5)
plot(x.coords,y.coords)

# Parameters for creating the mesh
lattice_spacing <- 0.25
extension_points <- 8

## We use the same mesh as in Regular Lattice case.
mesh <- create.regular.mesh(x.coords, y.coords, lattice_spacing, extension_points)

## Check how the mesh and observations look like
plot_mesh_and_obs(x.coords, y.coords, mesh, -5, 15, -5, 10)

## Matern Parameters for simulating the gaussian responses
range <- 0.5 #fixed in scenario 5
alpha <- 2

# Calculate Matern Covariance
M <- get.matern_cov(d=as.matrix(dist(cbind(x.coords,y.coords))), range = range, alpha = alpha)
#View(M)

################################
##### ESTIMATION STEPS #########
################################

x.coordinates <- rep(x.coords,rep(reps,length(x.coords))) #x-coordinate location of response
y.coordinates <- rep(y.coords,rep(reps,length(y.coords))) #y-coordinate location of response
range.values <- range #Range of Matern
alpha.values <- 2 #Alpha, nu+1
family <- gaussian(link=identity) #Response family
initial.values <- c(b0 = 2, phi0 = 1, s0 = 1) #Inital values of fixed effects, dispersion parameter, and marginal variance)
uid <-rep(1:length(x.coords), each = reps) ## unique location IDs

#Fix parameters for CorrHLfit
nu <- 1
rho <- sqrt(8*nu)/range.values

fixed_params <- list(nu = nu)

# Number of simulations.
number_of_iterations <- 2
# how many iterations to perform, for paper we used up to 200
iter <- seq(1,number_of_iterations) 

# Create empty vectors to store values in
alpha_vector<- c()
range_vector <- c()
beta_vector <- c()
s0_vector <- c()
phi_vector <- c()
iter_vector <- c()
time_vector <- c()
for (l in (1:length(iter))){
  message(sprintf("Simulation number: %i", l))
  v <- rnorm(N.obs)
  v2 <- rep(v,rep(reps,length(v)))
  sim.re <- t(chol(M))%*%v #Random Effects
  Y <- initial.values[1] + rep(sim.re,rep(reps,length(sim.re))) + rnorm(N.obs*reps)
  for (i in (1:length(alpha.values))){
    for (j in (1:length(range.values))){
      alpha_vector <- c(alpha_vector, alpha.values[i])
      range_vector <- c(range_vector, range.values[j])
      iter_vector <- c(iter_vector, iter[l])
      start_time <- proc.time()
      est <- estimate.fun(y.response = Y, x.coord = x.coordinates, y.coord = y.coordinates, 
                          range = range.values[j], 
                          alpha = alpha.values[i], 
                          family = family, 
                          init = initial.values, 
                          mesh = mesh,
                          lattice_spacing = lattice_spacing,
                          uid = uid)
      end.time <- proc.time() - start_time
      time_vector <- c(time_vector, end.time[3])
      
      # Fixed Effect
      est.beta <- est$beta
      beta_vector <- c(beta_vector, est.beta)
      # Mar. Variance
      est.s0 <- est$mar.var
      s0_vector <- c(s0_vector, est.s0)
      # Dispersion Parameter
      est.phi <- est$phi
      phi_vector <- c(phi_vector, est.phi)
      # Create Z matrix
      Z.mat <- create.Z(x.coordinates, y.coordinates, mesh, range, alpha, uid)
      # Random Effect
      est.re.ran <- Z.mat%*%(est$ran.eff)
      
      #Plot simulated random  effects against estimated random effects
      # uncomment following lines if you want plots 
      # file_name <- paste("alpha-", alpha.values[i], 
      #                   "_range-", range.values[j], 
      #                   "_iteration-", iter[l],
      #                   ".pdf", sep="")
      # pdf(file_name)
      # plot(v2,est.re.ran) 
      # dev.off()
      
    }
  }
}


results_data <- data.frame(iter_vector,
                           beta_vector,
                           s0_vector,
                           phi_vector,
                           time_vector)


# This gives the estimate of fixed effect, random effect and dispersion parameter 
# for range in line 82. For different range, change line 82 and rerun code from that line.
colMeans(results_data)
# Standard deviation of results
sapply(results_data, sd, na.rm = TRUE)



################################
##### spaMM Estimation##########
################################


# Estimate
alpha_vector<- c()
range_vector <- c()
beta_vector <- c()
s0_vector <- c()
phi_vector <- c()
iter_vector <- c()
time_vector <- c()
for (l in (1:length(iter))){
  message(sprintf("Simulation number: %i", l))
  v <- rnorm(N.obs)
  v2 <- rep(v,rep(reps,length(v)))
  sim.re <- t(chol(M))%*%v #Random Effects
  Y <- initial.values[1] + rep(sim.re,rep(reps,length(sim.re))) + rnorm(N.obs*reps)
  
  dataset <- data.frame(Y, x.coordinates, y.coordinates)
  
  start_time <- proc.time()
  model.spamm <- corrHLfit(Y ~ 1 + Matern(1|x.coordinates + y.coordinates),
                           data = dataset,
                           family = gaussian(),
                           rand.family = gaussian(),
                           #ranFix = fixed_params,
                           #init.corrHLfit = list(rho = 0.5, nu = 1),
                           #lower = list(rho = 0.001, nu = 0.0001),
                           #upper = list(rho = 10, nu = 3),
                           HLmethod="ML")
  end.time <- proc.time() - start_time
  
  time_vector <- c(time_vector, end.time[3])
  iter_vector <- c(iter_vector, iter[l])
  
  #nu <- model.spamm$corrPars[1]$'1'[1] #nu
  #alpha_vector <- c(alpha_vector, nu)
  
  rho <- model.spamm$corrPars[1]$'1'[2] #rho
  range_vector <- c(range_vector, rho)
  
  # Fixed Effect
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

#Uncomment the following lines to save the results in a csv

results_spaMM <- data.frame(iter_vector,
                                                                       beta_vector,
                                                                       s0_vector,
                                                                       phi_vector,
                                                                       as.vector(range_vector),
                                                                       time_vector)
colMeans(results_spaMM)
# Standard deviation of results
sapply(results_spaMM, sd, na.rm = TRUE)
