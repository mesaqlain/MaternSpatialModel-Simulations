# Simulation Scenario 4
#########################################
# This script contains the codes to run simulation scenario 4 in the paper 
# Spatial Matérn Models on a Lattice using the Hierarchical Likelihood Framework
# 100 observation locations are simulated before estimations are performed.
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
lattice_spacing <- 1
extension_points <- 8

## We use the same mesh as in Regular Lattice case.
mesh <- create.regular.mesh(x.coords, y.coords, lattice_spacing, extension_points)
dim(mesh) ##594 x 2
## Check how the mesh and observations look like
plot_mesh_and_obs(x.coords, y.coords, mesh, -10, 20, -10, 15)


## Matern Parameters for simulating the gaussian responses
range <- 2
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
family <- gaussian(link=identity) #Response family
initial.values <- c(b0 = 2, phi0 = 1, s0 = 1) #Inital values of fixed effects, dispersion parameter, and marginal variance)
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
  v <- rnorm(N.obs)
  #v2 <- rep(v,rep(3,length(v)))
  sim.re <- t(chol(M))%*%v #Random Effects
  Y <- initial.values[1] + rep(sim.re,rep(repeated_obs,length(sim.re))) + rnorm(N.obs*repeated_obs)
  for (i in (1:length(alpha.values))){
    for (j in (1:length(range.values))){
      alpha_vector <- c(alpha_vector, alpha.values[i])
      range_vector <- c(range_vector, range.values[j])
      iter_vector <- c(iter_vector, iter[l])
      est <- estimate.fun(y.response = Y, x.coord = x.coordinates, y.coord = y.coordinates, 
                          range = range.values[j], 
                          alpha = alpha.values[i], 
                          family = family, 
                          init = initial.values, 
                          mesh = mesh,
                          lattice_spacing = lattice_spacing,
                          uid = uid)
      # Fixed Effect
      est.beta <- est$beta
      beta_vector <- c(beta_vector, est.beta)
      # Mar. Variance
      est.s0 <- est$mar.var
      s0_vector <- c(s0_vector, est.s0)
      # Dispersion Parameter
      est.phi <- est$phi
      phi_vector <- c(phi_vector, est.phi)

    }
  }
}

#Uncomment the following lines to save the results in a csv

results_data <- data.frame(iter_vector,
                           beta_vector,
                           s0_vector,
                           phi_vector)


# This gives the estimate of fixed effect, random effect and dispersion parameter 
# for range in line 52. For different range, change line 52 and rerun code from that line.
colMeans(results_data)
# Standard deviation of results
sapply(results_data, sd, na.rm = TRUE)

# uncomment to export to csv
#write.csv(results_data, 
#         "results_data.csv", row.names = F, quote = F)
