# Simulation Scenarios 1-3
#########################################
# This script contains the codes to run simulations 1-3 in the paper 
# Spatial Matérn Models on a Lattice using the Hierarchical Likelihood Framework
# A regular lattice is created which corresponds to the observation locations; this can be
# modified to be irregular and(Scenario 2) to be shifted (Scenarios 2 and 3)
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

## First, create a regular lattice that forms our observation locations

## Parameters for the lattice
xstart <- 1
ystart <- 2
xsize <- 10
ysize <- 10
latspace <- 1 # space between each lattice poitns
# For scenario 1, xshift and yshift are 0
# For scenario 2, xshift is 0.2 units and yshift is 0
# For scenario 3, xshift and y shift are 0.5 each
xshift <- 0 #How many units to shift in the x direction )
yshift <- 0 #How many units to shift in the y direction
# For scenario 1, removepoints_number is 0, for scenarios 2 and 3 it is 20.
removepoints_number <- 0 #How many points to remove randomly.

## Create our observation matrix.
lattice_obs <- create.regular_lattice(xstart, ystart, xsize, ysize, latspace)

## Randomly remove some points
set.seed(123) # this seed used for simulations in the paper
if (removepoints_number == 0){
  lattice_obs <- lattice_obs
} else {
  lattice_obs <- lattice_obs[-sample(1:nrow(lattice_obs), removepoints_number), ]
  
}

## Plot to check how the observations look like.
# plot(lattice_obs)

## Parameters irregular_obs to create a mesh around our observation location:
x.coords <- lattice_obs[,1] + xshift # x coordinates
y.coords <- lattice_obs[,2] + yshift # y coordinates
N.obs <- dim(lattice_obs)[1] #Total number of observations

# Parameters for creating the mesh
lattice_spacing <- 1 # space between each lattice points of the mesh grid
extension_points <- 8 # this is to eliminate boundary effects

## We use the same mesh as in Regular Lattice case.
mesh <- create.regular.mesh(lattice_obs[,1], lattice_obs[,2], lattice_spacing, extension_points)

## Check how the mesh and observations look like
# xlim= c(-10,20) ylim = c(-10,20) these are the last 4 numbers in the argument,
# change accordingly for boundaries for other meshes
plot_mesh_and_obs(x.coords, y.coords, mesh, -10, 20, -10, 20)


## Matern Parameters for simulating the gaussian responses
range <- 2 # in the paper this is varied between 0.5, 1, and 2
alpha <- 2 # works for alpha = 2 and 3

# Calculate Matern Covariance
M <- get.matern_cov(d=as.matrix(dist(cbind(x.coords,y.coords))), range = range, alpha = alpha)
#View(M)



# So now we have simulated observation locations, observations values based on Matern covariance, and an example mesh.
# The observation locations coincide directly on the mesh. Next step is the estimations.


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

# Create empty vectors to store values in
alpha_vector<- c()
range_vector <- c()
beta_vector <- c()
s0_vector <- c()
phi_vector <- c()
iter_vector <- c()
for (l in (1:length(iter))){
  message(sprintf("Simulation number: %i", l))
  v <- rnorm(N.obs)
  v2 <- rep(v,rep(repeated_obs,length(v)))
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
                      phi_vector)


# This gives the estimate of fixed effect, random effect and dispersion parameter 
# for range in line 82. For different range, change line 82 and rerun code from that line.
colMeans(results_data)
# Standard deviation of results
sapply(results_data, sd, na.rm = TRUE)

# uncomment to export to csv
#write.csv(results_data, 
#         "results_data.csv", row.names = F, quote = F)
