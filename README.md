# MaternSpatialModel-Simulations
R Scripts for the simulation examples included in original research manuscript by Saqlain, Rönnegård, Alam, and Skarin

Scripts to perform the simulations in Scenarios 1-5 in the paper are given here.

## R packages used:

- **Matrix**
- **lgcp** (for circulant() function)
- **spaMM**

## Scripts with functions:

- **Functions_CreateZ.R**: Source file where the functions to create Z with 3 nearest corners are
- **Functions_CreateZ_4corners.R**: Source file where the functions to create Z with 4 nearest corners are
- **Functions_CreateQ.R**: Source file where the functions to create Q are
- **Functions_Estimation.R**: Source file where the functions to estimate the parameters are
- **Functions_Estimation_Binomial.R** Source file where the functions to estimate the parameters are for binomial response
- **Functions_LogL.R**: Source file where the log likelihood function for normal response is
- **Functions_LogL_Binomial.R**: Source file where the log likelihood function for binomial response is
- **Check_Q_and_Matern_functions_V4.R**: Source file where the functions to create regular lattice is
- **HL_Correction.R**: HL Correction file

## Scripts to run simulations:
- **Simulations_Scenarios_1-3.R**: Scenarios 1-3  
- **Simulations_Scenarios_4.R**: Scenario 4
- **Simulations_Scenarios_4_Binom.R**: Scenario 4 for binomial response
- **Simulations_Scenarios_5.R**: Scenario 5

## Directions:
- Copy the source function scripts in the working directory and set working directory in the simulations scripts.

**Simulations 1-3:**
- Set number_of_iterations = desired number of iterations (200 to follow the paper).
- For scenario 1 set: xshift = 0, yshift = 0, removepoints_number = 0
- For scenario 2 set: xshift = 0.2, yshift = 0, removepoints_number = 20
- For scenario 3 set: xshift = 0.5, yshift = 0.5, removepoints_number = 20
- Lattice parameters for the observation locations may be changed in line 35 to 39
- Matern parameters range and alpha may be changed in lines 82 and 83
- Means and standard deviation of the parameters beta, tau, and phi from all the simulations will be the output
- For Scenario 3 using the Z construction with 4 parameters, uncomment line 20 and comment out line 19.

**Simulations 4 (normal response):**
- Set number_of_iterations = desired number of iterations (200 to follow the paper).
- Means and standard deviation of the parameters beta, tau, and phi from all the simulations will be the output

**Simulations 4 (binomial response):**
- Set number_of_iterations = desired number of iterations (200 to follow the paper).
- Lattice spacing may be changed in line 43
- Matern parameters range and alpha may be changed in lines 54 and 55
- Means and standard deviation of the parameters beta and tau from all the simulations will be the output for both HL method and spaMM method

**Simulations 5:**
- Set number_of_iterations = desired number of iterations (200 to follow the paper).
- For 5(a) keep line 40, N.obs = 1000 fixed. Run different simulations by changing reps in line 29 (number of observations in each location)
- For 5(b) change line 40, N.obs = 500, 1000, 1500, and 2000 in the different simulations. Keep reps = 3 fixed in line 29
- Means and standard deviation of the parameters beta, tau, phi and time taken from all the simulations will be the output for both HL method and spaMM method

## Contact: 
- Murshid Saqlain (msq@du.se) for any questions related to the code.

## References:
- Lindgren, F., Rue, H. and Lindström, J., 2011. An explicit link between Gaussian fields and Gaussian Markov random fields: the stochastic partial differential equation approach. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 73(4), pp.423-498.
