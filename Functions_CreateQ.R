#Version 21st of September 2020 from MSQ
#########################################
# This script contains the func
# Packages needed
# library(lgcp) #For circulant() function
# library(Matrix) #For sparse matrix
# Version 4 Changes: 
#     The kappa function, my.inv.range(), is now a function of lattice_spacing, as we want to scale the 
#     as we want to rescale the range.
#     create.Q_matrix() is also a function of lattice_spacing.
# VERSION 5 CHANGES:
#     Added a diagonal identity matrix in the bottom right of Q.
#     Line 438 has been updated to create the bottom right matrix properly.

##################################################################################
## FOLLOWING CODES FROM LINDGREN
##################################################################################

#' Kappa Function
#'
#' This function returns the value of the scaling parameter Kappa in the SPDE representation.
#' @param range The parameter 'range' in the SPDE.
#' @param alpha The parameter 'alpha' in the SPDE (1, 2, or 3).
#' @return The value of Kappa.
#' @examples
#' my.inv.range(1,1)
#' my.inv.range(10,2)
my.inv.range <- function(range, alpha, lattice_spacing) {
  nu <- my.nu(alpha)
  return(2 * sqrt(2 * nu) / (range/lattice_spacing))
}


#' Kappa Function NOT SCALED
#'
#' This function returns the value of the scaling parameter Kappa in the SPDE representation.
#' @param range The parameter 'range' in the SPDE.
#' @param alpha The parameter 'alpha' in the SPDE (1, 2, or 3).
#' @return The value of Kappa.
#' @examples
#' my.inv.range(1,1)
#' my.inv.range(10,2)
my.inv.range.unscaled <- function(range, alpha) {
  nu <- my.nu(alpha)
  return(2 * sqrt(2 * nu) / (range))
}




#' Range Function
#'
#' This function returns the value of the range parameter in the SPDE representation.
#' @param kappa The parameter 'kappa' in the Matern Function.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3.
#' @return The value of range, which is denoted by rho in Lindgren and Rue 2011.
#' @examples
#' my.range(10,1)
#' my.range(1,2)
my.range <- function(kappa, alpha) {
  nu <- my.nu(alpha)
  return(2 * sqrt(2 * nu) / kappa)
}



#' Nu Function
#'
#' This function returns the value of the range parameter in the SPDE representation.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3.
#' @return The parameter nu in the matern function.
#' @examples
#' my.nu(2)
my.nu <- function(alpha) {
  nu <- alpha - 1.0
  return(nu)
}



#' Coefficients of the Precision Matrix Function 
#'
#' This function returns a matrix with the coefficients for the neighbors in the Precision Matrix.
#' Element [1,1] is the weight of the observation location, elements [1,2] and [2,1] are weights for the 
#' first neighbors, and so on. These weights are taken from Section 2.2.1 of Lindgren and Rue 2011 for a regular
#' two-dimensional lattice.
#' @param kappa The parameter 'kappa' in the Matern function.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3.
#' @param nh The parameter for the defining the dimensions of the matrix, usually for however many neighbors. Defaults to 4.
#' @return A matrix of coefficients for the neighbors in a Precision Matrix relative to an observation location.
#' @examples
#' my.make.coofs(10,2,4)
my.make.coofs <- function(kappa, alpha, nh = 4) {
  a <- 4 + Re(kappa^2)
  b <- Im(kappa^2)
  
  if (!is.complex(kappa)) {
    coofs <- matrix(0, nh + 1, nh + 1)
    
    if (alpha == 1.0) {
      coofs[2, 1] <- -1.0
      coofs[1, 1:2] <- c(a, -1.0)
    }
    else if (alpha == 2.0) {
      coofs[3, 1] <- 1
      coofs[2, 1:2] <- c(-2 * a, 2)
      coofs[1, 1:3] <- c(4 + a^2, -2 * a, 1)
    }
    else if (alpha == 3.0) {
      coofs[4, 1] <- -1
      coofs[3, 1:2] <- c(3 * a, -3)
      coofs[2, 1:3] <- c(-3 * (a^2 + 3), 6 * a, -3)
      coofs[1, 1:4] <- c(a * (a^2 + 12), -3 * (a^2 + 3), 3 * a, -1)
    }
    else if (alpha == 4.0) {
      b1 <- (a^2 + 6)^2 + 12 * a^2
      b2 <- -4 * a * (a^2 + 9)
      b3 <- 2 * (3 * a^2 + 8)
      b4 <- -4 * a
      b5 <- -12 * a
      b6 <- 12 * (a^2 + 2)
      
      coofs[5, 1] <- 1
      coofs[4, 1:2] <- c(b4, 4)
      coofs[3, 1:3] <- c(b3, b5, 6)
      coofs[2, 1:4] <- c(b2, b6, b5, 4)
      coofs[1, 1:5] <- c(b1, b2, b3, b4, 1)
    }
    else {
      stop(paste("unknown kappa ", kappa))
    }
  }
  else {
    coofs <- matrix(0, nh + 1, nh + 1)
    if (alpha == 2.0) {
      coofs[3, 1] <- 1
      coofs[2, 1:2] <- c(-2 * a, 2)
      coofs[1, 1:3] <- c(4 + a^2, -2 * a, 1)
      
      coofs[1, 1] <- coofs[1, 1] + b^2
      return(coofs)
    }
    else if (alpha == 3.0) {
      stop(paste("bad alpha ", alpha))
      
      coofs[4, 1] <- -1
      coofs[3, 1:2] <- c(3 * a, -3)
      coofs[2, 1:3] <- c(-3 * (a^2 + 3), 6 * a, -3)
      coofs[1, 1:4] <- c(a * (a^2 + 12), -3 * (a^2 + 3), 3 * a, -1)
      
      coofs2 <- matrix(0, nh + 1, nh + 1)
      coofs2[2, 1] <- -1.0
      coofs2[1, 1:2] <- c(a, -1.0)
      
      coofs <- coofs + b^2 * coofs2
      
      return(coofs)
    }
    else if (alpha == 4.0) {
      b1 <- (a^2 + 6)^2 + 12 * a^2
      b2 <- -4 * a * (a^2 + 9)
      b3 <- 2 * (3 * a^2 + 8)
      b4 <- -4 * a
      b5 <- -12 * a
      b6 <- 12 * (a^2 + 2)
      
      coofs[5, 1] <- 1
      coofs[4, 1:2] <- c(b4, 4)
      coofs[3, 1:3] <- c(b3, b5, 6)
      coofs[2, 1:4] <- c(b2, b6, b5, 4)
      coofs[1, 1:5] <- c(b1, b2, b3, b4, 1)
      
      coofs2 <- matrix(0, nh + 1, nh + 1)
      coofs2[3, 1] <- 1
      coofs2[2, 1:2] <- c(-2 * a, 2)
      coofs2[1, 1:3] <- c(4 + a^2, -2 * a, 1)
      
      coofs <- coofs + b^2 * coofs2
      
      return(coofs)
    }
    else {
      stop(paste("unknown complex kappa ", kappa))
    }
  }
  
  return(coofs)
}

##################################################################################
## END CODES FROM LINDGREN
##################################################################################



##################################################################################
## START OF CODES BY MSQ
##################################################################################



#' Get X Lattice Size Function
#'
#' This function returns the number of points on the x-axis of an existing mesh. It is used in a number of 
#' other functions.
#' @param mesh Mesh created by the create.regular.mesh() function which is a matrix of coordinates of the mesh.
#' @return Number of points on the x-axis from a mesh
#' @examples
#' get.x_lattice_size(example_mesh)
get.x_lattice_size <- function(mesh) {
  return( length(unique(mesh[,1])) )
}



#' Get Y Lattice Size Function
#'
#' This function returns the number of points on the y-axis of an existing mesh. It is used in a number of 
#' other functions.
#' @param mesh Mesh created by the create.regular.mesh() function which is a matrix of coordinates of the mesh.
#' @return Number of points on the y-axis from a mesh
#' @examples
#' get.y_lattice_size(example_mesh)
get.y_lattice_size <- function(mesh) {
  return( length(unique(mesh[,2])) )
}



#' Get Neighbours Function
#'
#' This function a vector of 0s and 1s that denote which neighbors are considered. 
#' For alpha=1, only the first neighbor is considered.
#' For alpha=2, the first and second and diagonal neighbors are considered.
#' For alpha=3, the first, second, third, diagonal, and off-diagonal neighbors are considered.
#' This function is used when creating blocks for the precision matrix.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3.
#' @return A vector of 0s and 1s that denote which neighbors are considered.
#' @examples
#' get.neighbours(1)
#' get.neighbours(2)
get.neighbours <- function(alpha){
  values <- c(rep(0,6)) 
  #Point itself, first neighbour, second neighbour, third neighbour, diagonal neighbour, off-diagonal neighbour
  if (alpha == 1.0){
    values[c(1,2)] <- c(1,1)
  }
  else if (alpha == 2.0){
    values[c(1,2,3,5)] <- c(1,1,1,1)
  }
  else if (alpha == 3.0){
    values[c(1,2,3,4,5,6)] <- c(1,1,1,1,1,1)
  }
  else{
    stop(paste("Unknown alpha: ", alpha))
  }
  return(values)
}


#' Create First Block Function
#'
#' This function creates the base vector of the circulant matrix that goes on the diagonal of the Precision Matrix.
#' The size of the block depends on the number of points on the x-axis of the lattice.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3.
#' @param mesh Mesh created by the create.regular.mesh() function which is a matrix of coordinates of the mesh.
#' @param coefs Coefficients of neighbours on a precision matrix, created by the my.make.coofs() function.
#' @return The base vector for the first block on the precision matrix.
#' @examples
#' create.base_first_block(1, example_mesh, example_coefs)
create.base_first_block <- function(alpha, mesh, coefs) {
  neighbors <- get.neighbours(alpha)
  x_lattice_size <- get.x_lattice_size(mesh)
  base_vector <- rep(0,x_lattice_size)
  base_vector[c(1, 2, 3, 4, x_lattice_size-2, x_lattice_size-1, x_lattice_size)] <- c(neighbors[1]*coefs[1,1], #Point itself
                                                                                      neighbors[2]*coefs[1,2], #First neighbour  
                                                                                      neighbors[3]*coefs[1,3], #Second neighbour
                                                                                      neighbors[4]*coefs[1,4], #Third neighbour
                                                                                      neighbors[4]*coefs[1,4], #Third neighbour
                                                                                      neighbors[3]*coefs[1,3], #Second neighbour
                                                                                      neighbors[2]*coefs[1,2]) #First neighbour   
  return(base_vector)
}



#' Create Second Block Function
#'
#' This function creates the base vector of the circulant matrix that goes on the first row 
#' and second column of the block circulant Precision Matrix. 
#' The size of the block depends on the number of points on the x-axis of the lattice.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3.
#' @param mesh Mesh created by the create.regular.mesh() function which is a matrix of coordinates of the mesh.
#' @param coefs Coefficients of neighbours on a precision matrix, created by the my.make.coofs() function.
#' @return The base vector for the second block on the precision matrix.
#' @examples
#' create.base_second_block(1, example_mesh, example_coefs)
create.base_second_block <- function(alpha, mesh, coefs) {
  neighbors <- get.neighbours(alpha)
  x_lattice_size <- get.x_lattice_size(mesh)
  base_vector <- rep(0,x_lattice_size)
  base_vector[c(1, 2, 3, x_lattice_size-1,  x_lattice_size)] <- c(neighbors[2]*coefs[1,2], #First neighbour
                                                                  neighbors[5]*coefs[2,2], #Diagonal neighbour
                                                                  neighbors[6]*coefs[2,3], #Off-Diagonal neighbour#Diagonal neighbour
                                                                  neighbors[6]*coefs[2,3], #Off-Diagonal neighbour#Diagonal neighbour
                                                                  neighbors[5]*coefs[2,2]) #DIagonal neighbour
  return(base_vector)
}



#' Create Third Block Function
#'
#' This function creates the base vector of the circulant matrix that goes on the first row 
#' and third column of the block circulant Precision Matrix. 
#' The size of the block depends on the number of points on the x-axis of the lattice.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3.
#' @param mesh Mesh created by the create.regular.mesh() function which is a matrix of coordinates of the mesh.
#' @param coefs Coefficients of neighbours on a precision matrix, created by the my.make.coofs() function.
#' @return The base vector for the third block on the precision matrix.
#' @examples
#' create.base_third_block(1, example_mesh, example_coefs)
create.base_third_block <- function(alpha, mesh, coefs) {
  neighbors <- get.neighbours(alpha)
  x_lattice_size <- get.x_lattice_size(mesh)
  base_vector <- rep(0,x_lattice_size)
  base_vector[c(1, 2, x_lattice_size)] <- c(neighbors[3]*coefs[1,3], #Second neighbour
                                            neighbors[6]*coefs[2,3], #Off-Diagonal neighbour
                                            neighbors[6]*coefs[2,3]) #Off-DIagonal neighbour
  return(base_vector)
}



#' Create Fourth Block Function
#'
#' This function creates the base vector of the circulant matrix that goes on the first row 
#' and fourth column of the block circulant Precision Matrix. 
#' The size of the block depends on the number of points on the x-axis of the lattice.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3.
#' @param mesh Mesh created by the create.regular.mesh() function which is a matrix of coordinates of the mesh.
#' @param coefs Coefficients of neighbours on a precision matrix, created by the my.make.coofs() function.
#' @return The base vector for the fourth block on the precision matrix.
#' @examples
#' create.base_fourth_block(1, example_mesh, example_coefs)
create.base_fourth_block <- function(alpha, mesh, coefs) {
  neighbors <- get.neighbours(alpha)
  x_lattice_size <- get.x_lattice_size(mesh)
  base_vector <- rep(0,x_lattice_size)
  base_vector[1] <- neighbors[4]*coefs[1,4]
  return(base_vector)
}




#' Create Zero Block Function
#'
#' This function creates a vector of 0s that will populate most of the Precision Matrix. 
#' These blocks of 0s will start from the 5th column, 
#' until the (n-4)th column where n is the number of points on the y-axis.
#' @param mesh Mesh created by the create.regular.mesh() function which is a matrix of coordinates of the mesh.
#' @return A vector of 0s.
#' @examples
#' create.base_fourth_block(1, example_mesh, example_coefs)
create.base_zero_block <- function(mesh) {
  x_lattice_size <- get.x_lattice_size(mesh)
  base_vector <- rep(0,x_lattice_size)
  return(base_vector)
}



#' Find Number of Zero Blocks Function
#'
#' Most of the blocks in the Precision Matrix is blocks of 0. 
#' This function finds how many of these blocks of 0s are needed on each row of the Block Circulant Precision Matrix.
#' These blocks of 0s will start from the 4th column, 
#' until the (n-4)th column where n is the number of points on the y-axis.
#' @param mesh Mesh created by the create.regular.mesh() function which is a matrix of coordinates of the mesh.
#' @return Number of 0 blocks needed.
#' @examples
#' create.base_fourth_block(1, example_mesh, example_coefs)
find.number_zero_blocks <- function(mesh){
  y_lattice_size <- get.y_lattice_size(mesh)
  if(y_lattice_size - 7 < 0){  
    return(0)
  }
  else{
    return(y_lattice_size - 7)
  }
}



#' Precision Matrix Function
#'
#' Function to create a sparse Precision Matrix for the regular lattice mesh overlayed on our observations
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3.
#' @param range The parameter 'range' in the SPDE.
#' @param mesh Mesh created by the create.regular.mesh() function which is a matrix of coordinates of the mesh.
#' @return A sparse precision matrix of the regular lattice mesh.
#' @examples
#' create.Q_matrix(1, 10, example_mesh)
create.Q_matrix <- function(alpha, range, mesh, lattice_spacing){
  kappa <- my.inv.range(range,alpha, lattice_spacing)
  #kappa <- my.inv.range.unscaled(range,alpha)
  coefs <- my.make.coofs(kappa, alpha, 4)
  base_first_block <- create.base_first_block(alpha, mesh, coefs)
  base_second_block <- create.base_second_block(alpha, mesh, coefs)
  base_third_block <- create.base_third_block(alpha, mesh, coefs)
  base_fourth_block <- create.base_fourth_block(alpha, mesh, coefs)
  base_zero_block <- create.base_zero_block(mesh)
  y_lattice_size <- get.y_lattice_size(mesh)
  base_vector <- c(base_first_block, base_second_block, base_third_block, base_fourth_block,
                   rep(base_zero_block,y_lattice_size-7),
                   base_fourth_block,base_third_block, base_second_block)
  #For the following matrix, the columns must be the base of each of the block matrices
  matrix_of_bases <- matrix(base_vector, ncol = y_lattice_size)
  mar.var <- my.marginal.variance(range, alpha, lattice_spacing)
  #mar.var <- 1
  Q_matrix <- mar.var*Matrix(circulant(matrix_of_bases),sparse = TRUE)
  return(Q_matrix)
}


#' Augmented Precision Matrix Function
#'
#' Function to create a sparse Precision Matrix for the regular lattice mesh overlayed on our observations
#' This augmented version has an identity matrix on the bottom right with dimensions the size of the number of observations.
#' @param Q Precision Matrix created with create.Q_matrix.
#' @param N The number of observations.
#' @return A sparse precision matrix of the regular lattice mesh.
#' @examples
#' create.Q_matrix(1, 10, example_mesh)
create.Q_matrix_augmented <- function(Q, id){
  unique_ids <- unique(id)
  bot_right_matrix <- sparseMatrix(i = rep(1:length(unique_ids), each=1), 
                                   j = seq(1,length(unique_ids)), 
                                   x = 1)
  Q.full <- rbind(cbind(Q, Matrix(0, nrow = nrow(Q), ncol = ncol(bot_right_matrix))),
                  cbind(Matrix(0, nrow = nrow(bot_right_matrix), ncol = ncol(Q)), bot_right_matrix))
  return(Q.full)
}
#create.Q_matrix_augmented <- function(Q, N){
#  Q.full <- rbind(cbind(Q, Matrix(0, nrow = nrow(Q), ncol=N)),
#                  cbind(Matrix(0, nrow = N, ncol = ncol(Q)), diag(N)))
#  return(Q.full)
#}
