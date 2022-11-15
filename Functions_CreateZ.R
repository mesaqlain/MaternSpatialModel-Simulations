#Version 13th of August 2020 from MSQ
#########################################
# This script contains the functions necessary to create the Z matrix for our mesh. 
# VERSION 4 Updates -
#   Z function is much more streamlined - 
#   (This is terribly vague information, for future reference, make more detailed
#   comments on what has been updated in this version and why)
# VERSION 5 Updates -
#   Changed get.wghts() function, lines 158-162 to "delete" the point farthest away by assigning it a weight of 0.
#   Change create.Z() function to add rho1 and rho2 and calculate the Matern correlation. 
#   Range and Alpha are now its arguments as well.
#   Added a get.ranks() function that is used in get.wghts() function.
# VERSION 6 Updates -
#   Changed the Create.Z() function so that the X_i values in the augmented Z are in the right columns.
#   Comment out lines 181-183 if you want the 4 locations version for Z.
#########################################


#' Create Regular Mesh Function
#'
#' This function returns a coordinates matrix of a regular mesh created around observation locations.
#' @param x Vector of x-coordinates of the observation locations.
#' @param y Vector of y-coordinates of the observation locations.
#' @param lattice_spacing Space between each point on the regular lattice mesh.
#' @param extension_points The number of points/neighbours beyond the extreme points to remove boundary effects. Defaults to 5.
#' @return A coordinates matrix of a regular lattice.
#' @examples
#' create.regular.mesh(c(1,2,1,2),c(1,1,2,2),0.5,5)
create.regular.mesh <- function(x, y, lattice_spacing, extension_points = 5){
  max_x <- max(x) #Max x-value in the dataset
  min_x <- min(x) #Min x-value in the dataset
  max_y <- max(y) #Max y-value in the dataset
  min_y <- min(y) #Min y-value in the dataset
  x_range <- max_x - min_x #Difference between max and min x value
  y_range <- max_y - min_y #Difference between max and min y value
  
  starting_x_value <- min_x - extension_points*lattice_spacing
  ending_x_value <- max_x + extension_points*lattice_spacing
  starting_y_value <- min_y - extension_points*lattice_spacing
  ending_y_value <- max_y+ extension_points*lattice_spacing
  x_lattice_size <- round( (ending_x_value-starting_x_value)/lattice_spacing ) +1  #Number of points on the x-axis of the mesh - 1.
  y_lattice_size <- round( (ending_y_value-starting_y_value)/lattice_spacing ) +1  #Number of points on the x-axis of the mesh - 1.
  #Update end points to avoid rounding errors
  ending_x_value <- starting_x_value + (x_lattice_size - 1)*lattice_spacing
  ending_y_value <- starting_y_value + (y_lattice_size - 1)*lattice_spacing
  xcoords <- rep( seq(starting_x_value, ending_x_value, length.out = x_lattice_size), times=y_lattice_size)
  ycoords <- rep( seq(starting_y_value, ending_y_value, length.out = y_lattice_size), each=x_lattice_size)  
  #Create a coordinates matrix
  coords.matrix <- cbind(xcoords,ycoords)
  return(coords.matrix)
}



#' Plot Mesh and Observation Locations Function
#'
#' This function plots the mesh and observation locations.
#' @param x Vector of x-coordinates of the observation locations.
#' @param y Vector of y-coordinates of the observation locations.
#' @param mesh Coordinates matrix of the mesh, created by create.regular.mesh() function.
#' @return A plot of the regular lattice and the observations (in red).
#' @examples
#' plot_mesh_and_obs(c(1,2,1,2),c(1,1,2,2),example_mesh)
plot_mesh_and_obs <- function(x,y,mesh, x_lower, x_upper, y_lower, y_upper){
  plot(mesh, pch=".", cex = 3, xlab = "X Coordinate", ylab = "Y Coordinate", 
       asp = 1,
       xlim = c(x_lower, x_upper),
       ylim = c(y_lower, y_upper)) 
  points(x,y, col="red")
}




#' Get Lattice Spacing Function
#'
#' This function extracts the lattice spacing from an existing regular lattice mesh. 
#' It's used in other functions for finding the indices.
#' @param mesh Coordinates matrix of the mesh, created by create.regular.mesh() function.
#' @return The lattice spacing between each points on the mesh.
#' @examples
#' get.lattice_spacing(example_mesh)
get.lattice_spacing <- function(mesh) {
  return( unique(mesh)[2] - unique(mesh)[1] )
}


#' Get Lattice Dimensions Function
#'
#' This function extracts the dimensions of an existing regular lattice mesh. 
#' It's used in other functions, such as when creating Z.
#' @param mesh Coordinates matrix of the mesh, created by create.regular.mesh() function.
#' @return Value of the x and y dimensions of the regular lattice mesh.
#' @examples
#' get.lattice_spacing(example_mesh)
get.lattice_dims <- function(mesh) {
  return( dim(mesh)[1] )
}


#' Get Rank Function
#'
#' This function returns the ranking of values in ascending order for a vector. 
#' It's used in other functions, such as when creating Z.
#' Ties are ranked randomly.
#' @param x Vector of values.
#' @return The ranks of the values in the vector in vector form.
#' @examples
#' get.ranks(c(1,10,4,5,10))
get.ranks <- function(x){
  return(rank(x,ties.method = "random"))
}




#' Get Bottom Left Index Function
#'
#' This function helps identify where an observation location falls on the regular lattice mesh, i.e. inside which
#' square of the mesh the observation position is included. 
#' This is done by identifying the index of the bottom left point in the mesh coordinates matrix.
#' @param x Vector of x-coordinates of the observation locations.
#' @param y Vector of y-coordinates of the observation locations.
#' @param mesh Coordinates matrix of the mesh, created by create.regular.mesh() function.
#' @return Matrix of indices of the bottom left corner of each square that encompasses the observation locations.
#' @examples
#' get.indices(c(1,2,1,2),c(1,1,2,2),example_mesh)
get.indices <- function(x, y, mesh) {
  lattice_spacing <- get.lattice_spacing(mesh)
  #standard_mesh <- cbind(mesh[,1]-min(mesh[,1]), mesh[,2]-min(mesh[,2]))/lattice_spacing
  standard_x <- (x-min(mesh[,1]))/lattice_spacing ##Absolute
  standard_y <- (y-min(mesh[,2]))/lattice_spacing
  indx_BL <- cbind(floor(standard_x), floor(standard_y)) + 1
  return(indx_BL)
}



#' Indices on Z Function
#'
#' This function helps identify where all the observation location fall on the regular lattice mesh, i.e. inside which
#' square of the mesh the observation position is included. 
#' It gets the indices in Z corresponding to the square where the observed position is included.
#' It is used when creating the Z function to assign weights to the correct points on our mesh.
#' @param mesh Coordinates matrix of the mesh, created by create.regular.mesh() function.
#' @param index.Matrix Matrix of indices created by the get.indices() function.
#' @return Matrix of indices of the four corners for each observation location.
#' @examples
#' get.indices(example_mesh,example_indices)
index_Matrix2Vector <- function(mesh, index.Matrix) {
  BL <- length(unique(mesh[,1]))*(index.Matrix[,2]-1) + index.Matrix[,1] 
  BR <- length(unique(mesh[,1]))*(index.Matrix[,2]-1) + index.Matrix[,1] + 1 
  TL <- length(unique(mesh[,1]))*index.Matrix[,2] + index.Matrix[,1] 
  TR <- length(unique(mesh[,1]))*index.Matrix[,2] + index.Matrix[,1] + 1
  return(cbind(BL,BR,TL,TR))
}



#' Get Weights on Z Function
#'
#' This function caculates the weight of each of the four corners of a square on the regular lattice
#' that surrounds an observation location. More weight is given closer the corner is to the observation location.
#' These weights are used along with the indices of these points, when creating the Z function.
#' @param x Vector of x-coordinates of the observation locations.
#' @param y Vector of y-coordinates of the observation locations.
#' @param mesh Coordinates matrix of the mesh, created by create.regular.mesh() function.
#' @return Matrix of indices of the four corners for each observation location.
#' @examples
#' get.wghts(x_vector, y_vector, example_mesh)
get.wghts <- function(x, y, mesh) {
  lattice_spacing <- get.lattice_spacing(mesh)
  N.obs <- length(x) # No. of observations
  standard_x <- (x-min(mesh[,1]))/lattice_spacing
  standard_y <- (y-min(mesh[,2]))/lattice_spacing
  dist_BL <- sqrt( (standard_x - floor(standard_x))^2 + (standard_y - floor(standard_y))^2 )
  dist_BR <- sqrt( (standard_x - (floor(standard_x)+1))^2 + (standard_y - floor(standard_y))^2 )
  dist_TL <- sqrt( (standard_x - floor(standard_x))^2 + (standard_y - (floor(standard_y)+1))^2 )
  dist_TR <- sqrt( (standard_x - (floor(standard_x)+1))^2 + (standard_y - (floor(standard_y)+1))^2 )
  Distances <- cbind(dist_BL, dist_BR, dist_TL, dist_TR)
  small.number = 1e-16
  Distances[Distances == 0] <- small.number #Avoids Inf for observations on a mesh corner 
  max_indices <- max.col((Distances))
  Weights <- 1/Distances
  for (i in 1:N.obs){
    Weights[i,max_indices[i]] <- 0
  }
  Sum <- rowSums(Weights)
  return(Weights/Sum)
}





#' Construct Z Matrix Function
#'
#' This function creates the Z matrix for the regular lattice mesh surrounding the observation area.
#' @param x Vector of x-coordinates of the observation locations.
#' @param y Vector of y-coordinates of the observation locations.
#' @param lattice_spacing Space between each point on the regular lattice mesh.
#' @param extension_points The number of points/neighbours beyond the extreme points to remove boundary effects. Defaults to 5.
#' @return A sparse Matrix of Z for the regular lattice mesh.
#' @examples
#' create.Z(example_mesh, example_indices, example_weights, N = 100)
create.Z <- function(x, y, mesh, range, alpha, id) {
  N <- length(x)
  lattice_spacing <- get.lattice_spacing(mesh)
  rho1 <- get.matern_cov(lattice_spacing, range, alpha)
  rho2 <- get.matern_cov(lattice_spacing*sqrt(2),range, alpha)
  indx.matrix <- get.indices(x, y, mesh) #Get indices of Bottom Left corner of each square that encompasses observations
  indx.vec <- index_Matrix2Vector(mesh, indx.matrix) #Get the index of each four corners
  wghts <- get.wghts(x, y, mesh) #Get weights for each of the four corners surrounding obsevations
  
  
  ranks <- apply(wghts,1,get.ranks)
  x_i.values <- c()
  for(i in 1:N){
    w1 <- wghts[i,which(ranks[,i] == 4)] #Max weight, closest point
    w2 <- wghts[i,which(ranks[,i] == 3)] #Second closest point
    w3 <- wghts[i,which(ranks[,i] == 2)] #Third closest point
    x_i <- sqrt(1- (((w1^2) + (w2^2) + (w3^2)) + 2*rho1*((w1*w2) + (w1*w3)) + 2*rho2*(w2*w3)))
    x_i.values <- c(x_i.values,x_i)
  }
  
  
  dims.lattice <- get.lattice_dims(mesh) #To make sure Z is of correct dimensions.
  Z <- sparseMatrix(i=rep(1:N, each=4), j=as.numeric(t(indx.vec)), x=as.numeric(t(wghts)), dims = c(N,dims.lattice))
  X.Vals <- sparseMatrix(i=rep(1:N, each=1), j=id, x = as.numeric(x_i.values))
  Z <- cbind(Z, X.Vals)
  return(Z)
}



#' Check Z Matrix Function
#'
#' This function performs a check to invesgitate whether our Z has been constructed reasonably. We expect each row
#' to have up to 4 non-zero values. In the special case where the observations fall on a regular lattice mesh exactly,
#' we expect to see only one-non-zero value on each row.
#' @param Z.matrix A Z Matrix created by the create.Z() function.
#' @return A sparse Matrix of Z for the regular lattice mesh.
#' @examples A matrix of indices of non-zero values on the Z matrix and the respective values on those indices.
#' check.Z_matrixValues(example_Z)
check.Z_matrixValues <- function(Z.matrix){
  nonzero_indices <- which(Z.matrix!=0, arr.ind = T) #Check where the non-zero values are on each row.
  Z.value <- c() #Create an empty vector.
  for (i in 1:dim(nonzero_indices)[1]){
    Z.value <- c(Z.value,Z.matrix[nonzero_indices[,1][i],nonzero_indices[,2][i]])
  }
  check_z <- as.data.frame(nonzero_indices) #Create dataframe with only non-zero indices.
  check_z$val <- round(Z.value,2) #Get the values on each of those indices on a new column, rounded to 2 decimal places.
  check_z <- check_z[order(check_z$row),]
  print("Row sums of Z:")
  print(rowSums(Z.matrix)) #Rowsums should be 1.
  print("Dimensions of Z:")
  print(dim(Z.matrix))
  print("Non-zero values:")
  return(check_z)
  print(check_z)
}
