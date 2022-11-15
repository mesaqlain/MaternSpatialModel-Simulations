# Version 21st of September 2020 from MSQ
#########################################
# This script contains the functions to create a regular lattice and 
# It also contains a function for the marginal variance of the Matern function and the Matern Covariance itself.
# Finally, it contains a function to plot the Matern covariance against Q inverse.
# Version 4 changes: 
#     my.marginal.variance() is now a function of lattice_spacing.
# Version 4: check.Q_matrix() has not been updated from Version 3.
#########################################



#' Create Coordinates Matrix Function
#'
#' This function returns a coordinates matrix from x and y vectors.
#' @param x Vector of x-coordinates of the observation locations.
#' @param y Vector of y-coordinates of the observation locations.
#' @return A coordinates matrix of x and y.
#' @examples
#' create.coords_mat(c(1,1,1,2,2,2),c(1,2,3,1,2,3))
create.coords_mat <- function(x,y){
  coords.mat <- matrix(c(x,y), ncol = 2, byrow = F)
  return(coords.mat)
}



#' Create Distance Matrix Function
#'
#' This function returns a distance  matrix for given coordinates. Potentially huge memory usage.
#' This function is used in check.Q_matrix() function.
#' @param coords.mat Coordinates matrix. Can be created by create.coords_mat() function if needed.
#' @return A distance matrix for given coordinates.
#' @examples
#' create.coords_mat(example_coordinates_matrix)
create.dist_matrix <- function(coords.mat){
  dist.matrix <- as.matrix(dist(coords.mat))
  return(dist.matrix)
}



#' Create Regular Lattice Function
#'
#' This function creates a regular lattice. The indexing starts from the bottom left most point and goes row by row.
#' This function is used in check.Q_matrix() function.
#' @param starting_x_value The value of x-coordinate of the bottom leftmost point on the lattice.
#' @param starting_y_value The value of y-coordinate of the bottom leftmost point on the lattice.
#' @param x_lattice_size The number of points in the x-direction.
#' @param y_lattice_size The number of points in the y-direction.
#' @param lattice_spacing The spacing between each point.
#' @return A matrix of coordinates for a regular lattice.
#' @examples
#' create.regular_lattice(0,0,20,20,1)
create.regular_lattice <- function(starting_x_value = 0,
                                   starting_y_value = 0,
                                   x_lattice_size = 10, 
                                   y_lattice_size = 10, 
                                   lattice_spacing = 0.5){
  ending_x_value <- starting_x_value + (x_lattice_size-1)*lattice_spacing
  ending_y_value <- starting_y_value + (y_lattice_size-1)*lattice_spacing
  xcoords <- rep( seq(starting_x_value, ending_x_value, length.out = x_lattice_size), times=y_lattice_size)
  ycoords <- rep( seq(starting_y_value, ending_y_value, length.out = y_lattice_size), each=x_lattice_size)  
  coords.matrix <- matrix(c(xcoords,ycoords),ncol=2,byrow=F)
  return(coords.matrix)
}



#' Marginal Variance of Matern Function
#'
#' This function returns the value of the marginal variance of the Matern function, 
#' defined in Lindgren and Rue 2011 pp 427.
#' This function is used in the get.matern_cov() function.
#' @param range The parameter 'range' in the Matern function.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3.
#' @return The value of marginal variance of the Matern function.
#' @examples
#' my.marginal.variance(10,2)
my.marginal.variance <- function(range, alpha, lattice_spacing){
  nu <- my.nu(alpha)
  kappa <- my.inv.range(range,alpha, lattice_spacing)
  m.var <- gamma(nu)/((gamma(nu+1))*(4*pi)*(kappa^(2*nu)))
  return(m.var)
}



#' Matern Covariance Function
#'
#' This function returns the Matern Covariance for two-dimensional data.
#' @param d Distance between two points.
#' @param range The parameter 'range' in the Matern function.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3. alpha = nu + d/2.
#' @return Matern covariance.
#' @examples
#' get.matern_cov(c(1,2,3),0.5,2)
get.matern_cov <-function(d, range, alpha) {
  nu <- my.nu(alpha)
  #kappa <- my.inv.range(range,alpha, lattice_spacing)
  kappa <- 2 * sqrt(2 * nu) / range
  #m.var <- my.marginal.variance(range,alpha)
  m.var <- 1
  d <- pmax(d, 1.0e-8) #To avoid NaN on the diagonal
  c <- (m.var/((2^(nu-1))*gamma(nu))) * ((kappa*d)^nu) * besselK(kappa * d, nu = nu)  
  return(c)
}



#' Function to Plot Q inverse against Matern Covariance for a regular lattice
#'
#' This function plots the Matern Covariance - get.matern_cov() - for a regular lattice against the 
#' inverse of the precision matrix found by using the create.Q_matrix() function.
#' @param range The parameter 'range' in the Matern function.
#' @param alpha The parameter 'alpha' in the SPDE. Values = 1, 2, or 3. alpha = nu + d/2.
#' @param starting_x_value The value of x-coordinate of the bottom leftmost point on the lattice.
#' @param starting_y_value The value of y-coordinate of the bottom leftmost point on the lattice.
#' @param x_lattice_size The number of points in the x-direction.
#' @param y_lattice_size The number of points in the y-direction.
#' @param lattice_spacing The spacing between each point.
#' @return A pdf file in the working directory of Matern Covariance vs Q inverse. 
#' @examples
#' check.Q_matrix(1,1,0,0,20,20,0.5)
check.Q_matrix <- function(range, alpha, starting_x_value, starting_y_value, x_lattice_size, y_lattice_size, lattice_spacing){
  
  #Create an example lattice
  example_lattice <- create.regular_lattice(starting_x_value,starting_y_value,
                                            x_lattice_size,y_lattice_size,
                                            lattice_spacing)
  
  #Create Q matrix for our coordinates
  example_Q.mat <- create.Q_matrix(alpha,range,example_lattice)
  
  #Create distance matrix
  dist_matrix <- (create.dist_matrix(example_lattice))
  #Create distance vector
  dist_vector <- as.vector(dist_matrix)
  
  #Matern Covariance
  example_maternCov <- get.matern_cov(d = dist_vector, 
                                      range = range, 
                                      alpha = alpha)
  
  #Turn that into vector to plot against distance
  Matern_vector <- as.vector(example_maternCov)
  
  #Create a file name and start saving the plots.
  file_name <- paste("alpha-", alpha, "_range-", range, "_lattice_spacing-", lattice_spacing, ".pdf", sep="")
  pdf(file_name)
  
  #Create a plot of Matern covariance against distance
  plot(dist_vector,Matern_vector, 
       xlab = "Distance",
       ylab = "Correlation",
       main = paste("Alpha-",alpha," Range-",range," Lat Space-",lattice_spacing))
  
  #Find inverse of Q
  Q.inv <- solve(example_Q.mat)
  #Divide by the max to normalize it.
  Q.inv <- Q.inv/max(Q.inv)
  #Take an index somewhere in the middle of Q
  ref_index <- round((dim(Q.inv)[1])/2)
  #Plot the points of the Q.inv on the graph
  points(dist_matrix[ref_index,],Q.inv[ref_index,], pch = 3, lty=3, lwd=1, col = "blue")
  #End plotting to save to pdf
  dev.off()
}
