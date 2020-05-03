#######################################################################
#' @title  Local Geary Statistics
#' @description  The function to compute univariate and bivariate local geary statistics
#' @param x A vector of numerical values
#' @param y A vector of numerical values, only to be used in the bivariate case
#' @param weights Weights structure from spdep, must be style "B"
#' @return local A vector of local geary statistics
#' @export

local_geary <- function(x,y = NULL, weights){
  #checking validity of weights
  check_weights(weights)
  W <- convert_matrix(weights)
  x2 <- x^2
  if (is.null(y)){
    lg <- x2 - 2*x*(W%*%x) + W%*%x2
  } else {
    y2 <- y^2
    lg <- (x2 - 2*x*(W%*%x) + W%*%x2) + (y2 - 2*y*(W%*%y) + W%*%y2)
  }
  return(lg)
}

#######################################################################
#' @title Local Geary Simulations
#' @description Function to compute reference distributions of local geary statistics
#' @param x A vector of numerical values
#' @param y A vector of numerical values, only to be used in the bivariate case
#' @param weights Weights structure from spdep, must be style "B"
#' @param permutations Number of permutations, the default is 999
#' @return local.sims Reference distributions of local geary statistics for each location in matrix form
#' @export


local_geary_sims <- function(x,y=NULL,weights,permutations = 999){

  #checking parameters
  check_permutations(permutations)
  check_weights(weights)

  #converting weights
  W <- convert_matrix(weights)

  #making id variable
  n <- nrow(W)
  id  <- 1:n

  #setting up data structures to store reference statistics
  local.sims  <- matrix(NA, nrow = n, ncol=permutations)
  x.sample <- matrix(NA, nrow = n, ncol = permutations)
  if(!is.null(y)){
    y.sample <- matrix(NA, nrow = n, ncol = permutations)
  }

  #Assigning sample values for x and y
  for(i in 1:n){
    sample.indices <- sample(id[-i], permutations, replace = TRUE)
    x.sample[i,] <- x[sample.indices]
    if(!is.null(y)){
      y.sample[i,] <- y[sample.indices]
    }
  }

  #Computing reference statistc values for x
  x2 <- x^2
  x2.sample <- x.sample ^ 2
  local.sims <- x2 - 2*x*W%*%x.sample + W%*%x2.sample

  #Adding the reference statistic component for y in the bivariate case
  if (!is.null(y)){
    y2 <- y^2
    y2.sample <- y.sample ^ 2
    local.sims  <- local.sims + (y2 - 2*y*W%*%y.sample + W%*%y2.sample)
  }
  return(local.sims)
}

#######################################################################
#' @title Local Geary P-values
#' @description Function to compute the p-value of the observed local geary statistics under the
#' conditional randomizaton approach
#' @param x A vector of numerical values
#' @param y A vector of numerical values, only to be used in the bivariate case
#' @param weights Weights structure from spdep, must be style "B"
#' @param permutations Number of permutations, the default is 999
#' @return pvalue A vector of p-values for each locations local geary statistic
#' @export

local_geary_pvalue <- function(x,y=NULL,weights,permutations = 999){

  #computing observed statistics and reference distributions
  observed <- local_geary(x,y=y,weights)
  sims <- local_geary_sims(x,y=y,weights,permutations = permutations)

  #computing p-value
  pvalue <- get_pvalue(sims,observed,type = "two-sided")
  return(pvalue)
}
