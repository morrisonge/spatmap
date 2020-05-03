#######################################################################
#' @title  Local Moran Statistics
#' @description  The function to compute univariate and bivariate local moran statistics
#' @param x A vector of numerical values
#' @param y A vector of numerical values, only to be used in the bivariate case
#' @param weights Weights structure from spdep, must be style "B"
#' @return local A vector of local moran statistics
#' @export

local_moran <- function(x, y = NULL, weights){
  #checking validity of weights
  check_weights(weights)

  #converting weights
  W <- convert_matrix(weights)

  #setting y = x in univariate case
  if(is.null(y)) y = x

  #standardizing x and y
  xs <- standardize(x)
  ys <- standardize(y)

  #computing the local moran
  local  <- (xs*W%*%ys)
  return(local)
}

#######################################################################
#' @title Local Moran Simulations
#' @description Function to compute reference distributions of local moran statistics
#' @param x A vector of numerical values
#' @param y A vector of numerical values, only to be used in the bivariate case
#' @param weights Weights structure from spdep, must be style "B"
#' @param permutations Number of permutations, the default is 999
#' @return local.sims Reference distributions of local moran statistics for each location in matrix form
#' @export
local_moran_sims <- function(x, y = NULL, weights, permutations = 999){

  #checking parameters
  check_permutations(permutations)
  check_weights(weights)


  #converting spdpe weights to full size weights matrix
  W <- convert_matrix(weights)


  #setting x = y in univariate case
  if(is.null(y)) y = x

  #making id variable
  n   <- nrow(W)
  id  <- 1:n

  #place to store results
  local.sims  <- matrix(NA, nrow = n, ncol=permutations)
  y.sample = matrix(NA, nrow = n, ncol = permutations)

  # filling each row of the sample matrix
  for(i in 1:n){
    sample.indices <- sample(id[-i], permutations, replace = TRUE)
    y.sample[i,] <- y[sample.indices]
  }

  #standardizing the y sample values
  y.sample <- (y.sample - apply(y.sample, 1, mean))/apply(y.sample, 1, sd)

  #standardizing x
  xs <- standardize(x)

  #calculating the local statistics
  local.sims  <- (xs*W%*%y.sample)
  return(local.sims)
}

#######################################################################
#' @title Local Moran P-values
#' @description Function to compute the p-value of the observed local moran statistics under the
#' conditional randomizaton approach
#' @param x A vector of numerical values
#' @param y A vector of numerical values, only to be used in the bivariate case
#' @param weights Weights structure from spdep, must be style "B"
#' @param permutations Number of permutations, the default is 999
#' @return pvalue A vector of p-values for each locations local moran statistic
#' @export

local_moran_pvalue <- function(x, y = NULL, weights, permutations = 999){

  #computing observed statistics and reference distributions
  observed <- local_moran(x,y=y,weights)
  sims <- local_moran_sims(x,y=y,weights,permutations = permutations)

  #computing p-value
  pvalue <- get_pvalue(sims,observed,type = "two-sided")
  return(pvalue)

}
