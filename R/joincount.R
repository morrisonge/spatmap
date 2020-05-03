#######################################################################
#' @title  Local Join Count Statistics
#' @description  The function to compute univariate and bivariate local join count statistics
#' @param x A vector of numerical values
#' @param y A vector of numerical values, default is NULL
#' @param weights Weights structure from spdep, must be style "B"
#' @return local A vector of local joint count statistics
#' @export 

local_jc <- function(x,y = NULL,weights){
  
  #checking for valid weights
  check_weights(weights)
  if(!is.binary(x)){
    stop("x must be a binary variable")
  }
  
  #converting weights to full size binary spatial weights matrix
  W <- convert_matrix(weights)
  B <- binarize(W,threshold = .00001)
  
  #computing the join count statistic for the bivariate and univariate case
  if(!is.null(y)){
    if(!is.binary(y)){
      stop("x must be a binary variable")
    }
    x <- x * y
  }
  jc <- x*B%*%x
  return(jc)
}

#######################################################################
#' @title Local Join Count Simulations
#' @description Function to compute reference distributions of local join count statistics
#' @param x A vector of numerical values
#' @param z A vector of numerical values, only to be used in the bivariate case, default is NULL
#' @param weights Weights structure from spdep, must be style "B"
#' @param permutations Number of permutations, the default is 999
#' @return local.sims Reference distributions of local join count statistics for each location in matrix form
#' @export

local_jc_sims <- function(x,y = NULL,weights,permutations){
  
  #checking parameters
  check_permutations(permutations)
  check_weights(weights)
  
  #checking x
  if(!is.binary(x)){
    stop("x must be a binary variable")
  }
  
  #converting weights to full size binary spatial weights matrix
  W <- convert_matrix(weights)
  B <- binarize(W,threshold = .00001)
  
  #making an id variable
  n <- nrow(W)
  id <- 1:n
  
  #creating matrix to store reference statistics
  local.sims <- matrix(NA, nrow = n, ncol=permutations)
  x.sample = matrix(NA, nrow = n, ncol = permutations)
  
  #setting up for the bivariate version of the statistic and checking y
  if(!is.null(y)){
    if(!is.binary(y)){
      stop("y must be a binary variable")
    }
    x <- x * y
  }
  
  #filling each row of the sample matrix
  for(i in 1:n){
    sample.indices <- sample(id[-i], permutations, replace = TRUE)
    x.sample[i,] <- x[sample.indices]
  }
  
  #computing the join count statistic
  local.sims <- x.sample * B%*%x.sample
  return(local.sims)
}

#######################################################################
#' @title Local Join Count P-values
#' @description Function to compute the p-value of the observed local join count statistics under the
#' conditional randomizaton approach
#' @param x A vector of numerical values
#' @param z A vector of numerical values, only to be used in the bivariate case,default is NULL
#' @param weights Weights structure from spdep, must be style "B"
#' @param permutations Number of permutations, the default is 999
#' @return pvalue A vector of p-values for each locations local geary statistic
#' @export

local_jc_pvalue <- function(x,y=NULL,weights,permutations = 999){
  
  #computing observed statistics and reference distributions
  observed <- local_jc(x,y=y,weights)
  sims <- local_jc_sims(x,y=y,weights,permutations = permutations)
  
  #computing p-value
  pvalue <- get_pvalue(observed,sims,type = "one-sided")
  return(pvalue)
}

