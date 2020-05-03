#######################################################################
#' @title  Local G Statistics
#' @description  The function to compute G and G* statistics
#' @param x A vector of numerical values
#' @param weights Weights structure from spdep, must be style "B"
#' @param type A string, designates g or gstar
#' @return local A vector of local geary statistics
#' @export   
local_g <- function(x,weights, type = "g"){
  
  #checking type
  supported_types <- c("g","gstar")
  if (!(type %in% supported_types)){
    stop("type not supported: choose g or gstar")
  }
  
  #checking for valid weights
  check_weights(weights)
  
  #computing initial steps in the local g and gstar statistic
  n <- nrow(W)
  lag <- W%*%x
  sum <- sum(x)
  
  #excluding element i in the sum for each location in the case of the g statistic
  if (type == "g"){
    sum <- rep(NA,n)
    for (i in 1:n){
      sum <- sum(x[-i])
    }
  }
  g <- lag/sum
  return(g)
}


#######################################################################
#' @title Local G Simulations
#' @description Function to compute reference distributions of local G and G* statistics
#' @param x A vector of numerical values
#' @param weights Weights structure from spdep, must be style "B"
#' @param permutations Number of permutations, the default is 999
#' @param type designates the type of statistic g or gstar, g is the default
#' @return local.sims Reference distributions of local geary statistics for each location in matrix form
#' @export

local_g_sims <- function(x,weights, permutations, type = "g"){
  
  #checking type
  supported_types <- c("g","gstar")
  if (!(type %in% supported_types)){
    stop("type not supported: choose g or gstar")
  }
  
  #checking parameters
  check_permutations(permutations)
  check_weights(weights)
  
  #making id variable
  n <- nrow(W)
  id <- 1:n
  
  #creating structure to store the reference distributions for each location
  local.sims  <- matrix(NA, nrow = n, ncol=permutations)
  x.sample = matrix(NA, nrow = n, ncol = permutations)
  
  #creating random samples for 999 permutations for each location
  for(i in 1:n){
    sample.indices <- sample(id[-i], permutations, replace = TRUE)
    x.sample[i,] <- x[sample.indices]
  }
  
  # applying the G statsitc formula to the sample values
  lag <- W%*%x.sample
  sum <- sum(x)
  
  #computing sums excluding element i for the g statistic
  if (type == "g"){
    sum <- rep(NA,n)
    for (i in 1:n){
      sum <- sum(x[-i])
    }
  }
  local.sims <- lag / sum
  return(local.sims)
}

#######################################################################
#' @title Local G P-values
#' @description Function to compute the p-value of the observed local G and G* statistics under the
#' conditional randomizaton approach
#' @param x A vector of numerical values
#' @param weights Weights structure from spdep, must be style "B"
#' @param permutations Number of permutations, the default is 999
#' @param type 
#' @return pvalue A vector of p-values for each locations local geary statistic

local_g_pvalue <- function(x,weights,permutations = 999,type = "g"){
  
  #computing observed statistics and reference distributions
  observed <- local_g(x,weights,type = type)
  sims <- local_g_sims(x,weights,permutations = permutations,type = type)
  
  #computing p-value
  pvalue <- get_pvalue(observed,sims,type = "two-sided")
  return(pvalue)
}
