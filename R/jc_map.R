#######################################################################
#' @title  Local Join Count Map
#' @description  The function to make local join count cluster maps. 
#' @param polys An sf dataframe
#' @xname string, the name of the x variable, this variable must be contained in the sf dataframe
#' @yname string, the name of the y variable, this variable must be contained in the sf dataframe, default option
#' is NULL
#' @param weights weights structure from spdep, must be style "B"; default is set equal to NULL, and first
#' order queen contiguity weights are used to construct the map
#' @param alpha numeric, cut-off level of significance, must be between 0 and 1, the default is .05
#' @param permutations numeric, number of permutations the conditional randimization approach to significance, maximum is 99999,
#' default is 999



jc_map <- function(polys, 
                   xname,
                   yname = NULL, 
                   weights = NULL, 
                   permutations = 999, 
                   alpha = .05){
  
  
  # creating weights if left to default 
  if (is.null(weights)){
    weights <- default_weights(polys)
  }
  
  #checking the validity of parameters
  check_parameters(polys,permutations,alpha,weights)
  
  
  #extracting x variable from sf dataframe
  x <- get_var(xname,polys)
  if (!(is.binary(x))){
    stop("xname must designate a binary variable")
  }
  
  
  if(!is.null(yname)){
    y <- get_var(yname,polys)
    if (!(is.binary(y))){
      stop("yname must designate a binary variable")
    } 
  } else {
    y <- NULL
  }

    
  #computing p_values from the observed statistics and reference distributions
  p_value <- local_jc_pvalue(x,y=y,weights,permutations = permutations)
  
  
  #Assigning 1 to significant locations
  n <- nrow(W)
  unique_value <- rep(0,n)
  unique_value[p_value <= alpha] <- 1
  
  #Creating the correct palette
  classes <- c(0,1)
  colors <- c("white","blue")
  pal <- match_palette(unique_value, classes,colors)
  
  #Adding the values to the sf dataframe
  polys["unique_value"] <- unique_value
  
  #Making the map
  tm_shape(polys) +
    tm_fill("unique_value", palette = pal)
}



















