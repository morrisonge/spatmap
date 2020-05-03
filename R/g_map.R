#######################################################################
#' @title  Local G Map
#' @description  The function to make local G and G* cluster maps. 
#' @param polys An sf dataframe
#' @xname string, the name of the x variable, this variable must be contained in the sf dataframe
#' @param weights weights structure from spdep, must be style "B"; default is set equal to NULL, and first
#' order queen contiguity weights are used to construct the map
#' @param alpha numeric, cut-off level of significance, must be between 0 and 1, the default is .05
#' @param permutations numeric, number of permutations the conditional randimization approach to significance, maximum is 99999,
#' default is 999
#' @param type string, can be "g" or "gstar"
#' @export

g_map <- function(polys,xname,weights = NULL, permutations = 999, alpha = .05,type = "g"){
  
  # creating weights if left to default 
  if (is.null(weights)){
    weights <- default_weights(polys)
  }
  
  #checking validity of parameters
  check_parameters(polys,permutations,alpha,weights)
  
  #extracting x variable from sf dataframe
  x <- get_var(xname,polys)
  
  #computing p_values from the observed statistics and reference distributions
  pvalue <- local_g_pvalue(x,weights,permutations = permutations,type = type)
  
  #assigning cluster classifications
  n <- nrow(W)
  g_patterns <- rep(NA,n)
  g_patterns[g > mean(g)] <- "High"
  g_patterns[g < mean(g)] <- "Low"
  g_patterns[p_value > alpha] <- "Not Significant" 
  
  
  #Creating the correct color palette
  classes <- c("High", "Low", "Not Significant")
  colors <- c("#DE2D26", "#3182BD", "#D3D3D3")
  pal <- match_palette(g_patterns,classes,colors)
  
  #Adding the patterns to the sf dataframe
  polys["g_patterns"] <- g_patterns 
  
  #Making the map
  tm_shape(polys) +
    tm_fill("g_patterns", palette = pal)
}
