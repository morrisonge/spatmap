#######################################################################
#' @title  Local Moran Map
#' @description  The function to make local moran cluster maps.
#' @param polys An sf dataframe
#' @xname string, the name of the x variable, this variable must be contained in the sf dataframe
#' @yname string, the name of the y variable, this variable must be contained in the sf dataframe, default option
#' is NULL, can only be used for moran, geary, and join count
#' @param weights weights structure from spdep, must be style "B"; default is set equal to NULL, and first
#' order queen contiguity weights are used to construct the map
#' @param alpha numeric, cut-off level of significance, must be between 0 and 1, the default is .05
#' @param permutations numeric, number of permutations the conditional randimization approach to significance, maximum is 99999,
#' default is 999

moran_map <- function(polys,
                      xname,
                      yname = NULL,
                      weights = NULL,
                      permutations = 999,
                      alpha){


  # creating weights if left to default
  if (is.null(weights)){
    weights <- default_weights(polys)
  }

  #checking validity of parameters
  check_parameters(polys,permutations,alpha,weights)

  # extracting x variable from sf dataframe
  x <- get_var(xname,polys)

  # assigning y variable
  if (is.null(yname)){
    y <- x
  } else {
    y <- get_var(yname,polys)
  }

  #computing p-value from observed statistics and reference distribution
  p_value <- local_moran_pvalue(x,y=y,weights,permutations = permutations)

  # Standardizing the x and y vars
  z_x <- standardize(x)
  z_y <- standardize(y)
  W <- convert_matrix(weights)

  #Assigning classifications
  lisa_patterns <- as.character( interaction(z_x > 0, W%*%z_y > 0) )
  lisa_patterns <- patterns %>%
    str_replace_all("TRUE","High") %>%
    str_replace_all("FALSE","Low")
  lisa_patterns[which(p_value > alpha)] <- "Not Significant"

  # Adding the classifications to the data frame
  polys["lisa_patterns"] <- lisa_patterns

  # Constructing the correct palette based on presense of classifications
  patts <- c("High.High", "High.Low", "Low.High", "Low.Low", "Not Significant")
  colors <- c("#DE2D26","#FCBBA1","#C6DBEF", "#3182BD", "#D3D3D3")
  pal <- match_palette(lisa_patterns,patts,colors)


  # Making the map
  tm_shape(area) +
    tm_fill("lisa_patterns", palette = pal)

}
