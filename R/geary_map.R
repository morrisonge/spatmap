#######################################################################
#' @title  Local Geary Map
#' @description  The function to make local geary cluster maps.
#' @param polys An sf dataframe
#' @param xname string, the name of the x variable, this variable must be contained in the sf dataframe
#' @param yname string, the name of the y variable, this variable must be contained in the sf dataframe, default option
#' is NULL, can only be used for moran, geary, and join count
#' @param weights weights structure from spdep, must be style "B"; default is set equal to NULL, and first
#' order queen contiguity weights are used to construct the map
#' @param alpha numeric, cut level of significance, must be between 0 and 1, the default is .05
#' @param permutations numeric, number of permutations the conditional randimization approach to significance, maximum is 99999,
#' default is 999
geary_map <- function(polys,
                      xname,
                      yname = NULL,
                      weights = NULL,
                      permutations = 999,
                      alpha = .05){


  # creating weights if left to default
  if (is.null(weights)){
    weights <- default_weights(polys)
  }

  #checking validity of parameters
  check_parameters(polys,permutations,alpha,weights)

  # converting weights to full spatial matrix
  W <- convert_matrix(weights)

  #extracting x variable from sf dataframe
  x <- get_var(xname,polys)
  if (any(is.na(x))){
    stop("x variable cannot have na values")
  }

  if (!is.numeric(x)){
    stop("x must be numeric")
  }

  #extracting y variable from sf dataframe, and assigining NULL if yname is null
  if (!is.null(yname)){
    y <- get_var(yname,polys)
    if (any(is.na(y))){
      stop("y variable cannot have na values")
    }

    if (!is.numeric(y)){
      stop("y must be numeric")
    }
  } else {
    y <- NULL
  }

  #computing observed and reference statistics
  observed <- local_geary(x,y=y,weights)
  sims <- local_geary_sims(x,y=y,weights,permutations=permutations)

  #computing p_values from the observed statistics and reference distributions
  pvalue <- get_pvalue(sims,observed,type = "two-sided")

  #computing mean local geary values for each location based on the reference distributions
  n <- nrow(polys)
  mean_lg <- rep(NA, n)
  for (i in 1:n){
    mean_lg[i] <- mean(sims[i,])
  }


  if(is.null(y)){

    # standardizing x
    z <- standardize(x)
    W <- convert_matrix(weights)

    #classifying univariate local geary patterns
    lg_patterns <- rep(NA,n)
    lg_patterns[z >= 0 & W%*%z >=0] <- "High-High"
    lg_patterns[z <= 0 & W%*%z <=0] <- "Low-Low"
    lg_patterns[observed > mean_lg] <- "Negative"
    lg_patterns[observed < mean_lg & z*W%*%z < 0] <- "Other-Positive"
    lg_patterns[pvalue > alpha] <- "Not Significant"


    # creating the correct palette for univariate geary
    classes <- c("High-High", "Low-Low", "Negative", "Other-Positive", "Not Significant")
    colors <- c("#DE2D26","#FCBBA1","#9ECAE1", "#FEE5D9","#D3D3D3")
    pal <- match_palette(lg_patterns,classes,colors)

  } else {

    #classifying bivariate local geary patterns
    lg_patterns <- as.character(interaction(mean_lg > observed))
    lg_patterns <- lg_patterns %>%
      str_replace_all("TRUE","Positive") %>%
      str_replace_all("FALSE","Negative")
    lg_patterns[pvalue > alpha] <- "Not Significant"

    #creating the correct palette for the bivariate geary
    classes <- c("Positive", "Negative", "Not Significant")
    colors <- c("#C6DBEF","#D3D3D3","#3182BD")
    pal <- match_palette(lg_patterns,classes,colors)
  }

  # Adding the pattern classifications to the sf data frame
  polys["lg_patterns"] <- lg_patterns

  # Making the local geary map
  tm_shape(polys) +
    tm_fill("lg_patterns", palette = pal)
}
