#######################################################################
#' @title  Significance Map
#' @description  The function to make significance maps for a variety of local statistics. These
#' statistics include moran, geary, G, G*, and join count. There are multivariate options for
#' geary and join count
#' @param polys An sf dataframe
#' @param vnames string or vector of strings, the name or names of the variables, they must be contained in the sf dataframe
#' @param type string, the type of local statistic, options are: "moran", "geary", "g", "gstar", and "joincount"
#' @param weights  weights structure from rgeoda, the default option is NULL and in this case,
#' the weights will be first order queen contiguity
#' @param alpha numeric, cut level of significance, must be between 0 and 1, the default is .05
#' @param permutations numeric, number of permutations the conditional randimization approach to significance, maximum is 99999,
#' default is 999
#' @export
significance_map <- function(polys,
                             vnames,
                             type,
                             weights = NULL,
                             alpha = .05,
                             permutations = 999){

  #checking for supported types
  supported_types <- c("moran","geary","g","gstar","joincount")
  if (!(type %in% supported_types)){
    stop("type is not support, choose from: moran, geary, g , gstar, join_count")
  }

  if (length(vnames) > 1){
    ysupported_types <- c("geary","joincount")
    if(!(type %in% ysupported_types)){
      stop("moran, g, and gstar do not have multivariate options")
    }
  }

  #converting sf to geoda
  gda <- sf_to_geoda(polys)


  # creating weights if left to default
  if (is.null(weights)){
    weights <- queen_weights(gda)
  }

  #getting variables or variable from the sf dataframe
  if (length(vnames) == 1){
    x <- get_var(vnames,polys)
  } else {
    df <- polys %>% select(vnames)
    st_geometry(df) <- NULL
  }

  #computing moran
  if (type == "moran"){
    lisa <- rgeoda::local_moran(weights, x,perm = permutations)
  }

  #computing geary
  if (type == "geary"){
    if (length(vnames) == 1){
      lisa <- rgeoda::local_geary(weights, x,perm = permutations)
    } else {
      lisa <- rgeoda::local_multigeary(weights, df,perm = permutations)
    }
  }

  #computing g
  if (type == "g"){
    lisa <- rgeoda::local_g(weights, x,perm = permutations)
  }

  #computing gstar
  if (type == "gstar"){
    lisa <- rgeoda::local_gstar(weights, x,perm = permutations)
  }

  #computing join count
  if (type == "join_count"){
    if (length(vnames) == 1){
      lisa <- rgeoda::local_joincount(weights, x,perm = permutations)
    } else {
      lisa <- rgeoda::local_multijoincount(weights, df,perm = permutations)
    }
  }

  #computing pvalues
  pvalue <- lisa_pvalues(lisa)

  #creating breaks based on p-values
  target_p <- 1 / (1 + permutations)
  potential_brks <- c(.00001, .0001, .001, .01)
  brks <- potential_brks[which(potential_brks > target_p & potential_brks < alpha)]
  brks2 <- c(target_p, brks, alpha)
  labels <- c(as.character(brks2), "Not Significant")
  brks3 <- c(-Inf, brks2, 1)

  cuts <- cut(pvalue, breaks = brks3,labels = labels)

  # Adding the p-value significance breaks to the data frame
  polys["sig"] <- cuts


  # Constructing the correct palette
  pal <- rev(brewer.pal(length(labels), "Greens"))
  pal[length(pal)] <- "#D3D3D3"


  # Making the map
  tm_shape(polys) +
    tm_fill("sig", palette = pal)
}
