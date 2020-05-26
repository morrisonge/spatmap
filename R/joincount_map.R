#######################################################################
#' @title  Local Join Count Map
#' @description  The function to make local join count cluster maps.
#' @param polys An sf dataframe
#' @param xname string, the name of the x variable, this variable must be contained in the sf dataframe
#' @param weights geoda weights structure, can be genrated with function from rgeoda, the default option is 1st order
#' queen contiguity weights
#' @param alpha numeric, cut-off level of significance, must be between 0 and 1, the default is .05
#' @param permutations numeric, number of permutations the conditional randimization approach to significance, maximum is 99999,
#' default is 999
joincount_map <- function(polys,
                          xname,
                          weights = NULL,
                          permutations = 999,
                          alpha = .05){

  #converting sf to geoda
  gda <- sf_to_geoda(polys)

  # creating weights if left to default
  if (is.null(weights)){
    weights <- queen_weights(gda)
  }

  #extracting x variable from sf dataframe
  x <- get_var(xname,polys)

  #computing local moran lisa
  lisa <- rgeoda::local_joincount(weights, x,perm = permutations)

  #computing lisa, pvalues, clusters, labels, and colors
  clusters <- lisa_clusters(lisa,cutoff = alpha)
  labels <- lisa_labels(lisa)
  pvalue <- lisa_pvalues(lisa)
  colors <- lisa_colors(lisa)

  lisa_patterns <- labels[clusters+1]

  # Constructing the correct palette based on presense of classifications
  pal <- match_palette(lisa_patterns,labels,colors)

  #getting the labels present in the data
  labels <- labels[labels %in% lisa_patterns]

  #adding the cluster values to the data frame
  polys["lisa_clusters"] <- clusters

  #making the map
  tm_shape(polys) +
    tm_fill("lisa_clusters",labels = labels, palette = pal,style = "cat")
}
