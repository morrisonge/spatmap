


multi_joincount_map <- function(polys,
                            vnames,
                            weights = NULL,
                            permutations = 999,
                            alpha = .05){
  
  
  #converting sf to geoda
  gda <- sf_to_geoda(polys)
  
  # creating weights if left to default
  if (is.null(weights)){
    weights <- queen_weights(gda)
  }
  
  df <- polys %>% select(vnames)
  st_geometry(df) <- NULL
  
  #computing local moran lisa
  lisa <- rgeoda::local_multijoincount(weights, df,perm = permutations)
  
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













