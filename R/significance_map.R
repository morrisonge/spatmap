#######################################################################
#' @title  Significance Map
#' @description  The function to make significance maps for a variety of local statistics. These
#' statistics include moran, geary, G, G*, and join count. There are bivariate options for moran,
#' geary, and join count
#' @param polys An sf dataframe
#' @xname string, the name of the x variable, this variable must be contained in the sf dataframe
#' @yname string, the name of the y variable, this variable must be contained in the sf dataframe, default option
#' is NULL, can only be used for moran, geary, and join count
#' @param type string, the type of local statistic, options are: "moran", "geary", "g", "gstar", and "join_count"
#' @param weights weights structure from spdep, must be style "B"; default is set equal to NULL, and first
#' order queen contiguity weights are used to construct the map
#' @param alpha numeric, cut level of significance, must be between 0 and 1, the default is .05
#' @param permutations numeric, number of permutations the conditional randimization approach to significance, maximum is 99999,
#' default is 999
#' @export

significance_map <- function(polys,
                             xname,
                             yname = NULL,
                             type,
                             weights = NULL,
                             alpha = .05,
                             permutations = 999){
  
  
  check_parameters(polys,permutations,alpha,weights)
  
  
  #checking for supported types
  supported_types <- c("moran","geary","g","gstar","join_count")
  if (!(type %in% supported_types)){
    stop("type is not support, choose from: moran, geary, g , gstar, join_count")
  }
  
  if (!is.null(yname)){
    ysupported_types <- c("moran","geary","join_count")
    if(!(type %in% ysupported_types)){
      stop("g and gstar do not have bivariate options")
    }
  }
  
  
  # creating weights if left to default 
  if (is.null(weights)){
    weights <- default_weights(polys)
  }
  
 
  
  #extracting x variable from sf dataframe
  x <- get_var(vname,polys)
  if (!is.null(yname)){
    y <- get_var(poly,yname)
  }
  
  
  #computing moran p-values
  if (type == "moran"){
    pvalue <- local_moran_pvalue(x,y=y,weights,permutations = permutations)
  }
  
  #computing geary p-values
  if (type == "geary"){
    pvalue <- local_geary_pvalue(x,y=y,weights,permutations = permutations)
  }
  
  #computing g p-values
  if (type == "g"){
    pvalue <- local_g_pvalue(x,weights,permutations = permutations)
  }
  
  #computing gstar p-values
  if (type == "gstar"){
    pvalue <- local_g_pvalue(x,weights,permutations = permutations,type = "gstar")
  }
  
  #computing join count p-values
  if (type == "join_count"){
    pvalue <- local_jc_pvalue(x,y=y,weights,permutations = permutations)
  }
  
  #Creating breaks based on p-values
  target_p <- 1 / (1 + permutations)
  potential_brks <- c(.00001, .0001, .001, .01)
  brks <- potential_brks[which(potential_brks > target_p & potential_brks < alpha)]
  brks2 <- c(target_p, brks, alpha)
  labels <- c(as.character(brks2), "Not Significant")
  brks3 <- c(-inf, brks2, 1)
  
  cuts <- cut(pvalue, breaks = brks3,labels = labels)
  
  # Adding the p-value significance breaks to the data frame
  polys <- polys %>% mutate(sig = cuts)
  
  
  # Constructing the correct palette
  pal <- rev(brewer.pal(length(labels), "Greens"))
  pal[length(pal)] <- "#D3D3D3"
  
  
  # Making the map
  tm_shape(polys) +
    tm_fill("sig", palette = pal)
}
