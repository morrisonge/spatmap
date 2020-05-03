# mapping function




###################################################################################################
# Univariate Moran
###################################################################################################


moran_map <- function(polys,xname,yname = NULL, weights = NULL,permutations = 999, alpha){
  
  
  
  if (permutations > 99999){
    stop("maximum number of permutations is 99999")
  }
  
  if (alpha > 1 | alpha < 0){
    stop("alpha must be between 0 and 1")
  }
  
  
  if (any(class(polys) == "sf")){
    if (!any(class(st_geometry(polys)) == "sfc_MULTIPOLYGON")){
      stop("geomtry type must be sfc_MULTIPOLYGON for sf dataframes")
    } 
  } else {
    stop("sf is only supported geometry")
  }
  
  # creating weights if left to default
  if (is.null(weights)){
    sgbp <- st_queen(polys)
    nb <- as_nb_sgbp(sgbp)
    weights <- nb2listw(queen.nb,style = "B", zero.policy = TRUE)
  }
  
  #converting weight to full spatial matrix
  W <- convert_matrix(weights)
  
  # extracting x variable from sf dataframe
  x <- get_var(xname,polys)
  
  # 
  if (is.null(yname)){
    y <- x
  } else {
    y <- get_var(yname,polys)
  }
  
  
  #computing observed statistic and reference distributions
  observed <- local_moran(x,y= y,W)
  sims <- local_moran_sims(x,y = y,W,permutations)
  
  
  #computing p-value from observed statistics and reference distribution
  p_value <- get_p_value(sims,observed, type = "two-sided")
  
  # Standardizing the x and y vars
  z_x <- standardize(x)
  z_y <- standardize(y)
  
  #Assigning classifications
  lisa_patterns <- as.character( interaction(z_x > 0, W%*%z_y > 0) ) 
  lisa_patterns <- patterns %>% 
    str_replace_all("TRUE","High") %>% 
    str_replace_all("FALSE","Low")
  lisa_patterns[which(p_value > alpha)] <- "Not Significant"
  
  # Adding the classifications to the data frame
  area <- area %>% mutate(lisa_patterns)
  
  # Constructing the correct palette based on presense of classifications
  patts <- c("High.High", "High.Low", "Low.High", "Low.Low", "Not Significant")
  colors <- c("#DE2D26","#FCBBA1","#C6DBEF", "#3182BD", "#D3D3D3")
  pal <- match_palette(lisa_patterns,patts,colors)
  
  
  # Making the map
  tm_shape(area) +
    tm_fill("lisa_patterns", palette = pal)
  
}

###################################################################################################
# Univariate Geary
###################################################################################################


geary_map <- function(polys,xname, weights = NULL,permutations = 999, alpha = .05){
  
  
  if (permutations > 99999){
    stop("maximum number of permutations is 99999")
  }
  
  if (alpha > 1 | alpha < 0){
    stop("alpha must be between 0 and 1")
  }
  
  
  if (any(class(polys) == "sf")){
    if (!any(class(st_geometry(polys)) == "sfc_MULTIPOLYGON")){
      stop("geomtry type must be sfc_MULTIPOLYGON for sf dataframes")
    } 
  } else {
    stop("sf is only supported geometry")
  }
  
  
  
  
  
  # creating weights if left to default 
  if (is.null(weights)){
    sgbp <- st_queen(polys)
    nb <- as_nb_sgbp(sgbp)
    weights <- nb2listw(queen.nb,style = "B", zero.policy = TRUE)
  }
  
  # converting weights to full spatial matrix
  W <- convert_matrix(weights)
  
  #extracting x variable from sf dataframe
  x <- get_var(xname,polys)
  
  
  #computing observed statistics and reference distributions for each location
  observed <- local_geary(x,W)
  sims <- local_geary_sims(x,W,permutations)
  
  #computing p_values from the observed statistics and reference distributions
  p_value <- get_p_value(sims,observed,type = "two-sided")
  
  #computing mean local geary values for each location based on the reference distributions
  mean_lg <- rep(NA, n)
  for (i in 1:n){
    mean_lg[i] <- mean(sims[i,]) 
  }
  
  # standardizing x
  z <- standardize(x)
  
  #classifying local geary patterns
  lg_patterns <- rep(NA,n)
  lg_patterns[z >= 0 & W%*%z >=0] <- "High-High"
  lg_patterns[z <= 0 & W%*%z <=0] <- "Low-Low"
  lg_patterns[lg > mean_lg] <- "Negative"
  lg_patterns[lg < mean_lg & z*W%*%z < 0] <- "Other-Positive"
  lg_patterns[p_value > alpha] <- "Not Significant"
  
  
  # creating the correct palette for the map
  classes <- c("High-High", "Low-Low", "Negative", "Other-Positive", "Not Significant")
  colors <- c("#DE2D26","#FCBBA1","#9ECAE1", "#FEE5D9","#D3D3D3")
  pal <- match_palette(lg_patterns,classes,colors)
  
  
  polys <- polys %>% mutate(lg_patterns)
  
  # Making the local geary map
  tm_shape(polys) +
    tm_fill("lg_patterns", palette = pal)
}
  


###################################################################################################
# Getis-Ord
###################################################################################################


g_map <- function(polys,vname,weights = NULL, permutations = 999, alpha = .05,type = "g"){
  
  if (permutations > 99999){
    stop("maximum number of permutations is 99999")
  }
  
  if (alpha > 1 | alpha < 0){
    stop("alpha must be between 0 and 1")
  }
  
  if (any(class(polys) == "sf")){
    if (!any(class(st_geometry(polys)) == "sfc_MULTIPOLYGON")){
      stop("geomtry type must be sfc_MULTIPOLYGON for sf dataframes")
    } 
  } else {
    stop("sf is only supported geometry")
  }
  
  # creating weights if left to default 
  if (is.null(weights)){
    sgbp <- st_queen(polys)
    nb <- as_nb_sgbp(sgbp)
    weights <- nb2listw(queen.nb,style = "B", zero.policy = TRUE)
  }
  
  
  # converting weights to full spatial matrix
  W <- convert_matrix(weights)
  
  #extracting x variable from sf dataframe
  x <- get_var(vname,polys)
  
  
  #computing observed statistics and reference distributions for each location
  observed <- local_g(x,W,type = type)
  sims <- local_g_sims(x,W,permutations,type = type)
  
  
  #computing p_values from the observed statistics and reference distributions
  p_value <- get_p_value(sims,observed,type = "two-sided")
  
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
  
  polys <- polys %>% mutate(g_patterns)
  
  #Making the map
  tm_shape(polys) +
    tm_fill("g_patterns", palette = pal)
}




###################################################################################################
# Local Join Count
###################################################################################################


jc_map <- function(polys, vname, weights = NULL, permutations = 999, alpha = .05){
  
  if (permutations > 99999){
    stop("maximum number of permutations is 99999")
  }
  
  if (alpha > 1 | alpha < 0){
    stop("alpha must be between 0 and 1")
  }
  
  
  if (any(class(polys) == "sf")){
    if (!any(class(st_geometry(polys)) == "sfc_MULTIPOLYGON")){
      stop("geomtry type must be sfc_MULTIPOLYGON for sf dataframes")
    } 
  } else {
    stop("sf is only supported geometry")
  }
  
  # creating weights if left to default 
  if (is.null(weights)){
    sgbp <- st_queen(polys)
    nb <- as_nb_sgbp(sgbp)
    weights <- nb2listw(queen.nb,style = "B", zero.policy = TRUE)
  }
  
  # converting weights to full spatial matrix
  W <- convert_matrix(weights)
  
  #extracting x variable from sf dataframe
  x <- get_var(vname,polys)
  
  if (!(is.binary(x))){
    stop("vanme must designate a binary variable")
  }
  
  
  #computing observed statistics and reference distributions for each location
  observed <- local_jc(x,W)
  sims <- local_jc_sims(x,W,permutations)
  
  
  #computing p_values from the observed statistics and reference distributions
  p_value <- get_p_value(sims,observed,type = "one-sided")
  
  n <- nrow(W)
  unique_value <- rep(0,n)
  unique_value[p_value <= alpha] <- 1
  
  
  polys <- polys %>% mutate(unique_value)
  
  classes <- c(0,1)
  colors <- c("white","blue")
  pal <- match_palette(unique_value, classes,colors)
  
  polys <- polys %>% mutate(unique_value)
  
  #Making the map
  tm_shape(polys) +
    tm_fill("unique_value", palette = pal)
}




###################################################################################################
# Significance Map
###################################################################################################



significance_map <- function(polys,
                             xname, 
                             type,
                             weights = NULL,
                             alpha = .05,
                             permutations = 999){
  
  if (permutations > 99999){
    stop("maximum number of permutations is 99999")
  }
  
  if (alpha > 1 | alpha < 0){
    stop("alpha must be between 0 and 1")
  }
  
  
  # testing data type of polys
  if (any(class(polys) == "sf")){
    if (!any(class(st_geometry(polys)) == "sfc_MULTIPOLYGON")){
      stop("geomtry type must be sfc_MULTIPOLYGON for sf dataframes")
    } 
  } else {
    stop("sf is only supported geometry")
  }
  
  
  #checking for supported types
  supported_types <- c("moran","geary","g","gstar","join_count")
  if (!(type %in% supported_types)){
    stop("type is not support, choose from: moran, geary, g , gstar, join_count")
  }
 
  #making sure custom weights are supported
  if (!(is.null(weights))){
    if (weights[[1]] != "B"){
      stop("Weights must be style B from spdep")
    }
  }
  
  
  # creating weights if left to default 
  if (is.null(weights)){
    sgbp <- st_queen(polys)
    nb <- as_nb_sgbp(sgbp)
    weights <- nb2listw(queen.nb,style = "B", zero.policy = TRUE)
  }
  
  # converting weights to full spatial matrix
  W <- convert_matrix(weights)
  
  #extracting x variable from sf dataframe
  x <- get_var(vname,polys)
  
  
  #computing moran observed statistics and sims
  if (type == "moran"){
    observed <- local_moran(x,W)
    sims <- local_moran_sims(x,W,permutations)
  }
  
  #computing geary observed statistics and sims
  if (type == "geary"){
    observed <- local_geary(x,W)
    sims <- local_geary_sims(x,W,permutations)
  }
  
  #computing g observed statistics and sims
  if (type == "g"){
    observed <- local_g(x,W,type = "g")
    sims <- local_g_sims(x,W,permutations, type = "g")
  }
  
  #computing gstar observed statistics and sims
  if (type == "gstar"){
    observed <- local_g(x,W,type = "gstar")
    sims <- local_g_sims(x,W,permutations, type = "gstar")
  }
  
  #computing join count observed statistics and sims
  if (type == "join_count"){
    if (is.binary(x)){
      observed <- local_jc(x,W)
      sims <- local_jc_sims(x,W,permutations)
    } else {
      stop("join_count is only used with binary variables")
    }
  }
  
  #computing moran observed statistics and sims
  if (type == "join_count"){
    tail <- "one-sided"
  } else {
    tail <- "two-sided"
  }
  #Computing pvalues
  p_value <- get_p_value(sims,observed,type = tail)
  
  
  #Creating breaks based on p-values
  target_p <- 1 / (1 + permutations)
  potential_brks <- c(.00001, .0001, .001, .01)
  brks <- potential_brks[which(potential_brks > target_p & potential_brks < alpha)]
  brks2 <- c(target_p, brks, alpha)
  labels <- c(as.character(brks2), "Not Significant")
  brks3 <- c(-inf, brks2, 1)
  
  cuts <- cut(p_value, breaks = brks3,labels = labels)
  
  # Adding the p-value significance breaks to the data frame
  polys <- polys %>% mutate(sig = cuts)
  
  
  # Constructing the correct palette
  pal <- rev(brewer.pal(length(labels), "Greens"))
  pal[length(pal)] <- "#D3D3D3"
  
  
  # Making the map
  tm_shape(polys) +
    tm_fill("sig", palette = pal)
}



bivariate_significance_map <- function(polys,
                                       xname,
                                       yname,
                                       type,
                                       weights = NULL,
                                       alpha = .05,
                                       permutations = 999){
  
  if (permutations > 99999){
    stop("maximum number of permutations is 99999")
  }
  
  if (alpha > 1 | alpha < 0){
    stop("alpha must be between 0 and 1")
  }
  
  
  # testing data type of polys
  if (any(class(polys) == "sf")){
    if (!any(class(st_geometry(polys)) == "sfc_MULTIPOLYGON")){
      stop("geomtry type must be sfc_MULTIPOLYGON for sf dataframes")
    } 
  } else {
    stop("sf is only supported geometry")
  }
  
  #checking for supported types
  supported_types <- c("moran","geary","join_count")
  if (!(type %in% supported_types)){
    stop("type is not support, choose from: moran, geary, join_count")
  }
  
  
  
  #making sure custom weights are supported
  if (!(is.null(weights))){
    if (weights[[1]] != "B"){
      stop("Weights must be style B from spdep")
    }
  }
  
  
  # creating weights if left to default 
  if (is.null(weights)){
    sgbp <- st_queen(polys)
    nb <- as_nb_sgbp(sgbp)
    weights <- nb2listw(queen.nb,style = "B", zero.policy = TRUE)
  }
  
  
  # converting weights to full spatial matrix
  W <- convert_matrix(weights)
  
  #extracting x variable from sf dataframe
  x <- get_var(xname,polys)
  y <- get_var(yname,polys)
  
  
  
}
























