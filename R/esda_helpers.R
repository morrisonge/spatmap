# helper functions




#function used to extract variable from sf dataframe with the variable name
get_var <- function(vname,df) {
  v <- df[vname] %>% st_set_geometry(NULL)
  v <- unname(v[,1])
  return(v)
}




#function used to contruct the correct palette for a set of classification
match_palette <- function(patterns, classifications, colors){
  classes_present <- base::unique(patterns)
  mat <- matrix(c(classifications,colors), ncol = 2)
  logi <- classifications %in% classes_present
  pre_col <- matrix(mat[logi], ncol = 2)
  pal <- pre_col[,2]
  return(pal)
}





# functions used to make sure the permutations parameter is valid
check_permutations <- function(permutations){
  if (permutations > 99999){
    stop("maximum number of permutations is 99999")
  }
  if (permutations < 1){
    stop("permutations must be greater than 1")
  }
}




#function to assess to validity of parameters used within the mapping functions
check_parameters <- function(polys,permutations,alpha){

  #setting maximum number of permutations
  check_permutations(permutations)


  #making sure alpha is in correct range
  if (alpha > 1 | alpha < 0){
    stop("alpha must be between 0 and 1")
  }


  min_sig <- 1 / (1 + permutations)
  if (alpha < min_sig){
    stop("More permutations are needed to have this level of alpha")
  }


  # testing data type of polys
  if (any(class(polys) == "sf")){
    if (!any(class(st_geometry(polys)) == "sfc_MULTIPOLYGON")){
      stop("geomtry type must be sfc_MULTIPOLYGON for sf dataframes")
    }
  } else {
    stop("sf is only supported geometry")
  }


}










