# helper functions



#used to assess whether a variable is binary or not for local joint count statisttics
is.binary <- function(x){
  vals <- unique(x)
  len <- length(vals)
  if(len > 2){
    return(FALSE)
  } else {
    if(min(x) == 0 || max(x) == 1){
      if(len == 1){
        return(TRUE)
      } else{
        if (min(x) == 0 && max(x) == 1){
          return(TRUE)
        } else {
          return(FALSE)
        }
      }
    } else {
      return(FALSE)
    }
  }

}



#function used to extract variable from sf dataframe with the variable name
get_var <- function(vname,df) {
  v <- df[vname] %>% st_set_geometry(NULL)
  v <- unname(v[,1])
  return(v)
}



#function used to compute pvalues based reference statistics and observed statistics
get_pvalue <- function(mat, observed,type = "one-sided"){
  nperm <- ncol(mat)
  nlocs <- nrow(mat)
  p_value <- rep(NA,nlocs)
  for(i in 1:nlocs){
    num_greater <- length(which(mat[i,] >= observed[i]))
    p_value[i] <- (num_greater + 1) / (nperm + 1)
  }
  if (type == "two-sided"){
    p_value <- ifelse(p_value > .5, 1-p_value, p_value)
    p_value
  } else {
    p_value
  }
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


#function used to assess whether the weights are of the proper form
check_weights <- function(weights){

  #making sure custom weights are supported
  if (!(is.null(weights))){
    if (weights[[1]] != "B"){
      stop("Weights must be style B from spdep")
    }
  }

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
check_parameters <- function(polys,permutations,alpha,weights){

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

  check_weights(weights)

}




#creates neighbors like structure for sf polygons
st_queen <- function(a, b = a) st_relate(a, b, pattern = "F***T****")

#converts the sf neighbors structure to spdep
as_nb_sgbp <- function(x, ...) {
  attrs <- attributes(x)
  x <- lapply(x, function(i) { if(length(i) == 0L) 0L else i } )
  attributes(x) <- attrs
  class(x) <- "nb"
  return(x)
}

#function to convert the matrix into a full size spatial weights matrix
convert_matrix <- function(weights){
  W  <- as(weights, "symmetricMatrix")
  W  <- as.matrix(W/Matrix::rowSums(W))
  W[which(is.na(W))] <- 0
  return(W)
}


#function to create default weights in the mapping functions when custom weights are not used
default_weights <- function(polys){
  sgbp <- st_queen(polys)
  nb <- as_nb_sgbp(sgbp)
  weights <- nb2listw(nb,style = "B", zero.policy = TRUE)
  return(weights)
}





