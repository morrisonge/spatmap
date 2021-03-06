% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/joincount.R
\name{local_jc_pvalue}
\alias{local_jc_pvalue}
\title{Local Join Count P-values}
\usage{
local_jc_pvalue(x, y = NULL, weights, permutations = 999)
}
\arguments{
\item{x}{A vector of numerical values}

\item{weights}{Weights structure from spdep, must be style "B"}

\item{permutations}{Number of permutations, the default is 999}

\item{z}{A vector of numerical values, only to be used in the bivariate case,default is NULL}
}
\value{
pvalue A vector of p-values for each locations local geary statistic
}
\description{
Function to compute the p-value of the observed local join count statistics under the
conditional randomizaton approach
}
