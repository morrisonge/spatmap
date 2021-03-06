% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/g.R
\name{local_g_pvalue}
\alias{local_g_pvalue}
\title{Local G P-values}
\usage{
local_g_pvalue(x, weights, permutations = 999, type = "g")
}
\arguments{
\item{x}{A vector of numerical values}

\item{weights}{Weights structure from spdep, must be style "B"}

\item{permutations}{Number of permutations, the default is 999}

\item{type}{String, type of Local G statistic, options are "g" or "gstar"}
}
\value{
pvalue A vector of p-values for each locations local geary statistic
}
\description{
Function to compute the p-value of the observed local G and G* statistics under the
conditional randomizaton approach
}
