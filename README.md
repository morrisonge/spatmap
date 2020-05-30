# spatmap

### Author: Grant Morrison(morrisonge7@gmail.com)

**spatmap** is an R package for the visualization of local spatial statistics. The
mapping functions support a variety of local spatial statistics with cluster maps and 
significance maps. These include Local Moran, Local Geary, Local Join Count, Local G, and Local G*.
There are multivariate options for the Local Geary and Local Join Count statistics. The statistical
computation is done by the package rgeoda by Xun Li. The visualization of these statistics is built
off of **tmap** to give a range range of formatting options and interactivity.


## Installation

Both **spatmap** and **rgeoda** are not available through CRAN, so you must install them remotely. We use 
`install_github` from **remotes** to install both of these packages.

```r
remotes::install_github("lixun910/rgeoda")
remotes::install_github("morrisonge/spatmap")
```


## Goals and Objectives

- Visualization of Local Moran Cluster and Significance Maps

- Visualization of Univariate and Multivariate Local Geary Cluster and Significance Maps

- Visualization of Local G and G* Cluster and Significance Maps

- Visualization of Univariate and Multivariate Local Join Count Cluster and Significance Maps


## Data Description

The data used in the vignette for this package is from the GeoDa website and will be loaded from 
the **geodaData** package by Angela Li. The dataset is originall from a classic social science 
study in the 1800's, and is now used as apart of the GeoDa workbook tutorials. To get access 
to this data, we will need to install **geodaData** from github.

```r
remotes::install_github("spatialanalysis/geodaData")
```

The dataset is from a classic social science study by Andre-Michel Guerry on crime, suicide ,
and other “moral” statistics. It is available on the GeoDa website. The spatial resolution is the
prefecture-level, which is similar to the county-level. Each polygon has measures for different
“moral” statistics of the area. There is some degree of temporal resolution with averages of the
variables pertaining to different sets of years. These years range from 1815 to 1834. There are
23 variables in the dataset, and 85 observations. I will only be using a few of these variables
for the project vignette. 

Variables included:

- **Donatns**: Donations to the poor

- **Infants**: Population per illegitimate birth

- **doncat**: Categorical Variable based on **Donatns**. The top quintile is assigned 1 and the rest 0

- **doncat**: Categorical Variable based on **Infants**. The top quintile is assigned 1 and the rest 0


## Local Statistics Maps


### Local Moran Cluster Map

![Local Moran Map](/images/moran.png)

### Local Geary Cluster Map

![Local Geary Map](/images/geary.png)

### Multivariate Local Geary Cluster Map

![Multivariate Local Geary Map](/images/multi_geary.png)

### Local Join Count Cluster Map

![Local Join Count Map](/images/joincount.png)

### Multivariate Local Join Count Cluster Map

![Multivariate Local Join Count Map](/images/multi_joincount.png)

### Local G Cluster Map

![Local G Map](/images/g.png)

### Significance Map 

![Significance Map](/images/moran_sig.png)




## Tutorials 












