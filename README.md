# spatmap

**spatmap** is an R package for the visualization and computation of local spatial statistics. The
mapping functions support a variety of local spatial statistics with cluster maps and 
significance maps. These include Local Moran, Local Geary, Local Join Count, Local G, and Local G*.
There are bivariate options for the Local Moran, Local Geary, and Local Join Count statistics.
All of these mapping functions have many formatting options through the **tmap** package.

## Local Statistics Maps

### Local Geary Cluster Maps

![Local Geary Map](/images/geary.png)

### Local Moran Cluster Maps

![Local Moran Map]("https://github.com/morrisonge/spatmap/images/moran.png")

### Local Join Count Maps

### Local G and G* Cluster Maps

![Local G Map]("https://github.com/morrisonge/spatmap/images/g.png")

### Significance Maps for each local statistic

![Moran Cluster and Significance Map]("https://github.com/morrisonge/spatmap/images/moranandsig.png")

### Bivariate options for Geary, Moran, and Join Count

## Installation

I am working to make **spatmap** available through CRAN, but for now, you will need to 
remote install the package from github. This is done done with `remotes::install_github`

```r
remotes::install_github("morrisonge/spatmap")
```



## Tutorials 












