BoostLines package
================
Philipp Hunziker
March 20, 2017

Fast geometric operations on lines using [Boost Geometries](http://www.boost.org/doc/libs/1_63_0/libs/geometry/doc/html/index.html).

<span style="color:red">CAUTION</span>: This package is in development - its interface may change.

Installation
------------

You can install the package directly from this repository:

``` r
library(devtools)
install_github("hunzikp/BoostLines")
```

Usage
-----

### Boost/Unboost

``` r
library(BoostLines)
library(sp)

## Make some lines
set.seed(0)
N <- 250
lines.ls <- vector('list', N)
for (i in 1:N) {
  strt <- runif(1)
  range <- runif(1)
  coords <- matrix(runif(4, strt, strt+range), 2, 2)
  lines.ls[[i]] <- Lines(slinelist = list(Line(coords = coords)), ID = i)
}
test.sl <- SpatialLines(lines.ls)

## Make BoostLines object
test.bl <- BoostLines(test.sl)

## Cast back to SpatialLines
test_cast.sl <- unboost(test.bl) 
```

### Noding

``` r
## Node
out.ls <- bNode(test.bl)

## Unpack & plot result
noded.sl <- unboost(out.ls[[1]])
noded.ids <- out.ls[[2]]  # integer vector of original IDs
plot(noded.sl, col=sample(colors(), length(noded.sl), replace=TRUE))
```

![](README_files/figure-markdown_github/node-1.png)

### Intersection

``` r
## Intersect query
lmat <- bIntersects(test.bl)

## Compare with gIntersects
library(rgeos)
```

    ## rgeos version: 0.3-21, (SVN revision 540)
    ##  GEOS runtime version: 3.5.0-CAPI-1.9.0 r4084 
    ##  Linking to sp version: 1.2-3 
    ##  Polygon checking: TRUE

``` r
lmat_geos <- gIntersects(test.sl, byid=TRUE)
all(lmat == lmat_geos)
```

    ## [1] TRUE

### Minimal Distance

``` r
## Cal Distances
dmat <- bDistance(test.bl)

## Compare with gDistance
dmat_geos <- gDistance(test.sl, byid=TRUE)
summary(abs(as.vector(dmat - dmat_geos)))
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## 0.000e+00 1.102e-09 2.617e-09 3.023e-09 4.536e-09 1.310e-08
