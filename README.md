<!-- README.md is generated from README.Rmd. Please edit that file -->
iQRAT
==========

The R package *iQRAT* provides tools for Integrated Quantile Rank Test for group-wise joint effect of rare and common variants in sequencing study.

Installation
------------

You can install the **development** version from [Github](https://github.com/tianyingw/iQRAT).

``` r
# install.packages("devtools")
devtools::install_github("tianyingw/iQRAT")
```

Usage
-----
Suppose you have a sample dataset stored in a list with Y = SampleData$y, covriates C = SampleData$c, genetic data X = SampleData$x.

``` r
library(quantreg)
library(SKAT)
library(iQRAT)
data("SampleData")

# Step 1: fit null model
null.fit = Null_model(Y = SampleData$y, C = SampleData$c)
# Step 2: run the test, p value will return

# SKAT version iQRAT
test.iQRAT1 = iQRAT(X = SampleData$x, C = SampleData$c, v =  null.fit, method.type = "S")

# If you want to specify weights
w = dbeta(colMeans(SampleData$x)/2,0.5,0.5) # we use beta density as an example
test.iQRAT2 = iQRAT(X = SampleData$x, C = SampleData$c, v =  null.fit, method.type = "S", w = w)

# Burden version iQRAT
test.iQRAT3 = iQRAT(X = SampleData$x, C = SampleData$c, v =  null.fit, method.type = "B")



#Y is trait, X is genotype in a region, e.g., a gene, and C is covriates such as gender, race, etc.
iQRAT(Y, X, C) 
```

License
-------

This package is free and open source software, licensed under GPL (&gt;=2).
