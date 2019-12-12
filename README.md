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

``` r
library(iQRAT)

#Y is trait, X is genotype in a region, e.g., a gene, and C is covriates such as gender, race, etc.
iQRAT(Y, X, C) 
```

License
-------

This package is free and open source software, licensed under GPL (&gt;=2).
