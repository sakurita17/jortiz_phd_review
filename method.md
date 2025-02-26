Numerator relationship matrix
================
2025-02-24

# Objective

This exercise want to recreate the concept and computation around the
Numerator relationship matrix and all three methods to its estimation.

# Library

For the aim of this exercise, I will use the library “pedigreeTools”:
<https://github.com/Rpedigree/pedigreeTools>

``` r
library(package = "pedigreeTools") # R package to create and create different pedigree
                                   # features
```

# 1) Create a pedigree

``` r
# Ax = b

# 1) Create a pedigree

ped = data.frame(iid = 1:7,
                 fid = c(0, 0, 1, 1, 3, 1, 5), 
                 mid = c(0, 0, 0, 2, 4, 4, 6))

(ped2 = pedigree(sire  = ped$fid, 
                dam   = ped$mid,
                label = ped$iid))
```

    ##   sire  dam
    ## 1 <NA> <NA>
    ## 2 <NA> <NA>
    ## 3    1 <NA>
    ## 4    1    2
    ## 5    3    4
    ## 6    1    4
    ## 7    5    6

# 2) Tabular method

If both parents $s$ and $d$ of animal $i$ are known:

- The diagonal elements $(a)_{ii}$ correspond to:
  $a_{ii} = 1 + F_i = 1 + 0.5(a_{sd})$

- The off-diagonal element $a_{ji}$ correspond to:
  $a_{ji} = a_{ij} = 0.5(a_{js} + a_{jd})$

If only one parent $s$ or $d$ of animal $i$ is known:

- The diagonal elements $(a)_{ii}$ correspond to: $a_{ii} = 1$

- The off-diagonal element $a_{ji}$ correspond to:
  $a_{ji} = a_{ij} = 0.5(a_{js/jd})$

If both parents are unknown and are assumed unrelated:

- The diagonal elements $(a)_{ii}$ correspond to: $a_{ii} = 1$

- The off-diagonal element $a_{ji}$ correspond to: $a_{ji} = a_{ij} = 0$

``` r
(A = getA(ped2))
```

    ## 'as(<dtTMatrix>, "dtCMatrix")' is deprecated.
    ## Use 'as(., "CsparseMatrix")' instead.
    ## See help("Deprecated") and help("Matrix-deprecated").

    ## 7 x 7 sparse Matrix of class "dsCMatrix"
    ##       1    2     3      4       5       6       7
    ## 1 1.000 .    0.500 0.5000 0.50000 0.75000 0.62500
    ## 2 .     1.00 .     0.5000 0.25000 0.25000 0.25000
    ## 3 0.500 .    1.000 0.2500 0.62500 0.37500 0.50000
    ## 4 0.500 0.50 0.250 1.0000 0.62500 0.75000 0.68750
    ## 5 0.500 0.25 0.625 0.6250 1.12500 0.56250 0.84375
    ## 6 0.750 0.25 0.375 0.7500 0.56250 1.25000 0.90625
    ## 7 0.625 0.25 0.500 0.6875 0.84375 0.90625 1.28125

# 3) Thompson’s method

The relationship matrix can be expressed as $A = TDT^T$, where $T$ is a
lower triangular matrix $D$ is a diagonal matrix.

Consider the relationshio with mendelian sampling

The matrix $T$ traces the flow of genes from one generation to the
other; in other words, it accounts only for direct (parent-offspring)
relationships.

Rules for $i_{th}$:

- $t_{ii} = 1$
- $t_{ij} = 0.5(t_{sj} + t_{dj})$ \# If both parents (s and d) are known
- $t_{ij} = 0.5(t_{sj})$ \# If only want parent (s) is known
- $t_{ij} = 0$ \# If neither parent is known

``` r
(T = getT(ped2))
```

    ## 7 x 7 sparse Matrix of class "dtCMatrix"
    ##       1    2    3   4   5   6 7
    ## 1 1.000 .    .    .   .   .   .
    ## 2 .     1.00 .    .   .   .   .
    ## 3 0.500 .    1.00 .   .   .   .
    ## 4 0.500 0.50 .    1.0 .   .   .
    ## 5 0.500 0.25 0.50 0.5 1.0 .   .
    ## 6 0.750 0.25 .    0.5 .   1.0 .
    ## 7 0.625 0.25 0.25 0.5 0.5 0.5 1

The diagonal matrix $D$ is the variance and covariance matrix for
Mendelian sampling.

``` r
(D = getD(ped2))
```

    ##       1       2       3       4       5       6       7 
    ## 1.00000 1.00000 0.75000 0.50000 0.50000 0.50000 0.40625
