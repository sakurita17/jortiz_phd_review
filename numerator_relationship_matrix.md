Numerator relationship matrix
================
2025-02-24

# Introduction

Additive relationships are a measure of the proportion of genes, which are identical
by descent, which are expected to be shared by two animals. The idea of tracing 
paths to establish the relationship among animals was introduced by Wright (1921), 
however was Malecot (1948) how defined relationship based on probabilities of 
individual genes at a locus being identical by descent.

Based on the path analysis of phenotypic values for pedigree members, Wrigth (1922)
defined relationship coefficient {$R$} between individuals {$i$} and {$j$} as the
correlation between their genetic values:

{$R_{i,j} = Cor(g_i,g_j)$}
{$Cov(g_i,g_j)/sqr(Var(g_i)Var(g_j))$}

The Numerator relationship matrix $(A)$ describes the additive relationship among 
individuals. $(A)$ is a symmetric and its diagonal elements represent twice 
the probability that two gametes taken at random from animal {$i$} will carry 
identical alleles by descent, in other words, the inbreeding coefficient of 
animal {$i$} (Wright, 1922).

The off-diagonal elements, {$a_{ij}$}, equals the numerator of the coefficient 
of relationship between animals {$i$} and {$j$}. When multiplied by the additive
genetic variance {$\sigma_u$}, {$A_\sigma_u$} is the covariance among breeding
values. {$A$} is twice the coefficient of coancestry (kinship).

The inverse of the numerator relationship matrix is used for breeding values 
computation, however as the size of the matrix increase we start to face compu-
tational barrers, memory is a bottleneck for large pedigree but attemps to reduce
memory usage increased the computational time. 

# Objective

To understand the numerator relationship matrix and method for its computation


# Library

For the aim of this exercise, I will use the library “pedigreeTools”:
<https://github.com/Rpedigree/pedigreeTools>


# Packages 

library(package = "pedigreeTools") 


# 1) Create a pedigree

``` r
ped = data.frame(iid = 1:7,
                 fid = c(0, 0, 1, 1, 3, 1, 5), 
                 mid = c(0, 0, 0, 2, 4, 4, 6))

(ped2 = pedigree(sire  = ped$fid, 
                 dam   = ped$mid,
                 label = ped$iid))
```
                
# 2) Tabular method

Applying the tabular method require to order the pedigree ensuring parents appear
before their progeny.Thus, the matrix {$A$} is formed following the nextrecursive 
method: 

If both parents of {$i^th$} individual are known, say $s$ and $d$

- The diagonal element:
  $a_{ii} = 1 + F_i = 1 + 0.5(a_{sd})$

- The off-diagonal element:
  $a_{ij} = a_{ji} = 0.5(a_{js} + a_{jd})$

If only one parent $s$ or $d$ is known:

- The diagonal element: 
  $a_{ii} = 1$

- The off-diagonal element:
  $a_{ij} = a_{ij} = 0.5(a_{js/jd})$

If neither parent is known:

- The diagonal element: 
  $a_{ii} = 1$

- The off-diagonal element:
  $a_{ij} = a_{ji} = 0$

Then, from the example pedigree from above, A is:

(A = getA(ped2))


# 3) Thompson’s method

The relationship matrix can be expressed as $A = TDT^T$, where $T$ is a
lower triangular matrix $D$ is a diagonal matrix.

Consider the relationship with mendelian sampling

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
