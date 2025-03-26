Numerator Relationship Matrix
================

================ 2025-02-24

# Load R packages

``` r
library(package = "pedigreeTools") 
library(package = "Matrix")
```

# Introduction

The additive relationships measure the proportion of genes that are
identical by descendant and are expected to be shared by two animals
(Wright, 1921). However, Malecot (1948) was a question of how defined
relationship-based probabilities of individual genes at a locus being
identical by descent.

Based on the path analysis of the phenotypic values of the members of
the pedigree, Wright (1922) defined the relationship coefficient $R$
between individual $i$ and $j$ as the correlation between their genetic
value (Gregor notes).

$$
\begin{aligned}
R_{i,j} &= \text{Cor}(g_i,g_j) \\
        &= \frac{\text{Cov}(g_i,g_j)}{\sqrt{\text{Var}(g_i) \cdot \text{Var}(g_j)}}
\end{aligned}
$$

The numerator relationship matrix $A$ describes the expected additive
relationship between individuals, before observing any phenotypic data,
it serves as the prior covariance matrix of additive genetic
values.Thus, The additive relationship $a_{ij}$ is used as a measure of
the covariance of breeding values between relative (gg).

The off-diagonal elements, $a_{ij}$, equals the numerator of the
coefficient of relationship between animals $i$ and $j$. When multiplied
by the additive genetic variance $\sigma^{2}_a$, $A\sigma^{2}_a$ is the
covariance among breeding values. The numerator relationships are equal
twice the coancestry, and they express the proportion of additive
genetic variance that two individuals have in common.

The coefficient of coancestry of two individuals reflects the
probability that two gametes taken at random, one from each, carry
alleles that are identical by descent (=inbreedig coefficient of their
progeny should they be mated together)

The diagonals $(a_{ii})$ contain the coefficient of inbreeding. That is
the probability that two genes at any locus in an individual are
identical by descent.

For the aim of this exercise, I will use the library “pedigreeTools”:
<https://github.com/Rpedigree/pedigreeTools>

# 1) Create a pedigree

First we will create a pedigree

``` r
ped = data.frame(ID = 1:6,
                 fid = c(0, 0, 1, 1, 4, 5), 
                 mid = c(0, 0, 2, 0, 3, 2))

ped2 = pedigree(sire  = ped$fid, 
                 dam   = ped$mid,
                 label = ped$ID)
print(ped2)
```

    ##   sire  dam
    ## 1 <NA> <NA>
    ## 2 <NA> <NA>
    ## 3    1    2
    ## 4    1 <NA>
    ## 5    4    3
    ## 6    5    2

# Construction of the relationship matrix

Systematic recurrent rules that are based on the flow of genes from
generation to generation with individual animals being specified.

- Path coefficient method: Suitable for small pedigrees with few
  generations and little inbreeding.
- Recursively using the tabular method:
- xxx

# 2) Tabular method

Method described by Henderson (1976). Initially, the pedigree is ordered
chronologically so that parents precede offspring. Base parents are
considered unrelated and non-inbred. Then, working one row at a time,
compute elements of $A$ using the following relationships:

$$
\[
\begin{array}{|l|c|c|}
\hline
\textbf{Condition} & \textbf{Diagonals} & \textbf{Off-diagonals} \\ 
\hline
\multirow{2}{*}{\parbox{4cm}{If both parents are\\ known}} 
 & a_{ii} = 1 + F_i = 1 + 0.5(a_{sd}) & a_{ij \vee ji} = 0.5(a_{js} + a_{jd}) \\ 
\hline
\multirow{2}{*}{\parbox{4cm}{If only one parent\\ is known}} 
 & a_{ii} = 1 & a_{ij \vee ji} = 0.5(a_{js}) \\ 
\hline
\multirow{2}{*}{\parbox{4cm}{If neither parent\\ is known}} 
 & a_{ii} = 1 & a_{ij \vee ji} = 0 \\ 
\hline
\end{array}
\]
$$ We can manually compute $A$ following the recursive method proposed
by Henderson (1976). To show how this work I used the function created
by Nilforooshan, M. A. (2021) in
<https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2021.655638/full>

``` r
tabularA <- function(ped) {
    ped$p1 = apply(ped[,2:3], 1, FUN=min)
    ped$p2 = apply(ped[,2:3], 1, FUN=max)
    ped = ped[,c("ID","p1","p2")]
    A = diag(nrow(ped))
    for(i in which(ped$p2 >0))
    {
        p1 = ped[i,]$p1
        p2 = ped[i,]$p2
        if(p1==0) {
            A[1:(i-1),i] = A[1:(i-1),p2]/2
        } else {
            A[1:(i-1),i] = (A[1:(i-1),p1] + A[1:(i-1),p2])/2
            A[i,i] = 1 + A[p1,p2]/2
        }
        A[i,1:(i-1)] = A[1:(i-1),i]
    }
    return(A)
}
```

``` r
tabularA(ped)
```

    ##      [,1]  [,2]   [,3]   [,4]   [,5]   [,6]
    ## [1,] 1.00 0.000 0.5000 0.5000 0.5000 0.2500
    ## [2,] 0.00 1.000 0.5000 0.0000 0.2500 0.6250
    ## [3,] 0.50 0.500 1.0000 0.2500 0.6250 0.5625
    ## [4,] 0.50 0.000 0.2500 1.0000 0.6250 0.3125
    ## [5,] 0.50 0.250 0.6250 0.6250 1.1250 0.6875
    ## [6,] 0.25 0.625 0.5625 0.3125 0.6875 1.1250

# 3) Generalized cholesky decomposition

The numerator relationship matrix $A$ can be expressed as $TDT^T$
(Thompson, 1977), where $T$ is a lower triangular matrix and $D$ is a
diagonal matrix.

## Matrix T

The matrix $T$ traces the flow of genes from one generation to the
other. For example, for our small pedigree, we know that the expected
and variance of the breeding values are as follow:

$$
\[
\begin{array}{|c|p{8cm}|}
\hline
\textbf{Marginal} & \textbf{Conditional} \\
\hline
a_1 \sim N(0, \sigma^2_{a}) & - \\
\hline
a_2 \sim N(0, \sigma^2_{a}) & - \\
\hline
a_3 \sim N(0, \sigma^2_{a}) & \(a_3 \mid a_1,a_2 \sim N\left(\frac{1}{2} (a_1 + a_2), \frac{1}{2} \sigma^2_{a} \right)\) \\
\hline
a_4 \sim N(0, \sigma^2_{a}) & \(a_4 \mid a_1 \sim N\left(\frac{1}{2} a_1, \frac{1}{2} \sigma^2_{a} \right)\) \\
\hline
a_5 \sim N(0, \sigma^2_{a}) & \(a_5 \mid a_4,a_3 \sim N\left(\frac{1}{2} (a_4 + a_3), \frac{1}{2} \sigma^2_{a} \right)\) \\
\hline
a_6 \sim N(0, \sigma^2_{a}) & \(a_6 \mid a_5,a_2 \sim N\left(\frac{1}{2} (a_5 + a_2), \frac{1}{2} \sigma^2_{a} \right)\) \\
\hline
\end{array}
\]
$$

Now, we define our system of equations:

$$
\begin{aligned}
a_1 &= r_1 \\
a_2 &= r_2 \\
a_3 &= \frac{1}{2} a_1 + \frac{1}{2} a_2 + r_3 \\
a_4 &= \frac{1}{2} a_1 + r_4 \\
a_5 &= \frac{1}{2} a_4 + \frac{1}{2} a_3 + r_5 \\
    &= \frac{1}{2} \left(\frac{1}{2} a_1 + r_4 \right) + \frac{1}{2} \left(\frac{1}{2} a_1 + \frac{1}{2} a_2 + r_3 \right) + r_5 \\
    &= \frac{1}{4} a_1 + \frac{1}{2} r_4 + \frac{1}{4} a_1 + \frac{1}{4} a_2 + \frac{1}{2} r_3 + r_5 \\
    &= \frac{1}{2} a_1 + \frac{1}{4} a_2 + \frac{1}{2} r_3 + \frac{1}{2} r_4 + r_5 \\
a_6 &= \frac{1}{2} a_5 + \frac{1}{2} a_2 + r_6 \\
    &= \frac{1}{2} \left(\frac{1}{2} a_4 + \frac{1}{2} a_3 + r_5 \right) + \frac{1}{2} a_2 + r_6 \\
    &= \frac{1}{4} a_4 + \frac{1}{4} a_3 + \frac{1}{2} r_5 + \frac{1}{2} a_2 + r_6 \\
    &= \frac{1}{4} \left(\frac{1}{2} a_1 + r_4 \right) + \frac{1}{4} \left(\frac{1}{2} a_1 + \frac{1}{2} a_2 + r_3 \right) + \frac{1}{2} r_5 + \frac{1}{2} a_2 + r_6 \\
    &= \frac{1}{8} a_1 + \frac{1}{4} r_4 + \frac{1}{8} a_1 + \frac{1}{8} a_2 + \frac{1}{4} r_3 + \frac{1}{2} r_5 + \frac{1}{2} a_2 + r_6 \\
    &= \frac{1}{4} a_1 + \frac{5}{8} a_2 + \frac{1}{4} r_3 + \frac{1}{4} r_4 + \frac{1}{2} r_5 + r_6
\end{aligned}
$$

Matrix form:

$$
\begin{bmatrix}
a_1 \\
a_2 \\
a_3 \\
a_4 \\
a_5 \\
a_6 
\end{bmatrix}
\=
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
\frac{1}{2} & \frac{1}{2} & 1 & 0 & 0 & 0 \\
\frac{1}{2} & 0 & 0 & 1 & 0 & 0 \\
\frac{1}{2} & \frac{1}{4} & \frac{1}{2} & \frac{1}{2} & 1 & 0 \\
\frac{1}{8} & \frac{5}{8} & \frac{1}{4} & \frac{1}{4} & \frac{1}{2} & 1
\end{bmatrix}
\cdot
\begin{bmatrix}
r_1 \\
r_2 \\
r_3 \\
r_4 \\
r_5 \\
r_6 
\end{bmatrix}
$$

However, we can follow the next rules to get $T$ following the next
rules.

$$
\[
\begin{array}{|l|c|c|}
\hline
\textbf{Condition} & \textbf{Diagonals} & \textbf{Off-diagonals} \\ 
\hline
\multirow{2}{*}{If both parents are known} 
 & a_{ii} = 1 & t_{ij} = 0.5(t_{sj} + t_{dj}) \\ 
\hline
\multirow{2}{*}{If only one parent is known} 
 & a_{ii} = 1 & t_{ij} = 0.5(t_{sj}) \\ 
\hline
\multirow{2}{*}{If neither parent is known} 
 & a_{ii} = 1 & t_{ij} = 0 \\ 
\hline
\end{array}
\]
$$

After, applying the rules:

$$
\begin{aligned}
T_{2,1} &= 0 \\
T_{3,1} &= 0.5(T_{1,1} + T_{2,1}) = 0.5(1 + 0) = 0.5 \\
T_{3,2} &= 0.5(T_{1,2} + T_{2,2}) = 0.5(0 + 1) = 0.5 \\
T_{4,1} &= 0.5(T_{1,1}) = 0.5(1) = 0.5 \\
T_{4,2} &= 0.5(T_{2,1}) = 0.5(0) = 0 \\
T_{4,3} &= 0.5(T_{3,1}) = 0.5(0.5) = 0.25 \\
T_{5,1} &= 0.5(T_{4,1} + T_{3,1}) = 0.5(0.5 + 0.5) = 0.5 \\
T_{5,2} &= 0.5(T_{4,2} + T_{3,2}) = 0.5(0 + 0.5) = 0.25 \\
T_{5,3} &= 0.5(T_{4,3} + T_{3,3}) = 0.5(0 + 1) = 0.5 \\
T_{5,4} &= 0.5(T_{4,4} + T_{3,4}) = 0.5(1 + 0) = 0.5 \\
T_{6,1} &= 0.5(T_{5,1} + T_{2,1}) = 0.5(0.5 + 0) = 0.25 \\
T_{6,2} &= 0.5(T_{5,2} + T_{2,2}) = 0.5(0.25 + 1) = 0.625 \\
T_{6,3} &= 0.5(T_{5,3} + T_{2,3}) = 0.5(0.5 + 0) = 0.25 \\
T_{6,4} &= 0.5(T_{5,4} + T_{2,4}) = 0.5(0.5 + 0) = 0.25 \\
T_{6,5} &= 0.5(T_{5,5} + T_{2,5}) = 0.5(1 + 0) = 0.5 
\end{aligned}
$$

We use the function $getT$ from pedigreeTools to get the matrix $T$.
This function compute $T$ from its inverse.

``` r
T = getT(ped2)
```

    ## 'as(<dtTMatrix>, "dtCMatrix")' is deprecated.
    ## Use 'as(., "CsparseMatrix")' instead.
    ## See help("Deprecated") and help("Matrix-deprecated").

``` r
T
```

    ## 6 x 6 sparse Matrix of class "dtCMatrix"
    ##      1     2    3    4   5 6
    ## 1 1.00 .     .    .    .   .
    ## 2 .    1.000 .    .    .   .
    ## 3 0.50 0.500 1.00 .    .   .
    ## 4 0.50 .     .    1.00 .   .
    ## 5 0.50 0.250 0.50 0.50 1.0 .
    ## 6 0.25 0.625 0.25 0.25 0.5 1

## Matrix D

The diagonal matrix $D$ is the variance and covariance matrix for
Mendelian sampling. If there is no inbreeding, the values for each
diagonal are 1 if no parents are known, 3/4 if one of the parents is
known and 1/2 if both parents are known.

For example, from our small pedigree.. see the example ipad…

When inbreeding is taking into account, the diagonal elements of $D$ are
1 for pedigree founders, $0.5 - 0.25(F_s + Fd)$ if both parents of
animal $i$ are known and $0.75 - 0.25(F_s)$ if only one parent $s$ is
known.

For example,

\$\$ \$\$

or for simplicity use the function $getD$ from pedigreeTools

``` r
D = getD(ped2)
D
```

    ##       1       2       3       4       5       6 
    ## 1.00000 1.00000 0.50000 0.75000 0.50000 0.46875

Finally, we can manually compute $A$ as follow:

``` r
gCholA = T %*% diag(D) %*% t(T)
gCholA
```

    ## 6 x 6 Matrix of class "dgeMatrix"
    ##      1     2      3      4      5      6
    ## 1 1.00 0.000 0.5000 0.5000 0.5000 0.2500
    ## 2 0.00 1.000 0.5000 0.0000 0.2500 0.6250
    ## 3 0.50 0.500 1.0000 0.2500 0.6250 0.5625
    ## 4 0.50 0.000 0.2500 1.0000 0.6250 0.3125
    ## 5 0.50 0.250 0.6250 0.6250 1.1250 0.6875
    ## 6 0.25 0.625 0.5625 0.3125 0.6875 1.1250

# 4) Cholesky decomposition

As the recursive method for computing $A$ can be demanding, Henderson
(1975) proposed the Cholesky decomposition to build $A$. According to
theory a matrix, as $A$ is symmetric and positive defined, it can be
decomposed as:

$$
\begin{aligned}
A &= LL^T\\
\end{aligned}
$$

Where

$L$ is a lower triangular matrix, and $L'$ is its transpose

$L$ can be manually computed by following there some rules. However,
$L = T %*%sqrt(D)$, and this is the method implemented by the function
$getL$ from the R package $pedigreeTools$

``` r
L = getL(ped2)
L
```

    ## 6 x 6 sparse Matrix of class "dtCMatrix"
    ##   1 2         3         4         5         6
    ## 1 1 . 0.5000000 0.5000000 0.5000000 0.2500000
    ## 2 . 1 0.5000000 .         0.2500000 0.6250000
    ## 3 . . 0.7071068 .         0.3535534 0.1767767
    ## 4 . . .         0.8660254 0.4330127 0.2165064
    ## 5 . . .         .         0.7071068 0.3535534
    ## 6 . . .         .         .         0.6846532

``` r
cholA = crossprod(L)
cholA
```

    ## 6 x 6 sparse Matrix of class "dsCMatrix"
    ##      1     2      3      4      5      6
    ## 1 1.00 .     0.5000 0.5000 0.5000 0.2500
    ## 2 .    1.000 0.5000 .      0.2500 0.6250
    ## 3 0.50 0.500 1.0000 0.2500 0.6250 0.5625
    ## 4 0.50 .     0.2500 1.0000 0.6250 0.3125
    ## 5 0.50 0.250 0.6250 0.6250 1.1250 0.6875
    ## 6 0.25 0.625 0.5625 0.3125 0.6875 1.1250

Note that $LL^T = L^TL$ only if $L$ is orthogonal. Remember that a matix
is defined orthogonal if $LL^T = L^TL = I$

``` r
L %*% t(L)
```

    ## 6 x 6 sparse Matrix of class "dgCMatrix"
    ##           1         2         3         4         5         6
    ## 1 1.8125000 0.5312500 0.5745243 0.7036456 0.4419417 0.1711633
    ## 2 0.5312500 1.7031250 0.5524272 0.2435696 0.3977476 0.4279082
    ## 3 0.5745243 0.5524272 0.6562500 0.1913664 0.3125000 0.1210307
    ## 4 0.7036456 0.2435696 0.1913664 0.9843750 0.3827328 0.1482318
    ## 5 0.4419417 0.3977476 0.3125000 0.3827328 0.6250000 0.2420615
    ## 6 0.1711633 0.4279082 0.1210307 0.1482318 0.2420615 0.4687500

``` r
t(L) %*% L
```

    ## 6 x 6 sparse Matrix of class "dgCMatrix"
    ##      1     2      3      4      5      6
    ## 1 1.00 .     0.5000 0.5000 0.5000 0.2500
    ## 2 .    1.000 0.5000 .      0.2500 0.6250
    ## 3 0.50 0.500 1.0000 0.2500 0.6250 0.5625
    ## 4 0.50 .     0.2500 1.0000 0.6250 0.3125
    ## 5 0.50 0.250 0.6250 0.6250 1.1250 0.6875
    ## 6 0.25 0.625 0.5625 0.3125 0.6875 1.1250

# 5) Inverse of the numerator relationship matrix

## Using LL’ method

Computational, inverse A directly, is inefficient when we need to
compute a huge number of animals, that is why both Cholesky
decomposition and generalized Cholesky decomposition help us to solve
it.

Thus, if

$$
A = LL^T 
$$

and

$$LL^T = (L^T)^{-1}L^{-1}$$

Then,

$$A^{-1} = (L^T)^{-1}L^{-1} =  (LL^T)^{-1}$$

where, - $L^{-1}$ is the inverse of the triangular matrix $L$ - $L^{-T}$
is the inverse of the transpose of $L$, which is the transpose of L^{-1}

``` r
solve(t(L) %*% L)
```

    ## 6 x 6 sparse Matrix of class "dgCMatrix"
    ##               1          2    3          4             5         6
    ## 1  1.833333e+00  0.5000000 -1.0 -0.6666667  1.110223e-16  .       
    ## 2  5.000000e-01  2.0333333 -1.0  .          5.333333e-01 -1.066667
    ## 3 -1.000000e+00 -1.0000000  2.5  0.5000000 -1.000000e+00  .       
    ## 4 -6.666667e-01  .          0.5  1.8333333 -1.000000e+00  .       
    ## 5  1.110223e-16  0.5333333 -1.0 -1.0000000  2.533333e+00 -1.066667
    ## 6  .            -1.0666667  .    .         -1.066667e+00  2.133333

## Using factorization method

Finally, $A^{-1}$ should be from generalized Cholesky decomposition.

$$A^{-1} = (T^{-1})^TD^{-1}T^{-1}$$

The diagonal $D^{-1}$ inverse is easy to compute because is the
reciprocal of diagonal elements of $D$

``` r
invD = solve(diag(D))
invD
```

    ##      [,1] [,2] [,3]     [,4] [,5]     [,6]
    ## [1,]    1    0    0 0.000000    0 0.000000
    ## [2,]    0    1    0 0.000000    0 0.000000
    ## [3,]    0    0    2 0.000000    0 0.000000
    ## [4,]    0    0    0 1.333333    0 0.000000
    ## [5,]    0    0    0 0.000000    2 0.000000
    ## [6,]    0    0    0 0.000000    0 2.133333

``` r
getDInv(ped2)
```

    ##        1        2        3        4        5        6 
    ## 1.000000 1.000000 2.000000 1.333333 2.000000 2.133333

Now, $T^-1$ is a lower triangular matrix with $1$ in the diagonal the
only non-zeros elements $-0.5$ correspond to know parents of the animal
$i$. It can be derived as $I - M$, where $I$ is an identuty matrix of
the order of animal on the pedigree and $M$ is a matrix of the
contributions of gametes from parents to progeny (Kennedy, 19)

Thus, the matrix $I$

``` r
n = length(ped2@label)
I = diag(n)
I
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6]
    ## [1,]    1    0    0    0    0    0
    ## [2,]    0    1    0    0    0    0
    ## [3,]    0    0    1    0    0    0
    ## [4,]    0    0    0    1    0    0
    ## [5,]    0    0    0    0    1    0
    ## [6,]    0    0    0    0    0    1

the matrix $M$

``` r
m = matrix(0, n, n)
m
```

    ##      [,1] [,2] [,3] [,4] [,5] [,6]
    ## [1,]    0    0    0    0    0    0
    ## [2,]    0    0    0    0    0    0
    ## [3,]    0    0    0    0    0    0
    ## [4,]    0    0    0    0    0    0
    ## [5,]    0    0    0    0    0    0
    ## [6,]    0    0    0    0    0    0

Now, for simplicity we can use the function

``` r
invT = getTInv(ped2)
```

Now,

``` r
gCholAinv = t(invT) %*% invD %*% invT
gCholAinv
```

    ## 6 x 6 Matrix of class "dgeMatrix"
    ##            1          2    3          4          5         6
    ## 1  1.8333333  0.5000000 -1.0 -0.6666667  0.0000000  0.000000
    ## 2  0.5000000  2.0333333 -1.0  0.0000000  0.5333333 -1.066667
    ## 3 -1.0000000 -1.0000000  2.5  0.5000000 -1.0000000  0.000000
    ## 4 -0.6666667  0.0000000  0.5  1.8333333 -1.0000000  0.000000
    ## 5  0.0000000  0.5333333 -1.0 -1.0000000  2.5333333 -1.066667
    ## 6  0.0000000 -1.0666667  0.0  0.0000000 -1.0666667  2.133333

The function $getAinv$ compute $(L^T)^{-1}L^{-1$, note that in this case
$L=sqrt(T)$

``` r
(getAInv(ped2))
```

    ## 6 x 6 sparse Matrix of class "dsCMatrix"
    ##            1          2    3          4          5         6
    ## 1  1.8333333  0.5000000 -1.0 -0.6666667  .          .       
    ## 2  0.5000000  2.0333333 -1.0  .          0.5333333 -1.066667
    ## 3 -1.0000000 -1.0000000  2.5  0.5000000 -1.0000000  .       
    ## 4 -0.6666667  .          0.5  1.8333333 -1.0000000  .       
    ## 5  .          0.5333333 -1.0 -1.0000000  2.5333333 -1.066667
    ## 6  .         -1.0666667  .    .         -1.0666667  2.133333

# Key word

Mendelian sampling

Marginal probability

Conditionally probability

Cholesky decomposition

Matrix properties
