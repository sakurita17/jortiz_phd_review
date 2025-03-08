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
defined relationship coefficient $R$ between individuals $i$ and $j$ as the
correlation between their genetic values (Gorjanc notes):

$$
\begin{aligned}
R_{i,j} = Cor(g_i,g_j) \\
Cov(g_i,g_j)/sqr(Var(g_i)Var(g_j))
\end{aligned}
$$

The Numerator relationship matrix $(A)$ describes the additive relationship among 
individuals. $(A)$ is a symmetric and its diagonal elements represent twice 
the probability that two gametes taken at random from animal $i$ will carry 
identical alleles by descent, in other words, the inbreeding coefficient of 
animal $i$ (Wright, 1922).

The off-diagonal elements, $a_{ij}$, equals the numerator of the coefficient 
of relationship between animals $i$ and $j$. When multiplied by the additive
genetic variance $\sigma^2_{a}$, $A\sigma^2_{a}$ is the covariance among breeding
values. $A$ is twice the coefficient of coancestry (kinship).

# Packages 

```r
library(package = "pedigreeTools") 
library(package = "Matrix")
```

# Objective

To understand the numerator relationship matrix and method for its computation

# Library

For the aim of this exercise, I will use the library “pedigreeTools”:
<https://github.com/Rpedigree/pedigreeTools>


# 1) Create a pedigree

``` r
ped = data.frame(iid = 1:6,
                 fid = c(0, 0, 1, 1, 4, 5), 
                 mid = c(0, 0, 2, 0, 3, 2))

(ped2 = pedigree(sire  = ped$fid, 
                 dam   = ped$mid,
                 label = ped$iid))
```

Output:

$$
\begin{array}{|c|c|c|}
\hline
idd & sire & dam \\
\hline
1 & NA & NA \\
\hline
2 & NA & NA \\
\hline
3 & 1 & 2 \\
\hline
4 & 1 & NA \\
\hline
5 & 4 & 3 \\
\hline
6 & 5 & 2 \\
\hline
\end{array}
$$

                
# 2) Tabular method

Applying the tabular method require to order the pedigree ensuring parents appear
before their progeny.Thus, the matrix {$A$} is formed following the recursive 
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

  $a_{ij} = a_{ji} = 0.5(a_{js/jd})$

If neither parent is known:

- The diagonal element:

  $a_{ii} = 1$

- The off-diagonal element:

  $a_{ij} = a_{ji} = 0$

Then, from the example pedigree from above, A is:

```r
(A = getA(ped2))
```

The prediction of breeding values requires the inverse of the relationship matrix, 
$A^-1$. This could be obtained by setting up $A$ by the recursive method and inverting
it. 

# 3) Generalised Cholesky decomposition 

The relationship matrix can be expressed as $A = TDT^T$ (Thompson, 1977), where $T$ 
is a lower triangular matrix and $D$ is a diagonal matrix.

## Matrix T
The matrix $T$ is a lower triangle matrix  traces the flow of genes from one generation to the other. Thus, 
for our small pedigree, we know that the expected and variance of the breeding 
values are as follow:

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

Now, the system of equations is:

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

We can get $T$ following the next rules, 

If both parents $s$ and $d$ are known 
 $t_{ij} = 0.5(t_{sj} + t_{dj})$

If only want parent $s$ or $d$ is known
 $t_{ij} = 0.5(t_{sj})$

If neither parent is known
 $t_{ij} = 0$

Diagonal is 
 $t_{ii} = 1$
 
For instance, for the pedigree above:

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
We use the function getT from pedigreeTools to get the matrix T

``` r
(T = getT(ped2))
```
Outcome:

$$
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 & 0 & 0 \\
0.5 & 0.5 & 1 & 0 & 0 & 0 \\
0.5 & 0 & 0 & 1 & 0 & 0 \\
0.5 & 0.25 & 0.5 & 0.5 & 1 & 0 \\
0.25 & 0.625 & 0.25 & 0.25 & 0.5 & 1 
\end{bmatrix}
$$

## Matrix D
The diagonal matrix $D$ is the variance and covariance matrix for Mendelian sampling. 
The Mendelian sampling $(m)$ for an animal $i$ with breeding value $a_i$ and $a_s$
and $a_d$ as breeding values for its sire and dam, respectively, is:

$$
\begin{aligned}
m_i = a_i - 0.5(a_s + a_d)
\end{aligned}
$$

Diagonal elements of $D$ are 1 for pedigree founders, $0.5 - 0.25(F_s + Fd)$ if both
parents of animal $i$ are known and $0.75 - 0.25(F_s)$ if only one parent $s$ is
known.

Thus, the diagonal matrix for our example is:

$$
\begin{equation}
\begin{aligned}
d_{1,1} &= 1 \\

d_{2,2} &= 1 \\

d_{3,3} &= 0.5 - 0.25(F_1 + F_2) 
         = 0.5 - 0.25(0)
         = 0 \\

d_{4,4} &= 0.75 - 0.25(F_1) 
         = 0.75 - 0.25(0)
         = 0.75 \\

d_{5,5} &= 0.5 - 0.25(F_4 + F_3) 
         = 0.5 - 0.25(0 + 0)
         = 0.5 \\

d_{6,6} &= 0.5 - 0.25(F_5 + F_2) 
         = 0.5 - 0.25(0.125 + 0)
         = 0.46875

\end{aligned}
\end{equation}
$$

``` r
(D = getD(ped2))
```
$$
\begin{bmatrix}
1.00000 & 1.00000 & 0.50000 & 0.75000 & 0.50000 & 0.46875 
\end{bmatrix}
$$

Finally, we can compute $A$

```r
(A_GCD =  T %*% diag(D) %*% Matrix::t(T))
```
    
# Cholesky decomposition 

Using the concepts of Cholesky decomposition, $A$ can be written as $A$ = LL', 
where $L$ is a lower triangle matrix with non-zero diagonal (Henderson, 1975). 
Quaas (1976) presented a strategy for obtaining the diagonal elements of A while
computing A−1 without setting up the relationship matrix.

```r
(L = t(relfactor(ped2)))
```


```r
(A_chol = L %*% t(L)) 
```

# Inverse of the relationship matrix 

If $A = LL^T$, where $L$ is a lower triangular and $L^T$ is its transpose, then 
then the inverse of $A$ can be computed as follows:

$A^-1 = (LL^T)^-1$

Using the property of matrix inversion:

$(LL^T) = (L^T)^-1L^-1$

Thus, 

$A^-1 = (L^T)^-1L^-1$

where, 
- L^{-1} is the inverse of the triangular matrix $L$
- L^{-T} is the inverse of the transpose of $L$, which is the transpose of L^{-1}










# Definitions

Mendelian sampling

Marginal probability

Conditionally probability

Cholesky decomposition

Matrix properties

