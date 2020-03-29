
<!-- README.md is generated from README.Rmd. Please edit that file -->

# snpcombineR

The goal of `snpcombineR` is to provide a simple and fast alternative
for combining SNP-level summary statistics at the gene-level (as well as
non-coding region-level). All underlying functions for combining
SNP-level statistics are written in
[`Rcpp`](http://dirk.eddelbuettel.com/code/rcpp.html), primarily using
[`RcppArmadillo`](http://dirk.eddelbuettel.com/code/rcpp.armadillo.html).
The package is currently under active development, with a primary usage
for simulating GWAS summary statistics to evaluate different approaches
of combining SNP-level p-values. Documentation will continue to be
updated as the package is developed. The remainder of the `README` walks
through a basic example for generating simulated GWAS summary
statistics.

## Installation

You can install `snpcombineR` from github with:

``` r
# install.packages("devtools")
devtools::install_github("ryurko/snpcombineR")
```

## Simulate GWAS Example

This is a basic example which shows you how to simulate GWAS summary
statistics using `snpcombineR`. In the current version of the package,
the user is expected to provide a genotype matrix. For simplicity, we
generate a genotype matrix below using the same code from a
[`bigsnpr`](https://privefl.github.io/bigsnpr/articles/pruning-vs-clumping.html)
vignette by Florian Privé, which relies on the `bindata` package. This
is only to generate a fake genotype matrix in this example and not
required for using `snpcombineR` (but you should check out
[`bigsnpr`](https://privefl.github.io/bigsnpr/index.html)).

``` r
gen <- function(n, m) {
  I <- 1:m
  p <- I / (2 * m + 1)

  mat <- outer(I, I, FUN = function(i, j) {
    1 / (abs(i - j) + 1)
  })

  bindata::rmvbin(n, p, bincorr = mat) + 
    bindata::rmvbin(n, p, bincorr = mat)
}

set.seed(1389)
X <- gen(5000, 20)
```

We can use `snpcombineR` to compute the correlation matrix of this
genotype matrix, with a comparison below of runtime with the `cor`
function (timed using `microbenchmark`):

``` r
microbenchmark::microbenchmark(
  "snpcombineR" = snpcombineR::compute_cor_matrix(X),
  "cor" = cor(X), times = 100
)
#> Unit: milliseconds
#>         expr      min       lq     mean   median       uq      max neval
#>  snpcombineR 1.145421 1.236781 1.426040 1.346145 1.414958 7.431948   100
#>          cor 1.244035 1.346676 1.480173 1.460104 1.545193 2.501123   100
#>  cld
#>    a
#>    a
```

Using this genotype matrix and a vector of case-rate probabilities for
each individual we can simulate GWAS summary statistics using the
`simulate_gene_gwas_data` function. In this example we assume the global
null is true and that every individual has a probability of 0.5, but
`create_gwas_case_prob` can be used to generate probabilities for
non-null settings with “causal” SNPs. The resulting dataset contains a
row for each SNP with columns (in order) for: beta, standard error,
z-statistic, and p-value.

``` r
case_rate_probs <- rep(0.5, nrow(X))
snpcombineR::simulate_gene_gwas_data(X, case_rate_probs)
#>               [,1]       [,2]        [,3]        [,4]
#>  [1,] -0.430039059 0.13154955 -3.26902726 0.001079179
#>  [2,] -0.112522154 0.09246153 -1.21696183 0.223618738
#>  [3,] -0.111593780 0.07704715 -1.44838299 0.147509969
#>  [4,] -0.058235862 0.06678970 -0.87192883 0.383247207
#>  [5,]  0.029052101 0.06151489  0.47227751 0.636728720
#>  [6,] -0.056882835 0.05701076 -0.99775606 0.318397662
#>  [7,] -0.037538756 0.05404239 -0.69461681 0.487295474
#>  [8,] -0.073754118 0.05041122 -1.46304966 0.143453791
#>  [9,] -0.037235576 0.04840096 -0.76931479 0.441706458
#> [10,] -0.024325751 0.04718675 -0.51552079 0.606189142
#> [11,]  0.009486209 0.04533704  0.20923751 0.834262833
#> [12,] -0.017862366 0.04412671 -0.40479711 0.685626656
#> [13,]  0.024019198 0.04326932  0.55510918 0.578819986
#> [14,]  0.000682254 0.04187923  0.01629099 0.987002248
#> [15,]  0.001691482 0.04152409  0.04073496 0.967507193
#> [16,] -0.040222506 0.04098881 -0.98130449 0.326442611
#> [17,] -0.016399413 0.04032441 -0.40668700 0.684237887
#> [18,] -0.006129125 0.04027484 -0.15218245 0.879043034
#> [19,]  0.035309925 0.04039983  0.87401164 0.382111917
#> [20,] -0.053188821 0.04034612 -1.31831325 0.187398808
```

## Contact

Ron Yurko: <ryurko@andrew.cmu.edu>
