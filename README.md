# CAPmethods

A unified R package bundling the **Covariate-Assisted Principal (CAP) regression**
family of methods for relating subject-level covariates to the
covariance / heterogeneity structure of multivariate and longitudinal outcomes.
The computational kernels are implemented in C++ via
[RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo).

## Methods

| Function | Method | Notes |
|----------|--------|-------|
| `hdcap()` | High-dimensional CAP | optional covariance shrinkage; `cov.shrinkage = FALSE` reproduces classical CAP |
| `lcap()` | Longitudinal CAP | time-invariant projection; random intercept + slopes |
| `mcap()` | Multilevel CAP | cluster-varying loadings (von Mises–Fisher) |
| `coc()` | Covariance-on-covariance regression | |
| `capmediation()` | CAP mediation | |
| `hcap()` | High-dimensional-covariance CAP (HCAP) | `glmnet`-based pipeline |
| `cappcl()` | CAP clustering | covariance-matrix clustering |

Each method that supports bootstrap inference also exposes a `*_boot()` companion
(`hdcap_boot()`, `lcap_boot()`, `mcap_boot()`, `coc_boot()`,
`capmediation_boot()`, `cappcl_boot()`). The escape hatch
`cap_internal(method, fn)` returns any internal function of a method by name.

## Installation

```r
# install.packages("remotes")
remotes::install_github("zhaoyi1026/CAPmethods")
```

A C++ toolchain (and BLAS/LAPACK) is needed to compile the kernels at install
time. On macOS this means the Xcode command-line tools; on Windows, Rtools.

Imported packages: MASS, glmnet, nlme, brglm2, foreach, doParallel, Rcpp
(plus base `parallel`, `stats`, `utils`, `graphics`, `grDevices`).

## Quick start

```r
library(CAPmethods)

## HDCAP on a small synthetic example (cov.shrinkage = FALSE => classical CAP)
set.seed(1)
n <- 40; p <- 5
X <- cbind(1, rnorm(n))                       # n x q covariate matrix (+intercept)
b <- c(0, 0.9)                                 # covariate drives variance on dir 1
Y_list <- lapply(1:n, function(i) {
  d   <- exp(as.numeric(X[i, ] %*% b))
  Sig <- diag(p); Sig[1, 1] <- d
  MASS::mvrnorm(50, rep(0, p), Sig)            # T_i x p response
})

fit <- hdcap(Y_list, X, stop.crt = "nD", nD = 1, cov.shrinkage = FALSE)
fit$beta    # ~ (0, 0.9)
fit$gamma   # estimated loading direction
```

## Design

Every method's original R implementation is sourced into its **own private
environment** at load time, so the identically named internal helpers
(`capReg`, `cov.ls`, `gamma.solve`, `obj.func`, …) do not collide in a single
package namespace. All C++ kernels compile into one shared object; each carries a
per-method symbol prefix (`hd_`, `lcap_`, `mlcap_`, `coc_`, `med_`, `hdcov_`,
`cluster_`) and is registered in `src/init.c`.

## License

GPL-3.
