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
| `hcap()` | CAP with high-dimensional covariates | regularization and post-selection inference by multiple splitting |
| `cappcl()` | CAP clustering | covariance-matrix clustering |

Each method that supports bootstrap inference also exposes a `*_boot()` companion
(`hdcap_boot()`, `lcap_boot()`, `mcap_boot()`, `coc_boot()`,
`capmediation_boot()`, `cappcl_boot()`). The escape hatch
`cap_internal(method, fn)` returns any internal function of a method by name.

Every method ships a built-in synthetic-data generator (`hdcap_example()`,
`lcap_example()`, …) returning data in the exact shape the wrapper expects plus
the ground truth, so a full run is two lines.

## Examples & visualization

See **[Examples.md](Examples.md)** for a worked, plotted walkthrough of every
method — generate example data, fit, and a figure checking the estimates against
the truth. A taste (HDCAP: estimated loadings vs truth, and subject variance
scores rising with the covariate):

![HDCAP example](man/figures/hdcap.png)

## Shiny app

A companion **CAP Methods Explorer** Shiny app (in [`app/`](app/)) provides a
point-and-click interface to every method — run on built-in examples or your own
uploaded data, tune parameters, and download result tables and plots. It uses
this package as its engine. To run it locally, clone the repo and:

```r
install.packages(c("shiny", "bslib", "bsicons", "DT", "plotly",
                   "shinycssloaders", "markdown", "mvtnorm"))
shiny::runApp("app", launch.browser = TRUE)   # from the repository root
```

See **[app/README.md](app/README.md)** for full install and launch instructions.

## Python package

A companion Python package, **`capcov`** (in [`Python/`](Python/)), ports the CAP
family to Python (numpy/scipy). Implemented and verified against this R package:
HDCAP/CAP, CAP-CoC, LCAP, CAP-mediation, CAP-clustering (MCAP and CAP-HDcov are
planned). Install from the cloned repository:

```bash
git clone https://github.com/zhaoyi1026/CAPmethods.git
pip install ./CAPmethods/Python            # core (numpy, scipy)
pip install "./CAPmethods/Python[full]"    # + scikit-learn / statsmodels / pandas
```

Requires Python ≥ 3.7. Quick check:

```python
import capcov
from capcov import examples
d   = examples.hdcap_example()
fit = capcov.cap_reg(d["Y"], d["X"], stop_crt="nD", nD=2, cov_shrinkage=False)
```

See **[Python/USAGE.md](Python/USAGE.md)** for how to use every method, and
**[Python/README.md](Python/README.md)** for status and verification notes.

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

## References

- **`hdcap()` — high-dimensional CAP** (and the classical CAP it subsumes with
  `cov.shrinkage = FALSE`):
  - Zhao, Y., Caffo, B., Luo, X., & Alzheimer's Disease Neuroimaging Initiative
    (2021). Principal regression for high dimensional covariance matrices.
    *Electronic Journal of Statistics*, 15(2), 4192.
    <https://doi.org/10.1214/21-EJS1887>
  - Zhao, Y., Wang, B., Mostofsky, S. H., Caffo, B. S., & Luo, X. (2021). Covariate
    assisted principal regression for covariance matrix outcomes. *Biostatistics*,
    22(3), 629–645. <https://doi.org/10.1093/biostatistics/kxz057>
- **`lcap()` — longitudinal CAP**: Zhao, Y., Caffo, B. S., & Luo, X. (2024).
  Longitudinal regression of covariance matrix outcomes. *Biostatistics*, 25(2),
  385–401. <https://doi.org/10.1093/biostatistics/kxac045>
- **`coc()` — covariance-on-covariance regression**: Zhao, Y., & Zhao, Y. (2025).
  Covariance-on-covariance regression. *Biometrics*, 81(3), ujaf097.
  <https://doi.org/10.1093/biomtc/ujaf097>
- **`capmediation()` — mediation with a graph mediator**: Xu, Y., & Zhao, Y. (2025).
  Mediation analysis with graph mediator. *Biostatistics*, 26(1), kxaf004.
  <https://doi.org/10.1093/biostatistics/kxaf004>

`mcap()` (multilevel CAP), `hcap()` (CAP with high-dimensional covariates), and
`cappcl()` (CAP clustering) do not have an associated publication listed yet.

## License

GPL-3.
