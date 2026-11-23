# CAPmethods 0.2.0 (2026-11-23)

The method code is updated to the latest version of each method:
HDCAP V9, LCAP V7, MCAP V5, CAP-CoC V4, CAP-mediation V5, HCAP V3 and
CAP clustering V5.

## Changes in results

* `hdcap()`, `hdcap_boot()` (HDCAP V9):
  - The objective that compares the initial values of the direction search is
    the function the fit minimizes: without covariance shrinkage the negative
    log-likelihood; with shrinkage the shrunk scores of the fit, with the fit's
    common shrinkage weights. (Before, it applied shrinkage weights of another
    formula to the sample covariances, also when the fit used no shrinkage.)
  - New argument `beta0`: the initial value(s) of β — a vector of length
    `ncol(X)`, or a matrix with one initial value per column. The default
    `NULL` uses two, `0` and the previous default start
    (`set.seed(100); c(0, rnorm(q - 1, 1, 1))`, now drawn without changing the
    caller's random numbers). Each initial value of the direction is run from
    each initial value of β; every β fit (including each bootstrap draw) keeps
    the converged estimate with the lowest objective.
* `coc()`, `coc_d1()`, `coc_boot()` (CAP-CoC V4):
  - With outcome covariance shrinkage (the default), the objective that compares
    initial values uses the shrunk covariances the fit iterates on
    (`cov.sk.y` at the fitted variances) instead of the per-subject linear
    shrinkage.
  - New argument `beta0` (as for `hdcap()`; α starts at 0). Default: two initial
    values, `0` and the previous start.

## Unchanged results

* `lcap()`, `lcap_boot()`, `mcap()`, `mcap_boot()`, `capmediation()`,
  `capmediation_boot()`, `cappcl()`, `cappcl_boot()` and `hcap()` give the
  results of 0.1.0. Their objectives already were the functions their fits
  minimize, and their coefficients start from data-driven values.
* `hcap()` has a new argument `obj.scale` (default `TRUE`, the previous
  comparison of initial values at unit-length γ). `obj.scale = FALSE` compares
  them by the function the fit minimizes, at γ'Hγ = 1.

## Robustness

* A failed iteration in a C++ kernel (e.g. a diverging β iteration) is now an R
  error, caught by the methods, instead of an error that could end the R
  session; such an initial value or bootstrap draw is skipped.
* The searches over initial values no longer reuse the fit of the previous
  initial value when one fails, and stop with a clear error when none succeeds.
* `hcap()` with `stop.crt = "DfD"` stops when the next direction cannot be
  estimated (it could retry the same direction indefinitely), and works with
  one initial value (`ninitial = 1`); `coc()` also works with `ninitial = 1`.
* Every kernel call names the package DLL, so the package can be used in a
  session that has also loaded a method's stand-alone kernels.
* The random EM start of `cappcl()` and the default β starts of `hdcap()` and
  `coc()` are drawn as before (`set.seed(100)`) but no longer reset the caller's
  random-number stream.

## Python package

* `capcov` 0.2.0 (in `Python/`) follows these changes: for HDCAP/CAP and CAP-CoC
  the objective that compares initial values and several initial values of β
  (`beta0`; the default ones are R's numbers), and for every method a clear
  error when the first direction cannot be estimated, the end of the search when
  a later one cannot.

# CAPmethods 0.1.0 (2026-06-05)

* First version: HDCAP, LCAP, MCAP, CAP-CoC, CAP-mediation, HCAP and CAP
  clustering, with C++ kernels, example-data generators, a vignette and the
  companion Shiny app and Python package.
