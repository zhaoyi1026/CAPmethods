# User-facing wrappers. Each forwards to the underlying method function living in
# that method's private environment (see zzz.R). Arguments and return values are
# exactly those of the wrapped function; see the source under inst/method/ and the
# package vignette for the full signatures.

# Fetch a function `fn` from method `method`'s private environment.
.cap_get <- function(method, fn) {
  e <- get0(method, envir = .cap_envs, inherits = FALSE)
  if (is.null(e))
    stop("CAPmethods: unknown method '", method, "'.", call. = FALSE)
  get(fn, envir = e)
}

#' Fit a CAP-family method
#'
#' Thin wrappers that run each covariate-assisted principal (CAP) regression
#' method. Each forwards its arguments to the underlying estimator (sourced into
#' a private environment at package load, so the methods' identically named
#' internal helpers do not collide) and returns that estimator's result
#' unchanged.
#'
#' - `hdcap()` — high-dimensional CAP; `cov.shrinkage = FALSE` recovers classical CAP.
#' - `lcap()` — longitudinal CAP with a time-invariant projection.
#' - `mcap()` — multilevel CAP with cluster-varying loadings.
#' - `coc()` — covariance-on-covariance regression (multi-direction selection);
#'   `coc_d1()` estimates a single direction.
#' - `capmediation()` — covariance-mediator analysis.
#' - `hcap()` — high-dimensional-covariate CAP (sparse, `glmnet`-based).
#' - `cappcl()` — clustering of covariance patterns.
#'
#' Each method has a matching [cap_examples] generator producing ready-to-use
#' data; `vignette("CAPmethods")` shows a worked, plotted run of every method.
#'
#' @param ... arguments passed to the underlying method estimator, e.g.
#'   `stop.crt` (`"nD"` or `"DfD"`), `nD`, `DfD.thred`, `cov.shrinkage`. See
#'   `vignette("CAPmethods")` for each method's full argument list.
#' @return The fitted-model object from the underlying estimator: a list whose
#'   elements depend on the method — typically the loading direction(s) `gamma`,
#'   coefficients (`beta` and/or `alpha`), per-subject variance `score`s, and
#'   direction-selection diagnostics (`DfD`).
#' @seealso [cap_examples] for matching example data, [cap_boot] for bootstrap
#'   inference, and [cap_internal] to reach any internal function by name.
#' @examples
#' d   <- hdcap_example(n = 30, Ti = 30)
#' fit <- hdcap(d$Y_list, d$X, nD = 1, cov.shrinkage = FALSE)
#' round(as.numeric(fit$beta), 2)
#' @name cap_fit
NULL

#' @rdname cap_fit
#' @export
hdcap <- function(...) .cap_get("hdcap", "capReg")(...)

#' @rdname cap_fit
#' @export
lcap <- function(...) .cap_get("lcap", "capReg")(...)

#' @rdname cap_fit
#' @export
mcap <- function(...) .cap_get("mcap", "lcapReg")(...)

#' @rdname cap_fit
#' @export
coc <- function(...) .cap_get("coc", "COCReg")(...)

#' @rdname cap_fit
#' @export
coc_d1 <- function(...) .cap_get("coc", "COCReg.D1")(...)

#' @rdname cap_fit
#' @export
capmediation <- function(...) .cap_get("mediation", "CAPMediation")(...)

#' @rdname cap_fit
#' @export
hcap <- function(...) .cap_get("hdcov", "run_HCAP_pipeline")(...)

#' @rdname cap_fit
#' @export
cappcl <- function(...) .cap_get("cluster", "capPCL")(...)

#' Bootstrap inference for CAP-family coefficients
#'
#' Companion bootstrap estimators that add standard errors and/or confidence
#' intervals for a method's coefficients. Pass `sims` to set the number of
#' bootstrap replicates.
#'
#' - `hdcap_boot()`, `lcap_boot()`, `mcap_boot()`, `coc_boot()`,
#'   `capmediation_boot()`, `cappcl_boot()` correspond to the like-named fits in
#'   [cap_fit].
#'
#' @param ... arguments passed to the underlying bootstrap estimator, notably
#'   `sims` (number of replicates), the fitted loading direction, and
#'   `conf.level`.
#' @return A list with bootstrap standard errors and/or confidence intervals for
#'   the method's coefficients.
#' @seealso [cap_fit]
#' @name cap_boot
NULL

#' @rdname cap_boot
#' @export
hdcap_boot <- function(...) .cap_get("hdcap", "cap_beta_boot")(...)

#' @rdname cap_boot
#' @export
lcap_boot <- function(...) .cap_get("lcap", "cap_beta_boot")(...)

#' @rdname cap_boot
#' @export
mcap_boot <- function(...) .cap_get("mcap", "lcap.beta.boot")(...)

#' @rdname cap_boot
#' @export
coc_boot <- function(...) .cap_get("coc", "COCReg.coef.boot")(...)

#' @rdname cap_boot
#' @export
capmediation_boot <- function(...) .cap_get("mediation", "CAPMediation_boot")(...)

#' @rdname cap_boot
#' @export
cappcl_boot <- function(...) .cap_get("cluster", "capPCL_coef_boot")(...)

#' Reach an internal function of a CAP method by name
#'
#' Escape hatch for advanced use: every internal function of a method (not only
#' the exported wrappers) lives in that method's private environment. This
#' returns the requested function object.
#'
#' @param method one of `"hdcap"`, `"lcap"`, `"mcap"`, `"coc"`, `"mediation"`,
#'   `"hdcov"`, `"cluster"`.
#' @param fn name of the function to retrieve.
#' @return The requested function object.
#' @seealso [cap_fit]
#' @examples
#' f <- cap_internal("coc", "COCReg.D1")
#' is.function(f)
#' @export
cap_internal <- function(method, fn) .cap_get(method, fn)
