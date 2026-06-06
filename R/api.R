# User-facing wrappers. Each forwards to the underlying method function living in
# that method's private environment (see zzz.R). Arguments and return values are
# exactly those of the wrapped function; see the source under inst/method/ for the
# full signatures.

# Fetch a function `fn` from method `method`'s private environment.
.cap_get <- function(method, fn) {
  e <- get0(method, envir = .cap_envs, inherits = FALSE)
  if (is.null(e))
    stop("CAPmethods: unknown method '", method, "'.", call. = FALSE)
  get(fn, envir = e)
}

#' Reach an internal function of a CAP method by name.
#'
#' Escape hatch for advanced use: every internal function of a method (not just
#' the exported wrappers) lives in a private environment. This returns it.
#'
#' @param method one of "hdcap", "lcap", "mcap", "coc", "mediation", "hdcov",
#'   "cluster".
#' @param fn name of the function to retrieve.
#' @return The requested function object.
#' @export
cap_internal <- function(method, fn) .cap_get(method, fn)

# ---- HDCAP (high-dimensional CAP; cov.shrinkage=FALSE recovers classical CAP)
#' @export
hdcap <- function(...) .cap_get("hdcap", "capReg")(...)
#' @export
hdcap_boot <- function(...) .cap_get("hdcap", "cap_beta_boot")(...)

# ---- LCAP (longitudinal CAP, time-invariant projection)
#' @export
lcap <- function(...) .cap_get("lcap", "capReg")(...)
#' @export
lcap_boot <- function(...) .cap_get("lcap", "cap_beta_boot")(...)

# ---- MCAP (multilevel CAP, cluster-varying loadings)
#' @export
mcap <- function(...) .cap_get("mcap", "lcapReg")(...)
#' @export
mcap_boot <- function(...) .cap_get("mcap", "lcap.beta.boot")(...)

# ---- CAP-CoC (covariance-on-covariance regression)
#' @export
coc <- function(...) .cap_get("coc", "COCReg")(...)
#' @export
coc_d1 <- function(...) .cap_get("coc", "COCReg.D1")(...)   # single direction
#' @export
coc_boot <- function(...) .cap_get("coc", "COCReg.coef.boot")(...)

# ---- CAP mediation
#' @export
capmediation <- function(...) .cap_get("mediation", "CAPMediation")(...)
#' @export
capmediation_boot <- function(...) .cap_get("mediation", "CAPMediation_boot")(...)

# ---- HDcov / HCAP (high-dimensional covariance CAP pipeline)
#' @export
hcap <- function(...) .cap_get("hdcov", "run_HCAP_pipeline")(...)

# ---- CAP clustering
#' @export
cappcl <- function(...) .cap_get("cluster", "capPCL")(...)
#' @export
cappcl_boot <- function(...) .cap_get("cluster", "capPCL_coef_boot")(...)
