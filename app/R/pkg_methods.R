# =============================================================================
# Bridge to the installed CAPmethods package.
#
# Instead of sourcing each method's sibling `Vn/` folder and dyn.load()-ing its
# kernel, every plugin pulls its method code from the CAPmethods package
# (>= 0.2.0 bundles HDCAP V9, LCAP V7, MCAP V5, CAP-CoC V4, CAP-mediation V5,
# HCAP V3, CAP-clustering V5). The package keeps each method's functions in its
# own private environment (CAPmethods:::.cap_envs[[method]]) with the compiled
# kernels registered in the package DLL; every .Call names that DLL
# (PACKAGE = .<prefix>_vN_dll = "CAPmethods", a variable cloned with the
# functions below), so the kernels resolve regardless of closure.
#
# cap_pkg_env(method) returns a fresh environment holding a clone of that
# method's functions, with:
#   * each function's environment rebound to the clone, so internal helper calls
#     — and any plugin override (e.g. hdcov's registerDoParallel) — resolve here;
#   * parent = the CAPmethods namespace, so imports (ginv, lme, brmultinom, ...)
#     and the registered kernels resolve.
# This reproduces the old `sys.source(..., envir = .X_env)` semantics exactly,
# without needing the source folders or a local C++ build.
#
# method ids (CAPmethods side): hdcap, lcap, mcap, coc, mediation, hdcov, cluster
#   CAP        -> "hdcap" (classical CAP = HDCAP with cov.shrinkage = FALSE)
#   HDCAP      -> "hdcap"
#   LCAP       -> "lcap"
#   MCAP       -> "mcap"
#   CAP-CoC    -> "coc"
#   CAP-medi.  -> "mediation"
#   CAP-HDcov  -> "hdcov"
#   CAP-clust. -> "cluster"
# =============================================================================

CAP_PKG_MIN_VERSION <- "0.2.0"   # the method versions the plugins are written for

cap_pkg_env <- function(method) {
  if (!requireNamespace("CAPmethods", quietly = TRUE))
    stop("The 'CAPmethods' package is required but not installed. Install it from ",
         "CAP-Rpkg/ (e.g. remotes::install_github(\"zhaoyi1026/CAPmethods\") or ",
         "install.packages(<CAP-Rpkg path>, repos = NULL, type = \"source\")).",
         call. = FALSE)
  if (utils::packageVersion("CAPmethods") < CAP_PKG_MIN_VERSION)
    warning("CAPmethods ", as.character(utils::packageVersion("CAPmethods")),
            " is installed; the app expects >= ", CAP_PKG_MIN_VERSION,
            " (run_local.R reinstalls it).", call. = FALSE)
  ns   <- asNamespace("CAPmethods")
  envs <- get0(".cap_envs", envir = ns, inherits = FALSE)
  if (is.null(envs) || !exists(method, envir = envs, inherits = FALSE))
    stop("CAPmethods has no method environment '", method, "'.", call. = FALSE)
  src <- get(method, envir = envs, inherits = FALSE)
  out <- new.env(parent = ns)
  for (nm in ls(src, all.names = TRUE)) {
    obj <- get(nm, envir = src, inherits = FALSE)
    if (is.function(obj)) environment(obj) <- out
    assign(nm, obj, envir = out)
  }
  out
}
