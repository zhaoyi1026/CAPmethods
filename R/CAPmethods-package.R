#' @keywords internal
"_PACKAGE"

## Namespace directives (collected here; the method code in inst/method/ resolves
## these imported symbols through the package namespace, and the compiled kernels
## are registered via useDynLib).
#' @import stats
#' @importFrom Rcpp evalCpp
#' @importFrom utils tail
#' @importFrom MASS ginv
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom nlme lme lmeControl
#' @importFrom brglm2 brmultinom
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @useDynLib CAPmethods, .registration = TRUE
NULL
