// =============================================================================
// CAP-HDcov/V2 -- RcppArmadillo kernels for the CAP matrix-algebra core
// -----------------------------------------------------------------------------
// Compiled to a shared object (build_kernels.R) and loaded with dyn.load(); the
// R wrappers in HCAP_Code.R call these via .Call(). The math mirrors V1
// (precompute_sigma, projected scores v'Sigma_i v, and the weighted covariance
// accumulation sum_i w_i Sigma_i) exactly.
//
// NOTE: CAP-HDcov is dominated by the high-dimensional lasso/GLM fitting
// (cv.glmnet / glm), which stays in R (already compiled C/Fortran). These
// kernels accelerate only the p x p covariance algebra (negligible when p is
// small, helpful when the response dimension p is large). See CAP-HDcov/README.md.
//
// Symbols are hdcov_-prefixed to avoid colliding with the other CAP kernels.
//
// CAP-HDcov/V3: the kernels of V2, unchanged, except that each .Call entry point
// is wrapped in BEGIN_RCPP / END_RCPP (as HDCAP V8).  In V2 a C++ exception (e.g.
// a conversion error, or a dimension mismatch) crossed the extern "C" boundary
// uncaught and aborted the R process; now it is an ordinary R error, which the
// try()/tryCatch() calls of HCAP_Code.R catch.  HCAP_Code.R names the V3 kernel
// in every .Call (PACKAGE).
// =============================================================================
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// per-subject sample covariance cube + sizes (== precompute_sigma):
// Sigma_i = crossprod(scale(Y_i, center)) / T_i   (denominator T_i)
static List precompute_core(const List& Y) {
  const int n = Y.size();
  mat y0 = as<mat>(Y[0]);
  const uword p = y0.n_cols;
  cube Sigma(p, p, n);
  vec Tvec(n);
  for (int i = 0; i < n; i++) {
    mat Yi = as<mat>(Y[i]);
    double Ti = Yi.n_rows;
    Tvec[i] = Ti;
    mat Yc = Yi.each_row() - mean(Yi, 0);
    Sigma.slice(i) = (Yc.t() * Yc) / Ti;
  }
  return List::create(_["Sigma"] = Sigma, _["Tvec"] = Tvec);
}

// projected scores  v' Sigma_i v  for all i  (== vapply(..., t(v) %*% S %*% v))
static vec score_core(const cube& S, const vec& v) {
  vec o(S.n_slices);
  for (uword i = 0; i < S.n_slices; i++) o[i] = as_scalar(v.t() * S.slice(i) * v);
  return o;
}

// weighted accumulation  sum_i w_i Sigma_i  (== Reduce(`+`, Sigma[,,i] * w_i))
static mat accum_core(const cube& S, const vec& w) {
  mat A(S.n_rows, S.n_cols, fill::zeros);
  for (uword i = 0; i < S.n_slices; i++) A += w[i] * S.slice(i);
  return A;
}

// ---- extern "C" wrappers ----------------------------------------------------
extern "C" SEXP hdcov_precompute_cpp(SEXP Y_) {
BEGIN_RCPP
  return precompute_core(as<List>(Y_));
END_RCPP
}
extern "C" SEXP hdcov_score_cpp(SEXP S_, SEXP v_) {
BEGIN_RCPP
  return wrap(score_core(as<cube>(S_), as<vec>(v_)));
END_RCPP
}
extern "C" SEXP hdcov_accum_cpp(SEXP S_, SEXP w_) {
BEGIN_RCPP
  return wrap(accum_core(as<cube>(S_), as<vec>(w_)));
END_RCPP
}
