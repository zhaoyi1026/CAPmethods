// =============================================================================
// CAP-mediation/V4 -- RcppArmadillo kernels
// -----------------------------------------------------------------------------
// Compiled to a shared object (build_kernels.R) and loaded with dyn.load(); the
// R wrappers in CAPMediation.R call these via .Call(). The math mirrors
// CAP-mediation/V3 (projected scores, the theta-update accumulation,
// eigen.solve, and per-subject sample covariance) exactly; only the hot loops
// move to C++. The mixed-effects (lme) fit stays in R.
//
// Symbols are med_-prefixed to avoid colliding with the CAP/HDCAP/CoC kernels
// that may be loaded in the same R process.
// =============================================================================
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// projected scores  theta' S_i theta  for all i  (== apply(S, 3, ...))
static vec score_core(const cube& S, const vec& v) {
  vec o(S.n_slices);
  for (uword i = 0; i < S.n_slices; i++) o[i] = as_scalar(v.t() * S.slice(i) * v);
  return o;
}

// weighted accumulation  sum_i w_i S_i  (== the theta-update A matrix loop)
static mat accum_core(const cube& S, const vec& w) {
  mat A(S.n_rows, S.n_cols, fill::zeros);
  for (uword i = 0; i < S.n_slices; i++) A += w[i] * S.slice(i);
  return A;
}

// generalized smallest eigenvector (== eigen.solve(A, H))
static vec eigen_solve_core(const mat& A, const mat& H) {
  vec d; mat U;
  eig_sym(d, U, symmatu(H));                 // ascending
  mat Dis = diagmat(1.0 / sqrt(d));
  mat M = Dis * U.t() * A * U * Dis;
  M = symmatu(0.5 * (M + M.t()));
  vec mv; mat V;
  eig_sym(mv, V, M);                          // ascending; smallest = column 0
  return U * Dis * V.col(0);
}

// per-subject sample covariance cube (== cov(M[[i]]); denominator T_i - 1)
static cube cov_core(const List& M) {
  const int n = M.size();
  mat m0 = as<mat>(M[0]);
  const uword p = m0.n_cols;
  cube out(p, p, n);
  for (int i = 0; i < n; i++) {
    mat Mi = as<mat>(M[i]);
    double Ti = Mi.n_rows;
    mat mc = Mi.each_row() - mean(Mi, 0);
    out.slice(i) = (mc.t() * mc) / (Ti - 1.0);
  }
  return out;
}

// ---- extern "C" wrappers ----------------------------------------------------
extern "C" SEXP med_score_cpp(SEXP S_, SEXP v_) {
  return wrap(score_core(as<cube>(S_), as<vec>(v_)));
}
extern "C" SEXP med_accum_cpp(SEXP S_, SEXP w_) {
  return wrap(accum_core(as<cube>(S_), as<vec>(w_)));
}
extern "C" SEXP med_eigen_solve_cpp(SEXP A_, SEXP H_) {
  return wrap(eigen_solve_core(as<mat>(A_), as<mat>(H_)));
}
extern "C" SEXP med_cov_cpp(SEXP M_) {
  return wrap(cov_core(as<List>(M_)));
}
