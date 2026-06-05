// =============================================================================
// CAP-clustering/V4 -- RcppArmadillo kernels
// -----------------------------------------------------------------------------
// Compiled to a shared object (build_kernels.R) and loaded with dyn.load(); the
// R wrappers in CAP_Cluster.R call these via .Call(). The math mirrors V3
// (per-subject second-moment Smat, projected scores, the gamma-update weighted
// accumulation, and eigen.solve) exactly. The multinomial membership fit
// (brglm2::brmultinom) stays in R.
//
// Symbols are cluster_-prefixed to avoid colliding with the other CAP kernels.
// =============================================================================
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// per-subject SECOND MOMENT cube + sizes (== t(Y_i) %*% Y_i / T_i, NO centering)
static List smat_core(const List& Y) {
  const int n = Y.size();
  mat y0 = as<mat>(Y[0]);
  const uword p = y0.n_cols;
  cube S(p, p, n);
  vec Tvec(n);
  for (int i = 0; i < n; i++) {
    mat Yi = as<mat>(Y[i]);
    double Ti = Yi.n_rows;
    Tvec[i] = Ti;
    S.slice(i) = (Yi.t() * Yi) / Ti;
  }
  return List::create(_["Smat"] = S, _["Tvec"] = Tvec);
}

// projected scores  v' S_i v  for all i  (== apply(Smat, 3, ...))
static vec score_core(const cube& S, const vec& v) {
  vec o(S.n_slices);
  for (uword i = 0; i < S.n_slices; i++) o[i] = as_scalar(v.t() * S.slice(i) * v);
  return o;
}

// weighted accumulation  sum_i w_i S_i  (gamma-update A matrix, H matrix)
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

// ---- extern "C" wrappers ----------------------------------------------------
extern "C" SEXP cluster_smat_cpp(SEXP Y_) {
  return smat_core(as<List>(Y_));
}
extern "C" SEXP cluster_score_cpp(SEXP S_, SEXP v_) {
  return wrap(score_core(as<cube>(S_), as<vec>(v_)));
}
extern "C" SEXP cluster_accum_cpp(SEXP S_, SEXP w_) {
  return wrap(accum_core(as<cube>(S_), as<vec>(w_)));
}
extern "C" SEXP cluster_eigen_solve_cpp(SEXP A_, SEXP H_) {
  return wrap(eigen_solve_core(as<mat>(A_), as<mat>(H_)));
}
