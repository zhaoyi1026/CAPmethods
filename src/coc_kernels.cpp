// =============================================================================
// CAP-CoC/V3 -- RcppArmadillo kernels for Covariance-on-Covariance regression
// -----------------------------------------------------------------------------
// Compiled to a shared object (build_kernels.R) and loaded with dyn.load(); the
// R wrappers in COCReg.R call these via .Call(). The math mirrors CAP-CoC/V2
// (cov.ls, cov.sk.x, cov.sk.y, eigen.solve, and the projected-score / weighted-
// accumulation helpers) exactly; only the hot loops move to C++.
//
// Symbols are coc_-prefixed to avoid colliding with the CAP/HDCAP kernels that
// may be loaded in the same R process.
// =============================================================================
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// standardized squared Frobenius norm: norm.F.std(A)^2 = sum(A_ij^2) / nrow(A)
static inline double nf2(const mat& A) { return accu(A % A) / (double)A.n_rows; }

// ---- linear shrinkage of one data matrix (== cov.ls) ------------------------
static mat cov_ls_core(const mat& Yin) {
  const uword nT = Yin.n_rows, p = Yin.n_cols;
  mat Xc = Yin.each_row() - mean(Yin, 0);
  mat S  = (Xc.t() * Xc) / (double)nT;
  mat Ip = eye(p, p);
  double m = trace(S) / (double)p;
  double d2 = nf2(S - m * Ip);
  double normS2 = accu(S % S);
  double acc = 0.0;
  for (uword i = 0; i < nT; i++) {
    rowvec x = Xc.row(i);
    double xx = dot(x, x);
    double xSx = as_scalar(x * S * x.t());
    acc += (xx * xx - 2.0 * xSx + normS2) / (double)p;
  }
  double b2bar = (acc / (double)nT) / (double)nT;
  double b2 = std::min(b2bar, d2);
  double a2 = d2 - b2;
  return b2 * m * Ip / d2 + a2 * S / d2;
}

// ---- shrinkage covariance cube for the predictor X (== cov.sk.x) ------------
static cube cov_sk_x_core(const List& X) {
  const int n = X.size();
  mat x0 = as<mat>(X[0]);
  const uword p = x0.n_cols;
  cube Sx(p, p, n);
  vec nx(n), nuv(n);
  std::vector<mat> Xc(n);
  for (int i = 0; i < n; i++) {
    mat Xi = as<mat>(X[i]);
    nx[i] = Xi.n_rows;
    mat xc = Xi.each_row() - mean(Xi, 0);
    Xc[i] = xc;
    mat S = (xc.t() * xc) / nx[i];
    Sx.slice(i) = S;
    nuv[i] = trace(S) / (double)p;          // norm.F.std(I, Sx_i)
  }
  double nuhat = mean(nuv);
  vec tau2(n), omega2(n), eps2(n);
  for (int i = 0; i < n; i++) {
    mat S = Sx.slice(i);
    tau2[i] = nf2(S - nuhat * eye(p, p));
    double normS2 = accu(S % S);
    double otmp = 0.0;
    const uword nxi = (uword)nx[i];
    for (uword ss = 0; ss < nxi; ss++) {
      rowvec x = Xc[i].row(ss);
      double xx = dot(x, x);
      double xSx = as_scalar(x * S * x.t());
      otmp += (xx * xx - 2.0 * xSx + normS2) / (double)p;
    }
    otmp /= (nx[i] * nx[i]);
    omega2[i] = std::min(otmp, tau2[i]);
    eps2[i]   = tau2[i] - omega2[i];
  }
  double tauh = mean(tau2), omh = mean(omega2), eph = mean(eps2);
  cube out(p, p, n);
  mat Ip = eye(p, p);
  for (int i = 0; i < n; i++)
    out.slice(i) = (omh / tauh) * nuhat * Ip + (eph / tauh) * Sx.slice(i);
  return out;
}

// ---- shrinkage covariance cube for the outcome Y given gamma, kappa --------
static cube cov_sk_y_core(const List& Y, const vec& gamma, const vec& kappa) {
  const int n = Y.size();
  mat y0 = as<mat>(Y[0]);
  const uword q = y0.n_cols;
  cube Sy(q, q, n);
  vec ny(n);
  for (int i = 0; i < n; i++) {
    mat Yi = as<mat>(Y[i]);
    ny[i] = Yi.n_rows;
    mat yc = Yi.each_row() - mean(Yi, 0);
    Sy.slice(i) = (yc.t() * yc) / ny[i];
  }
  double gg = dot(gamma, gamma);
  double muhat = mean(kappa) / gg;
  vec d2(n), psi2(n), phi2(n);
  for (int i = 0; i < n; i++) {
    double gs = as_scalar(gamma.t() * Sy.slice(i) * gamma);
    d2[i] = std::pow(gs - muhat * gg, 2.0);
    psi2[i] = std::min(std::pow(gs - kappa[i], 2.0) / ny[i], d2[i]);
    phi2[i] = d2[i] - psi2[i];
  }
  double dh = mean(d2), ph = mean(phi2), ps = mean(psi2);
  cube out(q, q, n);
  mat Iq = eye(q, q);
  for (int i = 0; i < n; i++)
    out.slice(i) = (ps / dh) * muhat * Iq + (ph / dh) * Sy.slice(i);
  return out;
}

// ---- projected scores  v' S_i v  for all i (== apply(S,3,...)) --------------
static vec score_core(const cube& S, const vec& v) {
  vec o(S.n_slices);
  for (uword i = 0; i < S.n_slices; i++) o[i] = as_scalar(v.t() * S.slice(i) * v);
  return o;
}

// ---- weighted accumulation  sum_i w_i S_i (== A1/A2/Hx/Hy loops) -----------
static mat accum_core(const cube& S, const vec& w) {
  mat A(S.n_rows, S.n_cols, fill::zeros);
  for (uword i = 0; i < S.n_slices; i++) A += w[i] * S.slice(i);
  return A;
}

// ---- generalized smallest eigenvector (== eigen.solve(A,H)) -----------------
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
extern "C" SEXP coc_cov_ls_cpp(SEXP Y_) {
  return wrap(cov_ls_core(as<mat>(Y_)));
}
extern "C" SEXP coc_cov_sk_x_cpp(SEXP X_) {
  return wrap(cov_sk_x_core(as<List>(X_)));
}
extern "C" SEXP coc_cov_sk_y_cpp(SEXP Y_, SEXP gamma_, SEXP kappa_) {
  return wrap(cov_sk_y_core(as<List>(Y_), as<vec>(gamma_), as<vec>(kappa_)));
}
extern "C" SEXP coc_score_cpp(SEXP S_, SEXP v_) {
  return wrap(score_core(as<cube>(S_), as<vec>(v_)));
}
extern "C" SEXP coc_accum_cpp(SEXP S_, SEXP w_) {
  return wrap(accum_core(as<cube>(S_), as<vec>(w_)));
}
extern "C" SEXP coc_eigen_solve_cpp(SEXP A_, SEXP H_) {
  return wrap(eigen_solve_core(as<mat>(A_), as<mat>(H_)));
}
