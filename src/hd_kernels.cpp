// =============================================================================
// HDCAP/V7 -- RcppArmadillo kernels (acceleration of HDCAP/V6, no sparsity)
// -----------------------------------------------------------------------------
// Compiled to a shared object (build_kernels.R) and loaded with dyn.load();
// the R wrappers in CAP_HD.R call these via .Call(). The math mirrors HDCAP/V6
// (cov.ls linear shrinkage, cap_beta, cap_D1, obj.func) exactly; only the hot
// loops move to C++. With covariance shrinkage, every per-subject quantity
// reduces to the score gamma' Sigma_i gamma and gg = gamma'gamma, so the beta
// and objective kernels need only those (not the full Sigma).
// =============================================================================
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// ---- linear shrinkage covariance (Ledoit-Wolf style; matches cov.ls) --------
static mat cov_ls_core(const mat& Yin) {
  const uword nT = Yin.n_rows, p = Yin.n_cols;
  mat Xc = Yin.each_row() - mean(Yin, 0);     // demean
  mat S  = (Xc.t() * Xc) / (double)nT;        // cov * (n-1)/n
  mat Ip = eye(p, p);
  double m = trace(S) / (double)p;            // norm.F.std(S, I)
  mat A = S - m * Ip;
  double d2 = accu(A % A) / (double)p;
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

// ---- shrinkage weights (cap_beta / cap_D1 style: /Tvec, /(psi2+phi2)) -------
static void shrink_weights(const vec& score, const vec& Xbeta, const vec& Tvec,
                           double gg, double& rho1, double& rho2) {
  double mu = mean(exp(Xbeta)) / gg;
  double delta2 = mean(square(score - mu * gg));
  double psi2 = std::min(mean(square(score - exp(Xbeta)) / Tvec), delta2);
  double phi2 = delta2 - psi2;
  rho1 = psi2 * mu / (psi2 + phi2);
  rho2 = phi2 / (psi2 + phi2);
}

// ---- beta given gamma (scores) ; identical recursion to V6 cap_beta ---------
static vec cap_beta_core(const vec& score, double gg, const mat& X, const vec& Tvec,
                         bool shrink, vec beta, int max_itr, double tol,
                         double& rho1, double& rho2, int& iter) {
  rho1 = 0.0; rho2 = 1.0;
  if (shrink) shrink_weights(score, X * beta, Tvec, gg, rho1, rho2);
  int s = 0; double diff = 1e9;
  while (s <= max_itr && diff > tol) {
    s++;
    vec sk = shrink ? (rho1 * gg + rho2 * score) : score;   // gamma' Stmp_i gamma
    vec e  = exp(-(X * beta));
    vec w  = Tvec % sk % e;
    mat Q1 = X.t() * (X.each_col() % w);
    vec Q2 = X.t() * (Tvec % (1.0 - sk % e));
    vec beta_new = beta - pinv(Q1) * Q2;
    if (shrink) shrink_weights(score, X * beta_new, Tvec, gg, rho1, rho2);
    diff = max(abs(beta_new - beta));
    beta = beta_new;
  }
  iter = s;
  return beta;
}

// ---- generalized smallest eigenvector (matches gamma.solve) -----------------
static vec gamma_solve_core(const mat& A, const mat& H) {
  vec d; mat U;
  eig_sym(d, U, symmatu(H));
  mat Dis = diagmat(1.0 / sqrt(d));
  mat M = Dis * U.t() * A * U * Dis;
  M = symmatu(0.5 * (M + M.t()));
  vec mv; mat V;
  eig_sym(mv, V, M);
  return U * Dis * V.col(0);
}

// ---- one HDCAP direction (method "CAP"); identical recursion to V6 cap_D1 ----
// score0 (from the INITIAL gamma) is held fixed in the weight update, exactly
// as in V6; the beta update uses the current gamma's score.
static vec cap_d1_core(const cube& Sig, const mat& X, const vec& Tvec, const mat& H,
                       vec gamma, bool shrink, vec beta, int max_itr, double tol,
                       double& rho1, double& rho2, int& iter) {
  const uword n = Sig.n_slices, p = Sig.n_rows;
  vec score0(n);
  for (uword i = 0; i < n; i++) score0[i] = as_scalar(gamma.t() * Sig.slice(i) * gamma);
  double gg0 = dot(gamma, gamma);
  rho1 = 0.0; rho2 = 1.0;
  if (shrink) shrink_weights(score0, X * beta, Tvec, gg0, rho1, rho2);

  int s = 0; double diff = 1e9;
  while (s <= max_itr && diff > tol) {
    s++;
    double cur_gg = dot(gamma, gamma);
    vec cur_score(n);
    for (uword i = 0; i < n; i++) cur_score[i] = as_scalar(gamma.t() * Sig.slice(i) * gamma);
    vec sk = shrink ? (rho1 * cur_gg + rho2 * cur_score) : cur_score;
    // beta update
    vec e  = exp(-(X * beta));
    vec w  = Tvec % sk % e;
    mat Q1 = X.t() * (X.each_col() % w);
    vec Q2 = X.t() * (Tvec % (1.0 - sk % e));
    vec beta_new = beta - pinv(Q1) * Q2;
    // weight update (uses fixed score0 with current gg, as in V6)
    if (shrink) shrink_weights(score0, X * beta_new, Tvec, cur_gg, rho1, rho2);
    // gamma update
    vec e2 = exp(-(X * beta_new));
    vec ww = Tvec % e2;
    mat S1(p, p, fill::zeros);
    for (uword i = 0; i < n; i++) S1 += ww[i] * Sig.slice(i);
    if (shrink) S1 = rho1 * sum(ww) * eye(p, p) + rho2 * S1;
    vec gamma_new = gamma_solve_core(S1, H);
    diff = std::max(max(abs(gamma_new - gamma)), max(abs(beta_new - beta)));
    gamma = gamma_new; beta = beta_new;
  }
  iter = s;
  return gamma;   // converged, unnormalized
}

// ---- HDCAP objective (obj.func: uses ^2 and /delta2 for the weights) --------
static double obj_core(const vec& score, double gg, const mat& X, const vec& Tvec,
                       const vec& beta) {
  vec Xbeta = X * beta;
  double mu = mean(exp(Xbeta)) / gg;
  double delta2 = mean(square(score - mu * gg));
  double psi2 = std::min(mean(square(score - exp(Xbeta))), delta2);
  double phi2 = delta2 - psi2;
  double rho1 = psi2 * mu / delta2;
  double rho2 = phi2 / delta2;
  vec sk = rho1 * gg + rho2 * score;
  return sum(Tvec % (Xbeta + exp(-Xbeta) % sk));
}

static vec score_core(const cube& Sig, const vec& gamma) {
  vec s(Sig.n_slices);
  for (uword i = 0; i < Sig.n_slices; i++) s[i] = as_scalar(gamma.t() * Sig.slice(i) * gamma);
  return s;
}

// ---- extern "C" wrappers ----------------------------------------------------
// Symbols are hd_-prefixed to avoid colliding with the CAP/V6 kernel (which is
// loaded in the same R process and exports score_cpp / objfunc_cpp too).
extern "C" SEXP hd_covls_cpp(SEXP Y_) {
  return wrap(cov_ls_core(as<mat>(Y_)));
}

extern "C" SEXP hd_score_cpp(SEXP Sig_, SEXP gamma_) {
  return wrap(score_core(as<cube>(Sig_), as<vec>(gamma_)));
}

extern "C" SEXP hd_capbeta_cpp(SEXP score_, SEXP gg_, SEXP X_, SEXP Tvec_, SEXP shrink_,
                            SEXP beta0_, SEXP maxitr_, SEXP tol_) {
  vec score = as<vec>(score_); double gg = as<double>(gg_);
  mat X = as<mat>(X_); vec Tvec = as<vec>(Tvec_); bool sh = as<bool>(shrink_);
  vec beta0 = as<vec>(beta0_);
  double rho1, rho2; int iter;
  vec beta = cap_beta_core(score, gg, X, Tvec, sh, beta0, as<int>(maxitr_),
                           as<double>(tol_), rho1, rho2, iter);
  return wrap(List::create(_["beta"] = beta, _["rho1"] = rho1, _["rho2"] = rho2,
                           _["iter"] = iter));
}

extern "C" SEXP hd_capd1_cpp(SEXP Sig_, SEXP X_, SEXP Tvec_, SEXP H_, SEXP gamma0_,
                          SEXP shrink_, SEXP beta0_, SEXP maxitr_, SEXP tol_) {
  cube Sig = as<cube>(Sig_); mat X = as<mat>(X_); vec Tvec = as<vec>(Tvec_);
  mat H = as<mat>(H_); vec g0 = as<vec>(gamma0_); bool sh = as<bool>(shrink_);
  vec beta0 = as<vec>(beta0_);
  double rho1, rho2; int iter;
  vec g = cap_d1_core(Sig, X, Tvec, H, g0, sh, beta0, as<int>(maxitr_),
                      as<double>(tol_), rho1, rho2, iter);
  return wrap(List::create(_["gamma"] = g, _["rho1"] = rho1, _["rho2"] = rho2,
                           _["iter"] = iter));
}

extern "C" SEXP hd_objfunc_cpp(SEXP score_, SEXP gg_, SEXP X_, SEXP Tvec_, SEXP beta_) {
  return wrap(obj_core(as<vec>(score_), as<double>(gg_), as<mat>(X_),
                       as<vec>(Tvec_), as<vec>(beta_)));
}
