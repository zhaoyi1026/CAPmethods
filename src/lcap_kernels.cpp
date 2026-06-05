// =============================================================================
// LCAP/V6 -- RcppArmadillo kernels (acceleration of LCAP_gamma-invar/V5)
// -----------------------------------------------------------------------------
// Longitudinal CAP (gamma-invariant) with subject-level random effects. Compiled
// to a shared object (build_kernels.R) and loaded with dyn.load(); the R wrappers
// in Longitudinal_HDCAP.R call these via .Call(). The math mirrors V5 exactly:
//
//   * cov.ls         -- per visit Ledoit-Wolf linear shrinkage (same as HDCAP/V7)
//   * gamma.solve    -- generalized smallest eigenvector (same as HDCAP/V7)
//   * cap.cov_beta   -- the random-effects EM recursion (lcap_recar_cpp), the
//                       bootstrap hot loop; shrinkage-aware (cov.ls.const folded
//                       in analytically: a shrunk covariance's score is
//                       rho1*gamma'gamma + rho2*raw_score, so no p x p rebuild)
//   * cap_D1         -- one direction, NO-shrinkage path (lcap_capd1_cpp): joint
//                       gamma + random-effects update; gamma via gamma.solve
//   * obj.func       -- the longitudinal objective used to pick initializations
//
// Random covariances enter only through the scalar score gamma' Sigma_iv gamma
// (beta side) and through weighted slice sums (gamma side). With cov.beta.diag
// (the website's only mode) the q x q random-effect covariance is diagonal, so
// MASS::ginv reduces to reciprocal-or-zero with threshold sqrt(eps)*max(diag) --
// replicated exactly here. Symbols are lcap_-prefixed to avoid colliding with the
// CAP/V6 (score_cpp/...) and HDCAP/V7 (hd_*) kernels loaded in the same process.
// =============================================================================
#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// ---- MASS::ginv for a diagonal PSD matrix: reciprocal-or-zero -----------------
// MASS uses Positive <- d > max(tol*d[1], 0) with tol = sqrt(.Machine$double.eps).
static vec ginv_diag(const vec& d) {
  vec r(d.n_elem, fill::zeros);
  if (d.n_elem == 0) return r;
  double tol = std::sqrt(datum::eps) * d.max();
  for (uword i = 0; i < d.n_elem; i++) if (d[i] > tol) r[i] = 1.0 / d[i];
  return r;
}

// ---- linear shrinkage covariance (Ledoit-Wolf; matches cov.ls) ---------------
static mat cov_ls_core(const mat& Yin) {
  const uword nT = Yin.n_rows, p = Yin.n_cols;
  mat Xc = Yin.each_row() - mean(Yin, 0);
  mat S  = (Xc.t() * Xc) / (double)nT;
  mat Ip = eye(p, p);
  double m = trace(S) / (double)p;
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

// ---- generalized smallest eigenvector (matches gamma.solve) ------------------
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

// ---- per-visit scores gamma' Sigma_r gamma for a flat cube (p x p x M) -------
static vec score_core(const cube& Sig, const vec& gamma) {
  vec s(Sig.n_slices);
  for (uword r = 0; r < Sig.n_slices; r++)
    s[r] = as_scalar(gamma.t() * Sig.slice(r) * gamma);
  return s;
}

// ---- cov.ls.const shrinkage parameters (rho1, rho2) from RAW scores ----------
// Faithful to V5 cov.ls.const, INCLUDING its nested averaging (running totals
// divided by nV[i] inside the subject loop). Returns rho1, rho2.
static void covlsconst_rho(const vec& scoreRaw, double gg, const mat& X,
                           const vec& Tvec, const ivec& vcount,
                           const vec& beta0vec, const mat& beta1mat,
                           uword n, uword q, double& rho1, double& rho2) {
  // mu
  double mu = 0.0; uword idx = 0;
  for (uword i = 0; i < n; i++) {
    double acc = 0.0;
    for (int v = 0; v < vcount[i]; v++) {
      rowvec bi(q); bi[0] = beta0vec[i];
      for (uword j = 1; j < q; j++) bi[j] = beta1mat(i, j - 1);
      double u1 = dot(X.row(idx), bi);
      acc += std::exp(u1);
      idx++;
    }
    mu += acc;
    mu /= (double)vcount[i];           // nested averaging, as in V5
  }
  mu /= ((double)n * gg);
  // delta, psi, phi
  double hd = 0.0, hp = 0.0, hf = 0.0; idx = 0;
  for (uword i = 0; i < n; i++) {
    double sd = 0.0, sp = 0.0, sf = 0.0;
    for (int v = 0; v < vcount[i]; v++) {
      rowvec bi(q); bi[0] = beta0vec[i];
      for (uword j = 1; j < q; j++) bi[j] = beta1mat(i, j - 1);
      double u1 = dot(X.row(idx), bi);
      double sc = scoreRaw[idx];
      double dlt = (sc - mu * gg) * (sc - mu * gg);
      double psi = (sc - std::exp(u1)) * (sc - std::exp(u1)) / Tvec[idx];
      psi = std::min(psi, dlt);
      sd += dlt; sp += psi; sf += (dlt - psi);
      idx++;
    }
    hd += sd; hd /= (double)vcount[i];
    hp += sp; hp /= (double)vcount[i];
    hf += sf; hf /= (double)vcount[i];
  }
  hd /= (double)n; hp /= (double)n; hf /= (double)n;
  rho1 = hp * mu / hd;
  rho2 = hf / hd;
}

// ---- random-effects EM recursion (matches cap.cov_beta / cap_beta) -----------
// scoreWork: working scores used for the beta step (cov.ls-based init when
//            shrink, else raw sample scores). scoreRaw: raw sample scores (only
//            used when shrink, to recompute rho each iteration). gg = gamma'gamma.
// Updates the per-subject (beta0vec, beta1mat) and fixed (beta) + variance
// components (sigma2, Omega diag) until convergence. cov.beta.diag assumed.
static void re_recursion(vec scoreWork, const vec& scoreRaw, double gg,
                         const mat& X, const vec& Tvec, const ivec& vcount,
                         bool shrink, uword n, uword q,
                         vec& beta0vec, mat& beta1mat, vec& beta,
                         double& sigma2, vec& OmegaDiag,
                         int max_itr, double tol, int& iter) {
  // precompute subject row ranges
  ivec start(n); int acc = 0;
  for (uword i = 0; i < n; i++) { start[i] = acc; acc += vcount[i]; }

  int s = 0; double diff = 100.0;
  while (s <= max_itr && diff > tol) {
    s++;
    // diagonal random-effect covariance + its ginv
    vec dcov(q); dcov[0] = sigma2;
    for (uword j = 1; j < q; j++) dcov[j] = OmegaDiag[j - 1];
    vec gdiag = ginv_diag(dcov);

    mat bmat(n, q);
    for (uword i = 0; i < n; i++) {
      vec b_i(q); b_i[0] = beta0vec[i];
      for (uword j = 1; j < q; j++) b_i[j] = beta1mat(i, j - 1);
      vec pt1(q, fill::zeros); mat pt2(q, q, fill::zeros);
      for (int v = 0; v < vcount[i]; v++) {
        uword r = start[i] + v;
        rowvec xr = X.row(r);
        double u1 = dot(xr, b_i);
        double e = std::exp(-u1);
        double w = Tvec[r] / 2.0;
        pt1 += (xr.t() - scoreWork[r] * e * xr.t()) * w;
        pt2 += (scoreWork[r] * e * (xr.t() * xr)) * w;
      }
      pt1 += (gdiag % (b_i - beta));        // ginv(diag) * (b_i - beta_fixed)
      pt2.diag() += gdiag;
      bmat.row(i) = (b_i - solve(pt2, pt1)).t();
    }
    rowvec bmean = mean(bmat, 0);
    mat Bc = bmat.each_row() - bmean;
    mat bcov = (Bc.t() * Bc) / (double)n;   // cov * (n-1)/n

    vec beta_upd = bmean.t();
    double sigma2_upd = bcov(0, 0);
    vec Omega_upd(q - 1);
    for (uword j = 1; j < q; j++) Omega_upd[j - 1] = bcov(j, j);

    diff = max(abs(beta_upd - beta));

    // commit updates
    beta0vec = bmat.col(0);
    beta1mat = bmat.cols(1, q - 1);
    beta = beta_upd; sigma2 = sigma2_upd; OmegaDiag = Omega_upd;

    if (shrink) {
      double r1, r2;
      covlsconst_rho(scoreRaw, gg, X, Tvec, vcount, beta0vec, beta1mat,
                     n, q, r1, r2);
      scoreWork = r1 * gg + r2 * scoreRaw;
    }
  }
  iter = s;
}

// ---- one direction, NO-shrinkage (matches cap_D1 / cap.cov_D1 method="CAP") ---
// Joint gamma + random-effects update. Score is recomputed from the CURRENT
// gamma each iteration; H is fixed (Hraw). Returns the unnormalized gamma; the R
// wrapper normalizes, sign-fixes, and re-estimates beta (V5 discards these betas).
static vec capd1_core(const cube& Sig, const mat& Hraw, const mat& X,
                      const vec& Tvec, const ivec& vcount, uword n, uword q,
                      vec gamma, vec beta0vec, mat beta1mat, vec beta,
                      double sigma2, vec OmegaDiag, int max_itr, double tol,
                      int& iter) {
  ivec start(n); int acc = 0;
  for (uword i = 0; i < n; i++) { start[i] = acc; acc += vcount[i]; }
  const uword M = Sig.n_slices, p = Sig.n_rows;

  int s = 0; double diff = 100.0;
  while (s <= max_itr && diff > tol) {
    s++;
    vec score(M);
    for (uword r = 0; r < M; r++) score[r] = as_scalar(gamma.t() * Sig.slice(r) * gamma);

    vec dcov(q); dcov[0] = sigma2;
    for (uword j = 1; j < q; j++) dcov[j] = OmegaDiag[j - 1];
    vec gdiag = ginv_diag(dcov);

    mat bmat(n, q);
    for (uword i = 0; i < n; i++) {
      vec b_i(q); b_i[0] = beta0vec[i];
      for (uword j = 1; j < q; j++) b_i[j] = beta1mat(i, j - 1);
      vec pt1(q, fill::zeros); mat pt2(q, q, fill::zeros);
      for (int v = 0; v < vcount[i]; v++) {
        uword r = start[i] + v;
        rowvec xr = X.row(r);
        double u1 = dot(xr, b_i);
        double e = std::exp(-u1);
        double w = Tvec[r] / 2.0;
        pt1 += (xr.t() - score[r] * e * xr.t()) * w;
        pt2 += (score[r] * e * (xr.t() * xr)) * w;
      }
      pt1 += (gdiag % (b_i - beta));
      pt2.diag() += gdiag;
      bmat.row(i) = (b_i - solve(pt2, pt1)).t();
    }
    rowvec bmean = mean(bmat, 0);
    mat Bc = bmat.each_row() - bmean;
    mat bcov = (Bc.t() * Bc) / (double)n;
    vec beta_upd = bmean.t();
    double sigma2_upd = bcov(0, 0);
    vec Omega_upd(q - 1);
    for (uword j = 1; j < q; j++) Omega_upd[j - 1] = bcov(j, j);

    // gamma update: Amat = sum_r exp(-x_r'b_upd_i) * Sigma_r * T_r/2
    mat Amat(p, p, fill::zeros);
    for (uword i = 0; i < n; i++) {
      vec bu(q); bu[0] = bmat(i, 0);
      for (uword j = 1; j < q; j++) bu[j] = bmat(i, j);
      for (int v = 0; v < vcount[i]; v++) {
        uword r = start[i] + v;
        double u1 = dot(X.row(r), bu);
        Amat += std::exp(-u1) * Sig.slice(r) * (Tvec[r] / 2.0);
      }
    }
    vec gamma_upd = gamma_solve_core(Amat, Hraw);

    diff = max(abs(beta_upd - beta));     // V5: convergence on beta only

    beta0vec = bmat.col(0);
    beta1mat = bmat.cols(1, q - 1);
    beta = beta_upd; sigma2 = sigma2_upd; OmegaDiag = Omega_upd;
    gamma = gamma_upd;
  }
  iter = s;
  return gamma;
}

// ---- longitudinal objective (matches obj.func) -------------------------------
static double obj_core(const vec& score, const mat& X, const vec& Tvec,
                       const ivec& vcount, const vec& beta0vec, const mat& beta1mat,
                       double beta0, const vec& beta1, double sigma2,
                       const vec& OmegaDiag, uword n, uword q) {
  double ll1 = 0.0; uword idx = 0;
  for (uword i = 0; i < n; i++) {
    for (int v = 0; v < vcount[i]; v++) {
      rowvec bi(q); bi[0] = beta0vec[i];
      for (uword j = 1; j < q; j++) bi[j] = beta1mat(i, j - 1);
      double u1 = dot(X.row(idx), bi);
      ll1 += (u1 + score[idx] * std::exp(-u1)) * (Tvec[idx] / 2.0);
      idx++;
    }
  }
  double ll2 = 0.0;
  for (uword i = 0; i < n; i++)
    ll2 += std::log(sigma2) / 2.0 +
           (beta0vec[i] - beta0) * (beta0vec[i] - beta0) / (2.0 * sigma2);
  double detO = 1.0; for (uword j = 0; j < OmegaDiag.n_elem; j++) detO *= OmegaDiag[j];
  double logdet = std::log(detO);
  vec Oginv = ginv_diag(OmegaDiag);
  double ll3 = 0.0;
  for (uword i = 0; i < n; i++) {
    double quad = 0.0;
    for (uword j = 0; j < q - 1; j++) {
      double d = beta1mat(i, j) - beta1[j];
      quad += d * d * Oginv[j];
    }
    ll3 += logdet / 2.0 + quad / 2.0;
  }
  return ll1 + ll2 + ll3;
}

// ============================ extern "C" wrappers =============================
extern "C" SEXP lcap_covls_cpp(SEXP Y_) {
  return wrap(cov_ls_core(as<mat>(Y_)));
}

extern "C" SEXP lcap_gammasolve_cpp(SEXP A_, SEXP H_) {
  return wrap(gamma_solve_core(as<mat>(A_), as<mat>(H_)));
}

extern "C" SEXP lcap_score_cpp(SEXP Sig_, SEXP gamma_) {
  return wrap(score_core(as<cube>(Sig_), as<vec>(gamma_)));
}

extern "C" SEXP lcap_recar_cpp(SEXP scoreWork_, SEXP scoreRaw_, SEXP gg_, SEXP X_,
                               SEXP Tvec_, SEXP vcount_, SEXP shrink_, SEXP n_,
                               SEXP q_, SEXP beta0vec_, SEXP beta1mat_, SEXP beta_,
                               SEXP sigma2_, SEXP OmegaDiag_, SEXP maxitr_, SEXP tol_) {
  vec scoreWork = as<vec>(scoreWork_), scoreRaw = as<vec>(scoreRaw_);
  double gg = as<double>(gg_);
  mat X = as<mat>(X_); vec Tvec = as<vec>(Tvec_); ivec vcount = as<ivec>(vcount_);
  bool shrink = as<bool>(shrink_);
  uword n = as<uword>(n_), q = as<uword>(q_);
  vec beta0vec = as<vec>(beta0vec_); mat beta1mat = as<mat>(beta1mat_);
  vec beta = as<vec>(beta_); double sigma2 = as<double>(sigma2_);
  vec OmegaDiag = as<vec>(OmegaDiag_);
  int iter;
  re_recursion(scoreWork, scoreRaw, gg, X, Tvec, vcount, shrink, n, q,
               beta0vec, beta1mat, beta, sigma2, OmegaDiag,
               as<int>(maxitr_), as<double>(tol_), iter);
  return wrap(List::create(_["beta"] = beta, _["beta0vec"] = beta0vec,
                           _["beta1mat"] = beta1mat, _["sigma2"] = sigma2,
                           _["OmegaDiag"] = OmegaDiag, _["iter"] = iter));
}

extern "C" SEXP lcap_capd1_cpp(SEXP Sig_, SEXP H_, SEXP X_, SEXP Tvec_, SEXP vcount_,
                               SEXP n_, SEXP q_, SEXP gamma0_, SEXP beta0vec_,
                               SEXP beta1mat_, SEXP beta_, SEXP sigma2_,
                               SEXP OmegaDiag_, SEXP maxitr_, SEXP tol_) {
  cube Sig = as<cube>(Sig_); mat H = as<mat>(H_); mat X = as<mat>(X_);
  vec Tvec = as<vec>(Tvec_); ivec vcount = as<ivec>(vcount_);
  uword n = as<uword>(n_), q = as<uword>(q_);
  vec gamma = as<vec>(gamma0_); vec beta0vec = as<vec>(beta0vec_);
  mat beta1mat = as<mat>(beta1mat_); vec beta = as<vec>(beta_);
  double sigma2 = as<double>(sigma2_); vec OmegaDiag = as<vec>(OmegaDiag_);
  int iter;
  vec g = capd1_core(Sig, H, X, Tvec, vcount, n, q, gamma, beta0vec, beta1mat,
                     beta, sigma2, OmegaDiag, as<int>(maxitr_), as<double>(tol_), iter);
  return wrap(List::create(_["gamma"] = g, _["iter"] = iter));
}

extern "C" SEXP lcap_objfunc_cpp(SEXP score_, SEXP X_, SEXP Tvec_, SEXP vcount_,
                                 SEXP n_, SEXP q_, SEXP beta0vec_, SEXP beta1mat_,
                                 SEXP beta0_, SEXP beta1_, SEXP sigma2_, SEXP OmegaDiag_) {
  vec score = as<vec>(score_); mat X = as<mat>(X_); vec Tvec = as<vec>(Tvec_);
  ivec vcount = as<ivec>(vcount_); uword n = as<uword>(n_), q = as<uword>(q_);
  vec beta0vec = as<vec>(beta0vec_); mat beta1mat = as<mat>(beta1mat_);
  double beta0 = as<double>(beta0_); vec beta1 = as<vec>(beta1_);
  double sigma2 = as<double>(sigma2_); vec OmegaDiag = as<vec>(OmegaDiag_);
  return wrap(obj_core(score, X, Tvec, vcount, beta0vec, beta1mat, beta0, beta1,
                       sigma2, OmegaDiag, n, q));
}
