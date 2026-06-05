// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Moore-Penrose pseudo-inverse matching MASS::ginv (SVD, tol = sqrt(eps)*max(sv)).
// Used instead of inv_sympd so a singular Omega (random-slope variance collapsing
// to ~0) does not throw -- matching V3's ginv(Omega) behaviour exactly.
static arma::mat ginv_arma(const arma::mat& X) {
  arma::mat U, V; arma::vec s;
  arma::svd(U, s, V, X);
  double tol = (s.n_elem ? s.max() : 0.0) * std::sqrt(arma::datum::eps);
  arma::vec si(s.n_elem, arma::fill::zeros);
  for (arma::uword i = 0; i < s.n_elem; i++) if (s[i] > tol) si[i] = 1.0 / s[i];
  return V * arma::diagmat(si) * U.t();
}

/***********************************************
 * MLCAP C++ kernels for computational bottlenecks
 *
 * These functions replace the inner R for-loops
 * in the block coordinate descent algorithm.
 ***********************************************/

// =============================================
// Kernel 1: Compute scores = gamma_i' S_ij gamma_i
//   score[i,j] = gamma_rnd[i,]' %*% Y.cov[[i]][,,j] %*% gamma_rnd[i,]
// =============================================
// [[Rcpp::export]]
NumericMatrix compute_scores_cpp(List Y_cov, NumericMatrix gamma_rnd,
                                 IntegerVector nvec) {
  int m = Y_cov.size();
  int nmax = max(nvec);

  NumericMatrix score(m, nmax);
  std::fill(score.begin(), score.end(), NA_REAL);

  for (int i = 0; i < m; i++) {
    NumericVector arr = Y_cov[i];
    IntegerVector dims = arr.attr("dim");
    int p = dims[0];
    int ni = nvec[i];

    // Map gamma row to arma vec
    vec gi(p);
    for (int k = 0; k < p; k++) gi(k) = gamma_rnd(i, k);

    for (int j = 0; j < ni; j++) {
      // Extract S_ij from the 3D array (column-major: p x p x n_i)
      mat Sij(p, p);
      int offset = j * p * p;
      for (int c = 0; c < p; c++)
        for (int r = 0; r < p; r++)
          Sij(r, c) = arr[offset + c * p + r];

      score(i, j) = as_scalar(gi.t() * Sij * gi);
    }
  }

  return score;
}

// =============================================
// Kernel 2: Objective function
//   l1 = sum((fit + score*exp(-fit)) * Tmat/2)
//   l2 = sum((beta0_rnd - beta0)^2 / (2*sigma2) + log(sigma2)/2)
//   l3 = sum over i of: log(det(Omega))/2 + (beta2_rnd[i,]-beta2)' Omega^{-1} (beta2_rnd[i,]-beta2)/2
//   l4 = sum(-log(Cp) - gamma_rnd %*% gamma * kappa)
// =============================================
// [[Rcpp::export]]
double obj_func_cpp(Nullable<List> X1_,
                    Nullable<List> X2_,
                    NumericMatrix Tmat,
                    NumericMatrix score,
                    NumericMatrix gamma_rnd,
                    NumericVector gamma_pop,
                    double kappa,
                    NumericVector beta0_rnd,
                    Nullable<NumericVector> beta1_,
                    Nullable<NumericMatrix> beta2_rnd_,
                    double beta0,
                    Nullable<NumericVector> beta2_,
                    double sigma2,
                    Nullable<NumericMatrix> Omega_,
                    IntegerVector nvec) {
  int m = beta0_rnd.size();
  int p = gamma_pop.size();

  bool has_X1 = X1_.isNotNull() && beta1_.isNotNull();
  bool has_X2 = X2_.isNotNull() && beta2_rnd_.isNotNull();

  List X1, X2;
  NumericVector beta1;
  NumericMatrix beta2_rnd;
  NumericVector beta2;
  NumericMatrix Omega;

  if (has_X1) { X1 = as<List>(X1_); beta1 = as<NumericVector>(beta1_); }
  if (has_X2) {
    X2 = as<List>(X2_);
    beta2_rnd = as<NumericMatrix>(beta2_rnd_);
    beta2 = as<NumericVector>(beta2_);
    Omega = as<NumericMatrix>(Omega_);
  }

  // l1: conditional likelihood
  double l1 = 0.0;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < nvec[i]; j++) {
      double fit_ij = beta0_rnd[i];
      if (has_X1) {
        NumericMatrix X1i = as<NumericMatrix>(X1[i]);
        for (int k = 0; k < X1i.ncol(); k++)
          fit_ij += X1i(j, k) * beta1[k];
      }
      if (has_X2) {
        NumericMatrix X2i = as<NumericMatrix>(X2[i]);
        for (int k = 0; k < X2i.ncol(); k++)
          fit_ij += X2i(j, k) * beta2_rnd(i, k);
      }
      l1 += (fit_ij + score(i, j) * exp(-fit_ij)) * Tmat(i, j) / 2.0;
    }
  }

  // l2: random intercept
  double l2 = 0.0;
  for (int i = 0; i < m; i++) {
    double d = beta0_rnd[i] - beta0;
    l2 += d * d / (2.0 * sigma2) + log(sigma2) / 2.0;
  }

  // l3: random slope
  double l3 = 0.0;
  if (has_X2) {
    int q2 = beta2_rnd.ncol();
    mat Om(Omega.begin(), q2, q2, false);
    mat Om_inv = ginv_arma(Om);          // pseudo-inverse (matches V3 ginv; no throw)
    double log_det_Om = log(det(Om));

    for (int i = 0; i < m; i++) {
      vec d(q2);
      for (int k = 0; k < q2; k++)
        d(k) = beta2_rnd(i, k) - beta2[k];
      l3 += log_det_Om / 2.0 + as_scalar(d.t() * Om_inv * d) / 2.0;
    }
  }

  // l4: vMF likelihood
  // C_p(kappa) = kappa^(p/2-1) / ((2*pi)^(p/2) * I_{p/2-1}(kappa))
  // R::bessel_i(x, nu, 1.0) returns unscaled I_nu(x)
  double l4 = 0.0;
  double nu = (double)p / 2.0 - 1.0;
  double log_Cp = nu * log(kappa) - ((double)p / 2.0) * log(2.0 * M_PI) - log(R::bessel_i(kappa, nu, 1.0));
  for (int i = 0; i < m; i++) {
    double dot = 0.0;
    for (int k = 0; k < p; k++)
      dot += gamma_rnd(i, k) * gamma_pop[k];
    l4 += -log_Cp - dot * kappa;
  }

  return l1 + l2 + l3 + l4;
}

// =============================================
// Kernel 3: Newton-Raphson update for beta0_rnd
//   For each i: beta0_rnd[i] -= gradient / hessian
// =============================================
// [[Rcpp::export]]
NumericVector update_beta0_rnd_cpp(NumericMatrix score,
                                   Nullable<List> X1_,
                                   Nullable<List> X2_,
                                   NumericMatrix Tmat,
                                   NumericVector beta0_rnd,
                                   Nullable<NumericVector> beta1_,
                                   Nullable<NumericMatrix> beta2_rnd_,
                                   double beta0, double sigma2,
                                   IntegerVector nvec) {
  int m = beta0_rnd.size();
  bool has_X1 = X1_.isNotNull() && beta1_.isNotNull();
  bool has_X2 = X2_.isNotNull() && beta2_rnd_.isNotNull();

  List X1, X2;
  NumericVector beta1;
  NumericMatrix beta2_rnd;
  if (has_X1) { X1 = as<List>(X1_); beta1 = as<NumericVector>(beta1_); }
  if (has_X2) { X2 = as<List>(X2_); beta2_rnd = as<NumericMatrix>(beta2_rnd_); }

  NumericVector result(m);

  for (int i = 0; i < m; i++) {
    double pt1 = 0.0, pt2 = 0.0;

    for (int j = 0; j < nvec[i]; j++) {
      double fit_ij = beta0_rnd[i];
      if (has_X1) {
        NumericMatrix X1i = as<NumericMatrix>(X1[i]);
        for (int k = 0; k < X1i.ncol(); k++)
          fit_ij += X1i(j, k) * beta1[k];
      }
      if (has_X2) {
        NumericMatrix X2i = as<NumericMatrix>(X2[i]);
        for (int k = 0; k < X2i.ncol(); k++)
          fit_ij += X2i(j, k) * beta2_rnd(i, k);
      }

      double s_exp = score(i, j) * exp(-fit_ij);
      double T_half = Tmat(i, j) / 2.0;
      pt1 += (1.0 - s_exp) * T_half;
      pt2 += s_exp * T_half;
    }

    pt1 += (beta0_rnd[i] - beta0) / sigma2;
    pt2 += 1.0 / sigma2;

    result[i] = beta0_rnd[i] - pt1 / pt2;
  }

  return result;
}

// =============================================
// Kernel 4: Newton-Raphson update for beta1
// =============================================
// [[Rcpp::export]]
NumericVector update_beta1_cpp(NumericMatrix score,
                               List X1,
                               Nullable<List> X2_,
                               NumericMatrix Tmat,
                               NumericVector beta0_rnd,
                               NumericVector beta1,
                               Nullable<NumericMatrix> beta2_rnd_,
                               IntegerVector nvec) {
  int m = beta0_rnd.size();
  int q1 = beta1.size();
  bool has_X2 = X2_.isNotNull() && beta2_rnd_.isNotNull();

  List X2;
  NumericMatrix beta2_rnd;
  if (has_X2) { X2 = as<List>(X2_); beta2_rnd = as<NumericMatrix>(beta2_rnd_); }

  vec grad(q1, fill::zeros);
  mat hess(q1, q1, fill::zeros);

  for (int i = 0; i < m; i++) {
    NumericMatrix X1i = as<NumericMatrix>(X1[i]);
    int ni = nvec[i];

    for (int j = 0; j < ni; j++) {
      double fit_ij = beta0_rnd[i];
      for (int k = 0; k < q1; k++)
        fit_ij += X1i(j, k) * beta1[k];
      if (has_X2) {
        NumericMatrix X2i = as<NumericMatrix>(X2[i]);
        for (int k = 0; k < X2i.ncol(); k++)
          fit_ij += X2i(j, k) * beta2_rnd(i, k);
      }

      double s_exp = score(i, j) * exp(-fit_ij);
      double T_half = Tmat(i, j) / 2.0;

      // gradient
      double w1 = (1.0 - s_exp) * T_half;
      for (int k = 0; k < q1; k++)
        grad(k) += w1 * X1i(j, k);

      // hessian
      double w2 = s_exp * T_half;
      for (int k1 = 0; k1 < q1; k1++)
        for (int k2 = 0; k2 < q1; k2++)
          hess(k1, k2) += w2 * X1i(j, k1) * X1i(j, k2);
    }
  }

  vec beta1_v(q1);
  for (int k = 0; k < q1; k++) beta1_v(k) = beta1[k];

  vec result = beta1_v - pinv(hess) * grad;

  NumericVector out(q1);
  for (int k = 0; k < q1; k++) out[k] = result(k);
  return out;
}

// =============================================
// Kernel 5: Newton-Raphson update for beta2_rnd
// =============================================
// [[Rcpp::export]]
NumericMatrix update_beta2_rnd_cpp(NumericMatrix score,
                                   Nullable<List> X1_,
                                   List X2,
                                   NumericMatrix Tmat,
                                   NumericVector beta0_rnd,
                                   Nullable<NumericVector> beta1_,
                                   NumericMatrix beta2_rnd,
                                   NumericVector beta2,
                                   NumericMatrix Omega,
                                   IntegerVector nvec) {
  int m = beta0_rnd.size();
  int q2 = beta2.size();
  bool has_X1 = X1_.isNotNull() && beta1_.isNotNull();

  List X1;
  NumericVector beta1;
  if (has_X1) { X1 = as<List>(X1_); beta1 = as<NumericVector>(beta1_); }

  mat Om(Omega.begin(), q2, q2, false);
  mat Om_inv = pinv(Om);

  NumericMatrix result(m, q2);

  for (int i = 0; i < m; i++) {
    NumericMatrix X2i = as<NumericMatrix>(X2[i]);
    int ni = nvec[i];

    vec grad(q2, fill::zeros);
    mat hess(q2, q2, fill::zeros);

    for (int j = 0; j < ni; j++) {
      double fit_ij = beta0_rnd[i];
      if (has_X1) {
        NumericMatrix X1i = as<NumericMatrix>(X1[i]);
        for (int k = 0; k < X1i.ncol(); k++)
          fit_ij += X1i(j, k) * beta1[k];
      }
      for (int k = 0; k < q2; k++)
        fit_ij += X2i(j, k) * beta2_rnd(i, k);

      double s_exp = score(i, j) * exp(-fit_ij);
      double T_half = Tmat(i, j) / 2.0;

      double w1 = (1.0 - s_exp) * T_half;
      for (int k = 0; k < q2; k++)
        grad(k) += w1 * X2i(j, k);

      double w2 = s_exp * T_half;
      for (int k1 = 0; k1 < q2; k1++)
        for (int k2 = 0; k2 < q2; k2++)
          hess(k1, k2) += w2 * X2i(j, k1) * X2i(j, k2);
    }

    // Add prior terms
    vec d(q2);
    for (int k = 0; k < q2; k++)
      d(k) = beta2_rnd(i, k) - beta2[k];
    grad += Om_inv * d;
    hess += Om_inv;

    vec b2i(q2);
    for (int k = 0; k < q2; k++) b2i(k) = beta2_rnd(i, k);

    vec upd = b2i - pinv(hess) * grad;
    for (int k = 0; k < q2; k++) result(i, k) = upd(k);
  }

  return result;
}

// =============================================
// Kernel 6: Compute A_i matrix for gamma solve
//   A_i = sum_j (T_ij/2) * exp(-mu_ij) * S_ij
// =============================================
// [[Rcpp::export]]
arma::mat compute_Ai_cpp(NumericVector Y_cov_i,
                          IntegerVector Y_cov_dim,
                          Nullable<NumericMatrix> X1i_,
                          Nullable<NumericMatrix> X2i_,
                          NumericVector Tmat_i,
                          double beta0_rnd_i,
                          Nullable<NumericVector> beta1_,
                          Nullable<NumericVector> beta2_rnd_i_,
                          int ni) {
  int p = Y_cov_dim[0];

  bool has_X1 = X1i_.isNotNull() && beta1_.isNotNull();
  bool has_X2 = X2i_.isNotNull() && beta2_rnd_i_.isNotNull();

  NumericMatrix X1i, X2i;
  NumericVector beta1, beta2_rnd_i;
  if (has_X1) { X1i = as<NumericMatrix>(X1i_); beta1 = as<NumericVector>(beta1_); }
  if (has_X2) { X2i = as<NumericMatrix>(X2i_); beta2_rnd_i = as<NumericVector>(beta2_rnd_i_); }

  mat A(p, p, fill::zeros);

  for (int j = 0; j < ni; j++) {
    double fit_ij = beta0_rnd_i;
    if (has_X1)
      for (int k = 0; k < X1i.ncol(); k++)
        fit_ij += X1i(j, k) * beta1[k];
    if (has_X2)
      for (int k = 0; k < X2i.ncol(); k++)
        fit_ij += X2i(j, k) * beta2_rnd_i[k];

    double w = (Tmat_i[j] / 2.0) * exp(-fit_ij);

    // Extract S_ij and accumulate
    int offset = j * p * p;
    for (int c = 0; c < p; c++)
      for (int r = 0; r < p; r++)
        A(r, c) += w * Y_cov_i[offset + c * p + r];
  }

  return A;
}

// =============================================
// Kernel 7: Compute H_i matrix
//   H_i = sum_j T_ij * S_ij / sum_j T_ij
// =============================================
// [[Rcpp::export]]
arma::mat compute_Hi_cpp(NumericVector Y_cov_i,
                          IntegerVector Y_cov_dim,
                          NumericVector Tmat_i,
                          int ni, std::string H_type) {
  int p = Y_cov_dim[0];

  if (H_type == "Identity") {
    return eye(p, p);
  }

  mat H(p, p, fill::zeros);
  double T_sum = 0.0;

  for (int j = 0; j < ni; j++) {
    double Tij = Tmat_i[j];
    T_sum += Tij;
    int offset = j * p * p;
    for (int c = 0; c < p; c++)
      for (int r = 0; r < p; r++)
        H(r, c) += Tij * Y_cov_i[offset + c * p + r];
  }

  H /= T_sum;
  return H;
}

// =============================================================================
// extern "C" wrappers for the precompiled-.so + dyn.load()/.Call() loading path
// (V4 originally used sourceCpp; the website builds this to mlcap_kernels.so via
// build_kernels.R and calls these from LCAP_varGamma.R). Symbols are mlcap_-
// prefixed to avoid clashes with the CAP/HDCAP/LCAP kernels in the same process.
// Only the two functions actually called from R are wrapped.
// =============================================================================
extern "C" SEXP mlcap_compute_scores_cpp(SEXP Ycov_, SEXP grnd_, SEXP nvec_) {
  try {
    return compute_scores_cpp(as<List>(Ycov_), as<NumericMatrix>(grnd_),
                              as<IntegerVector>(nvec_));
  } catch (std::exception& e) { Rf_error("%s", e.what()); }
  catch (...) { Rf_error("mlcap_compute_scores_cpp: unknown C++ exception"); }
  return R_NilValue;
}

extern "C" SEXP mlcap_obj_func_cpp(SEXP X1_, SEXP X2_, SEXP Tmat_, SEXP score_,
    SEXP grnd_, SEXP gpop_, SEXP kappa_, SEXP b0rnd_, SEXP beta1_, SEXP b2rnd_,
    SEXP beta0_, SEXP beta2_, SEXP sigma2_, SEXP Omega_, SEXP nvec_) {
  try {
    double r = obj_func_cpp(Nullable<List>(X1_), Nullable<List>(X2_),
        as<NumericMatrix>(Tmat_), as<NumericMatrix>(score_), as<NumericMatrix>(grnd_),
        as<NumericVector>(gpop_), as<double>(kappa_), as<NumericVector>(b0rnd_),
        Nullable<NumericVector>(beta1_), Nullable<NumericMatrix>(b2rnd_),
        as<double>(beta0_), Nullable<NumericVector>(beta2_), as<double>(sigma2_),
        Nullable<NumericMatrix>(Omega_), as<IntegerVector>(nvec_));
    return wrap(r);
  } catch (std::exception& e) { Rf_error("%s", e.what()); }
  catch (...) { Rf_error("mlcap_obj_func_cpp: unknown C++ exception"); }
  return R_NilValue;
}
