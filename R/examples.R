# Built-in example-data generators, one per method (mirroring the demo data in
# the Shiny app). Each returns a list with the inputs the corresponding wrapper
# expects plus a `truth` element holding the data-generating parameters, so the
# estimates can be checked. See the package vignette / Examples.md for full runs.

# symmetric PSD square root of  Phi %*% diag(eigs) %*% t(Phi)
.cap_rootSig <- function(Phi, eigs) {
  S <- Phi %*% diag(eigs, length(eigs)) %*% t(Phi); S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 1e-8))) %*% t(e$vectors)
}

#' Example data for the CAP methods
#'
#' Self-contained synthetic-data generators, one per method, returning data in
#' exactly the shape the matching wrapper expects, plus a `truth` element holding
#' the data-generating parameters so estimates can be checked. Pair each with its
#' [cap_fit] wrapper; `vignette("CAPmethods")` shows a worked run of each.
#'
#' @param n number of subjects.
#' @param m number of clusters (MCAP).
#' @param ni number of units per cluster (MCAP).
#' @param p response/mediator dimension; for CoC, the predictor (`X`) covariance
#'   dimension.
#' @param q number of covariates (HCAP); for CoC, the outcome (`Y`) covariance
#'   dimension.
#' @param Ti within-unit sample size (rows per response matrix).
#' @param nV number of visits per subject (LCAP).
#' @param Tx,Ty per-subject sample sizes for the predictor (`X`) and outcome
#'   (`Y`) covariance matrices (CoC).
#' @param kappa von Mises-Fisher concentration of the cluster loadings (MCAP).
#' @param seed RNG seed.
#' @return A named list of inputs (named as the matching wrapper expects) plus a
#'   `truth` list of the data-generating parameters.
#' @examples
#' d <- hdcap_example()
#' str(d, max.level = 1)
#' @name cap_examples
NULL

#' @rdname cap_examples
#' @export
hdcap_example <- function(n = 80L, p = 6L, Ti = 60L, seed = 1L) {
  set.seed(seed)
  X <- cbind(`(Intercept)` = 1, cov = rnorm(n))         # n x 2 covariate matrix
  Phi <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  if (Phi[which.max(abs(Phi[, 1])), 1] < 0) Phi[, 1] <- -Phi[, 1]
  beta <- c(0, 0.9)                                       # log-variance ~ X %*% beta
  Y <- lapply(seq_len(n), function(i) {
    eigs <- rep(0.4, p); eigs[1] <- exp(as.numeric(X[i, ] %*% beta))
    matrix(rnorm(Ti * p), Ti, p) %*% .cap_rootSig(Phi, eigs)
  })
  list(X = X, Y_list = Y,
       truth = list(gamma = Phi[, 1],
                    beta = matrix(beta, ncol = 1,
                                  dimnames = list(c("(Intercept)", "cov"), "D1"))))
}

#' @rdname cap_examples
#' @export
lcap_example <- function(n = 60L, p = 6L, nV = 6L, Ti = 40L, seed = 2024L) {
  set.seed(seed)
  tg <- seq(-0.5, 0.5, length.out = nV)
  Phi <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  if (Phi[which.max(abs(Phi[, 1])), 1] < 0) Phi[, 1] <- -Phi[, 1]
  be <- c(Intercept = 0, time = -0.5, dose = 0.6)
  Y <- X <- vector("list", n)
  for (i in seq_len(n)) {
    b0 <- be[["Intercept"]] + rnorm(1, 0, 0.3)
    bt <- be[["time"]]      + rnorm(1, 0, 0.2)
    bd <- be[["dose"]]      + rnorm(1, 0, 0.2)
    dose <- round(rnorm(nV), 2)
    Xi <- cbind(Intercept = 1, time = tg, dose = dose)
    Yi <- vector("list", nV)
    for (v in seq_len(nV)) {
      eigs <- rep(0.5, p); eigs[1] <- exp(b0 + bt * tg[v] + bd * dose[v])
      Yi[[v]] <- matrix(rnorm(Ti * p), Ti, p) %*% .cap_rootSig(Phi, eigs)
    }
    Y[[i]] <- Yi; X[[i]] <- Xi
  }
  names(Y) <- names(X) <- paste0("S", seq_len(n))
  list(Y = Y, X = X, cov_names = c("time", "dose"),
       truth = list(gamma = Phi[, 1],
                    beta = matrix(be, ncol = 1, dimnames = list(names(be), "D1"))))
}

#' @rdname cap_examples
#' @export
mcap_example <- function(m = 12L, ni = 15L, p = 5L, Ti = 80L, kappa = 150, seed = 2024L) {
  set.seed(seed)
  rvmf <- cap_internal("mcap", "rvmf")              # bundled vMF sampler
  gpop <- c(1, rep(0, p - 1))
  G  <- rvmf(m, gpop, kappa)
  Bg <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  Y <- X2 <- vector("list", m); cl_cos <- numeric(m)
  for (i in seq_len(m)) {
    gi <- G[i, ] / sqrt(sum(G[i, ]^2)); cl_cos[i] <- abs(sum(gi * gpop))
    B <- qr.Q(qr(cbind(gi, Bg)))[, seq_len(p)]; if (sum(B[, 1] * gi) < 0) B[, 1] <- -B[, 1]
    eps <- rnorm(1, 0, 0.2); th <- rnorm(1, 0, 0.2); dose <- round(rnorm(ni), 2)
    X2[[i]] <- matrix(dose, ni, 1, dimnames = list(NULL, "dose"))
    Yi <- vector("list", ni)
    for (j in seq_len(ni)) {
      eigs <- rep(0.2, p); eigs[1] <- exp(dose[j] * (1 + th) + eps)
      Yi[[j]] <- matrix(rnorm(Ti * p), Ti, p) %*% .cap_rootSig(B, eigs)
    }
    Y[[i]] <- Yi
  }
  names(Y) <- names(X2) <- paste0("Cl", seq_len(m))
  list(Y = Y, X2 = X2, cov_names = "dose",
       truth = list(gamma = gpop, kappa = kappa, cl_cos = cl_cos,
                    beta = matrix(c(0, 1), ncol = 1,
                                  dimnames = list(c("Intercept", "dose"), "D1"))))
}

#' @rdname cap_examples
#' @export
coc_example <- function(n = 150L, p = 10L, q = 5L, Tx = 150L, Ty = 150L, seed = 2024L) {
  # Case-1 simulation of Zhao et al. (Biometrics 2025): p-dim predictor X and
  # q-dim outcome Y covariances with TWO covariance-on-covariance pairs --
  # predictor directions {1, 3} drive outcome directions {2, 4} with alpha = (3, 2)
  # and covariate effects beta; the remaining directions carry covariate-free noise.
  stopifnot(p >= 3L, q >= 4L)
  set.seed(100); Gamma1 <- qr.Q(qr(matrix(runif(p * p), p, p)))      # predictor basis
  for (j in seq_len(p)) if (Gamma1[which.max(abs(Gamma1[, j])), j] < 0) Gamma1[, j] <- -Gamma1[, j]
  set.seed(500); Gamma2 <- qr.Q(qr(matrix(runif(q * q), q, q)))      # outcome basis
  for (j in seq_len(q)) if (Gamma2[which.max(abs(Gamma2[, j])), j] < 0) Gamma2[, j] <- -Gamma2[, j]
  x.eigen.m <- exp(seq(1, -2, length.out = p)); x.eigen.sd <- 0.5
  y.eigen.m <- exp(seq(1, -2, length.out = q)); y.eigen.sd <- 0.5
  x.idx <- c(1, 3); y.idx <- c(2, 4)
  alpha <- c(3, 2); beta <- cbind(c(1, -1), c(-1, 1))

  set.seed(seed)
  W <- cbind(Intercept = 1, group = rbinom(n, 1, 0.5))
  L1 <- matrix(NA_real_, n, p)
  for (j in seq_len(p)) L1[, j] <- exp(rnorm(n, log(x.eigen.m[j]), x.eigen.sd))
  L2 <- matrix(NA_real_, n, q)
  for (k in seq_len(q)) {
    f <- which(y.idx == k)
    L2[, k] <- if (length(f)) exp(alpha[f] * log(L1[, x.idx[f]]) + W %*% beta[, f])
               else exp(rnorm(n, log(y.eigen.m[k]), y.eigen.sd))
  }
  X <- Y <- vector("list", n)
  for (i in seq_len(n)) {
    X[[i]] <- MASS::mvrnorm(Tx, rep(0, p), Gamma1 %*% diag(L1[i, ]) %*% t(Gamma1))
    Y[[i]] <- MASS::mvrnorm(Ty, rep(0, q), Gamma2 %*% diag(L2[i, ]) %*% t(Gamma2))
  }
  names(X) <- names(Y) <- paste0("S", seq_len(n))
  Gt <- Gamma2[, y.idx, drop = FALSE]; colnames(Gt) <- c("P1", "P2")   # true outcome loadings
  Tt <- Gamma1[, x.idx, drop = FALSE]; colnames(Tt) <- c("P1", "P2")   # true predictor loadings
  list(Y = Y, X = X, W = W,
       truth = list(gamma = Gt, theta = Tt, alpha = alpha,
                    beta = matrix(beta, nrow = 2,
                                  dimnames = list(c("Intercept", "group"), c("P1", "P2")))))
}

#' @rdname cap_examples
#' @export
capmediation_example <- function(n = 60L, p = 10L, Ti = 40L, seed = 1L) {
  set.seed(seed)
  x <- rep(c(0, 1), length.out = n); w <- scale(rnorm(n))[, 1]
  X <- cbind(X = x, W1 = w)
  Phi <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  M <- vector("list", n); Y <- numeric(n)
  for (i in seq_len(n)) {
    lv <- 0 + 0.8 * x[i] + 0.4 * w[i] + rnorm(1, 0, 0.3)   # log M-variance on theta
    eigs <- rep(0.4, p); eigs[1] <- exp(lv)
    M[[i]] <- matrix(rnorm(Ti * p), Ti, p) %*% .cap_rootSig(Phi, eigs)
    Y[i] <- 0 + 0.5 * x[i] + 0.3 * w[i] + 0.6 * lv + rnorm(1, 0, 0.3)
  }
  list(X = X, M = M, Y = Y, truth = list(theta = Phi[, 1], alpha_x = 0.8, beta = 0.6))
}

#' @rdname cap_examples
#' @export
hcap_example <- function(n = 60L, p = 5L, q = 40L, Ti = 50L, seed = 2023L) {
  stopifnot(p == 5L, q >= 35L)
  set.seed(seed)
  s2 <- c(10, 20, 30); s3 <- c(15, 25, 35)
  b0 <- matrix(0, q, p); b0[s2, 2] <- c(2, 2, -2); b0[s3, 3] <- c(1, -1, 1)
  b0[1, ] <- sample(c(-10:-1, 1:10), p, replace = TRUE)
  phi <- matrix(c( 0.447,  0.447,  0.447,  0.447,  0.447,
                   0.447, -0.862,  0.138,  0.138,  0.138,
                   0.447,  0.138, -0.862,  0.138,  0.138,
                   0.447,  0.138,  0.138, -0.862,  0.138,
                   0.447,  0.138,  0.138,  0.138, -0.862), nrow = p)
  xmat <- cbind(1, matrix(rnorm(n * (q - 1)), n))
  lp <- matrix(NA_real_, n, p)
  idx1 <- colSums(b0[-1, , drop = FALSE]) != 0; idx0 <- !idx1
  lp[, idx1] <- exp(xmat %*% b0[, idx1, drop = FALSE])
  lp[, idx0] <- exp(sapply(b0[1, idx0], function(x) rnorm(n, x, 0.5)))
  Y <- lapply(seq_len(n), function(i) {
    Sig <- phi %*% diag(lp[i, ]) %*% t(phi); Sig <- (Sig + t(Sig)) / 2
    MASS::mvrnorm(Ti, rep(0, p), Sig)
  })
  list(X = xmat, Y_list = Y, truth = list(gamma = phi[, 2:3], beta = b0))
}

#' @rdname cap_examples
#' @export
cappcl_example <- function(n = 80L, p = 6L, Ti = 100L, seed = 1L) {
  set.seed(seed)
  W <- cbind(1, age = scale(rnorm(n))[, 1])              # membership covariates
  cl <- 1 + rbinom(n, 1, plogis(W %*% c(0, 2)))          # membership driven by W
  X <- cbind(1, group = rep(c(0, 1), length.out = n))    # variance covariate
  beta_true <- cbind(C1 = c(0.0, 0.8), C2 = c(0.5, -0.6))# cluster-specific variance models
  Phi <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  if (Phi[which.max(abs(Phi[, 1])), 1] < 0) Phi[, 1] <- -Phi[, 1]
  base <- seq(0.8, 0.1, length.out = p)
  Y <- lapply(seq_len(n), function(i) {
    eigs <- base; eigs[1] <- exp(sum(X[i, ] * beta_true[, cl[i]]))
    matrix(rnorm(Ti * p), Ti, p) %*% .cap_rootSig(Phi, eigs)
  })
  list(Y = Y, X = X, W = W,
       truth = list(gamma = Phi[, 1], cluster = cl, beta = beta_true))
}
