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

# Multivariate-normal sampler, bit-identical to mvtnorm::rmvnorm(method = "eigen")
# in modern mvtnorm (which fills the standard-normal matrix BY ROW). Lets the MCAP
# example reproduce the Shiny app's data exactly without an mvtnorm dependency.
.cap_rmvnorm <- function(n, mean, sigma) {
  ev <- eigen(sigma, symmetric = TRUE)
  R  <- ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0)))
  X  <- matrix(stats::rnorm(n * ncol(sigma)), nrow = n, byrow = TRUE) %*% R
  sweep(X, 2L, mean, "+")
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
hdcap_example <- function(n = 100L, Ti = 100L, seed = 100L) {
  # The manuscript HD-shrinkage simulation (210309/eg), identical to the Shiny app:
  # p = 20, common eigenbasis Gamma, one binary covariate (group). TWO directions
  # satisfy the CAP model -- basis cols 2 and 4, with group effects -1 and +1 --
  # while every other direction has a random covariate-free log-variance. With
  # covariance shrinkage on, nD = 2 recovers BOTH covariate-driven directions.
  p <- 20L; beta.sd <- 0.5
  # common orthonormal eigenbasis (runif(p) recycled into p x p, as in the paper)
  set.seed(seed)
  Gamma <- qr.Q(qr(matrix(stats::runif(p), p, p)))
  for (j in 1:p) if (Gamma[which.max(abs(Gamma[, j])), j] < 0) Gamma[, j] <- -Gamma[, j]
  beta1.vec <- c(seq(5, 1, length.out = 10), seq(0.5, -1, length.out = p - 10))
  beta2.vec <- c(c(0, -1, 0, 1, 0), rep(0, p - 5))
  beta.mat  <- rbind(beta1.vec, beta2.vec)
  cov_dirs  <- which(beta2.vec != 0)                    # = c(2, 4)

  set.seed(seed)
  X <- cbind(Intercept = 1, group = stats::rbinom(n, size = 1, prob = 0.5))
  Sigma <- array(NA, c(p, p, n))
  for (i in 1:n) {
    Delta <- matrix(0, p, p)
    for (j in 1:ncol(beta.mat)) {
      Delta[j, j] <- if (beta.mat[2, j] == 0)
        exp(stats::rnorm(1, mean = beta.mat[1, j], sd = beta.sd))
      else exp(as.numeric(X[i, ] %*% beta.mat[, j]))
    }
    Sigma[, , i] <- Gamma %*% Delta %*% t(Gamma)
  }
  set.seed(seed)
  Y <- lapply(1:n, function(i) .cap_rmvnorm(Ti, rep(0, p), Sigma[, , i]))

  # truth = the covariate-driven directions (the components satisfying the model)
  gamma_true <- Gamma[, cov_dirs, drop = FALSE]; colnames(gamma_true) <- paste0("D", seq_along(cov_dirs))
  beta_true  <- beta.mat[, cov_dirs, drop = FALSE]
  rownames(beta_true) <- c("Intercept", "group"); colnames(beta_true) <- paste0("D", seq_along(cov_dirs))
  list(X = X, Y_list = Y, truth = list(gamma = gamma_true, beta = beta_true))
}

#' @rdname cap_examples
#' @export
lcap_example <- function(n = 100L, nV = 5L, Ti = 100L, seed = 4L) {
  # The manuscript longitudinal simulation (p20_q3, case 1), identical to the Shiny
  # app: p = 20, common time-invariant eigenbasis Gamma, ~nV visits/subject with
  # ~Ti samples/visit, two within-subject covariates (x1, x2) and a subject random
  # intercept. TWO directions satisfy the CAP model -- basis cols 2 and 4 -- whose
  # log-variance depends on the covariates; nD = 2 with shrinkage recovers both.
  # (The estimator fits within-subject random slopes, so covariates are
  # time-varying; the manuscript's binary time-invariant x1 is rendered as one.)
  p <- 20L; nV.m <- nV; nT.m <- Ti; seedl <- seed
  set.seed(100)
  Gamma <- qr.Q(qr(matrix(stats::runif(p), p, p)))
  for (j in 1:p) if (Gamma[which.max(abs(Gamma[, j])), j] < 0) Gamma[, j] <- -Gamma[, j]
  beta.mat <- rbind(c(seq(3, -1, length.out = 5), seq(-1.5, -3, length.out = p - 5)),
                    c(c(0, -0.5, 0, 0.5, 0), rep(0, p - 5)),
                    c(c(0,  0.5, 0, -0.25, 0), rep(0, p - 5)))
  cov_dirs <- which(beta.mat[2, ] != 0)               # = c(2, 4)
  set.seed(100); nVvec <- round(stats::rnorm(n, nV.m, sd = 1))
  set.seed(100); nT <- matrix(NA_integer_, n, max(nVvec))
  for (i in 1:n) nT[i, 1:nVvec[i]] <- round(stats::rnorm(nVvec[i], mean = nT.m, sd = 5))
  beta0.mat <- matrix(NA, n, ncol(beta.mat))
  for (j in 1:ncol(beta.mat)) { set.seed(100); beta0.mat[, j] <- stats::rnorm(n, beta.mat[1, j], sd = 0.1) }
  set.seed(as.numeric("20201210") + seedl * 100)
  Y <- X <- vector("list", n)
  for (i in 1:n) {
    Xi <- cbind(Intercept = 1, x1 = stats::rnorm(nVvec[i], 0, 0.5), x2 = stats::rnorm(nVvec[i], 0, 0.5))
    Yi <- vector("list", nVvec[i])
    for (v in 1:nVvec[i]) {
      delta <- vapply(1:p, function(j) exp(sum(Xi[v, ] * c(beta0.mat[i, j], beta.mat[-1, j]))), numeric(1))
      Yi[[v]] <- .cap_rmvnorm(nT[i, v], rep(0, p), Gamma %*% diag(delta) %*% t(Gamma))
    }
    Y[[i]] <- Yi; X[[i]] <- Xi
  }
  names(Y) <- names(X) <- paste0("S", 1:n)
  gamma_true <- Gamma[, cov_dirs, drop = FALSE]; colnames(gamma_true) <- paste0("D", seq_along(cov_dirs))
  beta_true  <- beta.mat[, cov_dirs, drop = FALSE]
  rownames(beta_true) <- c("Intercept", "x1", "x2"); colnames(beta_true) <- paste0("D", seq_along(cov_dirs))
  list(Y = Y, X = X, cov_names = c("x1", "x2"),
       truth = list(gamma = gamma_true, beta = beta_true))
}

#' @rdname cap_examples
#' @export
mcap_example <- function(m = 20L, ni = 50L, Ti = 80L, kappa = 10, seed = 3L) {
  # The manuscript's gamma-varying simulation (p5_q4_2-1, case 1), identical to the
  # Shiny app demo: p = 5, q1 = 2 fixed covariates (X1), q2 = 1 random-slope
  # covariate (X2), TWO covariate-driven directions (population basis cols 2, 4),
  # cluster directions gamma_i ~ vMF(gamma, kappa). ni / Ti are the Poisson means
  # of the per-cluster unit count and per-unit sample size; `seed` is the
  # data-realisation seed (the other, structural seeds are fixed as in the paper).
  p <- 5L; q1 <- 2L; q2 <- 1L
  rv.idx <- c(2L, 4L); nvec.lambda <- ni; Tmat.lambda <- Ti; seedl <- seed
  rvmf <- cap_internal("mcap", "rvmf")              # bundled vMF sampler

  # population orthonormal basis (NB: runif(p) recycled into a p x p matrix, as in
  # the manuscript -- not runif(p*p)) + sign convention
  set.seed(100)
  Gamma <- qr.Q(qr(matrix(stats::runif(p), p, p)))
  for (j in 1:p) if (Gamma[which.max(abs(Gamma[, j])), j] < 0) Gamma[, j] <- -Gamma[, j]

  # rows: beta0 (intercept), beta1/beta2 (fixed X1), beta3 (random-slope X2 mean)
  beta.mat <- rbind(seq(5, -1, length.out = 5),
                    c(0, 1, 0, -1, 0),
                    c(0, -0.5, 0, 0.5, 0),
                    c(0, -0.5, 0, 0.5, 0))
  sigma <- 0.1; Omega <- matrix(0.1^2, 1, 1); beta.nr.sd <- 0.1
  beta.fix.idx <- 2:(q1 + 1); beta.rnd.idx <- (q1 + 2):(q1 + q2 + 1)

  set.seed(500); nvec <- stats::rpois(m, nvec.lambda)
  Tmat <- matrix(NA_integer_, m, max(nvec))
  for (i in 1:m) Tmat[i, 1:nvec[i]] <- stats::rpois(nvec[i], Tmat.lambda)

  set.seed(100); beta0.mat <- matrix(NA, m, p)
  for (kk in 1:p)
    beta0.mat[, kk] <- stats::rnorm(m, beta.mat[1, kk], if (kk %in% rv.idx) sigma else beta.nr.sd)
  set.seed(100); beta2.mat <- array(NA, c(m, q2, p))
  for (kk in 1:p) beta2.mat[, , kk] <- .cap_rmvnorm(m, beta.mat[beta.rnd.idx, kk], Omega)

  set.seed(100); gmat <- array(NA, c(m, p, length(rv.idx)))
  gmat[, , 1] <- rvmf(m, Gamma[, rv.idx[1]], kappa)
  for (ss in 2:length(rv.idx)) for (i in 1:m) {
    repeat {
      otmp <- rvmf(10000, Gamma[, rv.idx[ss]], kappa)
      oo <- apply(abs(otmp %*% gmat[i, , 1:(ss - 1), drop = FALSE]), 1, sum)
      idx <- which.min(oo)
      if (oo[idx] <= 1e-4) { gmat[i, , ss] <- otmp[idx, ]; break }
    }
  }
  Pi <- array(NA, c(p, p, m))
  for (i in 1:m) {
    Pt <- matrix(NA, p, p); Pt[, rv.idx] <- gmat[i, , ]; Pt[, -rv.idx] <- MASS::Null(gmat[i, , ])
    Pi[, , i] <- Pt
  }

  set.seed(as.numeric("20221206") + seedl * 100)
  X1.fix <- stats::rbinom(m, 1, 0.5)
  ids <- paste0("G", 1:m); x1_names <- paste0("X1_", 1:q1); x2_names <- paste0("X2_", 1:q2)
  Y <- X1 <- X2 <- vector("list", m); names(Y) <- names(X1) <- names(X2) <- ids
  cl_cos <- matrix(NA, m, length(rv.idx))
  for (i in 1:m) {
    cl_cos[i, ] <- vapply(seq_along(rv.idx),
      function(s) abs(sum(gmat[i, , s] * Gamma[, rv.idx[s]])), numeric(1))
    set.seed(1000 + i * 10)
    X1[[i]] <- cbind(rep(X1.fix[i], nvec[i]), stats::rnorm(nvec[i], 0, 0.5)); colnames(X1[[i]]) <- x1_names
    X2[[i]] <- cbind(stats::rnorm(nvec[i], 0, 0.5)); colnames(X2[[i]]) <- x2_names
    Yi <- vector("list", nvec[i])
    for (j in 1:nvec[i]) {
      delta <- numeric(p)
      for (kk in 1:p) {
        bt <- c(beta0.mat[i, kk], beta.mat[beta.fix.idx, kk], beta2.mat[i, , kk])
        delta[kk] <- exp(c(1, X1[[i]][j, ], X2[[i]][j, ]) %*% bt)
      }
      Sig <- Pi[, , i] %*% diag(delta) %*% t(Pi[, , i])
      set.seed(10000 + i * 100 + j * 10)
      Yi[[j]] <- .cap_rmvnorm(Tmat[i, j], rep(0, p), Sig)
    }
    Y[[i]] <- Yi
  }
  Gtrue <- Gamma[, rv.idx]; colnames(Gtrue) <- paste0("D", seq_along(rv.idx))
  Btrue <- beta.mat[, rv.idx]
  rownames(Btrue) <- c("Intercept", x1_names, x2_names); colnames(Btrue) <- paste0("D", seq_along(rv.idx))
  list(Y = Y, X1 = X1, X2 = X2, cov_names = c(x1_names, x2_names),
       truth = list(gamma = Gtrue, beta = Btrue, kappa = kappa, cl_cos = cl_cos))
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
capmediation_example <- function(n = 100L, p = 10L, Ti = 150L, seed = 2024L) {
  # Covariance-mediator example, identical to the Shiny app: a binary treatment
  # shifts the variance of a p-dim mediator M along one latent direction theta;
  # the outcome Y depends on that log-variance (mediator->outcome effect beta)
  # plus a direct treatment effect. Exposure->mediator-variance is the a-path
  # (alpha); the true indirect effect is IE = alpha * beta = 0.8 * 0.7 = 0.56.
  set.seed(seed)
  alpha0 <- 0.2; alpha_x <- 0.8; beta <- 0.7; gamma0 <- 0.1; gamma_x <- 0.4
  x <- rep(c(0, 1), length.out = n)
  X <- matrix(x, ncol = 1, dimnames = list(NULL, "treatment"))
  Phi <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
  if (Phi[which.max(abs(Phi[, 1])), 1] < 0) Phi[, 1] <- -Phi[, 1]
  base <- c(NA, exp(seq(log(1.2), log(0.3), length.out = p - 1)))   # background eigenvalues
  M <- vector("list", n); Y <- numeric(n)
  for (i in seq_len(n)) {
    lv <- alpha0 + alpha_x * x[i] + stats::rnorm(1, 0, 0.3)         # log mediator-variance along theta
    ev <- base; ev[1] <- exp(lv)
    M[[i]] <- matrix(stats::rnorm(Ti * p), Ti, p) %*% .cap_rootSig(Phi, ev)
    Y[i] <- gamma0 + gamma_x * x[i] + beta * lv + stats::rnorm(1, 0, 0.3)
  }
  names(M) <- paste0("S", seq_len(n))
  list(X = X, M = M, Y = Y,
       truth = list(theta = Phi[, 1], alpha = alpha_x, beta = beta,
                    gamma = gamma_x, IE = alpha_x * beta))
}

#' @rdname cap_examples
#' @export
hcap_example <- function(n = 100L, q = 200L, Ti = 100L, seed = 2023L) {
  # High-dimensional-covariate CAP example, identical to the Shiny app: p = 5
  # responses, q covariates, TWO covariate-driven directions (mediator basis cols
  # 2 and 3), each driven by a small SPARSE subset of the q covariates; the rest
  # carry a covariate-free baseline. The sparse estimator + multi-split inference
  # recover the directions and their true signal covariates.
  p <- 5L
  set.seed(seed)
  s2 <- c(10, 20, 30); s3 <- c(15, 25, 35)         # true signal covariates
  b0 <- matrix(0, q, p)
  b0[s2, 2] <- c(2, 2, -2); b0[s3, 3] <- c(1, -1, 1)
  b0[1, ] <- sample(c(-10:-1, 1:10), p, replace = TRUE)   # per-response intercept
  phi <- matrix(c(0.447, 0.447, 0.447, 0.447, 0.447,
                  0.447, -0.862, 0.138, 0.138, 0.138,
                  0.447, 0.138, -0.862, 0.138, 0.138,
                  0.447, 0.138, 0.138, -0.862, 0.138,
                  0.447, 0.138, 0.138, 0.138, -0.862), nrow = p)
  X <- cbind(1, matrix(stats::rnorm(n * (q - 1)), n))
  lp <- matrix(NA_real_, n, p); i1 <- colSums(b0[-1, ]) != 0
  lp[, i1]  <- exp(X %*% b0[, i1])
  lp[, !i1] <- exp(sapply(b0[1, !i1], function(z) stats::rnorm(n, z, 0.5)))
  Y <- lapply(seq_len(n), function(i)
    .cap_rmvnorm(Ti, rep(0, p), phi %*% diag(lp[i, ]) %*% t(phi)))
  names(Y) <- paste0("S", seq_len(n))
  cov_names <- c("Intercept", paste0("X", 1:(q - 1))); colnames(X) <- cov_names
  Gtrue <- phi[, c(2, 3)]; colnames(Gtrue) <- c("T2", "T3")
  Btrue <- b0[, c(2, 3)]; colnames(Btrue) <- c("T2", "T3"); rownames(Btrue) <- cov_names
  list(X = X, Y_list = Y, truth = list(gamma = Gtrue, beta = Btrue,
                                       signal = list(T2 = s2, T3 = s3)))
}

#' @rdname cap_examples
#' @export
cappcl_example <- function(n = 100L, p = 50L, Ti = 100L, seed = 1L) {
  # Covariance-clustering example, identical to the Shiny app: p = 50, K = 2
  # clusters, with TWO covariate-driven components -- D2 (well separated) and D4
  # (clusters differ only in covariate signs). Within a cluster the projected
  # log-variance is a log-linear model in X; cluster membership depends on W.
  set.seed(seed)
  Pi <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
  x1 <- stats::rbinom(n, 1, 0.5); x2 <- stats::rnorm(n)
  X <- cbind(Intercept = 1, x1 = x1, x2 = x2)               # q1 = 3
  w1 <- stats::rbinom(n, 1, 0.5)
  W <- cbind(Intercept = 1, w1 = w1)                        # q2 = 2
  cD2 <- 1 + stats::rbinom(n, 1, stats::plogis(W %*% c(0.5, -1)))
  cD4 <- 1 + stats::rbinom(n, 1, stats::plogis(W %*% c(-0.25, 0.5)))
  bD2 <- cbind(c(1, 1, -1), c(-1, -1, 1))                   # q1 x K (D2: well separated)
  bD4 <- cbind(c(0.5, 0.5, -0.5), c(0.5, -0.5, 0.5))        # D4: same baseline
  logmean <- seq(3, -1, length.out = p)
  Y <- vector("list", n)
  for (i in seq_len(n)) {
    lam <- exp(stats::rnorm(p, logmean, 0.2))
    lam[2] <- exp(sum(X[i, ] * bD2[, cD2[i]]))
    lam[4] <- exp(sum(X[i, ] * bD4[, cD4[i]]))
    Y[[i]] <- matrix(stats::rnorm(Ti * p), Ti, p) %*% .cap_rootSig(Pi, lam)
  }
  names(Y) <- paste0("S", seq_len(n))
  list(Y = Y, X = X, W = W,
       truth = list(gamma = cbind(D2 = Pi[, 2], D4 = Pi[, 4]),
                    cluster = cbind(D2 = cD2, D4 = cD4),
                    beta = list(D2 = bD2, D4 = bD4)))
}
