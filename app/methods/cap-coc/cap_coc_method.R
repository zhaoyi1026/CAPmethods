# =============================================================================
# CAP-CoC -- Covariance-on-Covariance regression  (plugin)
# -----------------------------------------------------------------------------
# Wraps the RcppArmadillo-accelerated implementation in ../CAP-CoC/V3/. Regresses
# the projected log-variance of an OUTCOME covariance (gamma' Sigma_i gamma) on
# the projected log-variance of a PREDICTOR covariance (theta' Delta_i theta),
# adjusting for covariates W. Both gamma and theta are learned from data.
#   log(gamma' Sigma_i gamma) = alpha * log(theta' Delta_i theta) + beta' W_i
# =============================================================================

# ---- 1. Load the method code from the CAPmethods package --------------------
# Uses the package's "coc" environment (COCReg(), COCReg.coef.boot(), ...); the
# coc_-prefixed kernel ships in the package DLL. See R/pkg_methods.R.
.coc_env <- cap_pkg_env("coc")

# ---- 2. Built-in example data (the manuscript's p=10, q=5 simulation) --------
# Faithful to the paper's case-1 design: two covariance-on-covariance pairs --
# predictor directions Gamma1[, c(1,3)] drive outcome directions Gamma2[, c(2,4)]
# with alpha = (3, 2) and covariate effects beta. All eigenvalues are random
# (heterogeneous) with decaying means; only the signal directions follow the CoC
# model. The app's default (nD = 1) recovers the leading pair cleanly.
coc_example <- function() {
  p <- 10L; q <- 5L; n <- 150L; nTx <- 150L; nTy <- 150L

  set.seed(100); G0 <- matrix(runif(p * p), p, p); Gamma1 <- qr.Q(qr(G0))
  for (j in 1:p) if (Gamma1[which.max(abs(Gamma1[, j])), j] < 0) Gamma1[, j] <- -Gamma1[, j]
  x.eigen.m <- exp(seq(1, -2, length.out = p)); x.eigen.sd <- 0.5
  set.seed(500); G0 <- matrix(runif(q * q), q, q); Gamma2 <- qr.Q(qr(G0))
  for (j in 1:q) if (Gamma2[which.max(abs(Gamma2[, j])), j] < 0) Gamma2[, j] <- -Gamma2[, j]
  y.eigen.m <- exp(seq(1, -2, length.out = q)); y.eigen.sd <- 0.5
  x.idx <- c(1, 3); y.idx <- c(2, 4)
  alpha <- c(3, 2); beta <- cbind(c(1, -1), c(-1, 1))

  set.seed(2024)
  W <- cbind(Intercept = 1, group = rbinom(n, 1, 0.5))
  L1 <- matrix(NA, n, p); for (j in 1:p) L1[, j] <- exp(rnorm(n, log(x.eigen.m[j]), x.eigen.sd))
  L2 <- matrix(NA, n, q)
  for (k in 1:q) {
    f <- which(y.idx == k)
    if (length(f)) L2[, k] <- exp(alpha[f] * log(L1[, x.idx[f]]) + W %*% beta[, f])
    else L2[, k] <- exp(rnorm(n, log(y.eigen.m[k]), y.eigen.sd))
  }
  X <- vector("list", n); Y <- vector("list", n)
  for (i in 1:n) {
    X[[i]] <- MASS::mvrnorm(nTx, rep(0, p), Gamma1 %*% diag(L1[i, ]) %*% t(Gamma1))
    Y[[i]] <- MASS::mvrnorm(nTy, rep(0, q), Gamma2 %*% diag(L2[i, ]) %*% t(Gamma2))
  }
  ids <- paste0("S", 1:n)
  names(Y) <- names(X) <- ids
  Gt <- Gamma2[, y.idx, drop = FALSE]; colnames(Gt) <- paste0("P", 1:2)
  Tt <- Gamma1[, x.idx, drop = FALSE]; colnames(Tt) <- paste0("P", 1:2)
  d <- list(Y = Y, X = X, W = W, ids = ids,
            y_names = paste0("Yv", 1:q), x_names = paste0("Xv", 1:p),
            w_names = colnames(W),
            Wdf = data.frame(id = ids, group = W[, "group"]),
            truth = list(gamma = Gt, theta = Tt, alpha = alpha,
                         beta = matrix(beta, nrow = 2, dimnames = list(c("Intercept", "group"), c("P1", "P2")))))
  d$preview_ui <- coc_preview(d)
  d
}

# ---- 3. Parser for uploaded data --------------------------------------------
# Y (outcome) and X (predictor): native list of matrices or long CSV. W: n x r
# covariate matrix / CSV (intercept auto-added, optional).
coc_parse <- function(files, opts) {
  yp <- coerce_Y_input(files$Y)
  xp <- coerce_Y_input(files$X)
  if (length(xp$Y) != length(yp$Y))
    stop("Outcome (Y) and predictor (X) must have the same number of subjects.")
  W <- coerce_X_input(files$W, ids = yp$ids, add_intercept = isTRUE(opts$add_intercept))
  d <- list(Y = yp$Y, X = xp$Y, W = W, ids = yp$ids,
            y_names = yp$var_names, x_names = xp$var_names, w_names = colnames(W),
            Wdf = data.frame(id = yp$ids, W[, setdiff(colnames(W), "Intercept"), drop = FALSE],
                             check.names = FALSE),
            truth = NULL)
  d$preview_ui <- coc_preview(d)
  d
}

# ---- 4. Preview / helpers ---------------------------------------------------
coc_preview <- function(d) {
  Txv <- vapply(d$X, nrow, integer(1)); Tyv <- vapply(d$Y, nrow, integer(1))
  p <- ncol(d$X[[1]]); q <- ncol(d$Y[[1]])
  w_cov <- setdiff(d$w_names, "Intercept")
  rng <- function(x) if (length(unique(x)) == 1) as.character(x[1]) else paste0(min(x), "–", max(x))
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", length(d$Y), theme = "primary"),
      bslib::value_box("Predictor dim (p)", p, theme = "secondary"),
      bslib::value_box("Outcome dim (q)", q, theme = "secondary"),
      bslib::value_box("Samples (X / Y)", paste0(rng(Txv), " / ", rng(Tyv)), theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           sprintf("Covariates W: %s.%s", paste(w_cov, collapse = ", "),
                   if (!is.null(d$truth)) "  Simulated data with known truth (two covariance-on-covariance pairs)." else "")))
}

coc_describe <- function(d)
  sprintf("%d subjects; predictor dim p = %d, outcome dim q = %d; covariates: %s.",
          length(d$Y), ncol(d$X[[1]]), ncol(d$Y[[1]]),
          paste(setdiff(d$w_names, "Intercept"), collapse = ", "))

coc_dfd_vec <- function(dfd, nD) {
  v <- if (is.list(dfd)) as.numeric(dfd$avg.level) else as.numeric(dfd)
  if (length(v) < nD) v <- c(v, rep(NA, nD - length(v)))
  v[seq_len(nD)]
}

# ---- 5. Run -----------------------------------------------------------------
coc_run <- function(d, params) {
  Y <- d$Y; X <- d$X; W <- d$W
  p <- ncol(X[[1]]); q <- ncol(Y[[1]])
  stop_crt <- params$stop.crt %||% "nD"
  nD <- if (stop_crt == "nD") (params$nD %||% 1) else NULL
  ninitial <- params$ninitial %||% 5
  fit <- .coc_env$COCReg(Y, X, W, stop.crt = stop_crt, nD = nD,
            DfD.thred = params$DfD.thred %||% 2,
            Hy = diag(q), Hx = diag(p),               # identity H (the recovering setting)
            cov.shrinkage.y = TRUE, cov.shrinkage.x = TRUE,
            burn.in = 200, max.itr = 1000, tol = 1e-4,
            ninitial = ninitial, seed = 100, score.return = TRUE, verbose = FALSE)
  # bootstrap inference for alpha / beta, per direction
  sims <- params$sims %||% 0
  inference <- NULL
  if (sims > 0) {
    inference <- lapply(seq_len(ncol(fit$gamma)), function(j)
      tryCatch(.coc_env$COCReg.coef.boot(Y, X, W, gamma = fit$gamma[, j], theta = fit$theta[, j],
                 cov.shrinkage.y = TRUE, cov.shrinkage.x = TRUE, boot = TRUE, sims = sims,
                 boot.ci.type = "se", conf.level = 0.95, verbose = FALSE),
               error = function(e) NULL))
  }
  list(fit = fit, inference = inference, W = W, ids = d$ids,
       y_names = d$y_names, x_names = d$x_names, w_names = d$w_names,
       truth = d$truth, params = params)
}

# match estimated directions to true pairs by combined gamma+theta cosine
coc_match <- function(Ghat, That, Gtrue, Ttrue) {
  J <- ncol(Ghat); K <- ncol(Gtrue); used <- integer(0); out <- list()
  for (j in seq_len(J)) {
    sc <- vapply(seq_len(K), function(k) if (k %in% used) NA_real_ else
      abs(sum(Ghat[, j] * Gtrue[, k])) + abs(sum(That[, j] * Ttrue[, k])), numeric(1))
    if (all(is.na(sc))) break
    k <- which.max(sc); used <- c(used, k)
    out[[length(out) + 1]] <- data.frame(est = j, k = k,
      gcos = abs(sum(Ghat[, j] * Gtrue[, k])), tcos = abs(sum(That[, j] * Ttrue[, k])))
  }
  do.call(rbind, out)
}

# ---- 6. Tables --------------------------------------------------------------
coc_summarize <- function(res) {
  fit <- res$fit; nD <- ncol(fit$gamma)
  dirs <- paste0("D", seq_len(nD)); out <- list()

  # alpha + beta coefficients (bootstrap inference if available, else point est.)
  coef_rows <- do.call(rbind, lapply(seq_len(nD), function(j) {
    inf <- if (!is.null(res$inference)) res$inference[[j]] else NULL
    if (!is.null(inf)) {
      a <- data.frame(Direction = dirs[j], Term = "alpha (log θ'Δθ)", inf$alpha,
                      check.names = FALSE, row.names = NULL)
      b <- data.frame(Direction = dirs[j], Term = paste0("beta: ", rownames(inf$beta)), inf$beta,
                      check.names = FALSE, row.names = NULL)
      rbind(a, b)
    } else {
      terms <- c("alpha (log θ'Δθ)", paste0("beta: ", res$w_names))
      data.frame(Direction = dirs[j], Term = terms,
                 Estimate = c(fit$alpha[j], fit$beta[, j]),
                 SE = NA, statistic = NA, pvalue = NA, LB = NA, UB = NA, check.names = FALSE)
    }
  }))
  out[["Coefficients (α, β)"]] <- coef_rows

  # gamma (outcome) and theta (predictor) loadings
  g <- as.data.frame(fit$gamma); names(g) <- dirs
  out[["γ outcome loadings"]] <- cbind(Variable = res$y_names, g)
  th <- as.data.frame(fit$theta); names(th) <- dirs
  out[["θ predictor loadings"]] <- cbind(Variable = res$x_names, th)

  # DfD (outcome / predictor) across dimensions
  out[["DfD across dimensions"]] <- data.frame(
    Dimension = dirs,
    `DfD outcome` = coc_dfd_vec(fit$DfD.y, nD),
    `DfD predictor` = coc_dfd_vec(fit$DfD.x, nD), check.names = FALSE)

  # truth comparison (example data)
  if (!is.null(res$truth)) {
    mt <- coc_match(fit$gamma, fit$theta, as.matrix(res$truth$gamma), as.matrix(res$truth$theta))
    rec <- do.call(rbind, lapply(seq_len(nrow(mt)), function(r) {
      k <- mt$k[r]; j <- mt$est[r]
      data.frame(`Est direction` = dirs[j], `Matched true pair` = colnames(res$truth$gamma)[k],
                 `γ cosine` = round(mt$gcos[r], 3), `θ cosine` = round(mt$tcos[r], 3),
                 `α true` = res$truth$alpha[k], `α est` = round(fit$alpha[j], 3),
                 check.names = FALSE)
    }))
    out[["Direction recovery"]] <- rec

    # alpha/beta truth vs estimate for matched pairs
    cmp <- do.call(rbind, lapply(seq_len(nrow(mt)), function(r) {
      k <- mt$k[r]; j <- mt$est[r]
      data.frame(`Est direction` = dirs[j], `Matched pair` = colnames(res$truth$gamma)[k],
                 Term = c("alpha", paste0("beta: ", rownames(res$truth$beta))),
                 True = c(res$truth$alpha[k], res$truth$beta[, k]),
                 Estimate = round(c(fit$alpha[j], fit$beta[, j]), 3), check.names = FALSE)
    }))
    out[["α, β: truth vs estimate"]] <- cmp
  }
  out
}

# ---- 7. Plots ---------------------------------------------------------------
coc_plots <- function(res) {
  fit <- res$fit; nD <- ncol(fit$gamma); dirs <- paste0("D", seq_len(nD)); plots <- list()

  mk_load <- function(M, vars, ttl) {
    g <- as.data.frame(M); names(g) <- dirs
    long <- do.call(rbind, lapply(seq_len(nD), function(j)
      data.frame(Variable = vars, Direction = dirs[j], Loading = g[, j])))
    long$Variable <- factor(long$Variable, levels = vars)
    plotly::plot_ly(long, x = ~Variable, y = ~Loading, color = ~Direction, type = "bar") |>
      plotly::layout(title = ttl, barmode = "group")
  }
  plots[["γ outcome loadings"]] <- list(plot = mk_load(fit$gamma, res$y_names, "Outcome projection loadings (γ)"))
  plots[["θ predictor loadings"]] <- list(plot = mk_load(fit$theta, res$x_names, "Predictor projection loadings (θ)"))

  # CoC fit: log(score.y) vs fitted = alpha*log(score.x) + W beta, per direction
  if (!is.null(fit$score.y) && !is.null(fit$score.x)) {
    sy <- as.matrix(fit$score.y); sx <- as.matrix(fit$score.x)
    W <- res$W
    sv <- do.call(rbind, lapply(seq_len(nD), function(j) {
      fitted <- fit$alpha[j] * log(sx[, j]) + as.numeric(W %*% fit$beta[, j])
      data.frame(id = res$ids, Direction = dirs[j],
                 logvar = log(sy[, j]), fitted = fitted)
    }))
    rng <- range(c(sv$fitted, sv$logvar), finite = TRUE)
    p_sc <- plotly::plot_ly(sv, x = ~fitted, y = ~logvar, color = ~Direction,
                            type = "scatter", mode = "markers",
                            marker = list(size = 7, opacity = 0.6),
                            text = ~id, hovertemplate = "%{text}<extra></extra>") |>
      plotly::add_lines(x = rng, y = rng, line = list(dash = "dot", color = "grey"),
                        showlegend = FALSE, inherit = FALSE) |>
      plotly::layout(title = "Outcome vs fitted log projected-variance",
                     xaxis = list(title = "α·log(θ'Δθ) + β'W"),
                     yaxis = list(title = "log(γ'Σγ)"))
    plots[["CoC fit"]] <- list(plot = p_sc, data = sv)
  }
  plots
}

# ---- 8. Register ------------------------------------------------------------
register_method(list(
  id = "cap-coc",
  name = "CAP-CoC",
  full_name = "Covariance-on-Covariance Regression",
  short = "Regress an outcome covariance's projected log-variance on a predictor covariance's projected log-variance (both directions learned), adjusting for covariates.",
  status = "ready",
  tags = c("covariance regression", "covariance-on-covariance", "two-way", "bootstrap"),
  paper = list(
    citation = "Zhao, Y., & Zhao, Y. (2025). Covariance-on-covariance regression. Biometrics, 81(3), ujaf097.",
    url = "https://doi.org/10.1093/biomtc/ujaf097"),
  explain = file.path(APP_DIR, "methods", "cap-coc", "explain.md"),
  example_note = paste("The manuscript's simulation (p=10 predictor, q=5 outcome,",
                       "150 subjects × 150 samples): two covariance-on-covariance pairs",
                       "with α=(3,2). The default (nD=1) recovers the leading pair."),
  x_intercept_option = TRUE,
  data_inputs = list(
    list(id = "Y", label = "Outcome (list of T×q matrices)",
         help = "An .rds/.RData holding a list of length n; element i is a T_iy × q outcome matrix. A long CSV (subject-id column + q outcome columns) is also accepted."),
    list(id = "X", label = "Predictor (list of T×p matrices)",
         help = "An .rds/.RData list; element i is a T_ix × p predictor matrix. A long CSV (subject-id + p predictor columns) is also accepted."),
    list(id = "W", label = "Covariates (n × r)",
         help = "An n × r covariate matrix in subject order, or a CSV whose first column is the subject id. Intercept auto-added (optional).")
  ),
  params = list(
    list(id = "stop.crt", label = "Number of directions chosen by", type = "select",
         choices = c("Fixed number (nD)" = "nD", "DfD criterion" = "DfD"), default = "nD",
         help = "Fix the count, or keep adding directions while DfD stays below the threshold. The leading direction is the most reliably recovered."),
    list(id = "nD", label = "Number of directions (nD)", type = "integer",
         default = 1, min = 1, max = 5,
         help = "Used when selecting by a fixed number. Each direction adds a bootstrap pass."),
    list(id = "DfD.thred", label = "DfD threshold", type = "numeric",
         default = 2, min = 1, max = 100, step = 0.5, help = "Used when selecting by the DfD criterion."),
    list(id = "ninitial", label = "# random initializations", type = "integer",
         default = 5, min = 1, max = 20, help = "Random restarts for the joint γ/θ search."),
    list(id = "sims", label = "Bootstrap replicates (α, β inference)", type = "integer",
         default = 200, min = 0, max = 2000, help = "Subject-level bootstrap; 0 = skip (point estimates only).")
  ),
  example = coc_example,
  export_example = function(d) list(Y = d$Y, X = d$X, W = d$Wdf),
  parse = coc_parse,
  describe_data = coc_describe,
  run = coc_run,
  summarize = coc_summarize,
  plots = coc_plots
))
