# =============================================================================
# HDCAP -- High-Dimensional Covariate Assisted Principal regression  (plugin)
# -----------------------------------------------------------------------------
# Wraps the RcppArmadillo-accelerated implementation in ../HDCAP/V7/ and
# registers it with the website. Mirrors the CAP plugin structure.
# =============================================================================

# ---- 1. Load the method code from the CAPmethods package --------------------
# Uses the package's "hdcap" environment (capReg(), cap_beta_boot(), ...); the
# compiled kernel ships in the package DLL. See R/pkg_methods.R.
.hdcap_env <- cap_pkg_env("hdcap")

# ---- 2. Built-in example data = the manuscript HD-shrinkage simulation -------
# Faithful translation of the paper's example (210309/eg): p = 20, n = 100,
# Tᵢ = 100, a common eigenbasis Gamma, one binary covariate "group". The
# log-variance is covariate-driven on TWO directions (basis cols 2 and 4, with
# group effects β = -1 and +1) -- the "components that satisfy the CAP model" --
# while every other direction has a random (covariate-free) log-variance. With
# covariance shrinkage on, nD = 2 recovers BOTH covariate-driven directions
# (γ cosine ≈ 0.99, β̂(group) ≈ ∓1); higher nD / DfD additionally pick up the
# high-variance covariate-free directions (β̂(group) ≈ 0).
hdcap_example <- function() {
  p <- 20L; n <- 100L; nT <- 100L; beta.sd <- 0.5

  # common orthonormal eigenbasis (manuscript seed; runif(p) recycled into p x p)
  set.seed(100)
  Gamma <- qr.Q(qr(matrix(stats::runif(p), p, p)))
  for (j in 1:p) if (Gamma[which.max(abs(Gamma[, j])), j] < 0) Gamma[, j] <- -Gamma[, j]

  # per-direction coefficients: row 1 = baseline (intercept) log-variance,
  # row 2 = group effect (nonzero only on directions 2 and 4).
  beta1.vec <- c(seq(5, 1, length.out = 10), seq(0.5, -1, length.out = p - 10))
  beta2.vec <- c(c(0, -1, 0, 1, 0), rep(0, p - 5))
  beta.mat  <- rbind(beta1.vec, beta2.vec)
  cov_dirs  <- which(beta2.vec != 0)                 # = c(2, 4)

  set.seed(100)
  X <- cbind(Intercept = 1, group = stats::rbinom(n, size = 1, prob = 0.5))

  # subject covariances: Sigma_i = Gamma diag(delta_i) Gamma'; delta covariate-
  # driven on cov_dirs, random elsewhere (drawn in-stream, as in the manuscript).
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
  set.seed(100)
  Y <- lapply(1:n, function(i) mvtnorm::rmvnorm(nT, rep(0, p), sigma = Sigma[, , i]))
  names(Y) <- paste0("S", 1:n)
  var_names <- paste0("V", 1:p)

  # truth = the covariate-driven directions (the components satisfying the model)
  gamma_true <- Gamma[, cov_dirs, drop = FALSE]
  colnames(gamma_true) <- paste0("D", seq_along(cov_dirs))
  beta_true <- beta.mat[, cov_dirs, drop = FALSE]
  rownames(beta_true) <- c("Intercept", "group")
  colnames(beta_true) <- paste0("D", seq_along(cov_dirs))

  d <- list(Y = Y, X = X, ids = names(Y), var_names = var_names,
            Xcov_df = data.frame(id = names(Y), group = X[, "group"]),
            truth = list(gamma = gamma_true, beta = beta_true))
  d$preview_ui <- hdcap_preview(d)
  d
}

# ---- 3. Parser for uploaded data -------------------------------------------
hdcap_parse <- function(files, opts) {
  yp <- coerce_Y_input(files$Y)
  X  <- coerce_X_input(files$X, ids = yp$ids,
                       add_intercept = isTRUE(opts$add_intercept))
  cov_df <- data.frame(id = yp$ids,
                       X[, setdiff(colnames(X), "Intercept"), drop = FALSE],
                       check.names = FALSE)
  d <- list(Y = yp$Y, X = X, ids = yp$ids, var_names = yp$var_names,
            Xcov_df = cov_df, truth = NULL)
  d$preview_ui <- hdcap_preview(d)
  d
}

# ---- 4. Data preview / helpers ---------------------------------------------
hdcap_preview <- function(d) {
  Tvec <- vapply(d$Y, nrow, integer(1))
  cov_names <- setdiff(colnames(d$X), "Intercept")
  p <- ncol(d$Y[[1]])
  hd_flag <- min(Tvec) - 5 < p
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", length(d$Y), theme = "primary"),
      bslib::value_box("Dimension (p)", p, theme = "secondary"),
      bslib::value_box("Obs / subject",
                       if (length(unique(Tvec)) == 1) Tvec[1]
                       else paste0(min(Tvec), "–", max(Tvec)),
                       theme = "secondary"),
      bslib::value_box("Covariates", length(cov_names), theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           sprintf("Covariates: %s.%s%s",
                   paste(cov_names, collapse = ", "),
                   if (hd_flag) "  High-dimensional (min Tᵢ − 5 < p): shrinkage auto-enabled." else "",
                   if (!is.null(d$truth)) "  Simulated data with known truth." else ""))
  )
}

hdcap_describe <- function(d) {
  sprintf("%d subjects, dimension p = %d, covariates: %s.",
          length(d$Y), ncol(d$Y[[1]]),
          paste(setdiff(colnames(d$X), "Intercept"), collapse = ", "))
}

# DfD avg.level vector (V7 returns DfD as a list even for one direction).
hdcap_dfd_vec <- function(fit) {
  if (is.list(fit$DfD)) as.numeric(fit$DfD$avg.level) else rep(1, ncol(fit$gamma))
}

# Greedy sign-aligned matching of estimated to true directions.
hdcap_match_dirs <- function(Ghat, Gtrue) {
  J <- ncol(Ghat); K <- ncol(Gtrue); used <- integer(0); out <- list()
  for (j in seq_len(J)) {
    ip <- vapply(seq_len(K), function(k)
      if (k %in% used) NA_real_ else sum(Ghat[, j] * Gtrue[, k]), numeric(1))
    if (all(is.na(ip))) break
    k <- which.max(abs(ip)); used <- c(used, k)
    out[[length(out) + 1]] <- data.frame(
      est = colnames(Ghat)[j] %||% paste0("D", j),
      true = colnames(Gtrue)[k] %||% paste0("D", k),
      k = k, sign = sign(ip[k]), cosine = abs(ip[k]))
  }
  do.call(rbind, out)
}

# ---- 5. Run -----------------------------------------------------------------
hdcap_run <- function(d, params) {
  Y <- d$Y; X <- d$X
  stop_crt <- params$stop.crt %||% "nD"
  shrink <- isTRUE(params$cov.shrinkage %||% TRUE)
  fit <- .hdcap_env$capReg(Y, X,
                stop.crt = stop_crt,
                nD = if (stop_crt == "nD") (params$nD %||% 2) else NULL,
                DfD.thred = params$DfD.thred %||% 2,
                method = "CAP",
                cov.shrinkage = shrink,
                ninitial = if (is.null(params$ninitial) || params$ninitial <= 0) NULL
                           else params$ninitial,
                score.return = TRUE, verbose = FALSE)

  # bootstrap inference for beta, per estimated direction
  sims <- params$sims %||% 0
  inference <- NULL
  if (sims > 0) {
    inference <- lapply(seq_len(ncol(fit$gamma)), function(j)
      tryCatch({
        bi <- .hdcap_env$cap_beta_boot(Y, X, gamma = fit$gamma[, j],
                            cov.shrinkage = shrink, sims = sims, verbose = FALSE)$Inference
        names(bi)[names(bi) == "Estiamte"] <- "Estimate"   # fix upstream typo
        bi
      }, error = function(e) NULL))
    names(inference) <- colnames(fit$gamma)
  }

  list(fit = fit, inference = inference, X = X, ids = d$ids,
       var_names = d$var_names, truth = d$truth, params = params)
}

# ---- 6. Tables --------------------------------------------------------------
hdcap_summarize <- function(res) {
  fit <- res$fit
  dirs <- colnames(fit$gamma)
  out <- list()

  # beta (bootstrap inference if available, else estimates)
  if (!is.null(res$inference)) {
    beta_rows <- do.call(rbind, lapply(seq_along(dirs), function(j) {
      inf <- res$inference[[j]]
      if (is.null(inf)) {
        data.frame(Direction = dirs[j], Term = rownames(fit$beta),
                   Estimate = fit$beta[, j], SE = NA, statistic = NA,
                   pvalue = NA, LB = NA, UB = NA, check.names = FALSE)
      } else {
        data.frame(Direction = dirs[j], Term = rownames(inf), inf,
                   check.names = FALSE, row.names = NULL)
      }
    }))
    out[["β (bootstrap inference)"]] <- beta_rows
  } else {
    bdf <- as.data.frame(fit$beta)
    out[["β (covariate effects)"]] <- cbind(Term = rownames(fit$beta), bdf)
  }

  # gamma loadings
  g <- as.data.frame(fit$gamma)
  out[["γ (loadings)"]] <- cbind(Variable = res$var_names %||% rownames(fit$gamma), g)

  # shrinkage weights per direction (fit$shrinkage is a list-matrix -> coerce)
  if (!is.null(fit$shrinkage)) {
    sk <- fit$shrinkage
    skdf <- as.data.frame(matrix(as.numeric(sk), nrow = nrow(sk)))
    names(skdf) <- colnames(sk)
    out[["Shrinkage weights"]] <- cbind(Direction = rownames(sk), skdf)
  }

  # orthogonality
  if (ncol(fit$gamma) > 1) {
    o <- as.data.frame(fit$orthogonality)
    out[["Orthogonality γ'γ"]] <- cbind(" " = rownames(fit$orthogonality), o)
  }

  # DfD across dimensions
  dfd_vec <- hdcap_dfd_vec(fit)
  out[["DfD across dimensions"]] <- data.frame(
    Dimension = dirs, `# directions (1..k)` = seq_along(dfd_vec),
    DfD = dfd_vec, check.names = FALSE)

  # truth comparison (example data)
  if (!is.null(res$truth)) {
    Gtrue <- as.matrix(res$truth$gamma)
    if (is.null(colnames(Gtrue))) colnames(Gtrue) <- paste0("D", seq_len(ncol(Gtrue)))
    mt <- hdcap_match_dirs(fit$gamma, Gtrue)
    cmp <- data.frame(Variable = res$var_names)
    for (r in seq_len(nrow(mt))) {
      cmp[[paste0(mt$true[r], "_true")]] <- Gtrue[, mt$k[r]]
      cmp[[paste0(mt$est[r], "_est")]]  <- fit$gamma[, mt$est[r]] * mt$sign[r]
    }
    out[["γ: truth vs estimate"]] <- cmp
    out[["Direction recovery"]] <- data.frame(
      Estimated = mt$est, `Matched true` = mt$true,
      `Cosine similarity` = round(mt$cosine, 4), check.names = FALSE)

    # beta: truth vs estimate (beta is invariant to gamma's sign, so no flip)
    if (!is.null(res$truth$beta)) {
      Btrue <- as.matrix(res$truth$beta)
      terms <- rownames(fit$beta)
      bcmp <- do.call(rbind, lapply(seq_len(nrow(mt)), function(r) {
        k <- mt$k[r]
        if (k > ncol(Btrue)) return(NULL)
        data.frame(Direction = mt$est[r], Term = terms,
                   True = Btrue[match(terms, rownames(Btrue)), k],
                   Estimate = fit$beta[, mt$est[r]],
                   check.names = FALSE, row.names = NULL)
      }))
      if (!is.null(bcmp)) out[["β: truth vs estimate"]] <- bcmp
    }
  }
  out
}

# ---- 7. Plots ---------------------------------------------------------------
hdcap_plots <- function(res) {
  fit <- res$fit
  plots <- list()
  dirs <- colnames(fit$gamma)
  vars <- res$var_names %||% rownames(fit$gamma)

  # (a) loadings
  g <- as.data.frame(fit$gamma)
  long <- do.call(rbind, lapply(seq_len(ncol(g)), function(j)
    data.frame(Variable = vars, Direction = colnames(g)[j], Loading = g[, j])))
  long$Variable <- factor(long$Variable, levels = vars)
  plots[["Loadings (γ)"]] <- list(plot = plotly::plot_ly(
    long, x = ~Variable, y = ~Loading, color = ~Direction, type = "bar") |>
      plotly::layout(title = "Projection loadings (γ)", barmode = "group"))

  # (b) subject scores vs fitted (all directions) + downloadable data
  if (!is.null(fit$score)) {
    sc_long <- do.call(rbind, lapply(seq_along(dirs), function(j)
      data.frame(Direction = dirs[j],
                 fitted = as.numeric(res$X %*% fit$beta[, j]),
                 logvar = log(pmax(fit$score[, j], 1e-12)))))
    rng <- range(c(sc_long$fitted, sc_long$logvar), finite = TRUE)
    p_sc <- plotly::plot_ly(sc_long, x = ~fitted, y = ~logvar, color = ~Direction,
                            type = "scatter", mode = "markers",
                            marker = list(size = 8, opacity = 0.7)) |>
      plotly::add_lines(x = rng, y = rng, line = list(dash = "dot", color = "grey"),
                        showlegend = FALSE, inherit = FALSE) |>
      plotly::layout(title = "Subject projected variance vs fitted",
                     xaxis = list(title = "x'β (fitted log-variance)"),
                     yaxis = list(title = "log(γ' Σ γ)"))
    score_csv <- data.frame(id = res$ids %||% seq_len(nrow(fit$score)))
    for (j in seq_along(dirs)) {
      score_csv[[paste0("score_", dirs[j])]] <- fit$score[, j]
      if (!is.null(fit$score.shrinkage))
        score_csv[[paste0("score_shrunk_", dirs[j])]] <- fit$score.shrinkage[, j]
      score_csv[[paste0("fitted_", dirs[j])]] <- as.numeric(res$X %*% fit$beta[, j])
    }
    plots[["Subject scores"]] <- list(plot = p_sc, data = score_csv)
  }

  # (c) DfD across dimensions
  dfd_vec <- hdcap_dfd_vec(fit)
  dfd_df <- data.frame(k = seq_along(dfd_vec), DfD = dfd_vec, Dimension = dirs)
  p_dfd <- plotly::plot_ly(dfd_df, x = ~k, y = ~DfD, type = "scatter",
                           mode = "lines+markers",
                           marker = list(size = 9, color = "#2c7fb8"),
                           line = list(color = "#2c7fb8"),
                           text = ~Dimension,
                           hovertemplate = "%{text}: %{y:.4f}<extra></extra>")
  if (identical(res$params$stop.crt, "DfD") && !is.null(res$params$DfD.thred)) {
    thr <- res$params$DfD.thred
    p_dfd <- plotly::add_lines(p_dfd, x = range(dfd_df$k), y = c(thr, thr),
                               line = list(dash = "dash", color = "firebrick"),
                               name = paste0("threshold = ", thr), inherit = FALSE)
  }
  p_dfd <- plotly::layout(p_dfd, title = "DfD across dimensions",
                          xaxis = list(title = "Number of directions (1..k)",
                                       tickvals = dfd_df$k, ticktext = dfd_df$Dimension),
                          yaxis = list(title = "DfD (avg. deviation from diagonality)"))
  plots[["DfD"]] <- list(plot = p_dfd)

  plots
}

# ---- 8. Register ------------------------------------------------------------
register_method(list(
  id = "hdcap",
  name = "HDCAP",
  full_name = "High-Dimensional Covariate Assisted Principal Regression",
  short = "CAP for high-dimensional covariance outcomes via linear covariance shrinkage.",
  status = "ready",
  tags = c("covariance regression", "high-dimensional", "shrinkage"),
  paper = list(
    citation = "Zhao, Y., Caffo, B., Luo, X., & Alzheimer's Disease Neuroimaging Initiative. (2021). Principal regression for high dimensional covariance matrices. Electronic journal of statistics, 15(2), 4192.",
    url = "https://doi.org/10.1214/21-EJS1887"),
  explain = file.path(APP_DIR, "methods", "hdcap", "explain.md"),
  example_note = paste("The manuscript's HD-shrinkage simulation (210309/eg): 100 subjects,",
                       "p=20, Tᵢ=100, one binary covariate (group). TWO directions satisfy the",
                       "CAP model (group effect β = -1 and +1); with shrinkage on, nD=2 recovers",
                       "both (γ cosine ≈ 0.99, β̂(group) ≈ ∓1). Higher nD / DfD also pick up the",
                       "high-variance covariate-free directions (β̂(group) ≈ 0)."),
  x_intercept_option = TRUE,
  data_inputs = list(
    list(id = "Y", label = "Response (list of matrices)",
         help = "An .rds/.RData file holding a list of length n; element i is a T_i x p matrix. A long CSV (subject-id column + p response columns) is also accepted."),
    list(id = "X", label = "Covariates (one row per subject)",
         help = "An n x q covariate matrix in subject order, or a CSV whose first column is the subject id.")
  ),
  params = list(
    list(id = "cov.shrinkage", label = "Covariance shrinkage", type = "checkbox",
         default = TRUE,
         help = "Linear shrinkage of each subject covariance (auto-enabled when min Tᵢ − 5 < p)."),
    list(id = "stop.crt", label = "Number of directions chosen by", type = "select",
         choices = c("Fixed number (nD)" = "nD", "DfD criterion" = "DfD"),
         default = "nD",
         help = "Fix the count, or keep adding directions while DfD stays below the threshold."),
    list(id = "nD", label = "Number of directions (nD)", type = "integer",
         default = 2, min = 1, max = 15,
         help = "Used when selecting by a fixed number. The example has two covariate-driven directions; nD=2 recovers both."),
    list(id = "DfD.thred", label = "DfD threshold", type = "numeric",
         default = 2, min = 1, max = 100, step = 0.5,
         help = "Used when selecting by the DfD criterion."),
    list(id = "ninitial", label = "# random initializations", type = "integer",
         default = 5, min = 1, max = 30,
         help = "More initializations = more robust optimization, slower."),
    list(id = "sims", label = "Bootstrap replicates (β inference)", type = "integer",
         default = 100, min = 0, max = 2000,
         help = "Subject-level bootstrap for β inference; 0 = skip (estimates only).")
  ),
  example = hdcap_example,
  export_example = function(d) list(Y = d$Y, X = d$Xcov_df),
  parse = hdcap_parse,
  describe_data = hdcap_describe,
  run = hdcap_run,
  summarize = hdcap_summarize,
  plots = hdcap_plots
))
