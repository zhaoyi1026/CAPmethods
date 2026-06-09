# =============================================================================
# CAP-HDcov (HCAP) -- CAP regression with High-Dimensional Covariates  (plugin)
# -----------------------------------------------------------------------------
# Wraps ../CAP-HDcov/V2/ (RcppArmadillo covariance kernel; HCAP_Code.R). Unlike
# the other CAP methods, the modelling challenge here is a LARGE number of
# COVARIATES q: a shared projection gamma summarizes each subject's response
# covariance, and the projected log-variance is regressed on the high-dim
# covariates with sparse (Gamma-lasso) estimation + multi-split selection
# inference (SPARE / SSHDI).
# =============================================================================

# ---- 1. Load the method code from the CAPmethods package --------------------
# Uses the package's "hdcov" environment (capReg(), SSHDI(), ...); the
# hdcov_-prefixed kernel ships in the package DLL. See R/pkg_methods.R.
.hdcov_env <- cap_pkg_env("hdcov")
# SSHDI's multi-split inference uses foreach %dopar%; its forked workers cannot
# see this private env / the package kernels, so every split would silently
# fail. Force the inference to run IN-PROCESS by making SSHDI's internal
# registerDoParallel() register the sequential backend instead. This override
# takes effect because cap_pkg_env() rebinds SSHDI's closure to .hdcov_env.
.hdcov_env$registerDoParallel <- function(...) foreach::registerDoSEQ()

# ---- 2. Built-in example data (high-dimensional covariates) ------------------
# n subjects, p responses, q covariates (q >> p). Two covariate-driven covariance
# directions; each is driven by a small sparse set of the q covariates. Mirrors
# the paper's simulate_data_hcap (q=200) -- the scale at which the sqrt(q)
# selection window works.
hdcov_example <- function() {
  set.seed(2023)
  n <- 100L; p <- 5L; q <- 200L; Ti <- 100L
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
  lp <- matrix(NA, n, p); i1 <- colSums(b0[-1, ]) != 0
  lp[, i1]  <- exp(X %*% b0[, i1])
  lp[, !i1] <- exp(sapply(b0[1, !i1], function(x) stats::rnorm(n, x, 0.5)))
  Y <- lapply(1:n, function(i)
    mvtnorm::rmvnorm(Ti, rep(0, p), phi %*% diag(lp[i, ]) %*% t(phi)))
  names(Y) <- paste0("S", 1:n)
  cov_names <- c("Intercept", paste0("X", 1:(q - 1)))
  colnames(X) <- cov_names
  var_names <- paste0("V", 1:p)
  # the two covariate-driven true directions and their sparse beta vectors
  Gtrue <- phi[, c(2, 3)]; colnames(Gtrue) <- c("T2", "T3")
  Btrue <- b0[, c(2, 3)]; colnames(Btrue) <- c("T2", "T3"); rownames(Btrue) <- cov_names
  d <- list(Y = Y, X = X, ids = names(Y), var_names = var_names, cov_names = cov_names,
            Xdf = data.frame(id = names(Y), X, check.names = FALSE),
            truth = list(gamma = Gtrue, beta = Btrue,
                         signal = list(T2 = s2, T3 = s3)))
  d$preview_ui <- hdcov_preview(d)
  d
}

# ---- 3. Parser for uploaded data (cross-sectional: reuse generic helpers) ----
hdcov_parse <- function(files, opts) {
  yp <- coerce_Y_input(files$Y)
  X  <- coerce_X_input(files$X, ids = yp$ids, add_intercept = isTRUE(opts$add_intercept))
  cov_names <- colnames(X) %||% paste0("X", seq_len(ncol(X)))
  colnames(X) <- cov_names
  d <- list(Y = yp$Y, X = X, ids = yp$ids, var_names = yp$var_names,
            cov_names = cov_names, Xdf = data.frame(id = yp$ids, X, check.names = FALSE),
            truth = NULL)
  d$preview_ui <- hdcov_preview(d)
  d
}

# ---- 4. Preview / helpers ---------------------------------------------------
hdcov_preview <- function(d) {
  Tvec <- vapply(d$Y, nrow, integer(1))
  q <- ncol(d$X); p <- ncol(d$Y[[1]])
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", length(d$Y), theme = "primary"),
      bslib::value_box("Responses (p)", p, theme = "secondary"),
      bslib::value_box("Covariates (q)", q, theme = "secondary"),
      bslib::value_box("Obs / subject",
                       if (length(unique(Tvec)) == 1) Tvec[1] else paste0(min(Tvec), "–", max(Tvec)),
                       theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           sprintf("High-dimensional covariates (q = %d).%s", q,
                   if (!is.null(d$truth)) "  Simulated data: true signal covariates are recovered by the selection inference." else "")))
}

hdcov_describe <- function(d)
  sprintf("%d subjects, p = %d responses, q = %d covariates.",
          length(d$Y), ncol(d$Y[[1]]), ncol(d$X))

# ---- 5. Run -----------------------------------------------------------------
hdcov_run <- function(d, params) {
  Y <- d$Y; X <- d$X
  stop_crt <- params$stop.crt %||% "nD"
  nD <- if (stop_crt == "nD") (params$nD %||% 1) else NULL
  ninitial <- max(2L, as.integer(params$ninitial %||% 2))   # ninitial=1 is unsupported upstream
  B <- as.integer(params$B %||% 20)
  foreach::registerDoSEQ()
  cap_res <- suppressWarnings(.hdcov_env$capReg(Y, X, stop.crt = stop_crt, nD = nD,
               DfD.thred = params$DfD.thred %||% 2, CAP.OC = TRUE,
               ninitial = ninitial, max.itr = 50L, tol = 1e-3, trace = FALSE))
  inference <- lapply(seq_len(ncol(cap_res$gamma)), function(j)
    tryCatch(suppressWarnings(.hdcov_env$SSHDI(Xmat = X, Ymat = Y, B = B,
             gamma_ini = cap_res$gamma[, j])), error = function(e) NULL))
  names(inference) <- colnames(cap_res$gamma) %||% paste0("D", seq_along(inference))
  list(fit = cap_res, inference = inference, X = X, ids = d$ids,
       var_names = d$var_names, cov_names = d$cov_names, truth = d$truth, params = params)
}

# ---- 6. Tables --------------------------------------------------------------
hdcov_summarize <- function(res) {
  fit <- res$fit
  dirs <- colnames(fit$gamma) %||% paste0("D", seq_len(ncol(fit$gamma)))
  cov_names <- res$cov_names %||% paste0("X", seq_len(nrow(fit$gamma)) )
  q <- length(res$cov_names)
  out <- list()

  # beta selection inference, one block per direction (sorted by selection freq)
  beta_rows <- do.call(rbind, lapply(seq_along(dirs), function(j) {
    inf <- res$inference[[j]]
    if (is.null(inf)) return(NULL)
    df <- data.frame(Direction = dirs[j], Covariate = res$cov_names,
                     Estimate = round(inf$ss.beta, 4), SD = round(inf$sd, 4),
                     `p-value` = signif(inf$p, 3),
                     `Sel.freq` = round(inf$sel.freq, 3),
                     Selected = ifelse(inf$sel.freq > 0.5, "✓", ""),
                     check.names = FALSE)
    df[order(-inf$sel.freq, inf$p), ]
  }))
  if (!is.null(beta_rows)) out[["β selection inference"]] <- beta_rows

  # selected covariates summary per direction
  sel_rows <- do.call(rbind, lapply(seq_along(dirs), function(j) {
    inf <- res$inference[[j]]; if (is.null(inf)) return(NULL)
    sel <- which(inf$sel.freq > 0.5)
    if (!length(sel)) return(data.frame(Direction = dirs[j], `# selected` = 0,
                                        Covariates = "(none)", check.names = FALSE))
    data.frame(Direction = dirs[j], `# selected` = length(sel),
               Covariates = paste(res$cov_names[sel], collapse = ", "), check.names = FALSE)
  }))
  if (!is.null(sel_rows)) out[["Selected covariates"]] <- sel_rows

  # gamma loadings
  g <- as.data.frame(fit$gamma); names(g) <- dirs
  out[["γ (loadings)"]] <- cbind(Variable = res$var_names, g)

  # orthogonality
  if (ncol(fit$gamma) > 1 && !is.null(fit$orthogonality)) {
    o <- as.data.frame(fit$orthogonality)
    out[["Orthogonality γ'γ"]] <- cbind(" " = dirs, o)
  }

  # truth comparison (example): which true signal covariates were recovered
  if (!is.null(res$truth)) {
    Gtrue <- as.matrix(res$truth$gamma)
    # match each estimated direction to its closest true direction
    rec <- do.call(rbind, lapply(seq_along(dirs), function(j) {
      ip <- abs(as.numeric(t(fit$gamma[, j]) %*% Gtrue))
      k <- which.max(ip); tname <- colnames(Gtrue)[k]
      inf <- res$inference[[j]]
      true_sig <- res$truth$signal[[tname]]
      sel <- if (is.null(inf)) integer(0) else which(inf$sel.freq > 0.5)
      data.frame(Direction = dirs[j], `Matched true` = tname,
                 `γ cosine` = round(ip[k], 4),
                 `True signal covariates` = paste(res$cov_names[true_sig], collapse = ", "),
                 `Recovered` = paste(intersect(res$cov_names[true_sig], res$cov_names[sel]), collapse = ", "),
                 `False positives` = length(setdiff(sel, c(1, true_sig))),
                 check.names = FALSE)
    }))
    out[["Direction & signal recovery"]] <- rec

    # beta truth vs estimate for the matched direction's true signal covariates
    bcmp <- do.call(rbind, lapply(seq_along(dirs), function(j) {
      ip <- abs(as.numeric(t(fit$gamma[, j]) %*% Gtrue)); k <- which.max(ip)
      tname <- colnames(Gtrue)[k]; inf <- res$inference[[j]]; if (is.null(inf)) return(NULL)
      sig <- res$truth$signal[[tname]]
      data.frame(Direction = dirs[j], Covariate = res$cov_names[sig],
                 True = res$truth$beta[sig, tname],
                 Estimate = round(inf$ss.beta[sig], 3),
                 `Sel.freq` = round(inf$sel.freq[sig], 3), check.names = FALSE)
    }))
    if (!is.null(bcmp)) out[["β: truth vs estimate (signal)"]] <- bcmp
  }
  out
}

# ---- 7. Plots ---------------------------------------------------------------
hdcov_plots <- function(res) {
  fit <- res$fit; plots <- list()
  dirs <- colnames(fit$gamma) %||% paste0("D", seq_len(ncol(fit$gamma)))
  vars <- res$var_names

  # (a) gamma loadings
  g <- as.data.frame(fit$gamma); names(g) <- dirs
  long <- do.call(rbind, lapply(seq_along(dirs), function(j)
    data.frame(Variable = vars, Direction = dirs[j], Loading = g[, j])))
  long$Variable <- factor(long$Variable, levels = vars)
  plots[["Loadings (γ)"]] <- list(plot = plotly::plot_ly(
    long, x = ~Variable, y = ~Loading, color = ~Direction, type = "bar") |>
      plotly::layout(title = "Projection loadings (γ)", barmode = "group"))

  # (b) selection frequency across covariates (per direction)
  seldf <- do.call(rbind, lapply(seq_along(dirs), function(j) {
    inf <- res$inference[[j]]; if (is.null(inf)) return(NULL)
    data.frame(cov = seq_along(inf$sel.freq), Covariate = res$cov_names,
               sel.freq = inf$sel.freq, Direction = dirs[j])
  }))
  if (!is.null(seldf)) {
    p_sel <- plotly::plot_ly(seldf, x = ~cov, y = ~sel.freq, color = ~Direction,
                             type = "bar", text = ~Covariate,
                             hovertemplate = "%{text}: %{y:.2f}<extra></extra>") |>
      plotly::layout(title = "Covariate selection frequency (multi-split)",
                     xaxis = list(title = "covariate index"),
                     yaxis = list(title = "selection frequency", range = c(0, 1)),
                     shapes = list(list(type = "line", x0 = 0, x1 = max(seldf$cov),
                                        y0 = 0.5, y1 = 0.5, line = list(dash = "dash", color = "firebrick"))))
    plots[["Selection frequency"]] <- list(plot = p_sel, data = seldf)
  }

  # (c) -log10(p) manhattan-style, with Bonferroni line
  pdf <- do.call(rbind, lapply(seq_along(dirs), function(j) {
    inf <- res$inference[[j]]; if (is.null(inf)) return(NULL)
    data.frame(cov = seq_along(inf$p), Covariate = res$cov_names,
               neglogp = -log10(pmax(inf$p, 1e-300)), Direction = dirs[j])
  }))
  if (!is.null(pdf)) {
    q <- length(res$cov_names); thr <- -log10(0.05 / q)
    p_man <- plotly::plot_ly(pdf, x = ~cov, y = ~neglogp, color = ~Direction,
                             type = "scatter", mode = "markers",
                             marker = list(size = 6), text = ~Covariate,
                             hovertemplate = "%{text}: -log10 p = %{y:.2f}<extra></extra>") |>
      plotly::layout(title = "Covariate significance (−log₁₀ p)",
                     xaxis = list(title = "covariate index"),
                     yaxis = list(title = "−log₁₀ p"),
                     shapes = list(list(type = "line", x0 = 0, x1 = q, y0 = thr, y1 = thr,
                                        line = list(dash = "dash", color = "firebrick"))))
    plots[["Significance (−log₁₀ p)"]] <- list(plot = p_man, data = pdf)
  }
  plots
}

# ---- 8. Register ------------------------------------------------------------
register_method(list(
  id = "cap-hdcov",
  name = "CAP-HDcov",
  full_name = "Covariate Assisted Principal Regression with High-Dimensional Covariates (HCAP)",
  short = "CAP where the number of covariates q is large: sparse (Gamma-lasso) projected-variance regression with multi-split selection inference.",
  status = "ready",
  tags = c("covariance regression", "high-dimensional covariates", "sparse", "selection inference"),
  paper = list(
    citation = "Covariate assisted principal regression with high-dimensional covariates (see AOAS2157.pdf).",
    url = NULL),
  explain = file.path(APP_DIR, "methods", "cap-hdcov", "explain.md"),
  example_note = paste("Simulated high-dimensional example: 100 subjects, p=5 responses,",
                       "q=200 covariates, Tᵢ=100, with covariate-driven directions whose",
                       "true signal covariates the selection inference recovers.",
                       "Runs in ~30–40 s (glmnet-bound)."),
  x_intercept_option = TRUE,
  data_inputs = list(
    list(id = "Y", label = "Response (list of matrices)",
         help = "An .rds/.RData holding a list of length n; element i is a T_i x p matrix. A long CSV (subject-id column + p response columns) is also accepted."),
    list(id = "X", label = "Covariates (n x q, high-dimensional)",
         help = "An n x q covariate matrix in subject order (q may be large), or a CSV whose first column is the subject id. Intercept auto-added (optional).")
  ),
  params = list(
    list(id = "stop.crt", label = "Number of directions chosen by", type = "select",
         choices = c("Fixed number (nD)" = "nD", "DfD criterion" = "DfD"), default = "nD",
         help = "Fix the count, or keep adding directions while DfD stays below the threshold."),
    list(id = "nD", label = "Number of directions (nD)", type = "integer",
         default = 1, min = 1, max = 5, help = "Used when selecting by a fixed number. Each direction adds an inference pass (slower)."),
    list(id = "DfD.thred", label = "DfD threshold", type = "numeric",
         default = 2, min = 1, max = 100, step = 0.5, help = "Used when selecting by the DfD criterion."),
    list(id = "ninitial", label = "# random initializations", type = "integer",
         default = 2, min = 2, max = 10, help = "Random restarts for direction estimation (≥2)."),
    list(id = "B", label = "Multi-split resamples (inference)", type = "integer",
         default = 20, min = 5, max = 200, help = "Sample-split + aggregation resamples for selection inference (more = stabler, slower).")
  ),
  example = hdcov_example,
  export_example = function(d) list(Y = d$Y, X = d$Xdf),
  parse = hdcov_parse,
  describe_data = hdcov_describe,
  run = hdcov_run,
  summarize = hdcov_summarize,
  plots = hdcov_plots
))
