# =============================================================================
# CAP -- Covariate Assisted Principal regression  (plugin)
# -----------------------------------------------------------------------------
# Wraps the original implementation in ../CAP/V6/ (RcppArmadillo-accelerated)
# and registers it with the website. This file is the template every other
# CAP-family method follows.
# =============================================================================

# ---- 1. Load the method code from the CAPmethods package --------------------
# Classical CAP = HDCAP with cov.shrinkage = FALSE, so this uses the package's
# "hdcap" environment (capReg(), cap_beta(), ...). See R/pkg_methods.R. The
# private env keeps these functions from colliding with other plugins.
# (Note: the CAP-C / common-PC variant is not in the package — multigroup was
# dropped — so only the optimization estimator "CAP" is offered here.)
.cap_env <- cap_pkg_env("hdcap")

# ---- 2. Built-in example data ----------------------------------------------
# Simulates the CAP generative model so the estimate can be checked against a
# known truth: log(gamma' Sigma_i gamma) = x_i' beta.
cap_example <- function() {
  set.seed(2024)
  n <- 80L            # subjects
  p <- 6L             # response variables (e.g. regions)
  Ti <- 120L          # observations per subject

  # covariates: intercept + binary group + continuous score
  group <- rep(c(0, 1), length.out = n)
  age   <- round(scale(stats::rnorm(n, 50, 10))[, 1], 2)
  Xcov  <- cbind(group = group, age = age)
  X     <- cbind(Intercept = 1, Xcov)

  # two covariate-driven directions with distinct effects
  beta_true <- cbind(D1 = c(Intercept = 0.0, group = 0.9, age = -0.5),
                     D2 = c(Intercept = 0.0, group = -0.7, age = 0.4))
  nD_true <- ncol(beta_true)

  # random orthonormal basis; the first nD_true columns are the true directions
  A <- matrix(stats::rnorm(p * p), p, p)
  Phi <- qr.Q(qr(A))
  for (k in seq_len(nD_true)) {           # sign convention: largest |loading| > 0
    if (Phi[which.max(abs(Phi[, k])), k] < 0) Phi[, k] <- -Phi[, k]
  }
  gamma_true <- Phi[, seq_len(nD_true), drop = FALSE]
  base_eigs <- c(0.6, 0.45, 0.3, 0.2)     # background directions: fixed variance

  Y <- vector("list", n)
  for (i in seq_len(n)) {
    eigs <- c(exp(as.numeric(X[i, ] %*% beta_true)),    # covariate-driven (D1,D2)
              base_eigs)[seq_len(p)]
    Sigma_i <- Phi %*% diag(eigs) %*% t(Phi)
    Sigma_i <- (Sigma_i + t(Sigma_i)) / 2
    ev <- eigen(Sigma_i, symmetric = TRUE)
    root <- ev$vectors %*% diag(sqrt(pmax(ev$values, 1e-8))) %*% t(ev$vectors)
    Y[[i]] <- matrix(stats::rnorm(Ti * p), Ti, p) %*% root
  }
  names(Y) <- paste0("S", seq_len(n))
  var_names <- paste0("V", seq_len(p))

  d <- list(
    Y = Y, X = X, ids = names(Y), var_names = var_names,
    Xcov_df = data.frame(id = names(Y), Xcov),
    truth = list(gamma = gamma_true, beta = beta_true)
  )
  d$preview_ui <- cap_preview(d)
  d
}

# ---- 3. Parser for uploaded data -------------------------------------------
# The response may be a native list of n matrices (each T_i x p, uploaded as
# .rds/.RData) or a long-format CSV; covariates may be an n x q matrix or a CSV.
cap_parse <- function(files, opts) {
  yp <- coerce_Y_input(files$Y)
  X  <- coerce_X_input(files$X, ids = yp$ids,
                       add_intercept = isTRUE(opts$add_intercept))
  cov_df <- data.frame(id = yp$ids,
                       X[, setdiff(colnames(X), "Intercept"), drop = FALSE],
                       check.names = FALSE)
  d <- list(Y = yp$Y, X = X, ids = yp$ids, var_names = yp$var_names,
            Xcov_df = cov_df, truth = NULL)
  d$preview_ui <- cap_preview(d)
  d
}

# ---- 4. Data preview UI -----------------------------------------------------
cap_preview <- function(d) {
  Tvec <- vapply(d$Y, nrow, integer(1))
  cov_names <- setdiff(colnames(d$X), "Intercept")
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", length(d$Y), theme = "primary"),
      bslib::value_box("Responses (p)", ncol(d$Y[[1]]), theme = "secondary"),
      bslib::value_box("Obs / subject",
                       if (length(unique(Tvec)) == 1) Tvec[1]
                       else paste0(min(Tvec), "–", max(Tvec)),
                       theme = "secondary"),
      bslib::value_box("Covariates", length(cov_names), theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           sprintf("Covariates: %s%s",
                   paste(cov_names, collapse = ", "),
                   if (!is.null(d$truth)) "  |  simulated data with known truth" else ""))
  )
}

cap_describe <- function(d) {
  sprintf("%d subjects, %d responses, covariates: %s.",
          length(d$Y), ncol(d$Y[[1]]),
          paste(setdiff(colnames(d$X), "Intercept"), collapse = ", "))
}

# ---- 5. Run -----------------------------------------------------------------
cap_run <- function(d, params) {
  Y <- d$Y; X <- d$X
  stop_crt <- params$stop.crt %||% "nD"
  fit <- .cap_env$capReg(Y, X,
                stop.crt = stop_crt,
                nD = if (stop_crt == "nD") params$nD else NULL,
                DfD.thred = params$DfD.thred %||% 2,
                OC = isTRUE(params$OC),
                ninitial = if (is.null(params$ninitial) || params$ninitial <= 0) NULL
                           else params$ninitial,
                score.return = TRUE)
  # asymptotic inference on every estimated direction's beta
  inference <- lapply(seq_len(ncol(fit$gamma)), function(j)
    tryCatch(.cap_env$cap_beta(Y, X, gamma = fit$gamma[, j], beta = fit$beta[, j],
                      method = "asmp"),
             error = function(e) NULL))
  names(inference) <- colnames(fit$gamma)

  list(fit = fit, inference = inference, X = X, ids = d$ids,
       var_names = d$var_names, truth = d$truth, params = params)
}

# Extract the DfD avg.level vector (one value per cumulative set of directions).
# capReg returns fit$DfD as a list(avg.level, sub.level) when >1 direction, or
# the scalar 1 when a single direction was estimated.
cap_dfd_vec <- function(fit) {
  if (is.list(fit$DfD)) as.numeric(fit$DfD$avg.level)
  else rep(1, ncol(fit$gamma))
}

# Match estimated directions to true directions by maximum absolute inner
# product (greedy), aligning signs. Returns a data.frame: est, true, cosine, sign.
cap_match_dirs <- function(Ghat, Gtrue) {
  J <- ncol(Ghat); K <- ncol(Gtrue)
  used <- integer(0); out <- list()
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

# ---- 6. Tables --------------------------------------------------------------
cap_summarize <- function(res) {
  fit <- res$fit
  dirs <- colnames(fit$gamma)
  out <- list()

  # beta with asymptotic inference -- one tidy block per estimated direction
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
  out[["β (covariate effects)"]] <- beta_rows

  # gamma loadings (all directions)
  g <- as.data.frame(fit$gamma)
  g <- cbind(Variable = res$var_names %||% rownames(fit$gamma), g)
  rownames(g) <- NULL
  out[["γ (loadings)"]] <- g

  # orthogonality (if >1 direction)
  if (ncol(fit$gamma) > 1) {
    o <- as.data.frame(fit$orthogonality)
    out[["Orthogonality γ'γ"]] <- cbind(" " = rownames(fit$orthogonality), o)
  }

  # DfD (deviation-from-diagonality) across dimensions: avg.level[k] is the DfD
  # of the joint set of the first k directions (1 = perfectly diagonal).
  dfd_vec <- cap_dfd_vec(fit)
  out[["DfD across dimensions"]] <- data.frame(
    Dimension = colnames(fit$gamma),
    `# directions (1..k)` = seq_along(dfd_vec),
    DfD = dfd_vec, check.names = FALSE)

  # truth comparison for the example data (matched & sign-aligned per direction)
  if (!is.null(res$truth)) {
    Gtrue <- as.matrix(res$truth$gamma)
    if (is.null(colnames(Gtrue)))
      colnames(Gtrue) <- paste0("D", seq_len(ncol(Gtrue)))
    mt <- cap_match_dirs(fit$gamma, Gtrue)
    cmp <- data.frame(Variable = res$var_names)
    for (r in seq_len(nrow(mt))) {
      tk <- Gtrue[, mt$k[r]]
      ek <- fit$gamma[, mt$est[r]] * mt$sign[r]
      cmp[[paste0(mt$true[r], "_true")]] <- tk
      cmp[[paste0(mt$est[r], "_est")]]  <- ek
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
cap_plots <- function(res) {
  fit <- res$fit
  plots <- list()

  # (a) loadings bar chart across directions
  g <- as.data.frame(fit$gamma)
  vars <- res$var_names %||% rownames(fit$gamma)
  long <- do.call(rbind, lapply(seq_len(ncol(g)), function(j) {
    data.frame(Variable = vars, Direction = colnames(g)[j], Loading = g[, j])
  }))
  long$Variable <- factor(long$Variable, levels = vars)
  p_load <- plotly::plot_ly(long, x = ~Variable, y = ~Loading, color = ~Direction,
                            type = "bar") |>
    plotly::layout(title = "Projection loadings (γ)",
                   barmode = "group", yaxis = list(title = "Loading"))
  plots[["Loadings (γ)"]] <- list(plot = p_load)

  # (b) projected variance (score) vs model-fitted log-variance, all directions
  if (!is.null(fit$score)) {
    dirs <- colnames(fit$gamma)
    sc_long <- do.call(rbind, lapply(seq_along(dirs), function(j) {
      data.frame(Direction = dirs[j],
                 fitted = as.numeric(res$X %*% fit$beta[, j]),
                 logvar = log(pmax(fit$score[, j], 1e-12)))
    }))
    rng <- range(c(sc_long$fitted, sc_long$logvar), finite = TRUE)
    p_sc <- plotly::plot_ly(sc_long, x = ~fitted, y = ~logvar,
                            color = ~Direction, type = "scatter",
                            mode = "markers",
                            marker = list(size = 8, opacity = 0.7)) |>
      plotly::add_lines(x = rng, y = rng, line = list(dash = "dot",
                        color = "grey"), showlegend = FALSE,
                        inherit = FALSE) |>
      plotly::layout(title = "Subject projected variance vs fitted",
                     xaxis = list(title = "x'β (fitted log-variance)"),
                     yaxis = list(title = "log(γ' Σ γ)"))
    # downloadable wide table: subject id + per-direction score and fitted value
    dirs <- colnames(fit$gamma)
    score_csv <- data.frame(id = res$ids %||% seq_len(nrow(fit$score)))
    for (j in seq_along(dirs)) {
      score_csv[[paste0("score_", dirs[j])]] <- fit$score[, j]
      score_csv[[paste0("fitted_", dirs[j])]] <- as.numeric(res$X %*% fit$beta[, j])
    }
    plots[["Subject scores"]] <- list(plot = p_sc, data = score_csv)
  }

  # (c) DfD across dimensions (with the threshold line when DfD selection used)
  dfd_vec <- cap_dfd_vec(fit)
  dfd_df <- data.frame(k = seq_along(dfd_vec), DfD = dfd_vec,
                       Dimension = colnames(fit$gamma))
  p_dfd <- plotly::plot_ly(dfd_df, x = ~k, y = ~DfD, type = "scatter",
                           mode = "lines+markers",
                           marker = list(size = 9, color = "#2c7fb8"),
                           line = list(color = "#2c7fb8"),
                           text = ~Dimension, hovertemplate = "%{text}: %{y:.4f}<extra></extra>")
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
  id = "cap",
  name = "CAP",
  full_name = "Covariate Assisted Principal (CAP) Regression",
  short = "Find linear projections whose variance is associated with covariates.",
  status = "ready",
  tags = c("covariance regression", "dimension reduction"),
  paper = list(
    citation = "Zhao, Y., Wang, B., Mostofsky, S. H., Caffo, B. S., & Luo, X. (2021). Covariate assisted principal regression for covariance matrix outcomes. Biostatistics, 22(3), 629-645.",
    url = "https://doi.org/10.1093/biostatistics/kxz057"),
  explain = file.path(APP_DIR, "methods", "cap", "explain.md"),
  example_note = paste("Simulated: 80 subjects, 6 responses, covariates group & age,",
                       "with TWO covariate-driven directions",
                       "(D1: β(group)=0.9, β(age)=-0.5; D2: β(group)=-0.7, β(age)=0.4).",
                       "Set nD=2 (or use DfD) to recover both."),
  x_intercept_option = TRUE,
  data_inputs = list(
    list(id = "Y", label = "Response (list of matrices)",
         help = "An .rds/.RData file holding a list of length n (one element per subject); element i is a T_i x p matrix (T_i = samples in subject i, p = data dimension, e.g. brain regions). A long CSV (subject-id column + p response columns) is also accepted."),
    list(id = "X", label = "Covariates (one row per subject)",
         help = "An n x q covariate matrix with rows in the same subject order as the response list, or a CSV whose first column is the subject id.")
  ),
  params = list(
    list(id = "stop.crt", label = "Number of directions chosen by", type = "select",
         choices = c("Fixed number (nD)" = "nD",
                     "DfD criterion" = "DfD"), default = "nD",
         help = "Fix the count, or keep adding directions while deviation-from-diagonality (DfD) stays below the threshold."),
    list(id = "nD", label = "Number of directions (nD)", type = "integer",
         default = 1, min = 1, max = 10,
         help = "Used when selecting by a fixed number."),
    list(id = "DfD.thred", label = "DfD threshold", type = "numeric",
         default = 2, min = 1, max = 50, step = 0.5,
         help = "Used when selecting by the DfD criterion (e.g. keep directions with DfD ≤ 2)."),
    list(id = "OC", label = "Orthogonal constraint (OC)", type = "checkbox",
         default = FALSE),
    list(id = "ninitial", label = "# random initializations", type = "integer",
         default = 5, min = 1, max = 30,
         help = "More initializations = more robust optimization, slower.")
  ),
  example = cap_example,
  export_example = function(d) list(
    Y = d$Y,           # native list of n matrices -> downloads as .rds
    X = d$Xcov_df      # covariate table          -> downloads as .csv
  ),
  parse = cap_parse,
  describe_data = cap_describe,
  run = cap_run,
  summarize = cap_summarize,
  plots = cap_plots
))
