# =============================================================================
# CAP-mediation -- CAP-based mediation analysis with a covariance/graph mediator
# -----------------------------------------------------------------------------
# Wraps the RcppArmadillo-accelerated implementation in ../CAP-mediation/V4/.
# The MEDIATOR is a subject-level covariance matrix M_i (e.g. brain connectivity);
# a projection theta summarizes it as a scalar log(theta' Sigma_i theta), and a
# hierarchical model links exposure -> mediator -> outcome:
#   mediator:  log(theta' Sigma_i theta) = alpha0 + x_i' alpha + b_i
#   outcome:   Y_i = gamma0 + x_i' gamma + beta * log(theta' Sigma_i theta) + e_i
#   indirect (mediation) effect:  IE = alpha[exposure] * beta
# =============================================================================

# ---- 1. Load the method code from the CAPmethods package --------------------
# Uses the package's "mediation" environment (CAPMediation(), CAPMediation_boot(),
# ...); the med_-prefixed kernel ships in the package DLL. See R/pkg_methods.R.
.med_env <- cap_pkg_env("mediation")

# ---- 2. Built-in example (the manuscript's single-treatment, p=10 setting) ----
# n subjects, single binary treatment (no extra covariates: q=0), p=10 mediator
# dimension. A treatment effect on the mediator (alpha) and a mediator effect on
# the outcome (beta) create a true indirect effect IE = alpha * beta.
med_example <- function() {
  set.seed(2024)
  n <- 100L; p <- 10L; Ti <- 150L
  alpha0 <- 0.2; alpha_x <- 0.8; beta <- 0.7; gamma0 <- 0.1; gamma_x <- 0.4
  x <- rep(c(0, 1), length.out = n)
  X <- matrix(x, ncol = 1, dimnames = list(NULL, "treatment"))
  Phi <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
  if (Phi[which.max(abs(Phi[, 1])), 1] < 0) Phi[, 1] <- -Phi[, 1]
  base <- c(NA, exp(seq(log(1.2), log(0.3), length.out = p - 1)))   # background eigenvalues
  M <- vector("list", n); Y <- numeric(n)
  for (i in 1:n) {
    lv <- alpha0 + alpha_x * x[i] + stats::rnorm(1, 0, 0.3)          # log mediator-variance along theta
    ev <- base; ev[1] <- exp(lv)
    S <- Phi %*% diag(ev) %*% t(Phi); S <- (S + t(S)) / 2
    eg <- eigen(S, symmetric = TRUE)
    rt <- eg$vectors %*% diag(sqrt(pmax(eg$values, 1e-8))) %*% t(eg$vectors)
    M[[i]] <- matrix(stats::rnorm(Ti * p), Ti, p) %*% rt
    Y[i] <- gamma0 + gamma_x * x[i] + beta * lv + stats::rnorm(1, 0, 0.3)
  }
  ids <- paste0("S", 1:n); names(M) <- ids
  d <- list(M = M, X = X, Y = Y, ids = ids,
            m_names = paste0("V", 1:p), x_names = "treatment",
            Xdf = data.frame(id = ids, treatment = x),
            Ydf = data.frame(id = ids, Y = Y),
            truth = list(theta = Phi[, 1], alpha = alpha_x, beta = beta,
                         gamma = gamma_x, IE = alpha_x * beta))
  d$preview_ui <- med_preview(d)
  d
}

# ---- 3. Parser for uploaded data --------------------------------------------
# M: mediator data (native list of T_i x p matrices, or long CSV id + p columns).
# X: n x nX exposure/covariate matrix (FIRST column = exposure of interest).
# Y: scalar outcome per subject (numeric vector, or CSV id + value).
med_parse <- function(files, opts) {
  mp <- coerce_Y_input(files$M)                       # list of T_i x p matrices
  X  <- coerce_X_input(files$X, ids = mp$ids, add_intercept = FALSE)  # no intercept (added internally)
  Y  <- med_coerce_Y(files$Y, ids = mp$ids)
  d <- list(M = mp$Y, X = X, Y = Y, ids = mp$ids,
            m_names = mp$var_names, x_names = colnames(X),
            Xdf = data.frame(id = mp$ids, X, check.names = FALSE),
            Ydf = data.frame(id = mp$ids, Y = Y), truth = NULL)
  d$preview_ui <- med_preview(d)
  d
}

# coerce the outcome into a numeric vector aligned to subject ids
med_coerce_Y <- function(obj, ids) {
  if (is.data.frame(obj)) {
    if (ncol(obj) >= 2 && all(as.character(ids) %in% as.character(obj[[1]]))) {
      v <- obj[[2]][match(as.character(ids), as.character(obj[[1]]))]
    } else v <- obj[[ncol(obj)]]
    y <- suppressWarnings(as.numeric(v))
  } else if (is.list(obj) && !is.null(names(obj))) {
    key <- intersect(c("Y", "y", "outcome"), names(obj))
    y <- as.numeric(obj[[if (length(key)) key[1] else 1]])
  } else y <- as.numeric(obj)
  if (length(y) != length(ids) || anyNA(y))
    stop("Outcome Y must give one numeric value per subject (n = ", length(ids), ").")
  y
}

# ---- 4. Preview / helpers ---------------------------------------------------
med_preview <- function(d) {
  Tv <- vapply(d$M, nrow, integer(1)); p <- ncol(d$M[[1]])
  exposure <- d$x_names[1]; covs <- d$x_names[-1]
  rng <- function(x) if (length(unique(x)) == 1) as.character(x[1]) else paste0(min(x), "–", max(x))
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", length(d$M), theme = "primary"),
      bslib::value_box("Mediator dim (p)", p, theme = "secondary"),
      bslib::value_box("Samples / subject", rng(Tv), theme = "secondary"),
      bslib::value_box("Exposure", exposure, theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           sprintf("Mediator = subject covariance matrix; exposure = %s%s.%s",
                   exposure, if (length(covs)) paste0("; covariates: ", paste(covs, collapse = ", ")) else "",
                   if (!is.null(d$truth)) "  Simulated data with a known indirect effect." else "")))
}

med_describe <- function(d)
  sprintf("%d subjects; mediator dim p = %d; exposure: %s.",
          length(d$M), ncol(d$M[[1]]), d$x_names[1])

med_dfd_vec <- function(dfd, nD) {
  v <- if (is.list(dfd)) as.numeric(dfd$avg.level) else as.numeric(dfd)
  if (length(v) < nD) v <- c(v, rep(NA, nD - length(v)))
  v[seq_len(nD)]
}

# ---- 5. Run -----------------------------------------------------------------
med_run <- function(d, params) {
  X <- d$X; M <- d$M; Y <- d$Y
  stop_crt <- params$stop.crt %||% "nD"
  nD <- if (stop_crt == "nD") (params$nD %||% 1) else NULL
  ninitial <- params$ninitial %||% 3
  fit <- .med_env$CAPMediation(X, M, Y, stop.crt = stop_crt, nD = nD,
            DfD.thred = params$DfD.thred %||% 2, ninitial = ninitial,
            seed = 100, score.return = TRUE, verbose = FALSE)
  sims <- params$sims %||% 0
  inference <- NULL
  if (sims > 0) {
    inference <- lapply(seq_len(ncol(fit$theta)), function(j)
      tryCatch(.med_env$CAPMediation_boot(X, M, Y, theta = fit$theta[, j],
                 boot = TRUE, sims = sims, boot.ci.type = "se", verbose = FALSE),
               error = function(e) NULL))
  }
  list(fit = fit, inference = inference, X = X, Y = Y, ids = d$ids,
       m_names = d$m_names, x_names = d$x_names, truth = d$truth, params = params)
}

# ---- 6. Tables --------------------------------------------------------------
med_summarize <- function(res) {
  fit <- res$fit; nD <- ncol(fit$theta); dirs <- paste0("C", seq_len(nD)); out <- list()
  eff_labels <- c(alpha = "α  exposure → mediator", beta = "β  mediator → outcome",
                  gamma = "γ  direct effect", IE = "IE  indirect (α·β)", DE = "DE  direct effect")

  # mediation effects (bootstrap inference if available, else point estimates)
  eff_rows <- do.call(rbind, lapply(seq_len(nD), function(j) {
    bi <- if (!is.null(res$inference)) res$inference[[j]] else NULL
    if (!is.null(bi)) {
      df <- as.data.frame(bi$coef)
      data.frame(Direction = dirs[j], Effect = eff_labels[rownames(bi$coef)], df,
                 check.names = FALSE, row.names = NULL)
    } else {
      cf <- fit$coef[, j]
      data.frame(Direction = dirs[j], Effect = eff_labels[rownames(fit$coef)],
                 Estimate = cf, SE = NA, statistics = NA, pvalue = NA, LB = NA, UB = NA,
                 check.names = FALSE, row.names = NULL)
    }
  }))
  out[["Mediation effects"]] <- eff_rows

  # other model coefficients (intercepts + covariate effects), if bootstrapped
  if (!is.null(res$inference) && !is.null(res$inference[[1]]$coef.other)) {
    oth <- do.call(rbind, lapply(seq_len(nD), function(j) {
      bi <- res$inference[[j]]; if (is.null(bi)) return(NULL)
      data.frame(Direction = dirs[j], Term = rownames(bi$coef.other),
                 as.data.frame(bi$coef.other), check.names = FALSE, row.names = NULL)
    }))
    if (!is.null(oth)) out[["Other coefficients"]] <- oth
  }

  # theta (mediator projection) loadings
  th <- as.data.frame(fit$theta); names(th) <- dirs
  out[["θ mediator loadings"]] <- cbind(Variable = res$m_names, th)

  # DfD + orthogonality
  out[["DfD across dimensions"]] <- data.frame(Dimension = dirs,
    DfD = med_dfd_vec(fit$DfD, nD), check.names = FALSE)
  if (nD > 1 && !is.null(fit$orthogonality)) {
    o <- as.data.frame(fit$orthogonality); names(o) <- dirs
    out[["Orthogonality θ'θ"]] <- cbind(" " = dirs, o)
  }

  # truth comparison (example data) -- match leading direction to true theta
  if (!is.null(res$truth)) {
    j <- which.max(vapply(seq_len(nD), function(k) abs(sum(fit$theta[, k] * res$truth$theta)), numeric(1)))
    cf <- fit$coef[, j]
    out[["Truth vs estimate"]] <- data.frame(
      Quantity = c("θ cosine", "α (exp→med)", "β (med→out)", "γ (direct)", "IE (α·β)"),
      True = c(1, res$truth$alpha, res$truth$beta, res$truth$gamma, res$truth$IE),
      Estimate = round(c(abs(sum(fit$theta[, j] * res$truth$theta)),
                         cf["alpha"], cf["beta"], cf["gamma"], cf["IE"]), 3),
      check.names = FALSE)
  }
  out
}

# ---- 7. Plots ---------------------------------------------------------------
med_plots <- function(res) {
  fit <- res$fit; nD <- ncol(fit$theta); dirs <- paste0("C", seq_len(nD)); plots <- list()

  # (a) theta loadings
  th <- as.data.frame(fit$theta); names(th) <- dirs
  long <- do.call(rbind, lapply(seq_len(nD), function(j)
    data.frame(Variable = res$m_names, Direction = dirs[j], Loading = th[, j])))
  long$Variable <- factor(long$Variable, levels = res$m_names)
  plots[["θ mediator loadings"]] <- list(plot = plotly::plot_ly(
    long, x = ~Variable, y = ~Loading, color = ~Direction, type = "bar") |>
      plotly::layout(title = "Mediator projection loadings (θ)", barmode = "group"))

  # (b) mediation effects with bootstrap 95% CI (leading direction)
  if (!is.null(res$inference) && !is.null(res$inference[[1]])) {
    bi <- res$inference[[1]]$coef
    keep <- intersect(c("alpha", "beta", "gamma", "IE"), rownames(bi))
    edf <- data.frame(Effect = factor(keep, levels = keep),
                      Estimate = bi[keep, "Estimate"], LB = bi[keep, "LB"], UB = bi[keep, "UB"])
    p_e <- plotly::plot_ly(edf, x = ~Effect, y = ~Estimate, type = "bar",
                           marker = list(color = "#2c7fb8"),
                           error_y = list(type = "data", symmetric = FALSE,
                                          array = ~(UB - Estimate), arrayminus = ~(Estimate - LB))) |>
      plotly::layout(title = "Mediation effects (D1) with bootstrap 95% CI",
                     yaxis = list(title = "estimate"),
                     shapes = list(list(type = "line", x0 = -0.5, x1 = length(keep) - 0.5,
                                        y0 = 0, y1 = 0, line = list(dash = "dot", color = "grey"))))
    plots[["Mediation effects"]] <- list(plot = p_e, data = edf)
  }

  # (c) mediator score vs outcome, colored by exposure (the mediator -> outcome path)
  if (!is.null(fit$score)) {
    sc <- log(pmax(as.matrix(fit$score)[, 1], 1e-12))
    expo <- res$X[, 1]
    sdf <- data.frame(id = res$ids, med_score = sc, Y = res$Y,
                      exposure = factor(expo))
    p_s <- plotly::plot_ly(sdf, x = ~med_score, y = ~Y, color = ~exposure,
                           type = "scatter", mode = "markers",
                           marker = list(size = 8, opacity = 0.7),
                           text = ~id, hovertemplate = "%{text}<extra></extra>") |>
      plotly::layout(title = "Outcome vs mediator score (log θ'Σθ), by exposure",
                     xaxis = list(title = "mediator score  log(θ'Σθ)"),
                     yaxis = list(title = "outcome Y"))
    plots[["Mediator → outcome"]] <- list(plot = p_s, data = sdf)
  }
  plots
}

# ---- 8. Register ------------------------------------------------------------
register_method(list(
  id = "cap-mediation",
  name = "CAP-mediation",
  full_name = "CAP-based Mediation Analysis with a Covariance Mediator",
  short = "Mediation analysis where the mediator is a subject-level covariance matrix; a learned projection θ summarizes it, giving exposure→mediator→outcome effects.",
  status = "ready",
  tags = c("covariance regression", "mediation", "graph mediator", "bootstrap"),
  paper = list(
    citation = "Xu, Y., & Zhao, Y. (2025). Mediation analysis with graph mediator. Biostatistics, 26(1), kxaf004.",
    url = "https://doi.org/10.1093/biostatistics/kxaf004"),
  explain = file.path(APP_DIR, "methods", "cap-mediation", "explain.md"),
  example_note = paste("The manuscript's single-treatment setting: 100 subjects,",
                       "p=10 mediator dimension, Tᵢ=150, a binary treatment with a",
                       "true indirect effect IE = α·β = 0.8·0.7 = 0.56."),
  data_inputs = list(
    list(id = "M", label = "Mediator (list of T×p covariance data)",
         help = "An .rds/.RData holding a list of length n; element i is a T_i × p matrix whose covariance is subject i's mediator. A long CSV (subject-id column + p columns) is also accepted."),
    list(id = "X", label = "Exposure / covariates (n × nX)",
         help = "An n × nX matrix in subject order (the FIRST column is the exposure of interest; any further columns are adjusted covariates), or a CSV whose first column is the subject id. No intercept column."),
    list(id = "Y", label = "Outcome (one value per subject)",
         help = "A numeric outcome per subject: a CSV with subject-id + value columns, or an .rds numeric vector in subject order.")
  ),
  params = list(
    list(id = "stop.crt", label = "Number of directions chosen by", type = "select",
         choices = c("Fixed number (nD)" = "nD", "DfD criterion" = "DfD"), default = "nD",
         help = "Fix the count, or keep adding mediator components while DfD stays below the threshold."),
    list(id = "nD", label = "Number of mediator components (nD)", type = "integer",
         default = 1, min = 1, max = 5, help = "Used when selecting by a fixed number. Each adds a bootstrap pass."),
    list(id = "DfD.thred", label = "DfD threshold", type = "numeric",
         default = 2, min = 1, max = 100, step = 0.5, help = "Used when selecting by the DfD criterion."),
    list(id = "ninitial", label = "# random initializations", type = "integer",
         default = 3, min = 1, max = 20, help = "Random restarts for the θ search."),
    list(id = "sims", label = "Bootstrap replicates (effect inference)", type = "integer",
         default = 200, min = 0, max = 2000, help = "Subject-level bootstrap for α/β/IE inference; 0 = skip.")
  ),
  example = med_example,
  export_example = function(d) list(M = d$M, X = d$Xdf, Y = d$Ydf),
  parse = med_parse,
  describe_data = med_describe,
  run = med_run,
  summarize = med_summarize,
  plots = med_plots
))
