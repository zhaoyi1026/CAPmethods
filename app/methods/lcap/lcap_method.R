# =============================================================================
# LCAP -- Longitudinal Covariate Assisted Principal regression (gamma-invariant)
# -----------------------------------------------------------------------------
# Wraps the implementation in ../LCAP_gamma-invar/V5/ and registers it with the
# website. Unlike CAP/HDCAP (cross-sectional), LCAP is LONGITUDINAL: each subject
# is observed over several VISITS, and the model uses a time-invariant projection
# gamma with subject-level RANDOM EFFECTS (random intercept + random slopes) on
# the log projected-variance. The gamma-varying multilevel sibling is "MCAP".
# =============================================================================

# ---- 1. Load the method code from the CAPmethods package --------------------
# Uses the package's "lcap" environment (capReg(), cap_beta_boot(), ...); the
# lcap_-prefixed kernel ships in the package DLL. See R/pkg_methods.R.
.lcap_env <- cap_pkg_env("lcap")

# ---- 2. Built-in example data = the manuscript LCAP simulation ---------------
# Faithful to the paper's gamma-invariant simulation (p20_q3, case 1): p = 20,
# a common time-invariant eigenbasis Gamma, n = 100 subjects each seen at ~5
# visits with ~100 samples/visit, two within-subject covariates (x1, x2), and a
# subject-level random intercept. TWO directions satisfy the CAP model -- basis
# cols 2 and 4 -- whose log-variance depends on the covariates (effects
# x1: -0.5/+0.5, x2: +0.5/-0.25 on dims 2/4); the rest carry only the random
# baseline. With shrinkage on, the default nD = 2 recovers BOTH covariate-driven
# directions. NB: the app's LCAP fits within-subject random slopes, so the
# covariates are time-varying here (the manuscript's binary time-invariant x1 is
# rendered as a within-subject covariate); the fixed beta and directions match.
lcap_example <- function() {
  p <- 20L; n <- 100L; nV.m <- 5; nT.m <- 100; seedl <- 4L

  # common orthonormal eigenbasis (runif(p) recycled into p x p, as in the paper)
  set.seed(100)
  Gamma <- qr.Q(qr(matrix(stats::runif(p), p, p)))
  for (j in 1:p) if (Gamma[which.max(abs(Gamma[, j])), j] < 0) Gamma[, j] <- -Gamma[, j]

  # per-direction coefficients: row 1 = baseline log-variance (random-intercept
  # mean), rows 2,3 = fixed effects of x1, x2 (nonzero only on dims 2 and 4).
  beta.mat <- rbind(c(seq(3, -1, length.out = 5), seq(-1.5, -3, length.out = p - 5)),
                    c(c(0, -0.5, 0, 0.5, 0), rep(0, p - 5)),
                    c(c(0,  0.5, 0, -0.25, 0), rep(0, p - 5)))
  cov_dirs <- which(beta.mat[2, ] != 0)               # = c(2, 4)

  # visit counts and per-visit sample sizes (manuscript seeds)
  set.seed(100); nV <- round(stats::rnorm(n, nV.m, sd = 1))
  set.seed(100); nT <- matrix(NA_integer_, n, max(nV))
  for (i in 1:n) nT[i, 1:nV[i]] <- round(stats::rnorm(nV[i], mean = nT.m, sd = 5))

  # subject random intercepts (same standard-normal draw per dimension, sd 0.1)
  beta0.mat <- matrix(NA, n, ncol(beta.mat))
  for (j in 1:ncol(beta.mat)) { set.seed(100); beta0.mat[, j] <- stats::rnorm(n, beta.mat[1, j], sd = 0.1) }

  set.seed(as.numeric("20201210") + seedl * 100)
  ids <- paste0("S", 1:n)
  Y <- X <- vector("list", n); names(Y) <- names(X) <- ids
  for (i in 1:n) {
    Xi <- cbind(Intercept = 1, x1 = stats::rnorm(nV[i], 0, 0.5), x2 = stats::rnorm(nV[i], 0, 0.5))
    Yi <- vector("list", nV[i])
    for (v in 1:nV[i]) {
      delta <- vapply(1:p, function(j)
        exp(sum(Xi[v, ] * c(beta0.mat[i, j], beta.mat[-1, j]))), numeric(1))
      Sig <- Gamma %*% diag(delta) %*% t(Gamma)
      Yi[[v]] <- mvtnorm::rmvnorm(nT[i, v], rep(0, p), sigma = Sig)
    }
    Y[[i]] <- Yi; X[[i]] <- Xi
  }
  var_names <- paste0("V", 1:p)

  gamma_true <- Gamma[, cov_dirs, drop = FALSE]; colnames(gamma_true) <- paste0("D", seq_along(cov_dirs))
  beta_true  <- beta.mat[, cov_dirs, drop = FALSE]
  rownames(beta_true) <- c("Intercept", "x1", "x2"); colnames(beta_true) <- paste0("D", seq_along(cov_dirs))

  d <- list(Y = Y, X = X, ids = ids, var_names = var_names,
            cov_names = c("x1", "x2"),
            truth = list(gamma = gamma_true, beta = beta_true))
  d$preview_ui <- lcap_preview(d)
  d
}

# ---- 3. Parsers for uploaded data -------------------------------------------
# Y: a native .rds nested list (subject -> visit -> T_iv x p matrix), OR a long
#    CSV whose first two columns are subject id and visit id, then p responses.
# X: a native .rds list (subject -> nV x q matrix), OR a long CSV whose first two
#    columns are subject id and visit id, then the covariates (one row / visit).
lcap_parse <- function(files, opts) {
  yp <- lcap_coerce_Y(files$Y)
  xp <- lcap_coerce_X(files$X, ids = yp$ids, nVvec = yp$nVvec,
                      add_intercept = isTRUE(opts$add_intercept))
  d <- list(Y = yp$Y, X = xp$X, ids = yp$ids, var_names = yp$var_names,
            cov_names = xp$cov_names, truth = NULL)
  d$preview_ui <- lcap_preview(d)
  d
}

# nested list Y[[i]][[v]] = T_iv x p  (validate or build from a long data.frame)
lcap_coerce_Y <- function(obj) {
  if (is.data.frame(obj)) return(lcap_Y_from_long(obj))

  # unwrap a named .RData container
  is_nested <- is.list(obj) && length(obj) > 0 && is.list(obj[[1]]) &&
    !is.matrix(obj[[1]]) && !is.data.frame(obj[[1]])
  if (is.list(obj) && !is.null(names(obj)) && !is_nested) {
    key <- intersect(c("Y", "Y_list", "Ylist", "response"), names(obj))
    if (length(key)) obj <- obj[[key[1]]]
    is_nested <- is.list(obj) && length(obj) > 0 && is.list(obj[[1]]) &&
      !is.matrix(obj[[1]]) && !is.data.frame(obj[[1]])
  }
  if (!is_nested)
    stop("Response must be a nested list (subject -> visit -> T x p matrix), or a long CSV with id + visit columns.")

  Y <- lapply(obj, function(subj) lapply(subj, function(m) {
    m <- as.matrix(m); storage.mode(m) <- "double"; m
  }))
  ps <- unlist(lapply(Y, function(subj) vapply(subj, ncol, integer(1))))
  if (length(unique(ps)) != 1L)
    stop("All visit matrices must have the same number of columns (p).")
  ids <- names(Y); if (is.null(ids) || any(ids == "")) ids <- paste0("S", seq_along(Y))
  names(Y) <- ids
  nVvec <- vapply(Y, length, integer(1))
  vn <- colnames(Y[[1]][[1]]); if (is.null(vn)) vn <- paste0("V", seq_len(ps[1]))
  list(Y = Y, ids = ids, nVvec = nVvec, var_names = vn)
}

lcap_Y_from_long <- function(df) {
  id <- as.character(df[[1]]); visit <- as.character(df[[2]])
  resp <- df[, -(1:2), drop = FALSE]
  resp[] <- lapply(resp, function(x) suppressWarnings(as.numeric(x)))
  if (any(vapply(resp, function(x) all(is.na(x)), logical(1))))
    stop("One or more response columns are non-numeric.")
  ids <- unique(id)
  Y <- lapply(ids, function(i) {
    vs <- unique(visit[id == i])
    lapply(vs, function(v) as.matrix(resp[id == i & visit == v, , drop = FALSE]))
  })
  names(Y) <- ids
  nVvec <- vapply(Y, length, integer(1))
  list(Y = Y, ids = ids, nVvec = nVvec, var_names = names(resp))
}

# list X[[i]] = nV x q  (validate or build from a long data.frame)
lcap_coerce_X <- function(obj, ids, nVvec, add_intercept = TRUE) {
  if (is.list(obj) && !is.data.frame(obj) && !is.matrix(obj) &&
      length(obj) > 0 && (is.matrix(obj[[1]]) || is.data.frame(obj[[1]]))) {
    X <- lapply(obj, function(m) { m <- as.matrix(m); storage.mode(m) <- "double"; m })
  } else if (is.list(obj) && !is.data.frame(obj) && !is.matrix(obj)) {
    key <- intersect(c("X", "covariates", "Xcov", "Z"), names(obj))
    if (!length(key)) stop("Could not find a covariate object (X) in the upload.")
    return(lcap_coerce_X(obj[[key[1]]], ids, nVvec, add_intercept))
  } else if (is.data.frame(obj)) {
    X <- lcap_X_from_long(obj, ids)
  } else {
    stop("Covariates must be a list of n matrices (one nV x q per subject), or a long CSV with id + visit columns.")
  }
  if (length(X) != length(ids))
    stop(sprintf("Covariates have %d subjects but the response has %d.", length(X), length(ids)))
  if (!all(vapply(seq_along(X), function(i) nrow(X[[i]]) == nVvec[i], logical(1))))
    stop("Each subject's covariate rows must equal that subject's number of visits.")
  has_int <- !is.null(colnames(X[[1]])) && "Intercept" %in% colnames(X[[1]])
  if (add_intercept && !has_int)
    X <- lapply(X, function(m) cbind(Intercept = 1, m))
  cov_names <- setdiff(colnames(X[[1]]) %||% paste0("X", seq_len(ncol(X[[1]]))), "Intercept")
  names(X) <- ids
  list(X = X, cov_names = cov_names)
}

lcap_X_from_long <- function(df, ids) {
  id <- as.character(df[[1]]); visit <- as.character(df[[2]])
  cov <- df[, -(1:2), drop = FALSE]
  cov[] <- lapply(cov, function(x) suppressWarnings(as.numeric(x)))
  if (anyNA(unlist(cov))) stop("Covariates contain non-numeric or missing values.")
  lapply(ids, function(i) {
    vs <- unique(visit[id == i])
    m <- as.matrix(cov[id == i, , drop = FALSE])
    rownames(m) <- NULL; m
  })
}

# ---- 4. Preview / helpers ---------------------------------------------------
lcap_preview <- function(d) {
  nVvec <- vapply(d$Y, length, integer(1))
  Tvec  <- unlist(lapply(d$Y, function(s) vapply(s, nrow, integer(1))))
  p <- ncol(d$Y[[1]][[1]])
  hd_flag <- min(Tvec) - 5 < p
  rng <- function(x) if (length(unique(x)) == 1) as.character(x[1])
                     else paste0(min(x), "–", max(x))
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", length(d$Y), theme = "primary"),
      bslib::value_box("Dimension (p)", p, theme = "secondary"),
      bslib::value_box("Visits / subject", rng(nVvec), theme = "secondary"),
      bslib::value_box("Samples / visit", rng(Tvec), theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           sprintf("Time-varying covariates: %s.%s%s",
                   paste(d$cov_names, collapse = ", "),
                   if (hd_flag) "  High-dimensional (min Tᵢᵥ − 5 < p): shrinkage auto-enabled." else "",
                   if (!is.null(d$truth)) "  Simulated data with known truth." else ""))
  )
}

lcap_describe <- function(d) {
  nVvec <- vapply(d$Y, length, integer(1))
  sprintf("%d subjects, %s visits each, dimension p = %d, covariates: %s.",
          length(d$Y),
          if (length(unique(nVvec)) == 1) nVvec[1] else paste0(min(nVvec), "-", max(nVvec)),
          ncol(d$Y[[1]][[1]]), paste(d$cov_names, collapse = ", "))
}

# DfD avg.level vector (returns a list of >1 directions, else scalar 1).
lcap_dfd_vec <- function(fit) {
  if (is.list(fit$DfD)) as.numeric(fit$DfD$avg.level) else rep(1, ncol(fit$gamma))
}

# Greedy sign-aligned matching of estimated to true directions.
lcap_match_dirs <- function(Ghat, Gtrue) {
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
lcap_run <- function(d, params) {
  Y <- d$Y; X <- d$X
  stop_crt <- params$stop.crt %||% "nD"
  shrink <- isTRUE(params$cov.shrinkage %||% TRUE)
  fit <- .lcap_env$capReg(Y, X,
              stop.crt = stop_crt,
              nD = if (stop_crt == "nD") (params$nD %||% 2) else NULL,
              DfD.thred = params$DfD.thred %||% 2,
              method = "CAP",
              cov.shrinkage = shrink,
              ninitial = if (is.null(params$ninitial) || params$ninitial <= 0) NULL
                         else params$ninitial,
              score.return = TRUE, verbose = FALSE)

  # bootstrap inference for the FIXED effects beta, per estimated direction
  sims <- params$sims %||% 0
  inference <- NULL
  if (sims > 0) {
    inference <- lapply(seq_len(ncol(fit$gamma)), function(j)
      tryCatch(
        .lcap_env$cap_beta_boot(Y, X, gamma = fit$gamma[, j],
                                cov.shrinkage = shrink, sims = sims,
                                verbose = FALSE)$Inference,
        error = function(e) NULL))
    names(inference) <- colnames(fit$gamma)
  }

  list(fit = fit, inference = inference, X = X, Y = Y, ids = d$ids,
       var_names = d$var_names, cov_names = d$cov_names,
       truth = d$truth, params = params)
}

# ---- 6. Tables --------------------------------------------------------------
lcap_summarize <- function(res) {
  fit <- res$fit
  dirs <- colnames(fit$gamma)
  terms <- rownames(fit$beta)
  out <- list()

  # fixed-effects beta (bootstrap inference if available, else estimates)
  if (!is.null(res$inference)) {
    beta_rows <- do.call(rbind, lapply(seq_along(dirs), function(j) {
      inf <- res$inference[[j]]
      if (is.null(inf) || is.null(inf$beta)) {
        data.frame(Direction = dirs[j], Term = terms, Estimate = fit$beta[, j],
                   SE = NA, statistic = NA, pvalue = NA, LB = NA, UB = NA,
                   check.names = FALSE)
      } else {
        data.frame(Direction = dirs[j], Term = rownames(inf$beta), inf$beta,
                   check.names = FALSE, row.names = NULL)
      }
    }))
    out[["β fixed effects (bootstrap)"]] <- beta_rows
  } else {
    bdf <- as.data.frame(fit$beta)
    out[["β fixed effects"]] <- cbind(Term = terms, bdf)
  }

  # random-effects variance components: random intercept sigma^2 + slope Omega diag
  re_rows <- do.call(rbind, lapply(seq_along(dirs), function(j) {
    s2 <- as.numeric(fit$beta0.sigma2[, j])
    Om <- as.matrix(fit$beta1.Omega[[j]])
    data.frame(Direction = dirs[j],
               Component = c("Intercept (random)", paste0(res$cov_names, " (random slope)")),
               Variance = c(s2, diag(Om)), check.names = FALSE)
  }))
  out[["Random-effect variances"]] <- re_rows

  # gamma loadings
  g <- as.data.frame(fit$gamma)
  out[["γ (loadings)"]] <- cbind(Variable = res$var_names %||% rownames(fit$gamma), g)

  # shrinkage weights per direction (only when shrinkage was used). fit$shrinkage
  # is a data.frame of list-columns (one cell per direction), so unlist per column.
  if (!is.null(fit$shrinkage)) {
    sk <- fit$shrinkage
    skdf <- as.data.frame(lapply(sk, function(col) as.numeric(unlist(col))))
    names(skdf) <- colnames(sk)
    out[["Shrinkage weights"]] <- cbind(Direction = rownames(sk) %||% dirs, skdf)
  }

  # orthogonality
  if (ncol(fit$gamma) > 1) {
    o <- as.data.frame(fit$orthogonality)
    out[["Orthogonality γ'γ"]] <- cbind(" " = rownames(fit$orthogonality), o)
  }

  # DfD across dimensions
  dfd_vec <- lcap_dfd_vec(fit)
  out[["DfD across dimensions"]] <- data.frame(
    Dimension = dirs, `# directions (1..k)` = seq_along(dfd_vec),
    DfD = dfd_vec, check.names = FALSE)

  # truth comparison (example data)
  if (!is.null(res$truth)) {
    Gtrue <- as.matrix(res$truth$gamma)
    if (is.null(colnames(Gtrue))) colnames(Gtrue) <- paste0("D", seq_len(ncol(Gtrue)))
    mt <- lcap_match_dirs(fit$gamma, Gtrue)
    cmp <- data.frame(Variable = res$var_names)
    for (r in seq_len(nrow(mt))) {
      cmp[[paste0(mt$true[r], "_true")]] <- Gtrue[, mt$k[r]]
      cmp[[paste0(mt$est[r], "_est")]]  <- fit$gamma[, mt$est[r]] * mt$sign[r]
    }
    out[["γ: truth vs estimate"]] <- cmp
    out[["Direction recovery"]] <- data.frame(
      Estimated = mt$est, `Matched true` = mt$true,
      `Cosine similarity` = round(mt$cosine, 4), check.names = FALSE)

    # fixed-effects truth vs estimate for matched directions
    # (beta is invariant to gamma's sign, so no flip)
    if (!is.null(res$truth$beta)) {
      Btrue <- as.matrix(res$truth$beta)
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
lcap_plots <- function(res) {
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
      plotly::layout(title = "Time-invariant projection loadings (γ)",
                     barmode = "group"))

  # (b) subject-visit scores vs fitted (flatten the n x maxV x nD score array)
  if (!is.null(fit$score)) {
    X <- res$X; ids <- res$ids %||% names(X) %||% seq_along(X)
    sv <- do.call(rbind, lapply(seq_along(dirs), function(j) {
      do.call(rbind, lapply(seq_along(X), function(i) {
        nv <- nrow(X[[i]])
        sc <- fit$score[i, seq_len(nv), j]
        ft <- as.numeric(X[[i]] %*% fit$beta[, j])
        data.frame(id = ids[i], visit = seq_len(nv), Direction = dirs[j],
                   score = sc, fitted = ft, logvar = log(pmax(sc, 1e-12)))
      }))
    }))
    rng <- range(c(sv$fitted, sv$logvar), finite = TRUE)
    p_sc <- plotly::plot_ly(sv, x = ~fitted, y = ~logvar, color = ~Direction,
                            type = "scatter", mode = "markers",
                            marker = list(size = 7, opacity = 0.6),
                            text = ~paste0(id, " v", visit),
                            hovertemplate = "%{text}<extra></extra>") |>
      plotly::add_lines(x = rng, y = rng, line = list(dash = "dot", color = "grey"),
                        showlegend = FALSE, inherit = FALSE) |>
      plotly::layout(title = "Subject-visit projected variance vs fitted",
                     xaxis = list(title = "x'β (fitted log-variance)"),
                     yaxis = list(title = "log(γ' Σ γ)"))
    plots[["Subject-visit scores"]] <- list(plot = p_sc, data = sv)
  }

  # (c) DfD across dimensions
  dfd_vec <- lcap_dfd_vec(fit)
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

# Long-format example tables (downloadable templates show the expected upload shape)
lcap_export_example <- function(d) {
  Ylong <- do.call(rbind, lapply(seq_along(d$Y), function(i) {
    do.call(rbind, lapply(seq_along(d$Y[[i]]), function(v) {
      m <- as.data.frame(d$Y[[i]][[v]]); names(m) <- d$var_names
      cbind(id = d$ids[i], visit = v, m)
    }))
  }))
  Xlong <- do.call(rbind, lapply(seq_along(d$X), function(i) {
    m <- as.data.frame(d$X[[i]][, setdiff(colnames(d$X[[i]]), "Intercept"), drop = FALSE])
    cbind(id = d$ids[i], visit = seq_len(nrow(m)), m)
  }))
  rownames(Ylong) <- rownames(Xlong) <- NULL
  list(Y = Ylong, X = Xlong)
}

# ---- 8. Register ------------------------------------------------------------
register_method(list(
  id = "lcap-invar",
  name = "LCAP",
  full_name = "Longitudinal Covariate Assisted Principal Regression (γ-invariant)",
  short = "CAP for longitudinal covariance outcomes: a time-invariant projection with subject-level random effects.",
  status = "ready",
  tags = c("covariance regression", "longitudinal", "random effects", "shrinkage"),
  paper = list(
    citation = "Zhao, Y., Caffo, B. S., & Luo, X. (2024). Longitudinal regression of covariance matrix outcomes. Biostatistics, 25(2), 385-401.",
    url = "https://doi.org/10.1093/biostatistics/kxac045"),
  explain = file.path(APP_DIR, "methods", "lcap", "explain.md"),
  example_note = paste("The manuscript's longitudinal simulation (p20_q3): 100 subjects,",
                       "~5 visits each, p=20, Tᵢᵥ≈100, two within-subject covariates and a",
                       "subject random intercept. TWO directions satisfy the CAP model",
                       "(basis cols 2,4). With shrinkage on, nD=2 recovers both (γ cosine",
                       "≈0.98/0.94); ~1 min."),
  x_intercept_option = TRUE,
  data_inputs = list(
    list(id = "Y", label = "Response (nested list or long CSV)",
         help = "An .rds/.RData holding a nested list: element i is itself a list over visits, each a T_iv x p matrix. A long CSV with columns [id, visit, V1..Vp] is also accepted."),
    list(id = "X", label = "Covariates (list or long CSV)",
         help = "An .rds list: element i is an nV_i x q matrix (one row per visit). A long CSV with columns [id, visit, cov1..covq] is also accepted. Covariates should vary within subject across visits.")
  ),
  params = list(
    list(id = "cov.shrinkage", label = "Covariance shrinkage", type = "checkbox",
         default = TRUE,
         help = "Linear shrinkage of each subject-visit covariance (recommended; auto-enabled when min Tᵢᵥ − 5 < p)."),
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
         help = "Subject-level bootstrap for fixed-effect β inference; 0 = skip (estimates only).")
  ),
  example = lcap_example,
  export_example = lcap_export_example,
  parse = lcap_parse,
  describe_data = lcap_describe,
  run = lcap_run,
  summarize = lcap_summarize,
  plots = lcap_plots
))
