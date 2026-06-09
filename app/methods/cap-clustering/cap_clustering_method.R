# =============================================================================
# CAP-clustering (CAPclust) -- Parsimonious Clustering of Covariance Matrices
# -----------------------------------------------------------------------------
# Wraps the RcppArmadillo-accelerated implementation in ../CAP-clustering/V4/.
# Subjects whose data are covariance matrices are clustered, sharing a common CAP
# projection gamma. Within cluster k the projected log-variance follows a
# log-linear model, and cluster membership depends on covariates:
#   within cluster k:  log(gamma' Sigma_i gamma) = x_i' beta_k
#   membership:        P(cluster k | w_i) proportional to exp(w_i' alpha_k)
# Estimated by EM (gamma via a generalized eigenproblem; alpha via penalized
# multinomial regression brglm2::brmultinom).
# =============================================================================

# ---- 1. Load the method code from the CAPmethods package --------------------
# Uses the package's "cluster" environment (capPCL(), capPCL_coef_boot(), ...);
# the cluster_-prefixed kernel ships in the package DLL. See R/pkg_methods.R.
.cluster_env <- cap_pkg_env("cluster")

# ---- 2. Built-in example (the manuscript's p=50, K=2 simulation) -------------
# Common eigenvectors Pi; two covariate-driven directions (D2, D4) whose
# eigenvalues follow a cluster-specific log-linear model; background eigenvalues
# decay exp(N(seq(3,-1), 0.2)). D2 has the well-separated clusters (the leading
# component the app recovers with nD=1); D4's clusters differ only in covariate
# signs (harder). x = (1, Bernoulli, N(0,1)); membership w = (1, Bernoulli).
clust_example <- function() {
  set.seed(1)
  n <- 100L; p <- 50L; Ti <- 100L
  Pi <- qr.Q(qr(matrix(stats::rnorm(p * p), p, p)))
  x1 <- stats::rbinom(n, 1, 0.5); x2 <- stats::rnorm(n)
  X <- cbind(Intercept = 1, x1 = x1, x2 = x2)              # q1 = 3
  w1 <- stats::rbinom(n, 1, 0.5)
  W <- cbind(Intercept = 1, w1 = w1)                       # q2 = 2
  cD2 <- 1 + stats::rbinom(n, 1, plogis(W %*% c(0.5, -1)))
  cD4 <- 1 + stats::rbinom(n, 1, plogis(W %*% c(-0.25, 0.5)))
  bD2 <- cbind(c(1, 1, -1), c(-1, -1, 1))                  # q1 x K (D2: well separated)
  bD4 <- cbind(c(0.5, 0.5, -0.5), c(0.5, -0.5, 0.5))       # D4: same baseline
  logmean <- seq(3, -1, length.out = p)
  Y <- vector("list", n)
  for (i in 1:n) {
    lam <- exp(stats::rnorm(p, logmean, 0.2))
    lam[2] <- exp(sum(X[i, ] * bD2[, cD2[i]]))
    lam[4] <- exp(sum(X[i, ] * bD4[, cD4[i]]))
    S <- Pi %*% diag(lam) %*% t(Pi); S <- (S + t(S)) / 2
    eg <- eigen(S, symmetric = TRUE)
    rt <- eg$vectors %*% diag(sqrt(pmax(eg$values, 1e-8))) %*% t(eg$vectors)
    Y[[i]] <- matrix(stats::rnorm(Ti * p), Ti, p) %*% rt
  }
  ids <- paste0("S", 1:n); names(Y) <- ids
  d <- list(Y = Y, X = X, W = W, ids = ids,
            var_names = paste0("V", 1:p),
            x_names = colnames(X), w_names = colnames(W),
            Xdf = data.frame(id = ids, x1 = x1, x2 = x2),
            Wdf = data.frame(id = ids, w1 = w1),
            truth = list(gamma = cbind(D2 = Pi[, 2], D4 = Pi[, 4]),
                         cluster = cbind(D2 = cD2, D4 = cD4),
                         beta = list(D2 = bD2, D4 = bD4)))
  d$preview_ui <- clust_preview(d)
  d
}

# ---- 3. Parser for uploaded data --------------------------------------------
# Y: native list of T_i x p matrices, or long CSV (id + p columns).
# X: variance-model covariates (n x q1); W: membership covariates (n x q2).
#    Both as n-row matrices / CSVs (subject id first col); intercepts auto-added.
clust_parse <- function(files, opts) {
  yp <- coerce_Y_input(files$Y)
  X  <- coerce_X_input(files$X, ids = yp$ids, add_intercept = TRUE)
  W  <- coerce_X_input(files$W, ids = yp$ids, add_intercept = TRUE)
  d <- list(Y = yp$Y, X = X, W = W, ids = yp$ids, var_names = yp$var_names,
            x_names = colnames(X), w_names = colnames(W),
            Xdf = data.frame(id = yp$ids, X[, setdiff(colnames(X), "Intercept"), drop = FALSE], check.names = FALSE),
            Wdf = data.frame(id = yp$ids, W[, setdiff(colnames(W), "Intercept"), drop = FALSE], check.names = FALSE),
            truth = NULL)
  d$preview_ui <- clust_preview(d)
  d
}

# ---- 4. Preview / helpers ---------------------------------------------------
clust_preview <- function(d) {
  Tv <- vapply(d$Y, nrow, integer(1)); p <- ncol(d$Y[[1]])
  rng <- function(x) if (length(unique(x)) == 1) as.character(x[1]) else paste0(min(x), "–", max(x))
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", length(d$Y), theme = "primary"),
      bslib::value_box("Dimension (p)", p, theme = "secondary"),
      bslib::value_box("Samples / subject", rng(Tv), theme = "secondary"),
      bslib::value_box("Membership covs", length(setdiff(d$w_names, "Intercept")), theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           sprintf("Variance-model covariates: %s.%s",
                   paste(setdiff(d$x_names, "Intercept"), collapse = ", "),
                   if (!is.null(d$truth)) "  Simulated data with known clusters (component D2 is well separated)." else "")))
}

clust_describe <- function(d)
  sprintf("%d subjects; dimension p = %d; variance covariates: %s; membership covariates: %s.",
          length(d$Y), ncol(d$Y[[1]]),
          paste(setdiff(d$x_names, "Intercept"), collapse = ", "),
          paste(setdiff(d$w_names, "Intercept"), collapse = ", "))

# best label-permutation agreement between two clusterings (K small)
clust_agree <- function(est, truth) {
  ke <- sort(unique(est)); kt <- sort(unique(truth))
  if (length(ke) != length(kt)) return(mean(est == truth))
  perms <- function(v) if (length(v) == 1) list(v) else
    do.call(c, lapply(seq_along(v), function(i)
      lapply(perms(v[-i]), function(p) c(v[i], p))))
  max(vapply(perms(kt), function(pm) {
    map <- setNames(pm, kt); mean(map[as.character(est)] == truth)
  }, numeric(1)))
}

clust_dfd_vec <- function(dfd, nD) {
  v <- if (is.list(dfd)) as.numeric(dfd$avg.level) else as.numeric(dfd)
  if (length(v) < nD) v <- c(v, rep(NA, nD - length(v)))
  v[seq_len(nD)]
}

# ---- 5. Run -----------------------------------------------------------------
clust_run <- function(d, params) {
  Y <- d$Y; X <- d$X; W <- d$W
  K <- as.integer(params$ncluster %||% 2)
  stop_crt <- params$stop.crt %||% "nD"
  nD <- if (stop_crt == "nD") (params$nD %||% 1) else NULL
  ninitial <- params$ninitial %||% 5
  fit <- suppressWarnings(.cluster_env$capPCL(Y, X, W, ncluster = K,
            stop.crt = stop_crt, nD = nD, DfD.thred = params$DfD.thred %||% 2,
            ninitial = ninitial, seed = 100, score.return = TRUE, verbose = FALSE))
  sims <- params$sims %||% 0
  inference <- NULL
  if (sims > 0) {
    inference <- lapply(seq_len(ncol(fit$gamma)), function(j)
      tryCatch(suppressWarnings(.cluster_env$capPCL_coef_boot(
                 Y, X, W, gamma = fit$gamma[, j], eta = fit$eta.est[[j]],
                 boot = TRUE, sims = sims, boot.ci.type = "se", verbose = FALSE)),
               error = function(e) NULL))
  }
  list(fit = fit, inference = inference, X = X, W = W, ids = d$ids, K = K,
       var_names = d$var_names, x_names = d$x_names, w_names = d$w_names,
       truth = d$truth, params = params)
}

# melt a (terms x 6 x K) bootstrap array (stats in dim 2, cluster in dim 3)
clust_melt_boot <- function(arr, term_names, kind) {
  K <- dim(arr)[3]
  do.call(rbind, lapply(seq_len(K), function(k)
    data.frame(Cluster = paste0("cluster", k), Term = term_names,
               Estimate = round(arr[, 1, k], 3), SE = round(arr[, 2, k], 3),
               pvalue = signif(arr[, 4, k], 3), LB = round(arr[, 5, k], 3),
               UB = round(arr[, 6, k], 3), check.names = FALSE, row.names = NULL)))
}

# ---- 6. Tables --------------------------------------------------------------
clust_summarize <- function(res) {
  fit <- res$fit; nD <- ncol(fit$gamma); dirs <- paste0("C", seq_len(nD)); out <- list()

  # cluster sizes per component
  sz <- do.call(rbind, lapply(seq_len(nD), function(j) {
    tab <- table(factor(fit$class[, j], levels = seq_len(res$K)))
    data.frame(Component = dirs[j], rbind(as.integer(tab)), check.names = FALSE)
  }))
  names(sz)[-1] <- paste0("cluster", seq_len(res$K))
  out[["Cluster sizes"]] <- sz

  # within-cluster variance-model beta (bootstrap if available, else point est.)
  if (!is.null(res$inference) && !is.null(res$inference[[1]])) {
    out[["β within-cluster (variance model)"]] <-
      cbind(Component = dirs[1], clust_melt_boot(res$inference[[1]]$beta, res$x_names, "beta"))
    out[["α membership (multinomial)"]] <-
      cbind(Component = dirs[1], clust_melt_boot(res$inference[[1]]$alpha, res$w_names, "alpha"))
  } else {
    b <- fit$beta[[1]]
    out[["β within-cluster (variance model)"]] <- data.frame(
      Cluster = rep(paste0("cluster", seq_len(ncol(b))), each = nrow(b)),
      Term = rep(res$x_names, ncol(b)), Estimate = round(as.numeric(b), 3), check.names = FALSE)
  }

  # gamma loadings
  g <- as.data.frame(fit$gamma); names(g) <- dirs
  out[["γ loadings"]] <- cbind(Variable = res$var_names, g)

  # DfD
  out[["DfD across dimensions"]] <- data.frame(Component = dirs,
    DfD = clust_dfd_vec(fit$DfD, nD), check.names = FALSE)

  # truth comparison (example data): match each component to closest true direction
  if (!is.null(res$truth)) {
    Gt <- as.matrix(res$truth$gamma)
    rec <- do.call(rbind, lapply(seq_len(nD), function(j) {
      ip <- abs(as.numeric(t(fit$gamma[, j]) %*% Gt)); k <- which.max(ip)
      tname <- colnames(Gt)[k]
      data.frame(Component = dirs[j], `Matched true` = tname,
                 `γ cosine` = round(ip[k], 3),
                 `cluster agreement` = round(clust_agree(fit$class[, j], res$truth$cluster[, k]), 3),
                 check.names = FALSE)
    }))
    out[["Direction & cluster recovery"]] <- rec
  }
  out
}

# ---- 7. Plots ---------------------------------------------------------------
clust_plots <- function(res) {
  fit <- res$fit; nD <- ncol(fit$gamma); dirs <- paste0("C", seq_len(nD)); plots <- list()

  # (a) gamma loadings
  g <- as.data.frame(fit$gamma); names(g) <- dirs
  long <- do.call(rbind, lapply(seq_len(nD), function(j)
    data.frame(Variable = res$var_names, Direction = dirs[j], Loading = g[, j])))
  long$Variable <- factor(long$Variable, levels = res$var_names)
  plots[["γ loadings"]] <- list(plot = plotly::plot_ly(
    long, x = ~Variable, y = ~Loading, color = ~Direction, type = "bar") |>
      plotly::layout(title = "Projection loadings (γ)", barmode = "group",
                     xaxis = list(showticklabels = FALSE)))

  # (b) projected log-variance score by estimated cluster (component 1)
  if (!is.null(fit$score)) {
    sc <- as.matrix(fit$score)[, 1]
    cl <- factor(fit$class[, 1])
    sdf <- data.frame(id = res$ids, score = log(pmax(sc, 1e-12)), cluster = cl)
    p_sc <- plotly::plot_ly(sdf, x = ~cluster, y = ~score, color = ~cluster,
                            type = "box", boxpoints = "all", jitter = 0.5, pointpos = 0,
                            text = ~id, hovertemplate = "%{text}<extra></extra>") |>
      plotly::layout(title = "Projected log-variance by estimated cluster (C1)",
                     yaxis = list(title = "log(γ' Σ γ)"), showlegend = FALSE)
    plots[["Score by cluster"]] <- list(plot = p_sc, data = sdf)
  }

  # (c) within-cluster beta with bootstrap CI (component 1)
  if (!is.null(res$inference) && !is.null(res$inference[[1]])) {
    arr <- res$inference[[1]]$beta; K <- dim(arr)[3]
    bdf <- do.call(rbind, lapply(seq_len(K), function(k)
      data.frame(Term = res$x_names, Cluster = paste0("cluster", k),
                 Estimate = arr[, 1, k], LB = arr[, 5, k], UB = arr[, 6, k])))
    bdf <- bdf[bdf$Term != "Intercept", , drop = FALSE]
    p_b <- plotly::plot_ly(bdf, x = ~Term, y = ~Estimate, color = ~Cluster, type = "bar",
                           error_y = list(type = "data", symmetric = FALSE,
                                          array = ~(UB - Estimate), arrayminus = ~(Estimate - LB))) |>
      plotly::layout(title = "Within-cluster variance coefficients β (C1) with 95% CI",
                     barmode = "group", yaxis = list(title = "estimate"))
    plots[["β by cluster"]] <- list(plot = p_b, data = bdf)
  }
  plots
}

# ---- 8. Register ------------------------------------------------------------
register_method(list(
  id = "cap-clustering",
  name = "CAP-clustering",
  full_name = "Parsimonious Clustering of Covariance Matrices (CAPclust)",
  short = "Model-based clustering of subjects whose data are covariance matrices, sharing a common CAP projection; clusters differ in their projected-variance model and membership depends on covariates.",
  status = "ready",
  tags = c("covariance regression", "clustering", "mixture model", "EM"),
  paper = list(
    citation = "Xu et al. (2026). Parsimonious clustering of covariance matrices. Biometrics (submitted).",
    url = NULL),
  explain = file.path(APP_DIR, "methods", "cap-clustering", "explain.md"),
  example_note = paste("The manuscript's simulation: 100 subjects, p=50, Tᵢ=100,",
                       "K=2 clusters, two covariate-driven components (D2 well separated,",
                       "D4 harder). The default (nD=1) recovers the leading component D2.",
                       "Runs in ~40 s (multinomial-EM bound)."),
  data_inputs = list(
    list(id = "Y", label = "Covariance data (list of T×p matrices)",
         help = "An .rds/.RData holding a list of length n; element i is a T_i × p matrix whose covariance is clustered. A long CSV (subject-id column + p columns) is also accepted."),
    list(id = "X", label = "Variance-model covariates (n × q1)",
         help = "Covariates for the within-cluster log-variance model, one row per subject (or a CSV with subject id first). Intercept added automatically."),
    list(id = "W", label = "Membership covariates (n × q2)",
         help = "Covariates that drive cluster membership (multinomial), one row per subject (or a CSV with subject id first). Intercept added automatically.")
  ),
  params = list(
    list(id = "ncluster", label = "Number of clusters (K)", type = "integer",
         default = 2, min = 2, max = 6, help = "Number of mixture components."),
    list(id = "stop.crt", label = "Number of components chosen by", type = "select",
         choices = c("Fixed number (nD)" = "nD", "DfD criterion" = "DfD"), default = "nD",
         help = "Fix the count, or keep adding components while DfD stays below the threshold. The leading component is the most reliably recovered."),
    list(id = "nD", label = "Number of components (nD)", type = "integer",
         default = 1, min = 1, max = 5, help = "Used when selecting by a fixed number. Each adds an EM + bootstrap pass (slower)."),
    list(id = "DfD.thred", label = "DfD threshold", type = "numeric",
         default = 2, min = 1, max = 100, step = 0.5, help = "Used when selecting by the DfD criterion."),
    list(id = "ninitial", label = "# random initializations", type = "integer",
         default = 5, min = 1, max = 20, help = "Random restarts for the γ / clustering EM."),
    list(id = "sims", label = "Bootstrap replicates (β, α inference)", type = "integer",
         default = 100, min = 0, max = 2000, help = "Subject-level bootstrap for coefficient inference; 0 = skip.")
  ),
  example = clust_example,
  export_example = function(d) list(Y = d$Y, X = d$Xdf, W = d$Wdf),
  parse = clust_parse,
  describe_data = clust_describe,
  run = clust_run,
  summarize = clust_summarize,
  plots = clust_plots
))
