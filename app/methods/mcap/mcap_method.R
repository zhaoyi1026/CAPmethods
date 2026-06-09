# =============================================================================
# MCAP -- Multilevel Covariate Assisted Principal regression (gamma-varying)
# -----------------------------------------------------------------------------
# Wraps the RcppArmadillo-accelerated implementation in ../LCAP_gamma-var/V4/ and
# registers it with the website. MCAP is the multilevel sibling of LCAP: data are
# nested in CLUSTERS, each cluster has its own projection gamma_i drawn from a
# von Mises-Fisher distribution around a population direction gamma (concentration
# kappa), and log(gamma_i' Sigma_ij gamma_i) follows a mixed model in the
# covariates with cluster-level random effects.
# =============================================================================

# ---- 1. Load the method code from the CAPmethods package --------------------
# Uses the package's "mcap" environment (lcapReg(), lcap.beta.boot(), and the
# bundled von Mises-Fisher sampler rvmf()); the mlcap_-prefixed kernel ships in
# the package DLL. See R/pkg_methods.R.
.mcap_env <- cap_pkg_env("mcap")

# ---- 2. Built-in example data = the manuscript simulation --------------------
# Faithful translation of the paper's gamma-varying simulation
# (p5_q4_2-1/gamma_unknown/case1): p=5 dimensions, q1=2 time-invariant FIXED
# covariates (X1), q2=1 time-varying RANDOM-slope covariate (X2), and TWO
# covariate-driven directions (the population basis columns rv.idx = c(2,4)).
# Cluster directions gamma_i ~ vMF(Gamma[,rv.idx], kappa = 10) are drawn mutually
# orthogonal; within a cluster
#   log(gamma_i' Sigma_ij gamma_i) = beta0_i + x1_ij' beta1 + x2_ij'(beta2_i)
# with a cluster random intercept (beta0_i) and random slope (beta2_i). Sizes
# m=20 clusters, n_i ~ Pois(50) units, T_ij ~ Pois(80) samples (a tractable point
# on the paper's grid). The fixed seed (seedl = 3) gives clean recovery of BOTH
# directions; nD=1 recovers the leading one, nD=2 / DfD recover both.
mcap_example <- function() {
  p <- 5L; q1 <- 2L; q2 <- 1L
  rv.idx <- c(2L, 4L)                                  # two covariate-driven dirs
  m <- 20L; nvec.lambda <- 50; Tmat.lambda <- 80; seedl <- 3L

  # population orthonormal basis (manuscript seed) + sign convention.
  # NB: the manuscript uses runif(p) recycled into a p x p matrix (not p*p) — keep
  # this exactly, or the basis (and the whole simulation) changes.
  set.seed(100)
  Gamma <- qr.Q(qr(matrix(stats::runif(p), p, p)))
  for (j in 1:p) if (Gamma[which.max(abs(Gamma[, j])), j] < 0) Gamma[, j] <- -Gamma[, j]

  # coefficient matrix: rows = beta0 (intercept), beta1/beta2 (fixed X1), beta3
  # (random-slope X2 mean); cols = dimensions.
  beta.mat <- rbind(seq(5, -1, length.out = 5),       # beta0  (random intercept mean)
                    c(0, 1, 0, -1, 0),                # beta1  (X1_1, fixed)
                    c(0, -0.5, 0, 0.5, 0),            # beta2  (X1_2, fixed)
                    c(0, -0.5, 0, 0.5, 0))            # beta3  (X2_1, random slope mean)
  sigma <- 0.1; Omega <- matrix(0.1^2, 1, 1); beta.nr.sd <- 0.1; kappa <- 10
  beta.fix.idx <- 2:(q1 + 1)                          # X1 rows -> 2:3
  beta.rnd.idx <- (q1 + 2):(q1 + q2 + 1)              # X2 row  -> 4

  set.seed(500); nvec <- stats::rpois(m, nvec.lambda)
  Tmat <- matrix(NA_integer_, m, max(nvec))
  for (i in 1:m) Tmat[i, 1:nvec[i]] <- stats::rpois(nvec[i], Tmat.lambda)

  # cluster random intercepts and random slopes
  set.seed(100); beta0.mat <- matrix(NA, m, p)
  for (kk in 1:p)
    beta0.mat[, kk] <- stats::rnorm(m, beta.mat[1, kk],
                                    if (kk %in% rv.idx) sigma else beta.nr.sd)
  set.seed(100); beta2.mat <- array(NA, c(m, q2, p))
  for (kk in 1:p) beta2.mat[, , kk] <- mvtnorm::rmvnorm(m, beta.mat[beta.rnd.idx, kk], sigma = Omega)

  # cluster-specific directions ~ vMF, orthogonalised across the two signals
  set.seed(100); gmat <- array(NA, c(m, p, length(rv.idx)))
  gmat[, , 1] <- .mcap_env$rvmf(m, Gamma[, rv.idx[1]], kappa)
  for (ss in 2:length(rv.idx)) for (i in 1:m) {
    repeat {
      otmp <- .mcap_env$rvmf(10000, Gamma[, rv.idx[ss]], kappa)
      oo <- apply(abs(otmp %*% gmat[i, , 1:(ss - 1), drop = FALSE]), 1, sum)
      idx <- which.min(oo)
      if (oo[idx] <= 1e-4) { gmat[i, , ss] <- otmp[idx, ]; break }
    }
  }
  Pi <- array(NA, c(p, p, m))
  for (i in 1:m) {
    Pt <- matrix(NA, p, p)
    Pt[, rv.idx]  <- gmat[i, , ]
    Pt[, -rv.idx] <- MASS::Null(gmat[i, , ])
    Pi[, , i] <- Pt
  }

  # generate covariates + responses
  set.seed(as.numeric("20221206") + seedl * 100)
  X1.fix <- stats::rbinom(m, 1, 0.5)                  # cluster-level binary covariate
  ids <- paste0("G", 1:m)
  x1_names <- paste0("X1_", 1:q1); x2_names <- paste0("X2_", 1:q2)
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
      Sig <- Pi[, , i] %*% diag(delta) %*% t(Pi[, , i])   # symmetric by construction
      set.seed(10000 + i * 100 + j * 10)
      Yi[[j]] <- mvtnorm::rmvnorm(Tmat[i, j], rep(0, p), sigma = Sig)
    }
    Y[[i]] <- Yi
  }

  var_names <- paste0("V", 1:p)
  Gtrue <- Gamma[, rv.idx]; colnames(Gtrue) <- paste0("D", seq_along(rv.idx))
  Btrue <- beta.mat[, rv.idx]
  rownames(Btrue) <- c("Intercept", x1_names, x2_names)
  colnames(Btrue) <- paste0("D", seq_along(rv.idx))
  d <- list(Y = Y, X1 = X1, X2 = X2, ids = ids, var_names = var_names,
            x1_names = x1_names, x2_names = x2_names,
            cov_names = c(x1_names, x2_names),
            truth = list(gamma = Gtrue, beta = Btrue, kappa = kappa, cl_cos = cl_cos))
  d$preview_ui <- mcap_preview(d)
  d
}

# ---- 3. Parsers for uploaded data -------------------------------------------
# Y: native .rds nested list (cluster -> unit -> T_ij x p matrix), OR long CSV
#    [cluster, unit, V1..Vp]. X1 (fixed-effect covariates) and X2 (random-effect
#    covariates) are each a list of n_i x q matrices, OR a long CSV
#    [cluster, unit, cov...] (no intercept column; intercept added internally).
mcap_parse <- function(files, opts) {
  yp <- mcap_coerce_Y(files$Y)
  x1 <- mcap_coerce_X(files$X1, ids = yp$ids, nivec = yp$nivec, prefix = "X1_")
  x2 <- mcap_coerce_X(files$X2, ids = yp$ids, nivec = yp$nivec, prefix = "X2_")
  d <- list(Y = yp$Y, X1 = x1$X, X2 = x2$X, ids = yp$ids, var_names = yp$var_names,
            x1_names = x1$cov_names, x2_names = x2$cov_names,
            cov_names = c(x1$cov_names, x2$cov_names), truth = NULL)
  d$preview_ui <- mcap_preview(d)
  d
}

mcap_coerce_Y <- function(obj) {
  if (is.data.frame(obj)) return(mcap_Y_from_long(obj))
  is_nested <- is.list(obj) && length(obj) > 0 && is.list(obj[[1]]) &&
    !is.matrix(obj[[1]]) && !is.data.frame(obj[[1]])
  if (is.list(obj) && !is.null(names(obj)) && !is_nested) {
    key <- intersect(c("Y", "Y_list", "response"), names(obj))
    if (length(key)) obj <- obj[[key[1]]]
    is_nested <- is.list(obj) && length(obj) > 0 && is.list(obj[[1]]) &&
      !is.matrix(obj[[1]]) && !is.data.frame(obj[[1]])
  }
  if (!is_nested)
    stop("Response must be a nested list (cluster -> unit -> T x p matrix), or a long CSV with cluster + unit columns.")
  Y <- lapply(obj, function(cl) lapply(cl, function(mm) {
    mm <- as.matrix(mm); storage.mode(mm) <- "double"; mm }))
  ids <- names(Y); if (is.null(ids) || any(ids == "")) ids <- paste0("Cl", seq_along(Y))
  names(Y) <- ids
  nivec <- vapply(Y, length, integer(1))
  vn <- colnames(Y[[1]][[1]]); if (is.null(vn)) vn <- paste0("V", seq_len(ncol(Y[[1]][[1]])))
  list(Y = Y, ids = ids, nivec = nivec, var_names = vn)
}

mcap_Y_from_long <- function(df) {
  cl <- as.character(df[[1]]); un <- as.character(df[[2]])
  resp <- df[, -(1:2), drop = FALSE]
  resp[] <- lapply(resp, function(x) suppressWarnings(as.numeric(x)))
  ids <- unique(cl)
  Y <- lapply(ids, function(i) {
    us <- unique(un[cl == i])
    lapply(us, function(u) as.matrix(resp[cl == i & un == u, , drop = FALSE]))
  })
  names(Y) <- ids
  list(Y = Y, ids = ids, nivec = vapply(Y, length, integer(1)), var_names = names(resp))
}

mcap_coerce_X <- function(obj, ids, nivec, prefix = "X") {
  if (is.list(obj) && !is.data.frame(obj) && !is.matrix(obj) && length(obj) > 0 &&
      (is.matrix(obj[[1]]) || is.data.frame(obj[[1]]))) {
    X <- lapply(obj, function(mm) { mm <- as.matrix(mm); storage.mode(mm) <- "double"; mm })
  } else if (is.list(obj) && !is.data.frame(obj) && !is.matrix(obj)) {
    key <- intersect(c("X1", "X2", "X", "covariates"), names(obj))
    if (!length(key)) stop("Could not find a covariate object in the upload.")
    return(mcap_coerce_X(obj[[key[1]]], ids, nivec, prefix))
  } else if (is.data.frame(obj)) {
    cl <- as.character(obj[[1]]); cov <- obj[, -(1:2), drop = FALSE]
    cov[] <- lapply(cov, function(x) suppressWarnings(as.numeric(x)))
    X <- lapply(ids, function(i) as.matrix(cov[cl == i, , drop = FALSE]))
  } else stop("Covariates must be a list of n_i x q matrices, or a long CSV with cluster + unit columns.")
  if (length(X) != length(ids)) stop("Covariate clusters do not match the response.")
  if (!all(vapply(seq_along(X), function(i) nrow(X[[i]]) == nivec[i], logical(1))))
    stop("Each cluster's covariate rows must equal that cluster's number of units.")
  cov_names <- colnames(X[[1]]) %||% paste0(prefix, seq_len(ncol(X[[1]])))
  for (i in seq_along(X)) colnames(X[[i]]) <- cov_names
  names(X) <- ids
  list(X = X, cov_names = cov_names)
}

# ---- 4. Preview / helpers ---------------------------------------------------
mcap_preview <- function(d) {
  nivec <- vapply(d$Y, length, integer(1))
  Tvec <- unlist(lapply(d$Y, function(cl) vapply(cl, nrow, integer(1))))
  p <- ncol(d$Y[[1]][[1]])
  rng <- function(x) if (length(unique(x)) == 1) as.character(x[1]) else paste0(min(x), "–", max(x))
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Clusters", length(d$Y), theme = "primary"),
      bslib::value_box("Dimension (p)", p, theme = "secondary"),
      bslib::value_box("Units / cluster", rng(nivec), theme = "secondary"),
      bslib::value_box("Samples / unit", rng(Tvec), theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           sprintf("Fixed-effect covariates (X1): %s | random-effect covariates (X2): %s.%s",
                   paste(d$x1_names %||% character(0), collapse = ", "),
                   paste(d$x2_names %||% character(0), collapse = ", "),
                   if (!is.null(d$truth)) "  Simulated data with known truth (cluster-varying γ; two directions)." else ""))
  )
}

mcap_describe <- function(d) {
  nivec <- vapply(d$Y, length, integer(1))
  sprintf("%d clusters, %s units each, dimension p = %d, fixed X1: %s, random X2: %s.",
          length(d$Y),
          if (length(unique(nivec)) == 1) nivec[1] else paste0(min(nivec), "-", max(nivec)),
          ncol(d$Y[[1]][[1]]),
          paste(d$x1_names %||% "—", collapse = ", "),
          paste(d$x2_names %||% "—", collapse = ", "))
}

mcap_dfd_vec <- function(fit) {
  d <- if (is.list(fit$DfD)) fit$DfD$DfD.avg else fit$DfD
  as.numeric(d)
}

# ---- 5. Run -----------------------------------------------------------------
mcap_run <- function(d, params) {
  stop_crt <- params$stop.crt %||% "nD"
  fit <- .mcap_env$lcapReg(
    d$Y, X1 = d$X1, X2 = d$X2, data.type = "Y",
    stop.crt = stop_crt,
    nD = if (stop_crt == "nD") (params$nD %||% 2) else NULL,
    DfD.thred = params$DfD.thred %||% 2,
    method = "CAP", H.type = "CAvgCov", Omega.diag = TRUE,
    ninitial = if (is.null(params$ninitial) || params$ninitial <= 0) NULL else params$ninitial,
    seed = 100, score.return = TRUE, verbose = FALSE)
  list(fit = fit, ids = d$ids, var_names = d$var_names,
       x1_names = d$x1_names, x2_names = d$x2_names, cov_names = d$cov_names,
       X1 = d$X1, X2 = d$X2, truth = d$truth, params = params)
}

# greedy sign-aligned match of estimated population directions to true
mcap_match_dirs <- function(Ghat, Gtrue) {
  J <- ncol(Ghat); K <- ncol(Gtrue); used <- integer(0); out <- list()
  for (j in seq_len(J)) {
    ip <- vapply(seq_len(K), function(k) if (k %in% used) NA_real_ else sum(Ghat[, j] * Gtrue[, k]), numeric(1))
    if (all(is.na(ip))) break
    k <- which.max(abs(ip)); used <- c(used, k)
    out[[length(out) + 1]] <- data.frame(est = colnames(Ghat)[j] %||% paste0("C", j),
      true = colnames(Gtrue)[k] %||% paste0("D", k), k = k, sign = sign(ip[k]), cosine = abs(ip[k]))
  }
  do.call(rbind, out)
}

# ---- 6. Tables --------------------------------------------------------------
mcap_summarize <- function(res) {
  fit <- res$fit
  dirs <- colnames(fit$gamma) %||% paste0("C", seq_len(ncol(fit$gamma)))
  terms <- rownames(fit$beta) %||% c("Intercept", res$cov_names)
  out <- list()

  # fixed-effects beta (population) per direction
  bdf <- as.data.frame(fit$beta); names(bdf) <- dirs
  out[["β fixed effects"]] <- cbind(Term = terms, bdf)

  # vMF concentration kappa per direction
  out[["κ (vMF concentration)"]] <- data.frame(Direction = dirs,
    kappa = as.numeric(fit$kappa), check.names = FALSE)

  # random-effect variances: random intercept sigma^2 + random-slope diag(Omega)
  re_rows <- do.call(rbind, lapply(seq_along(dirs), function(k) {
    s2 <- as.numeric(fit$beta0.sigma2[, k])
    Om <- if (!is.null(fit$beta2.Omega)) diag(as.matrix(fit$beta2.Omega[, , k])) else numeric(0)
    data.frame(Direction = dirs[k],
               Component = c("Intercept (random)", paste0(res$x2_names, " (random slope)")),
               Variance = c(s2, Om), check.names = FALSE)
  }))
  out[["Random-effect variances"]] <- re_rows

  # population gamma loadings
  g <- as.data.frame(fit$gamma); names(g) <- dirs
  out[["γ population loadings"]] <- cbind(Variable = res$var_names, g)

  # cluster direction spread: cosine of each cluster's gamma_i to population gamma
  cos_df <- data.frame(Cluster = res$ids)
  for (k in seq_along(dirs)) {
    gk <- fit$gamma[, k]
    cos_df[[dirs[k]]] <- vapply(seq_along(res$ids),
      function(i) abs(sum(fit$gamma.rnd[, k, i] * gk)), numeric(1))
  }
  out[["Cluster γᵢ cosine to population"]] <- cos_df

  # DfD across dimensions
  dfd <- mcap_dfd_vec(fit)
  out[["DfD across dimensions"]] <- data.frame(Dimension = dirs[seq_along(dfd)],
    `# directions (1..k)` = seq_along(dfd), DfD = dfd, check.names = FALSE)

  # orthogonality of population directions
  if (ncol(fit$gamma) > 1 && !is.null(fit$gamma.othogonality)) {
    o <- as.data.frame(fit$gamma.othogonality)
    out[["Orthogonality γ'γ"]] <- cbind(" " = dirs, o)
  }

  # truth comparison (example data)
  if (!is.null(res$truth)) {
    Gtrue <- as.matrix(res$truth$gamma)
    mt <- mcap_match_dirs(fit$gamma, Gtrue)
    cmp <- data.frame(Variable = res$var_names)
    for (r in seq_len(nrow(mt))) {
      cmp[[paste0(mt$true[r], "_true")]] <- Gtrue[, mt$k[r]]
      cmp[[paste0(mt$est[r], "_est")]]  <- fit$gamma[, mt$est[r]] * mt$sign[r]
    }
    out[["γ: truth vs estimate"]] <- cmp
    out[["Direction recovery"]] <- data.frame(Estimated = mt$est, `Matched true` = mt$true,
      `Cosine similarity` = round(mt$cosine, 4), check.names = FALSE)
    if (!is.null(res$truth$beta)) {
      Btrue <- as.matrix(res$truth$beta)
      bcmp <- do.call(rbind, lapply(seq_len(nrow(mt)), function(r) {
        k <- mt$k[r]; if (k > ncol(Btrue)) return(NULL)
        data.frame(Direction = mt$est[r], Term = terms,
                   True = Btrue[match(terms, rownames(Btrue)), k],
                   Estimate = fit$beta[, mt$est[r]], check.names = FALSE, row.names = NULL)
      }))
      if (!is.null(bcmp)) out[["β: truth vs estimate"]] <- bcmp
    }
  }
  out
}

# ---- 7. Plots ---------------------------------------------------------------
mcap_plots <- function(res) {
  fit <- res$fit; plots <- list()
  dirs <- colnames(fit$gamma) %||% paste0("C", seq_len(ncol(fit$gamma)))
  vars <- res$var_names

  # (a) population loadings
  g <- as.data.frame(fit$gamma); names(g) <- dirs
  long <- do.call(rbind, lapply(seq_along(dirs), function(k)
    data.frame(Variable = vars, Direction = dirs[k], Loading = g[, k])))
  long$Variable <- factor(long$Variable, levels = vars)
  plots[["Population loadings (γ)"]] <- list(plot = plotly::plot_ly(
    long, x = ~Variable, y = ~Loading, color = ~Direction, type = "bar") |>
      plotly::layout(title = "Population projection loadings (γ)", barmode = "group"))

  # (b) cluster direction spread (cosine to population) for D1
  cos1 <- vapply(seq_along(res$ids),
    function(i) abs(sum(fit$gamma.rnd[, 1, i] * fit$gamma[, 1])), numeric(1))
  cdf <- data.frame(Cluster = factor(res$ids, levels = res$ids), Cosine = cos1)
  plots[["Cluster γ spread (D1)"]] <- list(plot = plotly::plot_ly(
    cdf, x = ~Cluster, y = ~Cosine, type = "bar",
    marker = list(color = "#2c7fb8")) |>
      plotly::layout(title = "Cluster-specific γᵢ alignment to population γ (D1)",
                     yaxis = list(title = "|cos(γᵢ, γ)|", range = c(0, 1))),
    data = cdf)

  # (c) cluster-unit scores vs fitted (D1) -- fitted = cluster random intercept
  #     + X1 fixed effects + X2 cluster random slopes
  if (!is.null(fit$score)) {
    X1 <- res$X1; X2 <- res$X2
    bX1 <- if (length(res$x1_names)) fit$beta[res$x1_names, , drop = FALSE] else NULL
    sv <- do.call(rbind, lapply(seq_along(dirs), function(k) {
      do.call(rbind, lapply(seq_along(X2), function(i) {
        nu <- nrow(X2[[i]])
        b0 <- fit$beta0.rnd[i, k]
        b2 <- if (!is.null(fit$beta2.rnd)) fit$beta2.rnd[i, , k] else numeric(0)
        fx1 <- if (!is.null(bX1)) as.numeric(X1[[i]] %*% bX1[, k]) else 0
        ft <- b0 + fx1 + as.numeric(X2[[i]] %*% b2)
        sc <- fit$score[i, seq_len(nu), k]
        data.frame(Cluster = res$ids[i], unit = seq_len(nu), Direction = dirs[k],
                   score = sc, fitted = ft, logvar = log(pmax(sc, 1e-12)))
      }))
    }))
    rng <- range(c(sv$fitted, sv$logvar), finite = TRUE)
    p_sc <- plotly::plot_ly(sv, x = ~fitted, y = ~logvar, color = ~Direction,
                            type = "scatter", mode = "markers",
                            marker = list(size = 6, opacity = 0.5),
                            text = ~paste0(Cluster, " u", unit),
                            hovertemplate = "%{text}<extra></extra>") |>
      plotly::add_lines(x = rng, y = rng, line = list(dash = "dot", color = "grey"),
                        showlegend = FALSE, inherit = FALSE) |>
      plotly::layout(title = "Cluster-unit projected variance vs fitted",
                     xaxis = list(title = "fitted log-variance (β₀ᵢ + x'β₂ᵢ)"),
                     yaxis = list(title = "log(γᵢ' Σ γᵢ)"))
    plots[["Cluster-unit scores"]] <- list(plot = p_sc, data = sv)
  }

  # (d) DfD across dimensions
  dfd <- mcap_dfd_vec(fit)
  ddf <- data.frame(k = seq_along(dfd), DfD = dfd, Dimension = dirs[seq_along(dfd)])
  p_dfd <- plotly::plot_ly(ddf, x = ~k, y = ~DfD, type = "scatter", mode = "lines+markers",
                           marker = list(size = 9, color = "#2c7fb8"), line = list(color = "#2c7fb8"),
                           text = ~Dimension, hovertemplate = "%{text}: %{y:.4f}<extra></extra>")
  if (!is.null(res$params$DfD.thred)) {
    thr <- res$params$DfD.thred
    p_dfd <- plotly::add_lines(p_dfd, x = range(ddf$k), y = c(thr, thr),
                               line = list(dash = "dash", color = "firebrick"),
                               name = paste0("threshold = ", thr), inherit = FALSE)
  }
  p_dfd <- plotly::layout(p_dfd, title = "DfD across dimensions",
                          xaxis = list(title = "Number of directions (1..k)",
                                       tickvals = ddf$k, ticktext = ddf$Dimension),
                          yaxis = list(title = "DfD (avg. deviation from diagonality)"))
  plots[["DfD"]] <- list(plot = p_dfd)
  plots
}

# downloadable long-format example templates
mcap_export_example <- function(d) {
  Ylong <- do.call(rbind, lapply(seq_along(d$Y), function(i)
    do.call(rbind, lapply(seq_along(d$Y[[i]]), function(j) {
      mm <- as.data.frame(d$Y[[i]][[j]]); names(mm) <- d$var_names
      cbind(cluster = d$ids[i], unit = j, mm) }))))
  Xlong <- function(X) do.call(rbind, lapply(seq_along(X), function(i) {
    mm <- as.data.frame(X[[i]])
    cbind(cluster = d$ids[i], unit = seq_len(nrow(mm)), mm) }))
  X1long <- Xlong(d$X1); X2long <- Xlong(d$X2)
  rownames(Ylong) <- rownames(X1long) <- rownames(X2long) <- NULL
  list(Y = Ylong, X1 = X1long, X2 = X2long)
}

# ---- 8. Register ------------------------------------------------------------
register_method(list(
  id = "lcap-var",
  name = "MCAP",
  full_name = "Multilevel Covariate Assisted Principal Regression (cluster-varying γ)",
  short = "CAP for multilevel covariance outcomes: a cluster-varying projection drawn from a von Mises-Fisher distribution, with cluster-level random effects.",
  status = "ready",
  tags = c("covariance regression", "multilevel", "random effects", "von Mises-Fisher"),
  paper = list(
    citation = "Green, et al. Multilevel covariate-assisted principal regression (see green2026multilevel.pdf).",
    url = NULL),
  explain = file.path(APP_DIR, "methods", "mcap", "explain.md"),
  example_note = paste("The manuscript's γ-varying simulation (p5_q4_2-1, case 1):",
                       "20 clusters, nᵢ~Pois(50) units, Tᵢⱼ~Pois(80), p=5, with TWO",
                       "covariate-driven directions, cluster γᵢ ~ vMF(γ, κ=10), 2 fixed",
                       "covariates (X1) + 1 random-slope covariate (X2). Default nD=2",
                       "recovers both directions (γ cosine ≈0.99/0.94; ~2–3 min); nD=1",
                       "recovers the leading one faster."),
  data_inputs = list(
    list(id = "Y", label = "Response (nested list or long CSV)",
         help = "An .rds/.RData holding a nested list: element i is a list over cluster i's units, each a T_ij x p matrix. A long CSV with columns [cluster, unit, V1..Vp] is also accepted."),
    list(id = "X1", label = "Fixed-effect covariates (list or long CSV)",
         help = "An .rds list: element i is an n_i x q1 matrix (one row per unit) of covariates with FIXED (population) effects. A long CSV with columns [cluster, unit, cov...] is also accepted. No intercept column (β₀ is added internally)."),
    list(id = "X2", label = "Random-effect covariates (list or long CSV)",
         help = "An .rds list: element i is an n_i x q2 matrix (one row per unit) of covariates with cluster-level RANDOM slopes. A long CSV with columns [cluster, unit, cov...] is also accepted. No intercept column.")
  ),
  params = list(
    list(id = "stop.crt", label = "Number of directions chosen by", type = "select",
         choices = c("Fixed number (nD)" = "nD", "DfD criterion" = "DfD"), default = "nD",
         help = "Fix the count, or keep adding directions while deviation-from-diagonality (DfD) stays below the threshold."),
    list(id = "nD", label = "Number of directions (nD)", type = "integer",
         default = 2, min = 1, max = 5,
         help = "Used when selecting by a fixed number. The example has two covariate-driven directions: nD=1 recovers the leading one (faster), nD=2 recovers both."),
    list(id = "DfD.thred", label = "DfD threshold", type = "numeric",
         default = 2, min = 1, max = 100, step = 0.5,
         help = "Used when selecting by the DfD criterion: keep adding directions while the deviation-from-diagonality stays below this threshold."),
    list(id = "ninitial", label = "# random initializations", type = "integer",
         default = 8, min = 1, max = 30,
         help = "More initializations = more robust optimization (recommended ≥8 for clean recovery), slower.")
  ),
  example = mcap_example,
  export_example = mcap_export_example,
  parse = mcap_parse,
  describe_data = mcap_describe,
  run = mcap_run,
  summarize = mcap_summarize,
  plots = mcap_plots
))
