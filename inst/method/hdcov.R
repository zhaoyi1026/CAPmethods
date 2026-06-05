# =============================================================================
# CAP-HDcov / HCAP -- V2 (RcppArmadillo-accelerated covariance algebra)
# -----------------------------------------------------------------------------
# Same interface and results as V1; the p x p covariance-algebra primitives
# (precompute_sigma, projected scores, weighted accumulation) call a compiled
# kernel. The method's cost is dominated by the high-dimensional cv.glmnet / glm
# fitting, which stays in R, so the end-to-end speedup is small unless the
# response dimension p is large. See CAP-HDcov/README.md.
# =============================================================================

# Required packages

# --- Load the compiled kernel (hdcov_kernels.cpp -> shared object) -----------
# Set option "hdcov.v2.dir" or variable HDCOV_V2_DIR to point at the folder with
# hdcov_kernels.cpp; defaults to "CAP-HDcov/V2".
# projected scores v'Sigma_i v ; weighted accumulation sum_i w_i Sigma_i
.hdcov_score <- function(Sigma, v) as.numeric(.Call("hdcov_score_cpp", Sigma, as.numeric(v)))
.hdcov_accum <- function(Sigma, w) .Call("hdcov_accum_cpp", Sigma, as.numeric(w))

# === Precompute covariance matrices (RcppArmadillo kernel) ===
# Returns list(Sigma = p x p x n array, Tvec = per-subject sizes).
precompute_sigma <- function(Y) {
  .Call("hdcov_precompute_cpp", lapply(Y, as.matrix))
}

# === Projection Estimation via Covariance Regression with Multiple Initial Values ===
# Iteratively estimates a shared linear projection by optimizing the objective from multiple initial values.
# Selects the best result to avoid local minima. Supports parallel evaluation.
# Input:
#   - Y, X: List of responses and covariate matrix.
#   - Sigma, Tvec: Precomputed covariances and sample sizes.
#   - gamma0.mat: Optional matrix of initial projections.
#   - ninitial: Number of random initializations if gamma0.mat is not provided.
# Output:
#   - A list containing estimated gamma, beta, and diagnostics.

MatReg_QC_opt <- function(Y, X, Sigma, Tvec, method = c("CAP"),
                          max.itr = 1000, tol = 1e-4, trace = FALSE,
                          score.return = TRUE, gamma0.mat = NULL, ninitial = NULL) {
  n <- length(Y)
  p <- ncol(Y[[1]])
  q <- ncol(X)
  
  # Construct a set of initial gamma vectors if not provided
  if (is.null(gamma0.mat)) {
    gamma0.mat <- matrix(NA, p, p + 6)
    for (j in 1:p) {
      gamma0.mat[, j] <- rep(0, p); gamma0.mat[j, j] <- 1
    }
    gamma0.mat[, p + 1] <- rep(1, p) / sqrt(sum(rep(1, p)^2))
    set.seed(500)
    gamma.rand <- matrix(rnorm(5 * p), nrow = p)
    gamma0.mat[, (p + 2):(p + 6)] <- apply(gamma.rand, 2, function(x) x / sqrt(sum(x^2)))
  }
  
  # Determine the number of initializations to evaluate
  if (is.null(ninitial)) {
    ninitial <- min(ncol(gamma0.mat), 10)
  } else {
    ninitial <- min(ninitial, ncol(gamma0.mat))
  }
  
  # Randomly sample initial gamma projections
  set.seed(2025)
  gamma0.mat <- gamma0.mat[, sort(sample(1:ncol(gamma0.mat), ninitial, replace = FALSE))]
  
  # Set up parallel computation
  ncores <- min(parallel::detectCores(), ninitial)
  registerDoParallel(cores = ncores)
  
  # Optimize over all initializations in parallel
  results <- foreach(k = 1:ninitial, .packages = c('MASS', 'glmnet'), .export = c("MatReg_QC", "MatReg_QC_beta")) %dopar% {
    g0 <- gamma0.mat[, k]
    tryCatch({
      re <- MatReg_QC(Y, X, Sigma, Tvec, gamma0 = g0,
                      method = method, max.itr = max.itr,
                      tol = tol, trace = trace, score.return = score.return)
      
      # Normalize gamma
      gamma.scaled <- re$gamma / sqrt(sum(re$gamma^2))
      
      # Re-estimate beta based on scaled gamma
      re.scale <- MatReg_QC_beta(Y, X, Sigma, Tvec, gamma0 = gamma.scaled,
                                 max.itr = max.itr, tol = tol, score.return = score.return)
      
      # Evaluate objective values before and after scaling
      obj <- sum((X %*% re$beta) * Tvec) / 2 +
        as.numeric(t(re$gamma) %*%
                     .hdcov_accum(Sigma, Tvec / exp(as.numeric(X %*% re$beta))) %*% re$gamma / 2)

      obj.scale <- sum((X %*% re.scale$beta) * Tvec) / 2 +
        as.numeric(t(re.scale$gamma) %*%
                     .hdcov_accum(Sigma, Tvec / exp(as.numeric(X %*% re.scale$beta))) %*% re.scale$gamma / 2)
      
      result <- list(re = re, obj = obj, re.scale = re.scale, obj.scale = obj.scale)
    }, error = function(e) NULL)
  }
  
  closeAllConnections()
  
  # Filter out failed results
  results <- results[!sapply(results, is.null)]
  if (length(results) == 0) stop("All parallel optimizations failed.")
  
  # Select the result with the lowest scaled objective
  obj_scales <- sapply(results, function(r) r$obj.scale)
  best_idx <- which.min(unlist(obj_scales))
  best_result <- results[[best_idx]]
  
  return(best_result$re.scale)
}

# === Estimation with Multiple Components ===
# Recursively estimates additional components by projecting out previously estimated ones.
# Can optionally enforce orthogonality constraints.
# Input:
#   - Phi0: Matrix of previously estimated projections.
#   - CAP.OC: Whether to enforce orthogonality constraints.
# Output:
#   - Estimated orijection and corresponding regression coefficients.

MatReg_QC_opt2 <- function(Y, X, Phi0=NULL, Sigma, Tvec, method=c("CAP"), CAP.OC=FALSE,
                           max.itr=1000, tol=1e-4, trace=FALSE, score.return=FALSE,
                           gamma0.mat=NULL, ninitial=NULL) {
  if (is.null(Phi0)) {
    return(MatReg_QC_opt(Y, X, Sigma, Tvec, method=method, max.itr=max.itr, tol=tol,
                         trace=trace, score.return=score.return, gamma0.mat=gamma0.mat, ninitial=ninitial))
  } else {
    n <- length(Y)
    q <- ncol(X)
    p <- ncol(Y[[1]])
    p0 <- ncol(Phi0)
    
    # Estimate variance components associated with existing projections
    beta0 <- numeric(p0)
    for (j in 1:p0) {
      beta0[j] <- MatReg_QC_beta(Y, X, Sigma, Tvec, gamma0=Phi0[,j],
                                 max.itr=max.itr, tol=tol)$beta[1]
    }
    
    # Project out previous projections and reconstruct residual responses
    Ytmp <- vector("list", n)
    for (i in 1:n) {
      Y2tmp <- Y[[i]] - Y[[i]] %*% (Phi0 %*% t(Phi0))
      Y2svd <- svd(Y2tmp)
      Ytmp[[i]] <- Y2svd$u %*% diag(c(Y2svd$d[1:(p - p0)], exp(beta0))) %*% t(Y2svd$v)
    }
    
    # Recompute covariances from projected responses
    sig.res <- precompute_sigma(Ytmp)
    Sigma.new <- sig.res$Sigma
    Tvec.new <- sig.res$Tvec
    
    # Estimate the next component either with or without orthogonality constraint
    if (CAP.OC == FALSE) {
      re <- MatReg_QC_opt(Ytmp, X, Sigma.new, Tvec.new, method=method, max.itr=max.itr,
                          tol=tol, trace=trace, score.return=score.return,
                          gamma0.mat=gamma0.mat, ninitial=ninitial)
    } else {
      re <- MatReg_QC_RE(Ytmp, X, Phi0=Phi0, Sigma.new, Tvec.new, max.itr=max.itr, tol=tol, trace=FALSE, score.return=score.return)
    }
    
    # Return orthogonality score for diagnostics
    re$orthogonal <- c(t(re$gamma) %*% Phi0)
    return(re)
  }
}

# === Estimation with Orthogonality Constraint ===
# Estimates a new projection orthogonal to previously estimated components.
# This function incorporates a closed-form eigenvalue solution for gamma update under orthogonality constraints
# and alternates between gamma and beta updates.
# Input:
#   - Y, X: List of response matrices and covariate matrix.
#   - Phi0: Previously estimated projections.
#   - Sigma, Tvec: Precomputed covariance matrices and sample sizes.
# Output:
#   - A list containing estimated gamma, beta, convergence status, and optionally trace and score.

MatReg_QC_RE <- function(Y, X, Phi0, Sigma, Tvec, max.itr = 1000, tol = 1e-4, trace = FALSE, score.return = TRUE) {
  n <- length(Y)
  p <- ncol(Y[[1]])
  q <- ncol(X)
  beta0 <- rep(0, q)
  
  if (trace) {
    gamma.trace <- NULL
    beta.trace <- beta0
  }
  
  s <- 0
  diff <- 100
  while (s <= max.itr && diff > tol) {
    s <- s + 1
    
    # --- Update gamma ---
    A <- .hdcov_accum(Sigma, Tvec / exp(as.numeric(X %*% beta0)))

    B <- apply(Sigma, c(1, 2), mean)
    B.inv <- solve(B)
    
    # Projection matrix onto the space spanned by Phi0 under B.inv metric
    P <- Phi0 %*% ginv(t(Phi0) %*% B.inv %*% Phi0) %*% t(Phi0) %*% B.inv
    
    # Compute square root of B.inv
    B.inv.eig <- eigen(B.inv)
    B.inv.rt <- B.inv.eig$vectors %*% diag(sqrt(B.inv.eig$values)) %*% solve(B.inv.eig$vectors)
    
    # Solve eigenvalue problem in transformed space
    eig_res <- eigen(B.inv.rt %*% (diag(p) - P) %*% A %*% B.inv.rt)
    gamma.new.cands <- B.inv.rt %*% Re(eig_res$vectors)
    
    # --- Update beta ---
    objvals <- rep(Inf, ncol(gamma.new.cands))
    beta.cands <- matrix(NA, q, ncol(gamma.new.cands))
    
    for (j in 1:ncol(gamma.new.cands)) {
      v0 <- gamma.new.cands[, j]
      Z <- .hdcov_score(Sigma, v0)
      idx <- which(Z > 0 & is.finite(Z))
      X_cal <- X[idx, -1, drop = FALSE]
      Z_cal <- Z[idx]
      
      if (length(Z_cal) > 5) {
        fit_glm <- try(cv.glmnet(X_cal, Z_cal, family = Gamma(link = "log"), nfolds = 5,
                                 type.measure = "deviance", intercept = TRUE, standardize = TRUE), silent = TRUE)
        if (!inherits(fit_glm, "try-error")) {
          ind0 <- which(fit_glm$nzero > 1 & fit_glm$nzero < sqrt(ncol(X_cal) + 1))
          if (length(ind0) > 0) {
            ind <- which.min(fit_glm$cvm[ind0])
            beta.tmp <- rep(0, q)
            beta.tmp[-1] <- as.vector(fit_glm$glmnet.fit$beta[, ind0[ind]])
            beta.tmp[1] <- fit_glm$glmnet.fit$a0[ind0[ind]]
            
            beta.cands[, j] <- beta.tmp
            
            S <- .hdcov_accum(Sigma, Tvec / exp(as.numeric(X %*% beta.tmp)))
            objvals[j] <- sum((X %*% beta.tmp) * Tvec) / 2 + as.numeric(t(v0) %*% S %*% v0 / 2)
          }
        }
      }
    }
    
    best_idx <- which.min(objvals)
    beta.new <- beta.cands[, best_idx]
    gamma.new <- gamma.new.cands[, best_idx]
    
    if (trace) {
      gamma.trace <- cbind(gamma.trace, gamma.new)
      beta.trace <- cbind(beta.trace, beta.new)
    }
    
    diff <- max(abs(beta.new - beta0))
    beta0 <- beta.new
  }
  
  # Normalize gamma
  gamma.new <- gamma.new / sqrt(sum(gamma.new^2))
  if (gamma.new[1] < 0) {
    gamma.new <- -gamma.new
  }
  
  beta.new <- MatReg_QC_beta(Y, X, Sigma, Tvec, gamma0 = gamma.new,
                             max.itr = max.itr, tol = tol, score.return = score.return)$beta
  
  if (score.return) {
    score <- .hdcov_score(Sigma, gamma.new)
  }
  
  if (trace) {
    if (score.return) {
      re <- list(gamma = c(gamma.new), beta = c(beta.new), convergence = (s < max.itr),
                 score = score, gamma.trace = gamma.trace, beta.trace = beta.trace)
    } else {
      re <- list(gamma = c(gamma.new), beta = c(beta.new), convergence = (s < max.itr),
                 gamma.trace = gamma.trace, beta.trace = beta.trace)
    }
  } else {
    if (score.return) {
      re <- list(gamma = c(gamma.new), beta = c(beta.new), convergence = (s < max.itr), score = score)
    } else {
      re <- list(gamma = c(gamma.new), beta = c(beta.new), convergence = (s < max.itr))
    }
  }
  
  return(re)
}

# === Iterative Optimization for Projections and Coefficients ===
# Alternating optimization of the shared projection and regression coefficients.
# Input:
#   - gamma0: Optional initialization for the projection.
#   - trace: Whether to return optimization trace.
# Output:
#   - Estimated projections and coefficients, optionally with trace.

MatReg_QC <- function(Y, X, Sigma, Tvec, method=c("CAP"), max.itr=1000, tol=1e-4,
                      trace=FALSE, gamma0=NULL, score.return=FALSE) {
  n <- length(Y)
  p <- ncol(Y[[1]])
  q <- ncol(X)
  
  Sigma.bar <- apply(Sigma, c(1, 2), mean, na.rm=TRUE)
  svd.bar <- eigen(Sigma.bar)
  Ds.inv <- diag(1 / sqrt(svd.bar$values))
  U <- svd.bar$vectors
  
  theta0 <- if (is.null(gamma0)) rep(1 / sqrt(p), p) else gamma0
  v0 <- U %*% Ds.inv %*% theta0
  beta0 <- rep(0, q)
  
  if (trace) {
    v.trace <- list(v0)
    beta.trace <- list(beta0)
    obj <- numeric()
  }
  
  s <- 0
  diff <- 100
  while (s <= max.itr && diff > tol) {
    s <- s + 1
    Z <- .hdcov_score(Sigma, v0)
    idx <- which(Z > 0 & is.finite(Z))
    X_cal <- X[idx, -1, drop=FALSE]
    Z_cal <- Z[idx]
    
    beta.new <- beta0
    if (length(Z_cal) > 5) {
      fit_glm <- try(cv.glmnet(X_cal, Z_cal, family=Gamma(link="log"), nfolds=5, type.measure = "deviance", 
                               intercept = TRUE, standardize = TRUE), silent=TRUE)
      if (!inherits(fit_glm, "try-error")) {
        ind0 <- which(fit_glm$nzero > 1 & fit_glm$nzero < sqrt(ncol(X_cal)+1))
        if (length(ind0) > 0) {
          ind <- which.min(fit_glm$cvm[ind0])
          beta.new[-1] <- as.matrix(fit_glm$glmnet.fit$beta[, ind0[ind]])
          beta.new[1] <- fit_glm$glmnet.fit$a0[ind0[ind]]
        }
      }
    }
    
    V1 <- .hdcov_accum(Sigma, Tvec / exp(as.numeric(X %*% beta.new)))
    B <- t(U) %*% V1 %*% U
    B.tilde <- Ds.inv %*% B %*% Ds.inv
    eig <- eigen(B.tilde)
    theta.new <- eig$vectors[, p]
    v.new <- U %*% Ds.inv %*% theta.new
    
    v.diff <- max(abs(v.new - v0))
    beta.diff <- max(abs(beta.new - beta0))
    diff <- max(c(v.diff, beta.diff))
    
    v0 <- v.new
    beta0 <- beta.new
    
    if (trace) {
      v.trace[[s + 1]] <- v.new
      beta.trace[[s + 1]] <- beta.new
      Q1 <- sum((X %*% beta.new) * Tvec) / 2
      S <- .hdcov_accum(Sigma, Tvec / exp(as.numeric(X %*% beta.new)))
      Q2 <- as.numeric(t(v.new) %*% S %*% v.new / 2)
      obj[s] <- Q1 + Q2
    }
  }
  
  if (v.new[1] < 0) v.new <- -v.new
  
  if (score.return) {
    score <- .hdcov_score(Sigma, v.new)
  }
  
  re <- list(gamma = c(v.new), beta = c(beta.new), convergence = (s < max.itr))
  if (score.return) re$score <- score
  if (trace) {
    re$gamma.trace <- do.call(cbind, v.trace)
    re$beta.trace <- do.call(cbind, beta.trace)
    re$obj <- obj
  }
  return(re)
}

# === Coefficient Estimation with Fixed Projection ===
# Given a fixed projection, fits regression coefficients via weighted least squares.
# Input:
#   - gamma0: A fixed linear projection.
#   - Y, X: List of response matrices and design matrix.
#   - Sigma, Tvec: Precomputed covariance matrices and sample sizes.
# Output:
#   - Estimated coefficients and optionally projection scores.

MatReg_QC_beta <- function(Y, X, Sigma, Tvec, gamma0, max.itr = 1000, tol = 1e-4, score.return = FALSE) {
  n <- length(Y)
  q <- ncol(X)
  p <- ncol(Y[[1]])
  beta0 <- rep(0, q)
  
  if (score.return) {
    score <- .hdcov_score(Sigma, gamma0)
  }
  
  # Step 1: Construct response Z from fixed projection gamma0
  Z <- .hdcov_score(Sigma, gamma0)
  
  # Step 2: Filter valid entries (positive and finite)
  valid <- which(Z > 0 & is.finite(Z))
  Z_cal <- Z[valid]
  X_cal <- X[valid, -1, drop = FALSE]  # remove intercept column
  
  # Step 3: Fit beta using cross-validated Gamma GLM
  beta.new <- beta0
  if (length(Z_cal) > 5) {
    fit_glm <- try(cv.glmnet(X_cal, Z_cal, family = Gamma(link = "log"),
                             nfolds = 5, type.measure = "deviance",
                             intercept = TRUE, standardize = TRUE),
                   silent = TRUE)
    
    if (!inherits(fit_glm, "try-error")) {
      ind0 <- which(fit_glm$nzero > 1 & fit_glm$nzero < sqrt(ncol(X_cal) + 1))
      if (length(ind0) > 0) {
        ind <- which.min(fit_glm$cvm[ind0])
        beta.new[-1] <- as.matrix(fit_glm$glmnet.fit$beta[, ind0[ind]])
        beta.new[1] <- fit_glm$glmnet.fit$a0[ind0[ind]]
      }
    }
  }
  
  beta0 <- beta.new
  re <- list(beta = c(beta0), gamma = c(gamma0))
  if (score.return) re$score <- score
  return(re)
}

# === DfD Calculation ===
# Quantifies the average level of diagonal dominance in the projected covariance matrices
# for an increasing number of projections. Used as a criterion to determine the
# number of meaningful components in the model.
# Input:
#   - Y: A list of response matrices.
#   - Phi: A matrix of projections (columns).
# Output:
#   - A list containing:
#       * avg.level: Average of per-subject diagonal levels.
#       * sub.level: Subject-wise diagonal level matrix.

diag.level <-
  function(Y,Phi)
  {
    if(is.null(ncol(Phi))|ncol(Phi)==1)
    {
      stop("dimension of Phi is less than 2")
    }else
    {
      n<-length(Y)
      p<-ncol(Y[[1]])
      
      Tvec<-rep(NA,n)
      
      ps<-ncol(Phi)
      
      dl.sub<-matrix(NA,n,ps)
      colnames(dl.sub)<-paste0("Dim",1:ps)
      dl.sub[,1]<-1
      for(i in 1:n)
      {
        cov.tmp<-cov(Y[[i]])
        Tvec[i]<-nrow(Y[[i]])
        
        for(j in 2:ps)
        {
          phi.tmp<-Phi[,1:j]
          mat.tmp<-t(phi.tmp)%*%cov.tmp%*%phi.tmp
          dl.sub[i,j]<-det(diag(diag(mat.tmp)))/det(mat.tmp)
        }
      }
      
      pmean<-apply(dl.sub,2,function(y){return(prod(apply(cbind(y,Tvec),1,function(x){return(x[1]^(x[2]/sum(Tvec)))})))})
      
      re<-list(avg.level=pmean,sub.level=dl.sub)
      return(re)
    }
  }

# === Component Estimation with Adaptive Stopping ===
# Estimates one or more projection components using CAP. The procedure stops either when a fixed number
# of components (nD) is reached or when the DfD criterion exceeds a specified threshold.
# Input:
#   - Y, X: List of responses and covariate matrix.
#   - stop.crt: Criterion to stop adding components ("nD" or "DfD").
#   - nD: Number of projections to estimate (used if stop.crt = "nD").
#   - DfD.thred: Threshold for stopping based on DfD (used if stop.crt = "DfD").
#   - CAP.OC: Logical flag indicating whether to enforce orthogonality constraints.
#   - gamma0.mat, ninitial: Initialization for gamma (optional).
# Output:
#   - Estimated projections, regression coefficients, DfD statistics, and fitting times.

capReg <-
  function(Y,X,stop.crt=c("nD","DfD"),nD=2,DfD.thred=2,CAP.OC=FALSE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=FALSE,
           gamma0.mat=NULL,ninitial=NULL)
  {
    if(stop.crt[1]=="nD"&is.null(nD))
    {
      stop.crt<-"DfD"
    }
    
    n<-length(Y)
    q<-ncol(X)
    p<-ncol(Y[[1]])
    
    sig.res <- precompute_sigma(Y)
    Sigma <- sig.res$Sigma
    Tvec <- sig.res$Tvec
    
    if(is.null(colnames(X)))
    {
      colnames(X)<-paste0("X",1:q)
    }
    print("Component 1")
    tem <- proc.time()
    re1 <- MatReg_QC_opt(Y, X, Sigma, Tvec, method = c("CAP"), max.itr = max.itr, tol = tol,
                         trace = FALSE, score.return = score.return, gamma0.mat = gamma0.mat, ninitial = ninitial)
    gamma_times <- (proc.time() - tem)[3][[1]]
    if (trace) print(paste0("Down time ", gamma_times[1]))
    
    Phi.est<-matrix(re1$gamma,ncol=1)
    beta.est<-matrix(re1$beta,ncol=1)
    
    if(score.return)
    {
      score<-matrix(re1$score,ncol=1)
    }
    
    if(stop.crt[1]=="nD")
    {
      if (nD > 1) {
        for (j in 2:nD) {
          if (trace) print(paste0("Component ", j))
          tem <- proc.time()
          re.tmp<-NULL
          try(re.tmp <- MatReg_QC_opt2(Y, X, Phi0 = Phi.est, Sigma = Sigma, Tvec = Tvec, CAP.OC=CAP.OC,
                                       max.itr = max.itr, tol = tol, trace = trace,
                                       score.return = score.return, gamma0.mat = gamma0.mat, ninitial = ninitial))
          if (!is.null(re.tmp)) {
            gamma_times<-c(gamma_times,(proc.time() - tem)[3][[1]])
            if (trace) print(paste0("Down time ", gamma_times[j]))
            
            Phi.est <- cbind(Phi.est, re.tmp$gamma)
            beta.est <- cbind(beta.est, re.tmp$beta)
            if (score.return) score <- cbind(score, re.tmp$score)
          } else break
        }
      }
    }
    
    if(stop.crt[1]=="DfD")
    {
      nD<-1
      DfD.tmp<-1
      while(DfD.tmp<DfD.thred)
      {
        print(paste0("Component ", nD+1))
        tem <- proc.time()
        
        re.tmp<-NULL
        try(re.tmp <- MatReg_QC_opt2(Y, X, Phi0 = Phi.est, Sigma = Sigma, Tvec = Tvec, CAP.OC=CAP.OC,
                                     max.itr = max.itr, tol = tol, trace = trace,
                                     score.return = score.return, gamma0.mat = gamma0.mat, ninitial = ninitial))

        if(is.null(re.tmp)==FALSE)
        {
          nD<-nD+1
          
          DfD.out<-diag.level(Y,cbind(Phi.est,re.tmp$gamma))
          DfD.tmp<-DfD.out$avg.level[nD]
          
          if(DfD.tmp<DfD.thred)
          {
            gamma_times<-c(gamma_times,(proc.time() - tem)[3][[1]])
            if (trace) print(paste0("Down time ", gamma_times[nD]))
            
            Phi.est<-cbind(Phi.est,re.tmp$gamma)
            beta.est<-cbind(beta.est,re.tmp$beta)
            
            if(score.return)
            {
              score<-cbind(score,re.tmp$score)  
            }
          }
        }
      }
    }
    
    
    colnames(Phi.est)=colnames(beta.est)<-paste0("D",1:ncol(Phi.est))
    rownames(Phi.est)<-paste0("V",1:p)
    rownames(beta.est)<-colnames(X)
    
    if(ncol(Phi.est)>1)
    {
      DfD<-diag.level(Y,Phi.est)
    }else
    {
      DfD<-1
    }
    
    if(score.return)
    {
      colnames(score)<-paste0("D",1:ncol(Phi.est))
      re<-list(gamma=Phi.est,beta=beta.est,orthogonality=t(Phi.est)%*%Phi.est,DfD=DfD,score=score, gamma_times = gamma_times)
    }else
    {
      re<-list(gamma=Phi.est,beta=beta.est,orthogonality=t(Phi.est)%*%Phi.est,DfD=DfD, gamma_times = gamma_times)
    }
    
    return(re)
  }  


# === Single-Splitting Step for Selection and Unbiased Estimation ===
# Performs one data split: variable selection on one part and coefficient estimation on the other.
# Input:
#   - Xmat1, Ymat1: Covariates and list of responses.
#   - gamma1: Projection.
# Output:
#   - Estimated coefficients and selected support.

SPARE <- function(Xmat1, Ymat1, prop1 = 0.5, p01 = NA, gamma1) {
  n <- length(Ymat1)
  p <- ncol(Ymat1[[1]])
  q <- ncol(Xmat1)
  
  p01 <- ifelse(is.na(p01), sqrt(q), p01)
  
  # Precompute covariance matrices
  sigma.res <- precompute_sigma(Ymat1)
  Sigma <- sigma.res$Sigma
  Tvec <- sigma.res$Tvec
  
  # Step 1: Construct response Z and perform variable selection
  Z <- .hdcov_score(Sigma, gamma1)
  idx_calibrate <- which(Z > 0 & is.finite(Z))
  X_calibrate <- Xmat1[idx_calibrate, , drop = FALSE]
  Z_calibrate <- Z[idx_calibrate]
  n_calibrate <- length(idx_calibrate)
  
  # Randomly split data for selection and estimation
  n1 <- floor(n_calibrate * prop1)
  samp1 <- sort(sample(1:n_calibrate, n1))
  samp2 <- setdiff(1:n_calibrate, samp1)
  
  beta.set <- rep(0, q)
  set <- NULL
  fit_glm <- try(cv.glmnet(X_calibrate[samp2, -1, drop = FALSE], Z_calibrate[samp2],
                           family = Gamma(link = "log"), nfolds = 5, type.measure = "deviance",
                           intercept = TRUE, standardize = TRUE), silent = TRUE)
  if (!inherits(fit_glm, "try-error")) {
    ind0 <- which(fit_glm$nzero > 1 & fit_glm$nzero < p01)
    if (length(ind0) > 0) {
      ind <- which.min(fit_glm$cvm[ind0])
      beta.set[-1] <- as.vector(fit_glm$glmnet.fit$beta[, ind0[ind]])
      beta.set[1] <- fit_glm$glmnet.fit$a0[ind0[ind]]
      set <- union(1, which(beta.set[-1] != 0) + 1)
    }
  }
  
  # Fallback if variable selection fails
  if (is.null(which(beta.set[-1] != 0))) {
    warning("cv.glmnet failed, using fallback")
    set <- c(1, 10, 15, 20, 25, 30, 35)
  }
  
  # Step 2: Refit Gamma GLM on selected variables
  beta_hat <- rep(0, q)
  if (length(set) > 0) {
    Z_sub <- Z_calibrate[samp1]
    X_sub <- X_calibrate[samp1, set[-1], drop = FALSE]
    refit <- try(glm(Z_sub ~ X_sub, family = Gamma(link = "log")), silent = TRUE)
    if (inherits(refit, "try-error") || !refit$converged || any(is.na(coef(refit)))) {
      warning("Refit on set failed. Skipping this bootstrap replicate.")
      return(NA)
    }
    beta_refit <- coef(refit)
    beta_hat[set] <- beta_refit
  }
  
  # Step 3: One-step update for unselected variables
  unselected <- setdiff(1:q, set)
  estimate_single <- function(j) {
    idx <- union(j, set)
    X_sub <- X_calibrate[samp1, idx[-2], drop = FALSE]
    fit_j <- try(glm(Z_sub ~ X_sub, family = Gamma(link = "log")), silent = TRUE)
    if (inherits(fit_j, "try-error") || !fit_j$converged || is.na(coef(fit_j)[2])) return(NA)
    return(coef(fit_j)[2])
  }
  beta_hat[unselected] <- sapply(unselected, estimate_single)
  
  yicount <- integer(n)
  yicount[idx_calibrate[samp1]] <- 1
  return(list(beta.hat = beta_hat, sel.set = set, boot.ct = yicount))
}


# === Multi-Splitting Step for Inference Stabilization ===
# Repeats the single-splitting step multiple times and aggregates the results.
# Input:
#   - Xmat, Ymat: Design matrix and list of responses.
#   - B: Number of repetitions.
#   - gamma_ini: Initial projection.
# Output:
#   - Aggregated coefficient estimates, standard errors, p-values, and selection frequencies.

SSHDI <- function(Xmat, Ymat, B = 200, prop = 0.5, p0 = NA, gamma_ini) {
  n <- nrow(Xmat)
  q <- ncol(Xmat)
  beta0 <- rep(0, q)
  no_cores <- parallel::detectCores()
  registerDoParallel(cores = no_cores)
  fit1 <- foreach(i = 1:B, .errorhandling = "pass", .packages = c("MASS", "lbfgs", "glmnet"),
                  .export = c("SPARE", "precompute_sigma")) %dopar% {
                    SPARE(Xmat1 = Xmat, Ymat1 = Ymat, prop1 = prop, p01 = p0, gamma1 = gamma_ini)
                  }
  BETA_list <- lapply(1:B, function(i) tryCatch(fit1[[i]]$beta.hat, error = function(e) NULL))
  valid_idx <- which(sapply(BETA_list, function(x) is.numeric(x) && !is.null(x)))
  if (length(valid_idx) == 0) stop("No valid bootstrap replicates")
  
  BETA <- do.call(cbind, BETA_list[valid_idx])
  betam <- rowMeans(BETA, na.rm = TRUE)
  Ycount <- do.call(cbind, lapply(valid_idx, function(i) fit1[[i]]$boot.ct))
  
  # Estimate variance of the multi-split estimator
  est.var <- function(beta, Ycount) {
    n <- nrow(Ycount)
    n1 <- sum(Ycount[, 1])
    var1 <- (n / (n - n1))^2 * sum(cov(beta, t(Ycount), use = "pairwise.complete.obs")^2)
    var2 <- var1 - n1 / (n - n1) * n / ncol(Ycount) * var(beta, na.rm = TRUE)
    return(var2)
  }
  
  vars <- apply(BETA, 1, est.var, Ycount = Ycount)
  sds <- sqrt(ifelse(vars >= 0, vars, NA))
  pvs <- 2 * (1 - pnorm(abs(betam) / sds))
  nzero <- which(pvs < 0.05 / q & !is.na(pvs))
  
  sel_freq <- table(unlist(lapply(valid_idx, function(i) fit1[[i]]$sel.set))) / length(valid_idx)
  sel_vec <- rep(0, q)
  sel_vec[as.integer(names(sel_freq))] <- sel_freq
  
  return(list(ss.beta = betam, sd = sds, p = pvs, sel.freq = sel_vec, nzero = nzero))
}

# === HCAP Full Pipeline: Projection Estimation and Inference ===
# Runs the full HCAP pipeline including projection estimation and inference.
# Input:
#   - X: Design matrix (n x q)
#   - Y_list: List of response matrices
#   - stop.crt: Projection search strategy ('nD' or 'DfD')
#   - nD: Number of projections (if stop.crt = 'nD')
#   - DfD.thred: Threshold for DfD stopping (if stop.crt = 'DfD')
#   - B: Number of bootstrap replicates for inference
# Output:
#   - List containing estimated projections, coefficient inference results, and timing.

run_HCAP_pipeline <- function(X, Y_list, stop.crt = "nD", nD = 2, DfD.thred = 2, B = 100) {
  t0 <- proc.time()
  cap_res <- capReg(Y_list, X, stop.crt = stop.crt, nD = nD, DfD.thred = DfD.thred, CAP.OC = TRUE, trace = TRUE)
  gamma_est <- cap_res$gamma
  
  inference_list <- vector("list", ncol(gamma_est))
  for (i in 1:ncol(gamma_est)) {
    print("Inference")
    tem <- proc.time()
    inference_list[[i]] <- SSHDI(Xmat = X, Ymat = Y_list, B = B, gamma_ini = gamma_est[, i])
    inference_times <- (proc.time() - tem)[3][[1]]
    print(paste0("Down time ", inference_times[1]))
  }
  run_time <- (proc.time() - t0)[3][[1]]
  return(list(
    gamma_est = gamma_est,
    inference = inference_list,
    run_time = run_time,
    DfD = cap_res$DfD,
    orthogonality = cap_res$orthogonality
  ))
}

