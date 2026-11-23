#################################################
# Covariance-on-covariance regression
#
# V3: RcppArmadillo-accelerated update of V2. Same interface and outputs; the hot
# loops (cov.ls, cov.sk.x, cov.sk.y, eigen.solve) call a compiled kernel. See
# CAP-CoC/README.md.
#
# V4 (see CAP-CoC/README.md): the changes of HDCAP V8 and V9 that apply to CoC.
# (1) The objective that chooses among the initial values (obj.func, in
#     COCReg.D1.opt and COCReg.D1) is the function the fit minimizes,
#     mean_i (log(gamma'Sigma_i gamma) - alpha log(theta'Delta_i theta) - W_i'beta)^2
#     with gamma'Hy gamma = theta'Hx theta = 1, and the covariances of the fit:
#     with outcome shrinkage (cov.shrinkage.y, forced when min(T_i) - 5 < q)
#     Sigma_i = cov.sk.y(Y, gamma, kappa), the shrinkage of the fit at
#     kappa_i = exp(alpha log(theta'Delta_i theta) + W_i'beta) (the covariances
#     COCReg.coef and COCReg.D1.base iterate on); without, the sample covariance.
#     Up to V3 it used the per-subject linear shrinkage cov.ls (the starting
#     covariance of the iterations) with outcome shrinkage, a function the
#     iteration does not minimize.  Delta_i: cov.sk.x or the sample covariance,
#     as in the fit (unchanged).
# (2) Initial values of beta, argument beta0 of COCReg, COCReg.Dk, COCReg.D1.opt,
#     COCReg.D1, COCReg.D1.base, COCReg.coef, COCReg.coef.asmp and
#     COCReg.coef.boot: a vector of length ncol(W) (one initial value) or a matrix
#     with ncol(W) rows (one initial value per column); the initial alpha is 0.
#     beta0 = NULL (default) is two initial values, 0 and coc_beta0_v3() (the
#     V2-V3 start, set.seed(100); c(0, rnorm(r-1, 1, 1)), now drawn without
#     changing the caller's random numbers).  COCReg.D1.opt runs every initial
#     value of (gamma, theta) from every initial value of beta and keeps the
#     direction with the lowest objective; every (alpha, beta) fit (COCReg.coef:
#     the final coefficients of a direction, the objective of each candidate,
#     each bootstrap draw) keeps the converged estimate with the lowest objective.
#     With beta0 = coc_beta0_v3(ncol(W)) V4 runs the iterations of V3.
# (3) As HDCAP V8: the kernel no longer aborts R on a failed iteration (an R
#     error, caught by try()), one kernel per R version and architecture
#     (build_kernels.R), every .Call names the V4 kernel, the folder is found
#     from the file source() reads, and COCReg.D1.opt skips an initial value
#     whose coefficients or objective cannot be computed (V3 kept the coef.tmp of
#     the previous initial value) and stops with a clear error when none can.
#################################################


#################################################
# Load the compiled kernel (coc_kernels.cpp -> shared object), building it once
# if necessary.  The kernel file is named after the R version and architecture
# (coc_kernels_v4-R<major.minor>-<arch>.so, see build_kernels.R), so that
# machines sharing this folder with different R versions each load their own build.
# The folder: option "coc.v4.dir" or variable COC_V4_DIR if set, else the folder
# of this file when it is read with source(), else "CAP-CoC/V4", "V4" or ".".

# .coc_v4_dll: name of the loaded V4 kernel; every .Call() names it (PACKAGE), so
# that a V3 kernel loaded in the same session (same symbol names, but one that
# aborts R on a failed iteration) is never used by mistake
.coc_v4_dll <- "CAPmethods"   # CAPmethods: the kernels are registered in the package DLL

# The V2-V3 initial value of beta: 0 for the intercept, N(1,1) draws for the other
# coefficients under set.seed(100) (the start of HDCAP V5-V8, cap_beta0_v5).  V4
# draws it without changing the random number stream of the caller (V2-V3 left
# it at set.seed(100)).
coc_beta0_v3 <- function(r)
{
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  })
  set.seed(100)
  c(0, rnorm(r - 1, mean = 1, sd = 1))
}

# initial values of beta as an r x m matrix, one column per initial value:
#   beta0 = NULL (default): two, 0 and coc_beta0_v3(r)
#   a vector of length r: that one;  a matrix with r rows: its columns
.coc_beta0 <- function(r, beta0 = NULL)
{
  if (is.null(beta0)) return(cbind(rep(0, r), coc_beta0_v3(r)))
  B0 <- if (is.matrix(beta0)) beta0 else matrix(as.numeric(beta0), ncol = 1)
  if (!is.numeric(B0) || nrow(B0) != r || ncol(B0) < 1 || any(!is.finite(B0)))
    stop("beta0 must be a finite numeric vector of length ncol(W) = ", r,
         ", or a matrix with ", r, " rows (one initial value per column).")
  storage.mode(B0) <- "double"
  B0
}

# objective of a fit from its scores: mean_i (log sy_i - alpha log sx_i - W_i'beta)^2
# (obj.func with the covariances of the fit, whose scores COCReg.coef returns)
.coc_obj_scores <- function(sy, sx, W, alpha, beta)
{
  mean((log(sy) - alpha * log(sx) - c(W %*% beta))^2, na.rm = TRUE)
}

# objective of a direction, as COCReg.D1.opt compares its initial values: gamma and
# theta scaled to gamma'Hy gamma = theta'Hx theta = 1, (alpha, beta) from every
# initial value of beta (COCReg.coef), obj.func with the covariances of the fit
.coc_dir_obj <- function(Y, X, W, gamma, theta, Hy, Hx, Sy, Sx, cov.shrinkage.y, cov.shrinkage.x,
                         B0, max.itr, tol)
{
  gamma.unscale <- c(gamma) / sqrt(t(gamma) %*% Hy %*% gamma)[1, 1]
  theta.unscale <- c(theta) / sqrt(t(theta) %*% Hx %*% theta)[1, 1]
  coef.tmp <- COCReg.coef(Y, X, W, gamma = gamma.unscale, theta = theta.unscale,
                          cov.shrinkage.y = cov.shrinkage.y, cov.shrinkage.x = cov.shrinkage.x,
                          max.itr = max.itr, tol = tol, trace = FALSE, score.return = FALSE, beta0 = B0)
  obj.func(Y, X, W, Sigma = Sy, Delta = Sx, gamma = gamma.unscale, theta = theta.unscale,
           alpha = coef.tmp$alpha, beta = coef.tmp$beta, cov.shrinkage.y = cov.shrinkage.y)
}

#################################################
# standardized Frobenius norm
norm.F.std<-function(A1,A2=NULL)
{
  p<-nrow(A1)
  
  if(is.null(A2))
  {
    return(sqrt(sum(diag(A1%*%t(A1)))/p))
  }else
  {
    return(sum(diag(A1%*%t(A2)))/p)
  }
}
#################################################

#################################################
# objective function
# V4: argument cov.shrinkage.y.  TRUE: Sigma is replaced by the covariances the fit
# iterates on with outcome shrinkage, cov.sk.y(Y, gamma, kappa) at
# kappa_i = exp(alpha log(theta'Delta_i theta) + W_i'beta), so the value is the
# function the fit minimizes at (gamma, theta, alpha, beta) (COCReg.D1.opt,
# COCReg.D1).  FALSE (default): Sigma as given, as up to V3 (COCReg.D1.base
# monitors its iterations with the shrunk Sigma it holds).
obj.func<-function(Y,X,W,Sigma,Delta,gamma,theta,alpha,beta,cov.shrinkage.y=FALSE)
{
  # Y: list of outcome data
  # X: list of predictor data
  # W: covariates
  # Sigma: covariates of Y
  # Delta: covariates of X
  # gamma: linear projection on Y
  # theta: linear projection on X
  # alpha: model coefficient of log(theta'*Delta*theta)
  # beta: model coefficient of W

  # projected scores via the kernel, then vectorized objective
  sx<-as.numeric(.Call("coc_score_cpp",Delta,as.numeric(theta),PACKAGE=.coc_v4_dll))
  if(cov.shrinkage.y)
  {
    # covariance estimate of y of the fit (as COCReg.coef / COCReg.D1.base)
    kappa.tmp<-exp(alpha*log(sx)+c(W%*%beta))
    Sigma<-cov.sk.y(Y,gamma,kappa.tmp)
  }
  sy<-as.numeric(.Call("coc_score_cpp",Sigma,as.numeric(gamma),PACKAGE=.coc_v4_dll))
  lv<-(log(sy)-alpha*log(sx)-c(W%*%beta))^2
  re<-mean(lv,na.rm=TRUE)
  return(re)
}
#################################################

#################################################
# shrinkage estimator of the covariance matrices

# Delta
cov.sk.x<-function(X)
{
  # X: list of p dimensional data  (RcppArmadillo kernel)
  .Call("coc_cov_sk_x_cpp", lapply(X, as.matrix),PACKAGE=.coc_v4_dll)
}

# Sigma
cov.sk.y<-function(Y,gamma,kappa.vec)
{
  # Y: list of data; gamma: projection; kappa: eigenvalue means  (kernel)
  .Call("coc_cov_sk_y_cpp", lapply(Y, as.matrix), as.numeric(gamma), as.numeric(kappa.vec),PACKAGE=.coc_v4_dll)
}

# linear shrinkage of a n by p data matrix
cov.ls<-function(X)
{
  # X: n by p data matrix  (RcppArmadillo kernel)
  matrix(.Call("coc_cov_ls_cpp", as.matrix(X),PACKAGE=.coc_v4_dll), nrow = ncol(X))
}
#################################################

#################################################
# given gamma and theta, estimate alpha and beta

COCReg.coef<-function(Y,X,W,gamma,theta,cov.shrinkage.y=TRUE,cov.shrinkage.x=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,beta0=NULL)
{
  # Y: list of outcome data
  # X: list of predictor data
  # W: covariates
  # gamma: linear projection on Y
  # theta: linear projection on X
  # beta0 (V4): initial value(s) of beta (alpha starts at 0): a vector of length ncol(W), or a matrix with
  #        ncol(W) rows; NULL = two, 0 and coc_beta0_v3(ncol(W)).  Several: one fit per initial value, the
  #        converged estimate with the lowest objective (a non-converged one only when none converged)

  B0<-.coc_beta0(ncol(W),beta0)
  if(ncol(B0)>1)
  {
    best<-NULL
    for(jj in 1:ncol(B0))
    {
      re.jj<-NULL
      try(re.jj<-COCReg.coef(Y,X,W,gamma,theta,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,max.itr=max.itr,tol=tol,
                             trace=trace,score.return=TRUE,beta0=B0[,jj]))
      if(is.null(re.jj)) next
      if(any(!is.finite(c(re.jj$alpha,re.jj$beta)))) next
      # objective with the covariances of the fit (its scores: shrinkage at its alpha, beta)
      ob<-.coc_obj_scores(re.jj$score.y,re.jj$score.x,W,re.jj$alpha,re.jj$beta)
      if(!is.finite(ob)) ob<-Inf
      if(is.null(best)||(re.jj$convergence&&!best$re$convergence)||
         (re.jj$convergence==best$re$convergence&&ob<best$obj))
      {
        best<-list(re=re.jj,obj=ob)
      }
    }
    if(is.null(best))
    {
      stop("COCReg.coef: alpha and beta could not be estimated from any initial value of beta.")
    }
    re<-best$re
    if(!score.return)
    {
      re$score.y<-NULL
      re$score.x<-NULL
    }
    return(re)
  }

  n<-length(Y)
  p<-ncol(X[[1]])
  q<-ncol(Y[[1]])
  r<-ncol(W)
  
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-c("Intercept",paste0("W",1:(r-1)))
  }
  
  nxvec<-sapply(X,nrow)
  nyvec<-sapply(Y,nrow)
  
  if(min(nxvec)-5<p)
  {
    cov.shrinkage.x<-TRUE
  }
  if(min(nyvec)-5<q)
  {
    cov.shrinkage.y<-TRUE
  }
  
  # covariance estimate of X
  if(cov.shrinkage.x)
  {
    # shrinkage estimator
    Sx<-cov.sk.x(X)
  }else
  {
    # sample covariance matrix
    Sx<-array(NA,c(p,p,n))
    for(i in 1:n)
    {
      Xtmp<-scale(X[[i]],center=TRUE,scale=FALSE)
      Sx[,,i]<-t(Xtmp)%*%Xtmp/nxvec[i]
    }
  }
  
  # covariance estimate of Y
  Sy<-array(NA,c(q,q,n))
  for(i in 1:n)
  {
    if(cov.shrinkage.y)
    {
      Sy[,,i]<-cov.ls(Y[[i]])
    }else
    {
      Ytmp<-scale(Y[[i]],center=TRUE,scale=FALSE)
      Sy[,,i]<-t(Ytmp)%*%Ytmp/nyvec[i]
    }
  }
  
  # # H matrices
  # Hx<-matrix(0,p,p)
  # Hy<-matrix(0,q,q)
  # for(i in 1:n)
  # {
  #   Hx<-Hx+Sx[,,i]*nxvec[i]/sum(nxvec)
  #   Hy<-Hy+Sy[,,i]*nyvec[i]/sum(nyvec)
  # }
  
  score.x<-as.numeric(.Call("coc_score_cpp",Sx,as.numeric(theta),PACKAGE=.coc_v4_dll))
  score.y<-as.numeric(.Call("coc_score_cpp",Sy,as.numeric(gamma),PACKAGE=.coc_v4_dll))

  # V4: the initial value of beta given (V2-V3: set.seed(100); c(0,rnorm(r-1,mean=1,sd=1)),
  # i.e. coc_beta0_v3(r)); alpha starts at 0
  alpha0<-0
  beta0<-B0[,1]

  if(trace)
  {
    alpha.trace<-NULL
    beta.trace<-NULL

    obj<-NULL
  }
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update alpha
    alpha.new<-mean(log(score.y)*log(score.x)-W%*%beta0*log(score.x))/mean((log(score.x))^2)
    
    # update beta
    otmp<-rep(0,r)
    for(i in 1:n)
    {
      otmp<-otmp+(log(score.y[i])-alpha.new*log(score.x[i]))*W[i,]/n
    }
    beta.new<-ginv(t(W)%*%W/n)%*%otmp
    
    if(cov.shrinkage.y)
    {
      # update covariance estimate of y
      kappa.tmp<-rep(NA,n)
      for(i in 1:n)
      {
        kappa.tmp[i]<-exp(alpha.new*log(score.x[i])+t(beta.new)%*%W[i,])
      }
      Sy<-cov.sk.y(Y,gamma,kappa.tmp)
      score.y<-as.numeric(.Call("coc_score_cpp",Sy,as.numeric(gamma),PACKAGE=.coc_v4_dll))
    }
    
    diff<-max(c(abs(beta.new-beta0),abs(alpha.new-alpha0)))
    
    # update parameters
    beta0<-beta.new
    alpha0<-alpha.new
    
    if(trace)
    {
      beta.trace<-cbind(beta.trace,beta.new)
      alpha.trace<-c(alpha.trace,alpha.new)
      
      obj<-c(obj,obj.func(Y,X,W,Sy,Sx,gamma,theta,alpha.new,beta.new))
    }
    
    # print(diff)
  }
  
  # results
  alpha.est<-alpha.new
  names(alpha.est)<-"X"
  
  beta.est<-c(beta.new)
  names(beta.est)<-colnames(W)
  
  if(trace)
  {
    rownames(beta.trace)<-colnames(W)
    if(score.return)
    {
      re<-list(gamma=gamma,theta=theta,alpha=alpha.est,beta=beta.est,convergence=(s<max.itr),score.y=score.y,score.x=score.x,alpha.trace=alpha.trace,beta.trace=beta.trace,obj.trace=obj)
    }else
    {
      re<-list(gamma=gamma,theta=theta,alpha=alpha.est,beta=beta.est,convergence=(s<max.itr),alpha.trace=alpha.trace,beta.trace=beta.trace,obj.trace=obj)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma=gamma,theta=theta,alpha=alpha.est,beta=beta.est,convergence=(s<max.itr),score.y=score.y,score.x=score.x)
    }else
    {
      re<-list(gamma=gamma,theta=theta,alpha=alpha.est,beta=beta.est,convergence=(s<max.itr))
    }
  }
  
  return(re)
}

COCReg.coef.asmp<-function(Y,X,W,gamma,theta,conf.level=conf.level,cov.shrinkage.y=TRUE,cov.shrinkage.x=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,beta0=NULL)
{
  # Y: list of outcome data
  # X: list of predictor data
  # W: covariates
  # gamma: linear projection on Y
  # theta: linear projection on X
  
  n<-length(Y)
  p<-ncol(X[[1]])
  q<-ncol(Y[[1]])
  r<-ncol(W)
  
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-c("Intercept",paste0("W",1:(r-1)))
  }
  
  nxvec<-sapply(X,nrow)
  nyvec<-sapply(Y,nrow)
  
  if(min(nxvec)-5<p)
  {
    cov.shrinkage.x<-TRUE
  }
  if(min(nyvec)-5<q)
  {
    cov.shrinkage.y<-TRUE
  }
  
  # covariance estimate of X
  if(cov.shrinkage.x)
  {
    # shrinkage estimator
    Sx<-cov.sk.x(X)
  }else
  {
    # sample covariance matrix
    Sx<-array(NA,c(p,p,n))
    for(i in 1:n)
    {
      Xtmp<-scale(X[[i]],center=TRUE,scale=FALSE)
      Sx[,,i]<-t(Xtmp)%*%Xtmp/nxvec[i]
    }
  }
  
  # covariance estimate of Y
  Sy<-array(NA,c(q,q,n))
  for(i in 1:n)
  {
    if(cov.shrinkage.y)
    {
      Sy[,,i]<-cov.ls(Y[[i]])
    }else
    {
      Ytmp<-scale(Y[[i]],center=TRUE,scale=FALSE)
      Sy[,,i]<-t(Ytmp)%*%Ytmp/nyvec[i]
    }
  }
  
  # # H matrices
  # Hx<-matrix(0,p,p)
  # Hy<-matrix(0,q,q)
  # for(i in 1:n)
  # {
  #   Hx<-Hx+Sx[,,i]*nxvec[i]/sum(nxvec)
  #   Hy<-Hy+Sy[,,i]*nyvec[i]/sum(nyvec)
  # }
  
  score.x<-as.numeric(.Call("coc_score_cpp",Sx,as.numeric(theta),PACKAGE=.coc_v4_dll))
  score.y<-as.numeric(.Call("coc_score_cpp",Sy,as.numeric(gamma),PACKAGE=.coc_v4_dll))
  
  out.coef<-COCReg.coef(Y,X,W,gamma,theta,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,beta0=beta0)
  
  #-------------------
  # asymptotic inference
  Gx<-mean((log(score.x))^2,na.rm=TRUE)
  Qw<-t(W)%*%W/n
  Hxw<-apply(t(apply(cbind(log(score.x),W),1,function(x){return(x[1]*x[-1])})),2,mean,na.rm=TRUE)
  Fs.info<-rbind(cbind(Gx,t(Hxw)),cbind(Hxw,Qw))
  cov.asmp<-ginv(Fs.info)/sum(nyvec)
  colnames(cov.asmp)=rownames(cov.asmp)<-c("alpha",colnames(W))
  coef.se<-sqrt(diag(cov.asmp))
  #-------------------
  
  zvalue<-qnorm(1-(1-conf.level)/2)
  
  alpha.est<-out.coef$alpha
  alpha.se<-coef.se[1]
  alpha.stat<-alpha.est/alpha.se
  alpha.pv<-(1-pnorm(abs(alpha.stat)))*2
  alpha.ci<-c(alpha.est-zvalue*alpha.se,alpha.est+zvalue*alpha.se)
  
  beta.est<-out.coef$beta
  beta.se<-coef.se[-1]
  beta.stat<-beta.est/beta.se
  beta.pv<-(1-pnorm(abs(beta.stat)))*2
  beta.ci<-cbind(LB=beta.est-zvalue*beta.se,UB=beta.est+zvalue*beta.se)
  
  re.alpha<-data.frame(Estimate=alpha.est,SE=alpha.se,statistic=alpha.stat,pvalue=alpha.pv,LB=alpha.ci[1],UB=alpha.ci[2])
  rownames(re.alpha)<-"X"
  re.beta<-data.frame(Estimate=beta.est,SE=beta.se,statistic=beta.stat,pvalue=beta.pv,LB=beta.ci[,1],UB=beta.ci[,2])
  rownames(re.beta)<-colnames(W)
  
  re<-list(alpha=re.alpha,beta=re.beta,cov.asmp=cov.asmp)
  
  return(re)
}

COCReg.coef.lm<-function(Y,X,W,gamma,theta,conf.level=conf.level,cov.shrinkage.y=TRUE,cov.shrinkage.x=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE)
{
  # Y: list of outcome data
  # X: list of predictor data
  # W: covariates
  # gamma: linear projection on Y
  # theta: linear projection on X
  
  n<-length(Y)
  p<-ncol(X[[1]])
  q<-ncol(Y[[1]])
  r<-ncol(W)
  
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-c("Intercept",paste0("W",1:(r-1)))
  }
  
  nxvec<-sapply(X,nrow)
  nyvec<-sapply(Y,nrow)
  
  if(min(nxvec)-5<p)
  {
    cov.shrinkage.x<-TRUE
  }
  if(min(nyvec)-5<q)
  {
    cov.shrinkage.y<-TRUE
  }
  
  # covariance estimate of X
  if(cov.shrinkage.x)
  {
    # shrinkage estimator
    Sx<-cov.sk.x(X)
  }else
  {
    # sample covariance matrix
    Sx<-array(NA,c(p,p,n))
    for(i in 1:n)
    {
      Xtmp<-scale(X[[i]],center=TRUE,scale=FALSE)
      Sx[,,i]<-t(Xtmp)%*%Xtmp/nxvec[i]
    }
  }
  
  # covariance estimate of Y
  Sy<-array(NA,c(q,q,n))
  for(i in 1:n)
  {
    if(cov.shrinkage.y)
    {
      Sy[,,i]<-cov.ls(Y[[i]])
    }else
    {
      Ytmp<-scale(Y[[i]],center=TRUE,scale=FALSE)
      Sy[,,i]<-t(Ytmp)%*%Ytmp/nyvec[i]
    }
  }
  
  # # H matrices
  # Hx<-matrix(0,p,p)
  # Hy<-matrix(0,q,q)
  # for(i in 1:n)
  # {
  #   Hx<-Hx+Sx[,,i]*nxvec[i]/sum(nxvec)
  #   Hy<-Hy+Sy[,,i]*nyvec[i]/sum(nyvec)
  # }
  
  score.x<-as.numeric(.Call("coc_score_cpp",Sx,as.numeric(theta),PACKAGE=.coc_v4_dll))
  score.y<-as.numeric(.Call("coc_score_cpp",Sy,as.numeric(gamma),PACKAGE=.coc_v4_dll))
  
  dtmp<-data.frame(Y=log(score.y),X=log(score.x),W)
  fit<-lm(Y~0+.,data=dtmp)
  fit.confint<-confint(fit,level=conf.level)
  
  alpha.est<-fit$coefficients[1]
  alpha.se<-summary(fit)$coefficients["X",2]
  alpha.stat<-summary(fit)$coefficients["X",3]
  alpha.pv<-summary(fit)$coefficients["X",4]
  alpha.ci<-fit.confint["X",]
  
  beta.est<-fit$coefficients[-1]
  beta.se<-summary(fit)$coefficients[-1,2]
  beta.stat<-summary(fit)$coefficients[-1,3]
  beta.pv<-summary(fit)$coefficients[-1,4]
  beta.ci<-fit.confint[-1,]
  
  re.alpha<-data.frame(Estimate=alpha.est,SE=alpha.se,statistic=alpha.stat,pvalue=alpha.pv,LB=alpha.ci[1],UB=alpha.ci[2])
  rownames(re.alpha)<-"X"
  re.beta<-data.frame(Estimate=beta.est,SE=beta.se,statistic=beta.stat,pvalue=beta.pv,LB=beta.ci[,1],UB=beta.ci[,2])
  rownames(re.beta)<-colnames(W)
  
  re<-list(alpha=re.alpha,beta=re.beta,lm.out=fit)
  
  return(re)
}
#################################################

#################################################
# eigenvectors and eigenvalues of A with respect to H
# H positive definite and symmetric

eigen.solve<-function(A,H)
{
  # smallest generalized eigenvector of A w.r.t. H  (RcppArmadillo kernel)
  as.numeric(.Call("coc_eigen_solve_cpp", as.matrix(A), as.matrix(H),PACKAGE=.coc_v4_dll))
}
#################################################

#################################################
# 3 versions
# estimate gamma, theta, alpha, beta
# given theta: estimate gamma, alpha, beta
# given gamma: estimate theta, alpha, beta

COCReg.D1.base<-function(Y,X,W,gamma=NULL,theta=NULL,Hy=NULL,Hx=NULL,cov.shrinkage.y=TRUE,cov.shrinkage.x=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,gamma0=NULL,theta0=NULL,burn.in=1000,beta0=NULL)
{
  # Y: list of outcome data
  # X: list of predictor data
  # W: covariates
  
  n<-length(Y)
  p<-ncol(X[[1]])
  q<-ncol(Y[[1]])
  r<-ncol(W)
  
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-c("Intercept",paste0("W",1:(r-1)))
  }
  
  nxvec<-sapply(X,nrow)
  nyvec<-sapply(Y,nrow)
  
  if(min(nxvec)-5<p)
  {
    cov.shrinkage.x<-TRUE
  }
  if(min(nyvec)-5<q)
  {
    cov.shrinkage.y<-TRUE
  }
  
  # covariance estimate of X
  if(cov.shrinkage.x)
  {
    # shrinkage estimator
    Sx<-cov.sk.x(X)
  }else
  {
    # sample covariance matrix
    Sx<-array(NA,c(p,p,n))
    for(i in 1:n)
    {
      Xtmp<-scale(X[[i]],center=TRUE,scale=FALSE)
      Sx[,,i]<-t(Xtmp)%*%Xtmp/nxvec[i]
    }
  }
  
  # covariance estimate of Y
  Sy<-array(NA,c(q,q,n))
  for(i in 1:n)
  {
    if(cov.shrinkage.y)
    {
      Sy[,,i]<-cov.ls(Y[[i]])
    }else
    {
      Ytmp<-scale(Y[[i]],center=TRUE,scale=FALSE)
      Sy[,,i]<-t(Ytmp)%*%Ytmp/nyvec[i]
    }
    
    # print(i)
  }
  
  # H matrices
  if(is.null(Hx))
  {
    Hx<-matrix(0,p,p)
    for(i in 1:n)
    {
      Hx<-Hx+Sx[,,i]*nxvec[i]/sum(nxvec)
    }
  }
  if(is.null(Hy))
  {
    Hy<-matrix(0,q,q)
    for(i in 1:n)
    {
      Hy<-Hy+Sy[,,i]*nyvec[i]/sum(nyvec)
    }
  }
  
  # V4: initial value(s) of beta (r x m matrix).  Several: one run per initial value
  # of beta, the direction with the lowest objective, compared as in COCReg.D1.opt
  # (.coc_dir_obj: gamma'Hy gamma = theta'Hx theta = 1, alpha/beta from every initial
  # value, the covariances of the fit)
  B0<-.coc_beta0(r,beta0)
  if(ncol(B0)>1)
  {
    best<-NULL
    best.obj<-Inf
    for(jj in 1:ncol(B0))
    {
      re.jj<-NULL
      try(re.jj<-COCReg.D1.base(Y,X,W,gamma=gamma,theta=theta,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,
                                max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0=gamma0,theta0=theta0,burn.in=burn.in,
                                beta0=B0[,jj]))
      if(is.null(re.jj)==FALSE)
      {
        ob<-NA
        try(ob<-.coc_dir_obj(Y,X,W,re.jj$gamma,re.jj$theta,Hy,Hx,Sy,Sx,cov.shrinkage.y,cov.shrinkage.x,B0,max.itr,tol))
        if(is.finite(ob)&&ob<best.obj)
        {
          best<-re.jj
          best.obj<-ob
        }
      }
    }
    if(is.null(best))
    {
      stop("COCReg.D1.base: no initial value of beta gave an estimate.")
    }
    return(best)
  }

  # one initial value of beta (V2-V3: set.seed(100); c(0,rnorm(r-1,mean=1,sd=1)), i.e.
  # coc_beta0_v3(r)); alpha starts at 0
  alpha0<-0
  beta0<-B0[,1]
  if(is.null(gamma0)&is.null(gamma))
  {
    gamma0<-rep(1/sqrt(q),q) 
  }else
    if(is.null(gamma0))
    {
      gamma0<-gamma
    }
  if(is.null(theta0)&is.null(theta))
  {
    theta0<-rep(1/sqrt(p),p)
  }else
    if(is.null(theta0))
    {
      theta0<-theta
    }
  
  score.x<-as.numeric(.Call("coc_score_cpp",Sx,as.numeric(theta0),PACKAGE=.coc_v4_dll))
  score.y<-as.numeric(.Call("coc_score_cpp",Sy,as.numeric(gamma0),PACKAGE=.coc_v4_dll))
  
  if(trace)
  {
    alpha.trace<-NULL
    beta.trace<-NULL
    
    gamma.trace<-NULL
    theta.trace<-NULL
    
    obj<-NULL
  }
  
  obj0<-obj.func(Y,X,W,Sy,Sx,gamma0,theta0,alpha0,beta0)
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update alpha
    alpha.new<-mean(log(score.y)*log(score.x)-W%*%beta0*log(score.x))/mean((log(score.x))^2)
    
    # update beta
    otmp<-rep(0,r)
    for(i in 1:n)
    {
      otmp<-otmp+(log(score.y[i])-alpha.new*log(score.x[i]))*W[i,]/n
    }
    beta.new<-ginv(t(W)%*%W/n)%*%otmp
    
    if(is.null(theta))
    {
      # update theta
      Vtmp<-log(score.y)-c(W%*%beta.new)
      w2<-((alpha.new*log(score.x)-Vtmp)/score.x)*(2*alpha.new/n)
      A2<-.Call("coc_accum_cpp",Sx,as.numeric(w2),PACKAGE=.coc_v4_dll)
      theta.new<-c(eigen.solve(A2,Hx))
      # theta.new<-theta.new/sqrt(sum(theta.new^2))
      score.x<-as.numeric(.Call("coc_score_cpp",Sx,as.numeric(theta.new),PACKAGE=.coc_v4_dll))
    }else
    {
      theta.new<-theta
    }
    
    if(is.null(gamma))
    {
      # update gamma
      Utmp<-alpha.new*log(score.x)+c(W%*%beta.new)
      if(cov.shrinkage.y)
      {
        # update covariance estimate of y
        kappa.tmp<-rep(NA,n)
        for(i in 1:n)
        {
          kappa.tmp[i]<-exp(alpha.new*log(score.x[i])+t(beta.new)%*%W[i,])
        }
        Sy<-cov.sk.y(Y,gamma0,kappa.tmp)
        
        # Hy<-matrix(0,q,q)
        # for(i in 1:n)
        # {
        #   Hy<-Hy+Sy[,,i]*nyvec[i]/sum(nyvec)
        # }
      }
      w1<-((log(score.y)-Utmp)/score.y)*(2/n)
      A1<-.Call("coc_accum_cpp",Sy,as.numeric(w1),PACKAGE=.coc_v4_dll)
      gamma.new<-c(eigen.solve(A1,Hy))
      # gamma.new<-gamma.new/sqrt(sum(gamma.new^2))
      score.y<-as.numeric(.Call("coc_score_cpp",Sy,as.numeric(gamma.new),PACKAGE=.coc_v4_dll))
    }else
    {
      gamma.new<-gamma
    }
    
    diff1<-max(c(abs(gamma.new-gamma0),abs(theta.new-theta0)))
    diff2<-max(c(abs(beta.new-beta0),abs(alpha.new-alpha0)))
    diff<-max(diff1,diff2)
    
    obj.new<-obj.func(Y,X,W,Sy,Sx,gamma.new,theta.new,alpha.new,beta.new)
    
    if(obj.new>obj0&s>burn.in)
    {
      break
    }else
    {
      obj0<-obj.new
      
      # update parameters
      beta0<-beta.new
      alpha0<-alpha.new
      gamma0<-gamma.new
      theta0<-theta.new
      
      if(trace)
      {
        beta.trace<-cbind(beta.trace,beta.new)
        alpha.trace<-c(alpha.trace,alpha.new)
        
        gamma.trace<-cbind(gamma.trace,gamma.new)
        theta.trace<-cbind(theta.trace,theta.new)
        
        obj<-c(obj,obj.func(Y,X,W,Sy,Sx,gamma.new,theta.new,alpha.new,beta.new))
      }
    }
    
    # print(c(diff1,diff2,diff))
  }
  
  #---------------------------------------------
  # results
  gamma.est<-gamma0/sqrt(sum(gamma0^2))
  if(gamma.est[which.max(abs(gamma.est))]<0)
  {
    gamma.est<--gamma.est
  }
  theta.est<-theta0/sqrt(sum(theta0^2))
  if(theta.est[which.max(abs(theta.est))]<0)
  {
    theta.est<--theta.est
  }
  #---------------------------------------------
  
  otmp<-COCReg.coef(Y,X,W,gamma.est,theta.est,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,beta0=B0)

  if(trace)
  {
    rownames(beta.trace)<-colnames(W)
    if(score.return)
    {
      re<-list(gamma=gamma.est,theta=theta.est,alpha=otmp$alpha,beta=otmp$beta,convergence=(s<max.itr),score.y=otmp$score.y,score.x=otmp$score.x,
               alpha.trace=alpha.trace,beta.trace=beta.trace,gamma.trace=gamma.trace,theta.trace=theta.trace,obj.trace=obj)
    }else
    {
      re<-list(gamma=gamma.est,theta=theta.est,alpha=otmp$alpha,beta=otmp$beta,convergence=(s<max.itr),
               alpha.trace=alpha.trace,beta.trace=beta.trace,gamma.trace=gamma.trace,theta.trace=theta.trace,obj.trace=obj)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma=gamma.est,theta=theta.est,alpha=otmp$alpha,beta=otmp$beta,convergence=(s<max.itr),score.y=otmp$score.y,score.x=otmp$score.x)
    }else
    {
      re<-list(gamma=gamma.est,theta=theta.est,alpha=otmp$alpha,beta=otmp$beta,convergence=(s<max.itr))
    }
  }
  
  return(re)
}
COCReg.D1<-function(Y,X,W,gamma=NULL,theta=NULL,Hy=NULL,Hx=NULL,cov.shrinkage.y=TRUE,cov.shrinkage.x=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,gamma0=NULL,theta0=NULL,burn.in=1000,beta0=NULL)
{
  # Y: list of outcome data
  # X: list of predictor data
  # W: covariates
  # beta0 (V4): initial value(s) of beta of both COCReg.D1.base runs (each runs every initial value and keeps
  #        the direction with the lowest objective); NULL = two, 0 and coc_beta0_v3(ncol(W))

  if(is.null(gamma)&is.null(theta))
  {
    otmp<-COCReg.D1.base(Y,X,W,gamma=NULL,theta=NULL,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,
                         max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0=gamma0,theta0=theta0,burn.in=burn.in,beta0=beta0)
    # update theta estimate
    re<-COCReg.D1.base(Y,X,W,gamma=otmp$gamma,theta=NULL,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,
                       max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0=gamma0,theta0=theta0,burn.in=burn.in,beta0=beta0)
  }else
  {
    re<-COCReg.D1.base(Y,X,W,gamma=gamma,theta=theta,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,
                       max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0=gamma0,theta0=theta0,burn.in=burn.in,beta0=beta0)
  }
  return(re)
}

# try several initial value of gamma and theta and optimizer over the objective function
COCReg.D1.opt<-function(Y,X,W,gamma=NULL,theta=NULL,Hy=NULL,Hx=NULL,cov.shrinkage.y=TRUE,cov.shrinkage.x=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,burn.in=1000,
                        gamma0.mat=NULL,theta0.mat=NULL,ninitial=NULL,seed=100,beta0=NULL)
{
  # Y: list of outcome data
  # X: list of predictor data
  # W: covariates
  
  #----------------------------------
  # primary parameter setting/estimation
  n<-length(Y)
  p<-ncol(X[[1]])
  q<-ncol(Y[[1]])
  r<-ncol(W)
  
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-c("Intercept",paste0("W",1:(r-1)))
  }
  
  nxvec<-sapply(X,nrow)
  nyvec<-sapply(Y,nrow)
  
  if(min(nxvec)-5<p)
  {
    cov.shrinkage.x<-TRUE
  }
  if(min(nyvec)-5<q)
  {
    cov.shrinkage.y<-TRUE
  }
  
  # covariance estimate of X
  if(cov.shrinkage.x)
  {
    # shrinkage estimator
    Sx<-cov.sk.x(X)
  }else
  {
    # sample covariance matrix
    Sx<-array(NA,c(p,p,n))
    for(i in 1:n)
    {
      Xtmp<-scale(X[[i]],center=TRUE,scale=FALSE)
      Sx[,,i]<-t(Xtmp)%*%Xtmp/nxvec[i]
    }
  }
  Sx.bar<-apply(Sx,c(1,2),mean,na.rm=TRUE)
  Sx.bar.eigen<-eigen(Sx.bar)
  
  # covariance estimate of Y
  Sy<-array(NA,c(q,q,n))
  for(i in 1:n)
  {
    if(cov.shrinkage.y)
    {
      Sy[,,i]<-cov.ls(Y[[i]])
    }else
    {
      Ytmp<-scale(Y[[i]],center=TRUE,scale=FALSE)
      Sy[,,i]<-t(Ytmp)%*%Ytmp/nyvec[i]
    }
  }
  Sy.bar<-apply(Sy,c(1,2),mean,na.rm=TRUE)
  Sy.bar.eigen<-eigen(Sy.bar)

  # V3 fix: compute the H matrices here when not supplied. In V2 these were left
  # NULL and then used in the gamma/theta unscaling below (t(gamma)%*%Hy), which
  # crashed the default COCReg(Y,X,W) call. Defined as the sample-size-weighted
  # mean covariance, matching COCReg.D1.base.
  if(is.null(Hx))
  {
    Hx<-.Call("coc_accum_cpp",Sx,nxvec/sum(nxvec),PACKAGE=.coc_v4_dll)
  }
  if(is.null(Hy))
  {
    Hy<-.Call("coc_accum_cpp",Sy,nyvec/sum(nyvec),PACKAGE=.coc_v4_dll)
  }
  #----------------------------------

  #----------------------------------
  # decide the number of initial values going to try
  if(is.null(gamma0.mat)==FALSE&is.null(theta0.mat)==FALSE)
  {
    if(is.null(ninitial))
    {
      ninitial<-min(c(ncol(gamma0.mat),ncol(theta0.mat),10))
    }else
    {
      ninitial<-min(c(ncol(gamma0.mat),ncol(theta0.mat),ninitial))
    }
  }else
    if(is.null(gamma0.mat)==FALSE)
    {
      if(is.null(ninitial))
      {
        ninitial<-min(c(ncol(gamma0.mat),10))
      }else
      {
        ninitial<-min(c(ncol(gamma0.mat),ninitial))
      }
    }else
      if(is.null(theta0.mat)==FALSE)
      {
        if(is.null(ninitial))
        {
          ninitial<-min(c(ncol(theta0.mat),10))
        }else
        {
          ninitial<-min(c(ncol(theta0.mat),ninitial))
        }
      }else
      {
        if(is.null(ninitial))
        {
          ninitial<-min(p,q,10)
        }
      }
  
  # gamma initial
  if(is.null(gamma))
  {
    if(is.null(gamma0.mat))
    {
      # set.seed(seed)
      # gamma.tmp<-matrix(rnorm((max(q,ninitial)+1+5)*q,mean=0,sd=1),nrow=q)
      # gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      
      gamma0.mat<-Sy.bar.eigen$vectors[,1:ninitial,drop=FALSE]      # V4: drop=FALSE (ninitial = 1)
    }
    # set.seed(seed)
    # gamma0.mat<-matrix(gamma0.mat[,sort(sample(1:ncol(gamma0.mat),ninitial,replace=FALSE))],ncol=ninitial)
  }else
  {
    gamma0.mat<-matrix(rep(gamma,ninitial),nrow=q)
  }
  
  # theta initial
  if(is.null(theta))
  {
    if(is.null(theta0.mat))
    {
      # set.seed(seed)
      # theta.tmp<-matrix(rnorm((max(p,ninitial)+1+5)*p,mean=0,sd=1),nrow=p)
      # theta0.mat<-apply(theta.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      
      theta0.mat<-Sx.bar.eigen$vectors[,1:ninitial,drop=FALSE]      # V4: drop=FALSE (ninitial = 1)
    }
    # set.seed(seed)
    # theta0.mat<-matrix(theta0.mat[,sort(sample(1:ncol(theta0.mat),ninitial,replace=FALSE))],ncol=ninitial)
  }else
  {
    theta0.mat<-matrix(rep(theta,ninitial),nrow=p)
  }
  #----------------------------------
  
  #----------------------------------
  # try different initial values with the lowest objective function
  # V4: every initial value of (gamma, theta) from every initial value of beta (columns
  # of B0); candidate (kk-1)*m0+jj is the (gamma, theta) start kk with beta start jj
  B0<-.coc_beta0(r,beta0)
  m0<-ncol(B0)
  re.tmp<-vector("list",ninitial*m0)
  obj<-rep(NA,ninitial*m0)
  for(kk in 1:ninitial)
  {
    for(jj in 1:m0)
    {
      cc<-(kk-1)*m0+jj
      try(re.tmp[[cc]]<-COCReg.D1(Y,X,W,gamma=gamma,theta=theta,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,max.itr=max.itr,tol=tol,trace=trace,
                                  score.return=score.return,gamma0=gamma0.mat[,kk],theta0=theta0.mat[,kk],burn.in=burn.in,beta0=B0[,jj]))

      if(is.null(re.tmp[[cc]])==FALSE)
      {
        # objective at gamma'Hy gamma = theta'Hx theta = 1, alpha/beta from every initial value of
        # beta, the covariances of the fit (V4, .coc_dir_obj).  The coefficient fit and the objective
        # in one try() (V3 kept the coef.tmp of the previous initial value when this one's fit
        # failed): when they cannot be computed, obj[cc] stays NA and this candidate is not used.
        try(obj[cc]<-.coc_dir_obj(Y,X,W,re.tmp[[cc]]$gamma,re.tmp[[cc]]$theta,Hy,Hx,Sy,Sx,
                                  cov.shrinkage.y,cov.shrinkage.x,B0,max.itr,tol))
      }
    }
  }
  # V4: the candidates in increasing order of the objective, the first one whose final alpha
  # and beta (COCReg.coef from every initial value of beta) can be estimated; with one
  # initial value of beta these are the candidate's own (as V3)
  ord<-order(obj,na.last=NA)
  if(length(ord)==0)
  {
    stop("COCReg.D1.opt: no initial value of gamma and theta gave an estimate.")
  }
  re<-NULL
  for(opt.idx in ord)
  {
    if(m0>1)
    {
      fin<-NULL
      try(fin<-COCReg.coef(Y,X,W,gamma=re.tmp[[opt.idx]]$gamma,theta=re.tmp[[opt.idx]]$theta,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,
                           max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,beta0=B0))
      if(is.null(fin))
      {
        next
      }
      re<-re.tmp[[opt.idx]]
      re$alpha<-fin$alpha
      re$beta<-fin$beta
      if(score.return)
      {
        re$score.y<-fin$score.y
        re$score.x<-fin$score.x
      }
    }else
    {
      re<-re.tmp[[opt.idx]]
    }
    break
  }
  if(is.null(re))
  {
    stop("COCReg.D1.opt: alpha and beta could not be estimated for any initial value of gamma and theta.")
  }
  #----------------------------------
  
  return(re)
}
#################################################

#################################################
# second and higher direction
# Gamma0, Theta0
# option: can choose if remove the identified direction on Y or on X
COCReg.Dk<-function(Y,X,W,Gamma0=NULL,Theta0=NULL,remove.y=TRUE,remove.x=TRUE,Hy=NULL,Hx=NULL,cov.shrinkage.y=TRUE,cov.shrinkage.x=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,burn.in=1000,
                    gamma0.mat=NULL,theta0.mat=NULL,ninitial=NULL,seed=100,beta0=NULL)
{
  # Y: list of outcome data
  # X: list of predictor data
  # W: covariates
  # Gamma0: Y identified directions
  # Theta0: X identified directions
  # remove.y: if removed identified directions of Y
  # remove.x: if removed identified directions of X
  
  #----------------------------------
  # primary parameter setting/estimation
  n<-length(Y)
  p<-ncol(X[[1]])
  q<-ncol(Y[[1]])
  r<-ncol(W)
  
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-c("Intercept",paste0("W",1:(r-1)))
  }
  
  nxvec<-sapply(X,nrow)
  nyvec<-sapply(Y,nrow)
  
  if(min(nxvec)-5<p)
  {
    cov.shrinkage.x<-TRUE
  }
  if(min(nyvec)-5<q)
  {
    cov.shrinkage.y<-TRUE
  }
  #----------------------------------
  
  #----------------------------------
  if(is.null(Gamma0)&is.null(Theta0))
  {
    re<-COCReg.D1.opt(Y=Y,X=X,W=W,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,burn.in=burn.in,
                      gamma0.mat=gamma0.mat,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed,beta0=beta0)
    return(re)
  }else
  {
    if(is.null(Gamma0)==FALSE&remove.y==TRUE)
    {
      Ytmp<-vector("list",length=n)
      for(i in 1:n)
      {
        Ytmp[[i]]<-Y[[i]]-Y[[i]]%*%(Gamma0%*%t(Gamma0))
      }
    }else
    {
      Ytmp<-Y
    }
    if(is.null(Theta0)==FALSE&remove.x==TRUE)
    {
      Xtmp<-vector("list",length=n)
      for(i in 1:n)
      {
        Xtmp[[i]]<-X[[i]]-X[[i]]%*%(Theta0%*%t(Theta0))
      }
    }else
    {
      Xtmp<-X
    }
    
    re<-COCReg.D1.opt(Y=Ytmp,X=Xtmp,W=W,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,burn.in=burn.in,
                      gamma0.mat=gamma0.mat,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed,beta0=beta0)
    
    re$orthogonal.gamma<-c(t(re$gamma)%*%Gamma0)
    re$orthogonal.theta<-c(t(re$theta)%*%Theta0)
    
    return(re)
  }
  #----------------------------------
}
#################################################

#################################################
# level of diagonalization
diag.level<-function(Y,Gamma)
{
  # Y: list of outcome data
  # Gamma: common diagonalization matrix
  
  n<-length(Y)
  ps<-ncol(Gamma)
  
  if(is.null(ncol(Gamma))|ncol(Gamma)==1)
  {
    re<-list(avg.level=1,sub.level=rep(1,n))
  }else
  {
    nyvec<-sapply(Y,nrow)
    
    dl.sub<-matrix(NA,n,ps)
    colnames(dl.sub)<-paste0("C",1:ps)
    dl.sub[,1]<-1
    for(i in 1:n)
    {
      cov.tmp<-cov(Y[[i]])
      for(j in 2:ps)
      {
        gamma.tmp<-matrix(Gamma[,1:j],ncol=j)
        mat.tmp<-t(gamma.tmp)%*%cov.tmp%*%gamma.tmp
        dl.sub[i,j]<-det(diag(diag(mat.tmp)))/det(mat.tmp)
      }
    }
    pmean<-apply(dl.sub,2,function(y){return(prod(apply(cbind(y,nyvec),1,function(x){return(x[1]^(x[2]/sum(nyvec)))})))})
    re<-list(avg.level=pmean,sub.level=dl.sub)
  }
  
  return(re)
}
#################################################

#################################################
# finding first k directions
COCReg<-function(Y,X,W,stop.crt=c("nD","DfD.y","DfD"),nD=NULL,DfD.y.thred=2,DfD.thred=2,remove.y=TRUE,remove.x=TRUE,Hy=NULL,Hx=NULL,cov.shrinkage.y=TRUE,cov.shrinkage.x=TRUE,
                 max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,burn.in=1000,gamma0.mat=NULL,theta0.mat=NULL,ninitial=NULL,seed=100,verbose=TRUE,beta0=NULL)
{
  # Y: list of outcome data
  # X: list of predictor data
  # W: covariates
  # stop.crt: stopping criterion, nD=# of directions, DfD.y=DfD.y threshold
  # remove.y: if removed identified directions of Y
  # remove.x: if removed identified directions of X
  
  if(stop.crt[1]=="nD"&is.null(nD))
  {
    stop.crt<-"DfD"
  }
  
  #----------------------------------
  # primary parameter setting/estimation
  n<-length(Y)
  p<-ncol(X[[1]])
  q<-ncol(Y[[1]])
  r<-ncol(W)
  
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-c("Intercept",paste0("W",1:(r-1)))
  }
  
  nxvec<-sapply(X,nrow)
  nyvec<-sapply(Y,nrow)
  
  if(min(nxvec)-5<p)
  {
    cov.shrinkage.x<-TRUE
  }
  if(min(nyvec)-5<q)
  {
    cov.shrinkage.y<-TRUE
  }
  
  # covariance estimate of X
  if(cov.shrinkage.x)
  {
    # shrinkage estimator
    Sx<-cov.sk.x(X)
  }else
  {
    # sample covariance matrix
    Sx<-array(NA,c(p,p,n))
    for(i in 1:n)
    {
      Xtmp<-scale(X[[i]],center=TRUE,scale=FALSE)
      Sx[,,i]<-t(Xtmp)%*%Xtmp/nxvec[i]
    }
  }
  
  # covariance estimate of Y
  Sy<-array(NA,c(q,q,n))
  for(i in 1:n)
  {
    if(cov.shrinkage.y)
    {
      Sy[,,i]<-cov.ls(Y[[i]])
    }else
    {
      Ytmp<-scale(Y[[i]],center=TRUE,scale=FALSE)
      Sy[,,i]<-t(Ytmp)%*%Ytmp/nyvec[i]
    }
  }
  #----------------------------------
  
  #----------------------------------
  # First direction
  tm1<-system.time(re1<-COCReg.D1.opt(Y,X,W,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,max.itr=max.itr,tol=tol,trace=FALSE,
                                      score.return=score.return,burn.in=burn.in,gamma0.mat=gamma0.mat,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed,beta0=beta0))
  
  Gamma.est<-matrix(re1$gamma,ncol=1)
  Theta.est<-matrix(re1$theta,ncol=1)
  alpha.est<-c(re1$alpha)
  beta.est<-matrix(re1$beta,ncol=1)
  
  cp.time<-matrix(as.numeric(tm1[1:3]),ncol=1)
  
  if(score.return)
  {
    score.y<-matrix(re1$score.y,ncol=1)
    score.x<-matrix(re1$score.x,ncol=1)
  }
  
  if(verbose)
  {
    print(paste0("Component ",ncol(Gamma.est)))
  }
  #----------------------------------
  
  #----------------------------------
  if(stop.crt[1]=="nD")
  {
    if(nD>1)
    {
      for(j in 2:nD)
      {
        re.tmp<-NULL
        try(tm.tmp<-system.time(re.tmp<-COCReg.Dk(Y,X,W,Gamma0=Gamma.est,Theta0=Theta.est,remove.y=remove.y,remove.x=remove.x,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,
                                                  max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,burn.in=burn.in,gamma0.mat=gamma0.mat,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed,beta0=beta0)))
        if(is.null(re.tmp)==FALSE)
        {
          Gamma.est<-cbind(Gamma.est,re.tmp$gamma)
          Theta.est<-cbind(Theta.est,re.tmp$theta)
          alpha.est<-c(alpha.est,re.tmp$alpha)
          beta.est<-cbind(beta.est,re.tmp$beta)
          
          cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
          
          if(score.return)
          {
            score.y<-cbind(score.y,re.tmp$score.y)
            score.x<-cbind(score.x,re.tmp$score.x)
          }
          
          if(verbose)
          {
            print(paste0("Component ",ncol(Gamma.est)))
          }
        }else
        {
          break
        }
      }
    }
    
    colnames(Gamma.est)=colnames(Theta.est)=names(alpha.est)=colnames(beta.est)<-paste0("C",1:ncol(Gamma.est))
    rownames(beta.est)<-colnames(W)
    
    cp.time<-cbind(cp.time,apply(cp.time,1,sum))
    colnames(cp.time)<-c(paste0("C",1:ncol(Gamma.est)),"Total")
    rownames(cp.time)<-c("user","system","elapsed")
    
    DfD.y<-diag.level(Y,Gamma.est)
    DfD.x<-diag.level(X,Theta.est)
  }
  if(stop.crt[1]=="DfD.y")
  {
    nD<-1
    DfD.tmp<-1
    while(DfD.tmp<DfD.y.thred)
    {
      re.tmp<-NULL
      try(tm.tmp<-system.time(re.tmp<-COCReg.Dk(Y,X,W,Gamma0=Gamma.est,Theta0=Theta.est,remove.y=remove.y,remove.x=remove.x,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,
                                                max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,burn.in=burn.in,gamma0.mat=gamma0.mat,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed,beta0=beta0)))
      if(is.null(re.tmp)==FALSE)
      {
        nD<-nD+1
        
        DfD.y<-diag.level(Y,cbind(Gamma.est,re.tmp$gamma))
        DfD.tmp<-DfD.y$avg.level[nD]
        
        if(DfD.tmp<DfD.y.thred)
        {
          Gamma.est<-cbind(Gamma.est,re.tmp$gamma)
          Theta.est<-cbind(Theta.est,re.tmp$theta)
          alpha.est<-c(alpha.est,re.tmp$alpha)
          beta.est<-cbind(beta.est,re.tmp$beta)
          
          cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
          
          if(score.return)
          {
            score.y<-cbind(score.y,re.tmp$score.y)
            score.x<-cbind(score.x,re.tmp$score.x)
          }
          
          if(verbose)
          {
            print(paste0("Component ",ncol(Gamma.est)))
          }
        }else
        {
          break
        }
      }else
      {
        break
      }
    }
    
    colnames(Gamma.est)=colnames(Theta.est)=names(alpha.est)=colnames(beta.est)<-paste0("C",1:ncol(Gamma.est))
    rownames(beta.est)<-colnames(W)
    
    cp.time<-cbind(cp.time,apply(cp.time,1,sum))
    colnames(cp.time)<-c(paste0("C",1:ncol(Gamma.est)),"Total")
    rownames(cp.time)<-c("user","system","elapsed")
    
    DfD.y<-diag.level(Y,Gamma.est)
    DfD.x<-diag.level(X,Theta.est)
  }
  if(stop.crt[1]=="DfD")
  {
    nD<-1
    DfD.tmp<-1
    while(DfD.tmp<DfD.thred)
    {
      re.tmp<-NULL
      try(tm.tmp<-system.time(re.tmp<-COCReg.Dk(Y,X,W,Gamma0=Gamma.est,Theta0=Theta.est,remove.y=remove.y,remove.x=remove.x,Hy=Hy,Hx=Hx,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,
                                                max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,burn.in=burn.in,gamma0.mat=gamma0.mat,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed,beta0=beta0)))
      if(is.null(re.tmp)==FALSE)
      {
        nD<-nD+1
        
        DfD.y<-diag.level(Y,cbind(Gamma.est,re.tmp$gamma))
        DfD.x<-diag.level(X,cbind(Theta.est,re.tmp$theta))
        DfD.tmp<-max(DfD.y$avg.level[nD],DfD.x$avg.level[nD],na.rm=TRUE)
        
        if(DfD.tmp<DfD.thred&sum(is.na(c(DfD.y$avg.level,DfD.x$avg.level)))==0)
        {
          Gamma.est<-cbind(Gamma.est,re.tmp$gamma)
          Theta.est<-cbind(Theta.est,re.tmp$theta)
          alpha.est<-c(alpha.est,re.tmp$alpha)
          beta.est<-cbind(beta.est,re.tmp$beta)
          
          cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
          
          if(score.return)
          {
            score.y<-cbind(score.y,re.tmp$score.y)
            score.x<-cbind(score.x,re.tmp$score.x)
          }
          
          if(verbose)
          {
            print(paste0("Component ",ncol(Gamma.est)))
          }
        }else
        {
          break
        }
      }else
      {
        break
      }
    }
    
    colnames(Gamma.est)=colnames(Theta.est)=names(alpha.est)=colnames(beta.est)<-paste0("C",1:ncol(Gamma.est))
    rownames(beta.est)<-colnames(W)
    
    cp.time<-cbind(cp.time,apply(cp.time,1,sum))
    colnames(cp.time)<-c(paste0("C",1:ncol(Gamma.est)),"Total")
    rownames(cp.time)<-c("user","system","elapsed")
    
    DfD.y<-diag.level(Y,Gamma.est)
    DfD.x<-diag.level(X,Theta.est)
  }
  #----------------------------------
  
  #----------------------------------
  if(score.return)
  {
    colnames(score.y)=colnames(score.x)<-paste0("C",1:ncol(Gamma.est))
    re<-list(gamma=Gamma.est,theta=Theta.est,alpha=alpha.est,beta=beta.est,score.y=score.y,score.x=score.x,DfD.y=DfD.y,DfD.x=DfD.x,
             orthogonality.y=t(Gamma.est)%*%Gamma.est,orthogonality.x=t(Theta.est)%*%Theta.est,time=cp.time)
  }else
  {
    re<-list(gamma=Gamma.est,theta=Theta.est,alpha=alpha.est,beta=beta.est,DfD.y=DfD.y,DfD.x=DfD.x,
             orthogonality.y=t(Gamma.est)%*%Gamma.est,orthogonality.x=t(Theta.est)%*%Theta.est,time=cp.time)
  }
  
  return(re)
  #----------------------------------
}
#################################################

#################################################
# inference: bootstrap
COCReg.coef.boot<-function(Y,X,W,gamma=NULL,theta=NULL,cov.shrinkage.y=TRUE,cov.shrinkage.x=TRUE,boot=TRUE,sims=1000,boot.ci.type=c("se","perc"),conf.level=0.95,
                           max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,outlier.remove=FALSE,verbose=TRUE,beta0=NULL)
{
  # Y: list of outcome data
  # X: list of predictor data
  # W: covariates
  # gamma: linear projection on Y
  # theta: linear projection on X
  
  #-------------------------------
  n<-length(Y)
  p<-ncol(X[[1]])
  q<-ncol(Y[[1]])
  r<-ncol(W)
  
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-c("Intercept",paste0("W",1:(r-1)))
  }
  
  nxvec<-sapply(X,nrow)
  nyvec<-sapply(Y,nrow)
  #-------------------------------
  
  #-------------------------------
  if(boot)
  {
    if(is.null(gamma)==FALSE&is.null(theta)==FALSE)
    {
      #-------------------------------
      B0<-.coc_beta0(r,beta0)     # V4: initial value(s) of beta of every draw
      alpha.boot<-rep(NA,sims)
      beta.boot<-matrix(NA,r,sims)
      
      for(b in 1:sims)
      {
        set.seed(100+b)
        idx.tmp<-sample(1:n,n,replace=TRUE)
        
        Ytmp<-Y[idx.tmp]
        Xtmp<-X[idx.tmp]
        Wtmp<-matrix(W[idx.tmp,],ncol=ncol(W))
        colnames(Wtmp)<-colnames(W)
        
        re.tmp<-NULL
        try(re.tmp<-COCReg.coef(Y=Ytmp,X=Xtmp,W=Wtmp,gamma=gamma,theta=theta,cov.shrinkage.y=cov.shrinkage.y,cov.shrinkage.x=cov.shrinkage.x,max.itr=max.itr,tol=tol,trace=FALSE,score.return=score.return,beta0=B0))
        if(is.null(re.tmp)==FALSE)
        {
          alpha.boot[b]<-re.tmp$alpha
          beta.boot[,b]<-re.tmp$beta
        }
        
        if(verbose)
        {
          print(paste0("Bootstrap sample ",b))
        }
      }
      #-------------------------------
      
      #-------------------------------
      # remove outliers
      if(outlier.remove)
      {
        if(sum(is.na(alpha.boot))>0)
        {
          itmp.nna<-which(is.na(alpha.boot)==FALSE)
          dis.cook<-cooks.distance(lm(alpha.boot~1))
          cook.thred<-4/((length(itmp.nna)-2-2))
          itmp<-itmp.nna[which(dis.cook>cook.thred)]
          alpha.boot[itmp]<-NA
          
          # print(length(itmp))
          
          while(length(itmp)>0)
          {
            itmp.nna<-which(is.na(alpha.boot)==FALSE)
            dis.cook<-cooks.distance(lm(alpha.boot~1))
            cook.thred<-4/((length(itmp.nna)-2-2))
            itmp<-itmp.nna[which(dis.cook>cook.thred)]
            alpha.boot[itmp]<-NA
            
            # print(length(itmp))
          }
        }
        for(j in 1:r)
        {
          if(sum(is.na(beta.boot[j,]))>0)
          {
            itmp.nna<-which(is.na(beta.boot[j,])==FALSE)
            dis.cook<-cooks.distance(lm(beta.boot[j,]~1))
            cook.thred<-4/((length(itmp.nna)-2-2))
            itmp<-itmp.nna[which(dis.cook>cook.thred)]
            beta.boot[j,itmp]<-NA
            
            # print(length(itmp))
            
            while(length(itmp)>0)
            {
              itmp.nna<-which(is.na(beta.boot[j,])==FALSE)
              dis.cook<-cooks.distance(lm(beta.boot[j,]~1))
              cook.thred<-4/((length(itmp.nna)-2-2))
              itmp<-itmp.nna[which(dis.cook>cook.thred)]
              beta.boot[j,itmp]<-NA
              
              # print(length(itmp))
            }
          }
        }
      }
      #-------------------------------
      
      #-------------------------------
      alpha.est<-mean(alpha.boot,na.rm=TRUE)
      alpha.se<-sd(alpha.boot,na.rm=TRUE)
      alpha.stat<-alpha.est/alpha.se
      alpha.pv<-(1-pnorm(abs(alpha.stat)))*2
      
      beta.est<-apply(beta.boot,1,mean,na.rm=TRUE)
      beta.se<-apply(beta.boot,1,sd,na.rm=TRUE)
      beta.stat<-beta.est/beta.se
      beta.pv<-(1-pnorm(abs(beta.stat)))*2
      
      if(boot.ci.type[1]=="se")
      {
        zv<-qnorm(1-(1-conf.level)/2,mean=0,sd=1)
        alpha.ci<-c(alpha.est-zv*alpha.se,alpha.est+zv*alpha.se)
        beta.ci<-cbind(beta.est-zv*beta.se,beta.est+zv*beta.se)
      }
      if(boot.ci.type[1]=="perc")
      {
        alpha.ci<-quantile(alpha.boot,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
        beta.ci<-t(apply(beta.boot,1,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE))
      }
      
      re.alpha<-data.frame(Estimate=alpha.est,SE=alpha.se,statistic=alpha.stat,pvalue=alpha.pv,LB=alpha.ci[1],UB=alpha.ci[2])
      rownames(re.alpha)<-"X"
      re.beta<-data.frame(Estimate=beta.est,SE=beta.se,statistic=beta.stat,pvalue=beta.pv,LB=beta.ci[,1],UB=beta.ci[,2])
      rownames(re.beta)<-colnames(W)
      
      re<-list(alpha=re.alpha,beta=re.beta,alpha.boot=alpha.boot,beta.boot=beta.boot)
      return(re)
      #-------------------------------
    }else
    {
      stop("Error! Need gamma and theta value.")
    }
  }else
  {
    stop("Error!")
  }
}
#################################################


