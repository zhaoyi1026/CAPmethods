#################################################
# Multilevel Covariate Assisted Principal regression (MCAP, gamma-varying)
# linear shrinkage on the covariance matrix
#
# V4: RcppArmadillo-accelerated update of V3. Same interface and results as V3
# (bit-identical); the score and objective hot kernels are compiled. The kernel
# (mlcap_kernels.cpp) is built once to a shared object by build_kernels.R and
# loaded via dyn.load() + .Call() -- no compiler is needed to *run*, only to
# build. See LCAP_gamma-var/README.md.
#################################################


#################################################
# Locate this folder (holds mlcap_kernels.cpp / rvmf_function.R). Set option
# "mlcap.v4.dir" or variable MLCAP_V4_DIR to override; defaults to a search.

# von Mises-Fisher sampler: use Rfast if available, else the bundled fallback.

# Load the compiled kernel (build once if missing).
#################################################

###########################################
# Bias corrected confidence interval for bootstrap
BC.CI <- function(theta,sims,conf.level=0.95) 
{
  low <- (1 - conf.level)/2
  high <- 1 - low
  z.inv <- length(theta[theta < mean(theta)])/sims
  z <- qnorm(z.inv)
  U <- (sims - 1) * (mean(theta) - theta)
  top <- sum(U^3)
  under <- (1/6) * (sum(U^2))^{3/2}
  a <- top/under
  lower.inv <- pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
  lower2 <- lower <- quantile(theta, lower.inv)
  upper.inv <- pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
  upper2 <- upper <- quantile(theta, upper.inv)
  return(c(lower, upper))
}
###########################################

###########################################
# subspace similarity: Krzanowski (1979)
# library("Matrix")
space_sim<-function(M1,M2,thred=1e-2)
{
  k1<-length(which(abs(eigen(t(M1)%*%M1)$values)>thred))
  k2<-length(which(abs(eigen(t(M2)%*%M2)$values)>thred))
  # k1<-rankMatrix(M1)
  # k2<-rankMatrix(M2)
  
  M1.new<-matrix(M1[,1:k1],ncol=k1)
  M2.new<-matrix(M2[,1:k2],ncol=k2)
  
  S<-t(M1.new)%*%M2.new%*%t(M2.new)%*%M1.new
  
  re<-list(similarity=sum(diag(S))/min(k1,k2),k1=k1,k2=k2)
  
  return(re)
}
###########################################

#################################################
# helper: compute scores using C++ kernel
compute.scores<-function(Y.cov,gamma.rnd)
{
  nvec<-sapply(Y.cov,function(x){return(dim(x)[[3]])})
  .Call("mlcap_compute_scores_cpp",Y.cov,gamma.rnd,as.integer(nvec))
}

#################################################
# objective function (C++ accelerated)
obj.func<-function(X1=NULL,X2=NULL,Tmat,Sigma,gamma.rnd,gamma,kappa,beta0.rnd,beta1=NULL,beta2.rnd=NULL,beta0,beta2=NULL,sigma2,Omega=NULL)
{
  nvec<-as.integer(sapply(Sigma,function(x){return(dim(x)[[3]])}))
  score<-.Call("mlcap_compute_scores_cpp",Sigma,gamma.rnd,nvec)
  .Call("mlcap_obj_func_cpp",X1,X2,Tmat,score,gamma.rnd,gamma,kappa,beta0.rnd,beta1,beta2.rnd,beta0,beta2,sigma2,Omega,nvec)
}
#################################################

#################################################
# given gammas, estimate coefficients

# only has fixed effects
lcap.beta.cov.fix<-function(Y.cov,X1,Tmat,gamma.rnd,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE)
{
  # Y.cov: a list of covariance matrix of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # Tmat: m by nmax matrix of sample size
  # gamma.rnd: random gamma
  
  #-----------------------------
  m<-length(Y.cov)
  nvec<-sapply(Y.cov,dim)[3,]
  p<-sapply(Y.cov,dim)[1,1]
  Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
  
  if(is.null(names(Y.cov))==FALSE)
  {
    cluster.names<-names(Y.cov)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  q1<-ncol(X1[[1]])
  if(is.null(names(X1))==TRUE)
  {
    names(X1)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X1[[i]]))==TRUE)
    {
      X1.names<-paste0("X1",1:q1)
      colnames(X1[[i]])<-X1.names
    }else
    {
      X1.names<-colnames(X1[[i]])
    }
  }
  q<-q1+1
  #-----------------------------
  
  #-----------------------------
  score<-compute.scores(Y.cov,gamma.rnd)
  rownames(score)<-cluster.names
  colnames(score)<-paste0("unit",1:ncol(score))
  #-----------------------------
  
  #-----------------------------
  gamma.bar<-apply(gamma.rnd,2,mean,na.rm=TRUE)
  gamma.R<-sqrt(sum(gamma.bar^2))
  gamma.est<-gamma.bar/gamma.R
  kappa.est<-(gamma.R*(p-gamma.R^2))/(1-gamma.R^2)
  #-----------------------------
  
  #-----------------------------
  # Initial estimate: using a mixed effects model
  dtmp<-NULL
  for(i in 1:m)
  {
    dtmp<-rbind(dtmp,data.frame(ID=rep(i,nvec[i]),score=log(score[i,1:nvec[i]]),X1[[i]]))
  }
  colnames(dtmp)<-c("ID","score",X1.names)
  eval(parse(text=paste0("fit.tmp<-lme(score~",paste(X1.names,collapse="+"),",random=~1|ID,data=dtmp,control=lmeControl(opt='optim'))")))
  coef.fix<-fit.tmp$coefficients$fixed
  beta1.new<-coef.fix[-1]
  beta0.rnd.new<-coef.fix[1]+fit.tmp$coefficients$random$ID
  beta0.new<-mean(beta0.rnd.new,na.rm=TRUE)
  sigma2.new<-mean((beta0.rnd.new-beta0.new)^2,na.rm=TRUE)
  #-----------------------------
  
  #-----------------------------
  if(trace)
  {
    beta1.trace<-cbind(beta1.new)
    beta0.rnd.trace<-cbind(beta0.rnd.new)
    beta0.trace<-cbind(beta0.new)
    sigma2.trace<-cbind(sigma2.new)
    
    obj.trace<-obj.func(X1=X1,X2=NULL,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd,gamma=gamma.est,kappa=kappa.est,beta0.rnd=beta0.rnd.new,beta1=beta1.new,beta0=beta0.new,sigma2=sigma2.new)
  }
  #-----------------------------
  
  #-----------------------------
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update random beta0
    beta0.rnd.upd<-rep(NA,m)
    for(i in 1:m)
    {
      fit.tmp<-beta0.rnd.new[i]+c(X1[[i]]%*%beta1.new)
      
      pt1<-sum((1-score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2))+(beta0.rnd.new[i]-beta0.new)/sigma2.new
      pt2<-sum(score[i,1:nvec[i]]*exp(-fit.tmp)*Tmat[i,1:nvec[i]]/(2))+1/sigma2.new
      
      beta0.rnd.upd[i]<-beta0.rnd.new[i]-pt1/pt2
    }
    # update beta0
    beta0.upd<-mean(beta0.rnd.upd,na.rm=TRUE)
    # update beta0 variance
    sigma2.upd<-mean((beta0.rnd.upd-beta0.upd)^2,na.rm=TRUE)
    
    # update beta1
    pt1<-rep(0,q1)
    pt2<-matrix(0,q1,q1)
    for(i in 1:m)
    {
      fit.tmp<-beta0.rnd.upd[i]+c(X1[[i]]%*%beta1.new)
      
      tmp1<-(1-score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2)
      pt1<-pt1+apply(t(apply(cbind(tmp1,X1[[i]]),1,function(x){return(x[1]*x[-1])})),2,sum,na.rm=TRUE)
      
      tmp2<-score[i,1:nvec[i]]*exp(-fit.tmp)*Tmat[i,1:nvec[i]]/(2)
      tmp3<-matrix(0,q1,q1)
      colnames(tmp3)=rownames(tmp3)<-X1.names
      for(j in 1:nvec[i])
      {
        tmp3<-tmp3+tmp2[j]*(X1[[i]][j,]%*%t(X1[[i]][j,]))
      }
      pt2<-pt2+tmp3
    }
    beta1.upd<-c(beta1.new-ginv(pt2)%*%pt1)
    
    # calculate converge criterion
    diff<-max(abs(c(beta0.upd-beta0.new,beta1.upd-beta1.new)))
    
    # update parameters
    beta0.rnd.new<-beta0.rnd.upd
    beta0.new<-beta0.upd
    sigma2.new<-sigma2.upd
    beta1.new<-beta1.upd
    
    if(trace)
    {
      beta1.trace<-cbind(beta1.trace,beta1.new)
      beta0.rnd.trace<-cbind(beta0.rnd.trace,beta0.rnd.new)
      beta0.trace<-cbind(beta0.trace,beta0.new)
      sigma2.trace<-cbind(sigma2.trace,sigma2.new)
      
      obj.trace<-c(obj.trace,obj.func(X1=X1,X2=NULL,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd,gamma=gamma.est,kappa=kappa.est,beta0.rnd=beta0.rnd.new,beta1=beta1.new,beta0=beta0.new,sigma2=sigma2.new))
    }
    
    # print(c(diff,obj.func(X1=X1,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd,gamma=gamma.est,kappa=kappa.est,beta0.rnd=beta0.rnd.new,beta1=beta1.new,beta0=beta0.new,sigma2=sigma2.new)))
  }
  #-----------------------------
  
  beta.est<-c(beta0.new,beta1.new)
  names(beta.est)<-c("Intercept",X1.names)
  
  if(trace)
  {
    colnames(beta1.trace)=colnames(beta0.rnd.trace)=colnames(beta0.trace)=colnames(sigma2.trace)<-paste0("iteration",0:(ncol(beta1.trace)-1))
    
    beta.trace<-rbind(beta0.trace,beta1.trace)
    rownames(beta.trace)[1]<-"Intercept"
    
    if(score.return)
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),score=score,
               beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,obj=obj.trace)
    }else
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
               beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,obj=obj.trace)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),score=score)
    }else
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr))
    }
  }
  
  return(re)
}

# only has random effects
lcap.beta.cov.rnd<-function(Y.cov,X2,Tmat,gamma.rnd,Omega.diag=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE)
{
  # Y.cov: a list of covariance matrix of Y
  # X2: a list of covariates corresponding to random-effect beta
  # Tmat: m by nmax matrix of sample size
  # gamma.rnd: random gamma
  
  #-----------------------------
  m<-length(Y.cov)
  nvec<-sapply(Y.cov,dim)[3,]
  p<-sapply(Y.cov,dim)[1,1]
  Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
  
  if(is.null(names(Y.cov))==FALSE)
  {
    cluster.names<-names(Y.cov)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  q2<-ncol(X2[[1]])
  if(is.null(names(X2))==TRUE)
  {
    names(X2)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X2[[i]]))==TRUE)
    {
      X2.names<-paste0("X2",1:q2)
      colnames(X2[[i]])<-X2.names
    }else
    {
      X2.names<-colnames(X2[[i]])
    }
  }
  q<-q2+1
  #-----------------------------
  
  #-----------------------------
  score<-compute.scores(Y.cov,gamma.rnd)
  rownames(score)<-cluster.names
  colnames(score)<-paste0("unit",1:ncol(score))
  #-----------------------------

  #-----------------------------
  gamma.bar<-apply(gamma.rnd,2,mean,na.rm=TRUE)
  gamma.R<-sqrt(sum(gamma.bar^2))
  gamma.est<-gamma.bar/gamma.R
  kappa.est<-(gamma.R*(p-gamma.R^2))/(1-gamma.R^2)
  #-----------------------------

  #-----------------------------
  # Initial estimate: using a mixed effects model
  beta0.rnd.new<-rep(NA,m)
  beta2.rnd.new<-matrix(NA,m,q2)
  rownames(beta2.rnd.new)<-names(X2)
  colnames(beta2.rnd.new)<-X2.names
  for(i in 1:m)
  {
    dtmp<-data.frame(score=log(score[i,1:nvec[i]]),X2[[i]])
    fit.tmp<-lm(score~.,data=dtmp)
    beta0.rnd.new[i]<-fit.tmp$coefficients[1]
    beta2.rnd.new[i,]<-fit.tmp$coefficients[-1]
  }
  beta0.new<-mean(beta0.rnd.new,na.rm=TRUE)
  sigma2.new<-mean((beta0.rnd.new-beta0.new)^2,na.rm=TRUE)
  beta2.new<-apply(beta2.rnd.new,2,mean,na.rm=TRUE)
  Omega.new<-cov(beta2.rnd.new)*(m-1)/m
  if(Omega.diag==TRUE)
  {
    Omega.tmp<-matrix(0,q2,q2)
    colnames(Omega.tmp)=rownames(Omega.tmp)<-X2.names
    diag(Omega.tmp)<-diag(Omega.new)
    Omega.new<-Omega.tmp
  }
  #-----------------------------
  
  #-----------------------------
  if(trace)
  {
    beta0.rnd.trace<-cbind(beta0.rnd.new)
    beta0.trace<-cbind(beta0.new)
    sigma2.trace<-cbind(sigma2.new)
    beta2.rnd.trace<-cbind(beta2.rnd.new)
    beta2.trace<-cbind(beta2.new)
    Omega.trace<-cbind(Omega.new)
    
    obj.trace<-obj.func(X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd,gamma=gamma.est,kappa=kappa.est,
                        beta0.rnd=beta0.rnd.new,beta2.rnd=beta2.rnd.new,beta0=beta0.new,beta2=beta2.new,sigma2=sigma2.new,Omega=Omega.new)
  }
  #-----------------------------
  
  #-----------------------------
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update random beta0
    beta0.rnd.upd<-rep(NA,m)
    for(i in 1:m)
    {
      fit.tmp<-beta0.rnd.new[i]+X2[[i]]%*%beta2.rnd.new[i,]
      
      pt1<-sum((1-score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2))+(beta0.rnd.new[i]-beta0.new)/sigma2.new
      pt2<-sum((score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2))+1/sigma2.new
      
      beta0.rnd.upd[i]<-beta0.rnd.new[i]-pt1/pt2
    }
    # update beta0
    beta0.upd<-mean(beta0.rnd.upd,na.rm=TRUE)
    # update beta0 variance
    sigma2.upd<-mean((beta0.rnd.upd-beta0.upd)^2,na.rm=TRUE)
    
    # update beta2
    beta2.rnd.upd<-matrix(NA,m,q2)
    rownames(beta2.rnd.upd)<-names(X2)
    colnames(beta2.rnd.upd)<-X2.names
    for(i in 1:m)
    {
      fit.tmp<-beta0.rnd.upd[i]+X2[[i]]%*%beta2.rnd.new[i,]
      
      pt1<-rep(0,q2)
      pt2<-matrix(0,q2,q2)
      for(j in 1:nvec[i])
      {
        pt1<-pt1+((1-score[i,j]*exp(-fit.tmp[j]))*Tmat[i,j]/(2))*X2[[i]][j,]
        
        pt2<-pt2+(score[i,j]*exp(-fit.tmp[j])*Tmat[i,j]/(2))*(X2[[i]][j,]%*%t(X2[[i]][j,]))
      }
      pt1<-pt1+c(ginv(Omega.new)%*%(beta2.rnd.new[i,]-beta2.new))
      pt2<-pt2+ginv(Omega.new)
      
      beta2.rnd.upd[i,]<-beta2.rnd.new[i,]-ginv(pt2)%*%pt1
    }
    beta2.upd<-apply(beta2.rnd.upd,2,mean,na.rm=TRUE)
    Omega.upd<-cov(beta2.rnd.upd)*(m-1)/m
    if(Omega.diag==TRUE)
    {
      Omega.tmp<-matrix(0,q2,q2)
      colnames(Omega.tmp)=rownames(Omega.tmp)<-X2.names
      diag(Omega.tmp)<-diag(Omega.upd)
      Omega.upd<-Omega.tmp
    }
    
    # calculate converge criterion
    diff<-max(abs(c(beta0.upd-beta0.new,beta2.upd-beta2.new)))
    
    # update parameters
    beta0.rnd.new<-beta0.rnd.upd
    beta0.new<-beta0.upd
    sigma2.new<-sigma2.upd
    beta2.rnd.new<-beta2.rnd.upd
    beta2.new<-beta2.upd
    Omega.new<-Omega.upd
    
    if(trace)
    {
      beta0.rnd.trace<-cbind(beta0.rnd.trace,beta0.rnd.new)
      beta0.trace<-cbind(beta0.trace,beta0.new)
      sigma2.trace<-cbind(sigma2.trace,sigma2.new)
      beta2.rnd.trace<-cbind(beta2.rnd.trace,beta2.rnd.new)
      beta2.trace<-cbind(beta2.trace,beta2.new)
      Omega.trace<-cbind(Omega.trace,Omega.new)
      
      obj.trace<-c(obj.trace,obj.func(X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd,gamma=gamma.est,kappa=kappa.est,
                                      beta0.rnd=beta0.rnd.new,beta2.rnd=beta2.rnd.new,beta0=beta0.new,beta2=beta2.new,sigma2=sigma2.new,Omega=Omega.new))
    }
    
    # print(c(diff,obj.func(X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd,gamma=gamma.est,kappa=kappa.est,
    #                       beta0.rnd=beta0.rnd.new,beta2.rnd=beta2.rnd.new,beta0=beta0.new,beta2=beta2.new,sigma2=sigma2.new,Omega=Omega.new)))
  }
  #-----------------------------
  
  beta.est<-c(beta0.new,beta2.new)
  names(beta.est)<-c("Intercept",X2.names)
  
  if(trace)
  {
    colnames(beta2.trace)=colnames(beta0.rnd.trace)=colnames(beta0.trace)=colnames(sigma2.trace)<-paste0("iteration",0:(ncol(beta0.trace)-1))
    
    beta.trace<-rbind(beta0.trace,beta2.trace)
    rownames(beta.trace)[1]<-"Intercept"
    
    if(score.return)
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,beta2.rnd=beta2.rnd.new,beta2.Omega=Omega.new,convergence=(s<max.itr),score=score,
               beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,beta2.rnd.trace=beta2.rnd.trace,beta2.Omega.trace=Omega.trace,obj=obj.trace)
    }else
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,beta2.rnd=beta2.rnd.new,beta2.Omega=Omega.new,convergence=(s<max.itr),
               beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,beta2.rnd.trace=beta2.rnd.trace,beta2.Omega.trace=Omega.trace,obj=obj.trace)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,beta2.rnd=beta2.rnd.new,beta2.Omega=Omega.new,convergence=(s<max.itr),score=score)
    }else
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,beta2.rnd=beta2.rnd.new,beta2.Omega=Omega.new,convergence=(s<max.itr))
    }
  }
  
  return(re)
}

# has both fixed and random effects
lcap.beta.cov.mix<-function(Y.cov,X1,X2,Tmat,gamma.rnd,Omega.diag=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE)
{
  # Y.cov: a list of covariance matrix of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # Tmat: m by nmax matrix of sample size
  # gamma.rnd: random gamma
  
  #-----------------------------
  m<-length(Y.cov)
  nvec<-sapply(Y.cov,dim)[3,]
  p<-sapply(Y.cov,dim)[1,1]
  Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
  
  if(is.null(names(Y.cov))==FALSE)
  {
    cluster.names<-names(Y.cov)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  q1<-ncol(X1[[1]])
  if(is.null(names(X1))==TRUE)
  {
    names(X1)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X1[[i]]))==TRUE)
    {
      X1.names<-paste0("X1",1:q1)
      colnames(X1[[i]])<-X1.names
    }else
    {
      X1.names<-colnames(X1[[i]])
    }
  }
  q2<-ncol(X2[[1]])
  if(is.null(names(X2))==TRUE)
  {
    names(X2)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X2[[i]]))==TRUE)
    {
      X2.names<-paste0("X2",1:q2)
      colnames(X2[[i]])<-X2.names
    }else
    {
      X2.names<-colnames(X2[[i]])
    }
  }
  q<-q1+q2+1
  #-----------------------------
  
  #-----------------------------
  score<-compute.scores(Y.cov,gamma.rnd)
  rownames(score)<-cluster.names
  colnames(score)<-paste0("unit",1:ncol(score))
  #-----------------------------

  #-----------------------------
  gamma.bar<-apply(gamma.rnd,2,mean,na.rm=TRUE)
  gamma.R<-sqrt(sum(gamma.bar^2))
  gamma.est<-gamma.bar/gamma.R
  kappa.est<-(gamma.R*(p-gamma.R^2))/(1-gamma.R^2)
  #-----------------------------
  
  #-----------------------------
  # Initial estimate: using a mixed effects model
  dtmp<-NULL
  for(i in 1:m)
  {
    dtmp<-rbind(dtmp,data.frame(ID=rep(i,nvec[i]),score=log(score[i,1:nvec[i]]),X1[[i]],X2[[i]]))
  }
  colnames(dtmp)<-c("ID","score",X1.names,X2.names)
  fit.tmp<-NULL
  try(eval(parse(text=paste0("fit.tmp<-lme(score~",paste(c(X1.names,X2.names),collapse="+"),",random=~1+",paste(X2.names,collapse="+"),"|ID,data=dtmp,control=lmeControl(opt='optim'))"))),silent=TRUE)
  if(is.null(fit.tmp)==FALSE)
  {
    coef.fix<-fit.tmp$coefficients$fixed
    beta1.new<-coef.fix[2:(q1+1)]
    beta0.rnd.new<-coef.fix[1]+fit.tmp$coefficients$random$ID[,1]
    beta0.new<-mean(beta0.rnd.new,na.rm=TRUE)
    sigma2.new<-mean((beta0.rnd.new-beta0.new)^2,na.rm=TRUE)
    beta2.rnd.new<-matrix(NA,m,q2)
    colnames(beta2.rnd.new)<-X2.names
    rownames(beta2.rnd.new)<-names(X2)
    for(kk in 1:q2)
    {
      beta2.rnd.new[,kk]<-coef.fix[q1+1+kk]+fit.tmp$coefficients$random$ID[,kk+1]
    }
    beta2.new<-apply(beta2.rnd.new,2,mean,na.rm=TRUE)
    Omega.new<-cov(beta2.rnd.new)*(m-1)/m
    if(Omega.diag==TRUE)
    {
      Omega.tmp<-matrix(0,q2,q2)
      colnames(Omega.tmp)=rownames(Omega.tmp)<-X2.names
      diag(Omega.tmp)<-diag(Omega.new)
      Omega.new<-Omega.tmp
    }
  }else
  {
    dtmp<-NULL
    for(i in 1:m)
    {
      dtmp<-rbind(dtmp,data.frame(ID=rep(i,nvec[i]),score=log(score[i,1:nvec[i]]),X1[[i]],X2[[i]]))
    }
    colnames(dtmp)<-c("ID","score",X1.names,X2.names)
    # fixed effects
    eval(parse(text=paste0("fit.fix<-lm(score~",paste(c(X1.names),collapse="+"),",data=dtmp)")))
    beta1.new<-fit.fix$coefficients[-1]
    
    beta0.rnd.new<-rnorm(m,mean=0,sd=1)
    beta0.new<-mean(beta0.rnd.new,na.rm=TRUE)
    sigma2.new<-mean((beta0.rnd.new-beta0.new)^2,na.rm=TRUE)
    
    beta2.rnd.new<-matrix(rnorm(m*q2,mean=0,sd=1),m,q2)
    colnames(beta2.rnd.new)<-X2.names
    rownames(beta2.rnd.new)<-names(X2)
    beta2.new<-apply(beta2.rnd.new,2,mean,na.rm=TRUE)
    Omega.new<-cov(beta2.rnd.new)*(m-1)/m
    if(Omega.diag==TRUE)
    {
      Omega.tmp<-matrix(0,q2,q2)
      colnames(Omega.tmp)=rownames(Omega.tmp)<-X2.names
      diag(Omega.tmp)<-diag(Omega.new)
      Omega.new<-Omega.tmp
    }
  }
  #-----------------------------
  
  #-----------------------------
  if(trace)
  {
    beta0.rnd.trace<-cbind(beta0.rnd.new)
    beta0.trace<-cbind(beta0.new)
    sigma2.trace<-cbind(sigma2.new)
    beta1.trace<-cbind(beta1.new)
    beta2.rnd.trace<-cbind(beta2.rnd.new)
    beta2.trace<-cbind(beta2.new)
    Omega.trace<-cbind(Omega.new)
    
    obj.trace<-obj.func(X1=X1,X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd,gamma=gamma.est,kappa=kappa.est,
                        beta0.rnd=beta0.rnd.new,beta1=beta1.new,beta2.rnd=beta2.rnd.new,beta0=beta0.new,beta2=beta2.new,sigma2=sigma2.new,Omega=Omega.new)
  }
  #-----------------------------
  
  #-----------------------------
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update random beta0
    beta0.rnd.upd<-rep(NA,m)
    for(i in 1:m)
    {
      fit.tmp<-beta0.rnd.new[i]+X1[[i]]%*%beta1.new+X2[[i]]%*%beta2.rnd.new[i,]
      
      pt1<-sum((1-score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2))+(beta0.rnd.new[i]-beta0.new)/sigma2.new
      pt2<-sum(score[i,1:nvec[i]]*exp(-fit.tmp)*Tmat[i,1:nvec[i]]/(2))+1/sigma2.new
      
      beta0.rnd.upd[i]<-beta0.rnd.new[i]-pt1/pt2
    }
    # update beta0
    beta0.upd<-mean(beta0.rnd.upd,na.rm=TRUE)
    # update beta0 variance
    sigma2.upd<-mean((beta0.rnd.upd-beta0.upd)^2,na.rm=TRUE)
    
    # update beta1
    pt1<-rep(0,q1)
    pt2<-matrix(0,q1,q1)
    for(i in 1:m)
    {
      fit.tmp<-beta0.rnd.upd[i]+X1[[i]]%*%beta1.new+X2[[i]]%*%beta2.rnd.new[i,]
      
      tmp1<-(1-score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2)
      pt1<-pt1+apply(t(apply(cbind(tmp1,X1[[i]]),1,function(x){return(x[1]*x[-1])})),2,sum,na.rm=TRUE)
      
      tmp2<-score[i,1:nvec[i]]*exp(-fit.tmp)*Tmat[i,1:nvec[i]]/(2)
      tmp3<-matrix(0,q1,q1)
      colnames(tmp3)=rownames(tmp3)<-X1.names
      for(j in 1:nvec[i])
      {
        tmp3<-tmp3+tmp2[j]*(X1[[i]][j,]%*%t(X1[[i]][j,]))
      }
      pt2<-pt2+tmp3
    }
    beta1.upd<-c(beta1.new-ginv(pt2)%*%pt1)
    
    # update beta2
    beta2.rnd.upd<-matrix(NA,m,q2)
    rownames(beta2.rnd.upd)<-names(X2)
    colnames(beta2.rnd.upd)<-X2.names
    for(i in 1:m)
    {
      fit.tmp<-beta0.rnd.upd[i]+X1[[i]]%*%beta1.upd+X2[[i]]%*%beta2.rnd.new[i,]
      
      pt1<-rep(0,q2)
      pt2<-matrix(0,q2,q2)
      for(j in 1:nvec[i])
      {
        pt1<-pt1+((1-score[i,j]*exp(-fit.tmp[j]))*Tmat[i,j]/(2))*X2[[i]][j,]
        
        pt2<-pt2+(score[i,j]*exp(-fit.tmp[j])*Tmat[i,j]/(2))*(X2[[i]][j,]%*%t(X2[[i]][j,]))
      }
      pt1<-pt1+c(ginv(Omega.new)%*%(beta2.rnd.new[i,]-beta2.new))
      pt2<-pt2+ginv(Omega.new)
      
      beta2.rnd.upd[i,]<-beta2.rnd.new[i,]-ginv(pt2)%*%pt1
    }
    beta2.upd<-apply(beta2.rnd.upd,2,mean,na.rm=TRUE)
    Omega.upd<-cov(beta2.rnd.upd)*(m-1)/m
    if(Omega.diag==TRUE)
    {
      Omega.tmp<-matrix(0,q2,q2)
      colnames(Omega.tmp)=rownames(Omega.tmp)<-X2.names
      diag(Omega.tmp)<-diag(Omega.upd)
      Omega.upd<-Omega.tmp
    }
    
    # calculate converge criterion
    diff<-max(abs(c(beta0.upd-beta0.new,beta1.upd-beta1.new,beta2.upd-beta2.new)))
    
    # update parameters
    beta0.rnd.new<-beta0.rnd.upd
    beta0.new<-beta0.upd
    sigma2.new<-sigma2.upd
    beta1.new<-beta1.upd
    beta2.rnd.new<-beta2.rnd.upd
    beta2.new<-beta2.upd
    Omega.new<-Omega.upd
    
    if(trace)
    {
      beta0.rnd.trace<-cbind(beta0.rnd.trace,beta0.rnd.new)
      beta0.trace<-cbind(beta0.trace,beta0.new)
      sigma2.trace<-cbind(sigma2.trace,sigma2.new)
      beta1.trace<-cbind(beta1.trace,beta1.new)
      beta2.rnd.trace<-cbind(beta2.rnd.trace,beta2.rnd.new)
      beta2.trace<-cbind(beta2.trace,beta2.new)
      Omega.trace<-cbind(Omega.trace,Omega.new)
      
      obj.trace<-c(obj.trace,obj.func(X1=X1,X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd,gamma=gamma.est,kappa=kappa.est,
                                      beta0.rnd=beta0.rnd.new,beta1=beta1.new,beta2.rnd=beta2.rnd.new,beta0=beta0.new,beta2=beta2.new,sigma2=sigma2.new,Omega=Omega.new))
    }
    
    # print(c(diff,obj.func(X1=X1,X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd,gamma=gamma.est,kappa=kappa.est,
    #                       beta0.rnd=beta0.rnd.new,beta1=beta1.new,beta2.rnd=beta2.rnd.new,beta0=beta0.new,beta2=beta2.new,sigma2=sigma2.new,Omega=Omega.new)))
  }
  #-----------------------------
  
  beta.est<-c(beta0.new,beta1.new,beta2.new)
  names(beta.est)<-c("Intercept",X1.names,X2.names)
  
  if(trace)
  {
    colnames(beta2.trace)=colnames(beta1.trace)=colnames(beta0.rnd.trace)=colnames(beta0.trace)=colnames(sigma2.trace)<-paste0("iteration",0:(ncol(beta0.trace)-1))
    
    beta.trace<-rbind(beta0.trace,beta1.trace,beta2.trace)
    rownames(beta.trace)[1]<-"Intercept"
    
    if(score.return)
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,beta2.rnd=beta2.rnd.new,beta2.Omega=Omega.new,convergence=(s<max.itr),score=score,
               beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,beta2.rnd.trace=beta2.rnd.trace,beta2.Omega.trace=Omega.trace,obj=obj.trace)
    }else
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,beta2.rnd=beta2.rnd.new,beta2.Omega=Omega.new,convergence=(s<max.itr),
               beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,beta2.rnd.trace=beta2.rnd.trace,beta2.Omega.trace=Omega.trace,obj=obj.trace)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,beta2.rnd=beta2.rnd.new,beta2.Omega=Omega.new,convergence=(s<max.itr),score=score)
    }else
    {
      re<-list(gamma.rnd=gamma.rnd,beta=beta.est,beta0.rnd=beta0.rnd.new,beta0.sigma2=sigma2.new,beta2.rnd=beta2.rnd.new,beta2.Omega=Omega.new,convergence=(s<max.itr))
    }
  }
  
  return(re)
}

# compile the functions
lcap.beta.cov<-function(Y.cov,X1=NULL,X2=NULL,Tmat,gamma.rnd,Omega.diag=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE)
{
  # Y.cov: a list of covariance matrix of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # Tmat: m by nmax matrix of sample size
  # gamma.rnd: random gamma
  
  if(is.null(X1)==FALSE&is.null(X2)==FALSE)
  {
    re<-lcap.beta.cov.mix(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,gamma.rnd=gamma.rnd,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace)
  }else
    if(is.null(X1)==FALSE)
    {
      re<-lcap.beta.cov.fix(Y.cov=Y.cov,X1=X1,Tmat=Tmat,gamma.rnd=gamma.rnd,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace)
    }else
      if(is.null(X2)==FALSE)
      {
        re<-lcap.beta.cov.rnd(Y.cov=Y.cov,X2=X2,Tmat=Tmat,gamma.rnd=gamma.rnd,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace)
      }
  
  return(re)
}

lcap.beta<-function(Y,X1=NULL,X2=NULL,gamma.rnd,Omega.diag=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE)
{
  # Y: a list of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # gamma.rnd: random gamma
  
  #-----------------------------
  m<-length(Y)
  nvec<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  if(is.null(names(Y))==TRUE)
  {
    names(Y)<-paste0("G",1:m)
  }
  
  Y.cov<-vector("list",length=m)
  names(Y.cov)<-names(Y)
  Tmat<-matrix(NA,m,max(nvec))
  rownames(Tmat)<-names(Y)
  for(i in 1:m)
  {
    Y.cov[[i]]<-array(NA,c(p,p,nvec[i]))
    for(j in 1:nvec[i])
    {
      Tmat[i,j]<-nrow(Y[[i]][[j]])
      Y.cov[[i]][,,j]<-cov(Y[[i]][[j]])*(Tmat[i,j]-1)/Tmat[i,j]
    }
  }
  
  re<-lcap.beta.cov(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,gamma.rnd=gamma.rnd,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace)
  
  return(re)
}
#################################################

#################################################
# estimate both gamma and beta
# finding the first component
gamma.solve.QP<-function(A,H)
{
  p<-ncol(H)
  
  H.svd<-svd(H)
  H.d.sqrt<-diag(sqrt(H.svd$d))
  H.d.sqrt.inv<-diag(1/sqrt(H.svd$d))
  H.sqrt.inv<-H.svd$u%*%H.d.sqrt.inv%*%t(H.svd$v)
  
  #---------------------------------------------------
  # svd decomposition method
  
  eigen.tmp<-eigen(H.d.sqrt.inv%*%t(H.svd$u)%*%A%*%H.svd$u%*%H.d.sqrt.inv)
  eigen.tmp.vec<-Re(eigen.tmp$vectors)
  re<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,p]
  
  # obj<-rep(NA,ncol(eigen.tmp$vectors))
  # for(j in 1:ncol(eigen.tmp$vectors))
  # {
  #   otmp<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,j]
  #   obj[j]<-t(otmp)%*%A%*%otmp
  # }
  # re<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,which.min(obj)]
  #---------------------------------------------------
  
  #---------------------------------------------------
  # eigenvector of A with respect to H
  
  # xtmp<-H.sqrt.inv%*%Re(eigen(H.sqrt.inv%*%A%*%H.sqrt.inv)$vectors)
  # opt.idx<-which.min(diag(t(xtmp)%*%A%*%xtmp))
  # re<-xtmp[,opt.idx]
  # re<-xtmp[,opt.idx]/sqrt(sum((xtmp[,opt.idx])^2))
  
  # eigen.tmp<-eigen(H.sqrt.inv%*%A%*%H.sqrt.inv)
  # eigen.tmp.vec<-Re(eigen.tmp$vectors)
  # 
  # obj<-rep(NA,ncol(eigen.tmp$vectors))
  # for(j in 1:ncol(eigen.tmp$vectors))
  # {
  #   otmp<-H.sqrt.inv%*%eigen.tmp.vec[,j]
  #   obj[j]<-t(otmp)%*%A%*%otmp
  # }
  # re<-H.sqrt.inv%*%eigen.tmp.vec[,p]
  #---------------------------------------------------
  
  return(re)
}
gamma.solve<-function(A,kappa,gamma,H)
{
  # minimize t(gamma.est)%*%A%*%gamma.est-kappa*t(gamma)%*%gamma.est
  # such that t(gamma.est)%*%A%*%gamma.est=1
  
  # 1. find eigenvector of A with respect to H
  # 2. find the eigenvector minimizes t(gamma.est)%*%A%*%gamma.est-kappa*t(gamma)%*%gamma.est
  
  p<-ncol(H)
  
  H.svd<-svd(H)
  H.d.sqrt<-diag(sqrt(H.svd$d))
  H.d.sqrt.inv<-diag(1/sqrt(H.svd$d))
  H.sqrt.inv<-H.svd$u%*%H.d.sqrt.inv%*%t(H.svd$v)
  
  #---------------------------------------------------
  # eigenvector of A with respect to H
  
  xtmp<-H.sqrt.inv%*%Re(eigen(H.sqrt.inv%*%A%*%H.sqrt.inv)$vectors)
  for(ss in 1:ncol(xtmp))
  {
    if(xtmp[which.max(abs(xtmp[,ss])),ss]<0)
    {
      xtmp[,ss]<--xtmp[,ss]
    }
  }
  
  obj<-rep(NA,ncol(xtmp))
  for(ss in 1:ncol(xtmp))
  {
    obj[ss]<-c(t(xtmp[,ss])%*%A%*%xtmp[,ss])-kappa*c(t(gamma)%*%xtmp[,ss])
  }
  
  gamma.est<-xtmp[,which.min(obj)]
  
  re<-list(gamma=gamma.est,obj=obj[which.min(obj)],constraint=c(t(gamma.est)%*%H%*%gamma.est))
  
  return(re)
}

# only has fixed effects
lcap.cov.D1.fix<-function(Y.cov,X1,Tmat,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0=NULL,kappa0=NULL)
{
  # Y.cov: a list of covariance matrix of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # Tmat: m by nmax matrix of sample size
  
  #-----------------------------
  m<-length(Y.cov)
  nvec<-sapply(Y.cov,dim)[3,]
  p<-sapply(Y.cov,dim)[1,1]
  Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
  
  if(is.null(names(Y.cov))==FALSE)
  {
    cluster.names<-names(Y.cov)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  q1<-ncol(X1[[1]])
  if(is.null(names(X1))==TRUE)
  {
    names(X1)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X1[[i]]))==TRUE)
    {
      X1.names<-paste0("X1",1:q1)
      colnames(X1[[i]])<-X1.names
    }else
    {
      X1.names<-colnames(X1[[i]])
    }
  }
  q2<-0
  q<-q1+1
  #-----------------------------
  
  #-----------------------------
  H.mat<-array(NA,c(p,p,m))
  dimnames(H.mat)[[3]]<-cluster.names
  for(i in 1:m)
  {
    if(H.type[1]=="CAvgCov"|H.type[1]=="AvgCov")
    {
      Htmp<-matrix(0,p,p)
      for(j in 1:nvec[i])
      {
        Htmp<-Htmp+Y.cov[[i]][,,j]*Tmat[i,j]
      }
      Htmp<-Htmp/sum(Tmat[i,1:nvec[i]])
      H.mat[,,i]<-Htmp
    }
    if(H.type[1]=="Identity")
    {
      H.mat[,,i]<-diag(rep(1,p))
    }
  }
  if(H.type[1]=="AvgCov")
  {
    Htmp<-matrix(0,p,p)
    for(i in 1:m)
    {
      Htmp<-Htmp+H.mat[,,i]*sum(Tmat[i,1:nvec[i]])/sum(Tmat,na.rm=TRUE)
    }
    for(i in 1:m)
    {
      H.mat[,,i]<-Htmp
    }
  }
  #-----------------------------
  
  if(method[1]=="CAP")
  {
    
    #-----------------------------
    # set initial value of gamma0
    if(is.null(gamma0))
    {
      set.seed(100)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      gamma0<-gamma0.mat[,sample(ncol(gamma0.mat),1)]
    }
    if(is.null(kappa0))
    {
      kappa0<-1
    }
    # generate random gamma
    gamma.new<-gamma0
    kappa.new<-kappa0
    gamma.rnd.new<-rvmf(m,mu=gamma.new,k=kappa.new) 
    rownames(gamma.rnd.new)<-cluster.names
    
    # calculate score
    score<-compute.scores(Y.cov,gamma.rnd.new)
    rownames(score)<-cluster.names
    colnames(score)<-paste0("unit",1:ncol(score))
    
    # Initial estimate: using a mixed effects model
    dtmp<-NULL
    for(i in 1:m)
    {
      dtmp<-rbind(dtmp,data.frame(ID=rep(i,nvec[i]),score=log(score[i,1:nvec[i]]),X1[[i]]))
    }
    colnames(dtmp)<-c("ID","score",X1.names)
    eval(parse(text=paste0("fit.tmp<-lme(score~",paste(X1.names,collapse="+"),",random=~1|ID,data=dtmp,control=lmeControl(opt='optim'))")))
    coef.fix<-fit.tmp$coefficients$fixed
    beta1.new<-coef.fix[-1]
    beta0.rnd.new<-coef.fix[1]+fit.tmp$coefficients$random$ID
    beta0.new<-mean(beta0.rnd.new,na.rm=TRUE)
    sigma2.new<-mean((beta0.rnd.new-beta0.new)^2,na.rm=TRUE)
    
    obj0<-obj.func(X1=X1,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd.new,gamma=gamma.new,kappa=kappa.new,beta0.rnd=beta0.rnd.new,beta1=beta1.new,beta0=beta0.new,sigma2=sigma2.new)
    #-----------------------------
    
    #-----------------------------
    if(trace)
    {
      beta1.trace<-cbind(beta1.new)
      beta0.rnd.trace<-cbind(beta0.rnd.new)
      beta0.trace<-cbind(beta0.new)
      sigma2.trace<-cbind(sigma2.new)
      
      gamma.rnd.trace<-cbind(gamma.rnd.new)
      gamma.trace<-cbind(gamma.new)
      kappa.trace<-c(kappa.new)
      
      obj.trace<-c(obj0)
    }
    #-----------------------------
    
    #-----------------------------
    s<-0
    diff<-100
    while(s<=max.itr&diff>tol)
    {
      s<-s+1
      
      # calculate score
      score<-compute.scores(Y.cov,gamma.rnd.new)
      
      # update random beta0
      beta0.rnd.upd<-rep(NA,m)
      for(i in 1:m)
      {
        fit.tmp<-beta0.rnd.new[i]+c(X1[[i]]%*%beta1.new)
        
        pt1<-sum((1-score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2))+(beta0.rnd.new[i]-beta0.new)/sigma2.new
        pt2<-sum(score[i,1:nvec[i]]*exp(-fit.tmp)*Tmat[i,1:nvec[i]]/(2))+1/sigma2.new
        
        beta0.rnd.upd[i]<-beta0.rnd.new[i]-pt1/pt2
      }
      # update beta0
      beta0.upd<-mean(beta0.rnd.upd,na.rm=TRUE)
      # update beta0 variance
      sigma2.upd<-mean((beta0.rnd.upd-beta0.upd)^2,na.rm=TRUE)
      
      # update beta1
      pt1<-rep(0,q1)
      pt2<-matrix(0,q1,q1)
      for(i in 1:m)
      {
        fit.tmp<-beta0.rnd.upd[i]+c(X1[[i]]%*%beta1.new)
        
        tmp1<-(1-score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2)
        pt1<-pt1+apply(t(apply(cbind(tmp1,X1[[i]]),1,function(x){return(x[1]*x[-1])})),2,sum,na.rm=TRUE)
        
        tmp2<-score[i,1:nvec[i]]*exp(-fit.tmp)*Tmat[i,1:nvec[i]]/(2)
        tmp3<-matrix(0,q1,q1)
        colnames(tmp3)=rownames(tmp3)<-X1.names
        for(j in 1:nvec[i])
        {
          tmp3<-tmp3+tmp2[j]*(X1[[i]][j,]%*%t(X1[[i]][j,]))
        }
        pt2<-pt2+tmp3
      }
      beta1.upd<-c(beta1.new-ginv(pt2)%*%pt1)
      
      # update gamma's
      gamma.rnd.upd<-matrix(NA,m,p)
      rownames(gamma.rnd.upd)<-cluster.names
      gamma.rnd.upd.std<-gamma.rnd.upd
      for(i in 1:m)
      {
        Atmp<-matrix(0,p,p)
        fit.tmp<-beta0.rnd.upd[i]+c(X1[[i]]%*%beta1.upd)
        for(j in 1:nvec[i])
        {
          Atmp<-Atmp+(Tmat[i,j]/2)*exp(-fit.tmp[j])*Y.cov[[i]][,,j]
        }
        
        gamma.est.tmp.out<-NULL
        try(gamma.est.tmp.out<-gamma.solve(A=Atmp,kappa=kappa.new,gamma=gamma.new,H=H.mat[,,i]),silent=TRUE)
        if(is.null(gamma.est.tmp.out)==FALSE)
        {
          gamma.est.tmp<-gamma.est.tmp.out$gamma
        }else
        {
          gamma.est.tmp<-c(gamma.solve.QP(A=Atmp,H=H.mat[,,i]))
        }
        if(gamma.est.tmp[which.max(abs(gamma.est.tmp))]<0)
        {
          gamma.est.tmp<--gamma.est.tmp
        }
        gamma.rnd.upd[i,]<-gamma.est.tmp
        gamma.rnd.upd.std[i,]<-gamma.est.tmp/sqrt(sum(gamma.est.tmp^2))
      }
      gamma.rnd.upd.avg<-apply(gamma.rnd.upd.std,2,mean,na.rm=TRUE)
      gamma.rnd.upd.avg.norm<-sqrt(sum(gamma.rnd.upd.avg^2))
      gamma.upd<-gamma.rnd.upd.avg/gamma.rnd.upd.avg.norm           # unit norm update
      kappa.upd<-gamma.rnd.upd.avg.norm*(p-gamma.rnd.upd.avg.norm^2)/(1-gamma.rnd.upd.avg.norm^2)
      
      # calculate converge criterion
      diff<-max(abs(c(beta0.upd-beta0.new,beta1.upd-beta1.new)))
      
      # calculate objective function
      obj.upd<-obj.func(X1=X1,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd.upd,gamma=gamma.upd,kappa=kappa.upd,beta0.rnd=beta0.rnd.upd,beta1=beta1.upd,beta0=beta0.upd,sigma2=sigma2.upd)
      
      # if(obj.upd<obj0&obj.upd>0)
      # {
      # update parameters
      beta0.rnd.new<-beta0.rnd.upd
      beta0.new<-beta0.upd
      sigma2.new<-sigma2.upd
      beta1.new<-beta1.upd
      
      gamma.rnd.new<-gamma.rnd.upd.std                # unit norm update
      # gamma.rnd.new<-gamma.rnd.upd                      # gamma update (without standardizing)
      gamma.new<-gamma.upd
      kappa.new<-kappa.upd
      
      obj0<-obj.upd
      
      if(trace)
      {
        beta1.trace<-cbind(beta1.trace,beta1.new)
        beta0.rnd.trace<-cbind(beta0.rnd.trace,beta0.rnd.new)
        beta0.trace<-cbind(beta0.trace,beta0.new)
        sigma2.trace<-cbind(sigma2.trace,sigma2.new)
        
        gamma.rnd.trace<-cbind(gamma.rnd.trace,gamma.rnd.new)
        gamma.trace<-cbind(gamma.trace,gamma.new)
        kappa.trace<-c(kappa.trace,kappa.new)
        
        obj.trace<-c(obj.trace,obj.upd)
      }
      
      # print(c(diff,obj.upd))
      # }else
      # {
      #   break
      # }
    }
    #-----------------------------
    
    #-----------------------------
    # standardize gamma estimate
    for(i in 1:m)
    {
      gamma.rnd.new[i,]<-gamma.rnd.new[i,]/sqrt(sum(gamma.rnd.new[i,]^2))
      if(gamma.rnd.new[i,which.max(abs(gamma.rnd.new[i,]))]<0)
      {
        gamma.rnd.new[i,]<--gamma.rnd.new[i,]
      }
    }
    gamma.rnd.new.avg<-apply(gamma.rnd.new,2,mean,na.rm=TRUE)
    gamma.rnd.new.avg.norm<-sqrt(sum(gamma.rnd.new.avg^2))
    gamma.new<-gamma.rnd.new.avg/gamma.rnd.new.avg.norm
    if(gamma.new[which.max(abs(gamma.new))]<0)
    {
      gamma.new<--gamma.new
    }
    kappa.new<-gamma.rnd.new.avg.norm*(p-gamma.rnd.new.avg.norm^2)/(1-gamma.rnd.new.avg.norm^2)
    
    # reestimate beta
    beta.out<-lcap.beta.cov.fix(Y.cov=Y.cov,X1=X1,Tmat=Tmat,gamma.rnd=gamma.rnd.new,max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE)
    beta.est<-beta.out$beta
    
    # calculate score
    score<-compute.scores(Y.cov,gamma.rnd.new)
    rownames(score)<-cluster.names
    colnames(score)<-paste0("unit",1:ncol(score))
    #-----------------------------
    
    if(trace)
    {
      colnames(beta1.trace)=colnames(beta0.rnd.trace)=colnames(beta0.trace)=colnames(sigma2.trace)=colnames(gamma.trace)<-paste0("iteration",0:(ncol(beta0.trace)-1))
      
      beta.trace<-rbind(beta0.trace,beta1.trace)
      rownames(beta.trace)[1]<-"Intercept"
      
      if(score.return)
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,
                 convergence=(s<max.itr),score=score,
                 gamma.trace=gamma.trace,gamma.rnd.trace=gamma.rnd.trace,kappa.trace=kappa.trace,beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,obj=obj.trace)
      }else
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,
                 convergence=(s<max.itr),
                 gamma.trace=gamma.trace,gamma.rnd.trace=gamma.rnd.trace,kappa.trace=kappa.trace,beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,obj=obj.trace)
      }
    }else
    {
      if(score.return)
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,
                 convergence=(s<max.itr),score=score)
      }else
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,convergence=(s<max.itr))
      }
    }
    
    return(re)
  }
}
lcap.cov.D1.fix.opt<-function(Y.cov,X1,Tmat,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,kappa0=NULL,seed=100)
{
  # Y.cov: a list of covariance matrix of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # Tmat: m by nmax matrix of sample size
  
  #-----------------------------
  m<-length(Y.cov)
  nvec<-sapply(Y.cov,dim)[3,]
  p<-sapply(Y.cov,dim)[1,1]
  Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
  
  if(is.null(names(Y.cov))==FALSE)
  {
    cluster.names<-names(Y.cov)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  q1<-ncol(X1[[1]])
  if(is.null(names(X1))==TRUE)
  {
    names(X1)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X1[[i]]))==TRUE)
    {
      X1.names<-paste0("X1",1:q1)
      colnames(X1[[i]])<-X1.names
    }else
    {
      X1.names<-colnames(X1[[i]])
    }
  }
  q2<-0
  q<-q1+1
  #-----------------------------
  
  if(method[1]=="CAP")
  {
    # set initial values
    if(is.null(gamma0.mat))
    {
      #--------------------------------
      set.seed(seed)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      #--------------------------------
    }
    if(is.null(ninitial))
    {
      ninitial<-min(ncol(gamma0.mat),10)
    }else
    {
      if(ninitial>ncol(gamma0.mat))
      {
        ninitial<-ncol(gamma0.mat)
      }
    }
    set.seed(seed)
    gamma0.mat<-matrix(gamma0.mat[,sort(sample(1:ncol(gamma0.mat),ninitial,replace=FALSE))],ncol=ninitial)
    
    re.tmp<-vector("list",ncol(gamma0.mat))
    obj<-rep(NA,ncol(gamma0.mat))
    for(kk in 1:ncol(gamma0.mat))
    {
      try(re.tmp[[kk]]<-lcap.cov.D1.fix(Y.cov=Y.cov,X1=X1,Tmat=Tmat,method=method,H.type=H.type,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0=gamma0.mat[,kk],kappa0=kappa0))
      
      if(is.null(re.tmp[[kk]])==FALSE)
      {
        if(re.tmp[[kk]]$convergence==TRUE)
        {
          obj[kk]<-obj.func(X1=X1,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=re.tmp[[kk]]$gamma.rnd,gamma=re.tmp[[kk]]$gamma,kappa=re.tmp[[kk]]$kappa,
                            beta0.rnd=re.tmp[[kk]]$beta0.rnd,beta1=re.tmp[[kk]]$beta[2:(q1+1)],beta0=re.tmp[[kk]]$beta[1],sigma2=re.tmp[[kk]]$beta0.sigma2)
        }
      }
    }
    opt.idx<-which.min(obj)
    re<-re.tmp[[opt.idx]]
  }
  
  return(re)
}

# only has random effects
lcap.cov.D1.rnd<-function(Y.cov,X2,Tmat,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0=NULL,kappa0=NULL)
{
  # Y.cov: a list of covariance matrix of Y
  # X2: a list of covariates corresponding to random-effect beta
  # Tmat: m by nmax matrix of sample size
  
  #-----------------------------
  m<-length(Y.cov)
  nvec<-sapply(Y.cov,dim)[3,]
  p<-sapply(Y.cov,dim)[1,1]
  Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
  
  if(is.null(names(Y.cov))==FALSE)
  {
    cluster.names<-names(Y.cov)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  q1<-0
  q2<-ncol(X2[[1]])
  if(is.null(names(X2))==TRUE)
  {
    names(X2)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X2[[i]]))==TRUE)
    {
      X2.names<-paste0("X2",1:q2)
      colnames(X2[[i]])<-X2.names
    }else
    {
      X2.names<-colnames(X2[[i]])
    }
  }
  q<-q2+1
  #-----------------------------
  
  #-----------------------------
  H.mat<-array(NA,c(p,p,m))
  dimnames(H.mat)[[3]]<-cluster.names
  for(i in 1:m)
  {
    if(H.type[1]=="CAvgCov"|H.type[1]=="AvgCov")
    {
      Htmp<-matrix(0,p,p)
      for(j in 1:nvec[i])
      {
        Htmp<-Htmp+Y.cov[[i]][,,j]*Tmat[i,j]
      }
      Htmp<-Htmp/sum(Tmat[i,1:nvec[i]])
      H.mat[,,i]<-Htmp
    }
    if(H.type[1]=="Identity")
    {
      H.mat[,,i]<-diag(rep(1,p))
    }
  }
  if(H.type[1]=="AvgCov")
  {
    Htmp<-matrix(0,p,p)
    for(i in 1:m)
    {
      Htmp<-Htmp+H.mat[,,i]*sum(Tmat[i,1:nvec[i]])/sum(Tmat,na.rm=TRUE)
    }
    for(i in 1:m)
    {
      H.mat[,,i]<-Htmp
    }
  }
  #-----------------------------
  
  if(method[1]=="CAP")
  {
    
    #-----------------------------
    # set initial value of gamma0
    if(is.null(gamma0))
    {
      set.seed(100)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      gamma0<-gamma0.mat[,sample(ncol(gamma0.mat),1)]
    }
    if(is.null(kappa0))
    {
      kappa0<-1
    }
    # generate random gamma
    gamma.new<-gamma0
    kappa.new<-kappa0
    gamma.rnd.new<-rvmf(m,mu=gamma.new,k=kappa.new) 
    rownames(gamma.rnd.new)<-cluster.names

    # calculate score
    score<-compute.scores(Y.cov,gamma.rnd.new)
    rownames(score)<-cluster.names
    colnames(score)<-paste0("unit",1:ncol(score))

    # Initial estimate: using a mixed effects model
    beta0.rnd.new<-rep(NA,m)
    beta2.rnd.new<-matrix(NA,m,q2)
    rownames(beta2.rnd.new)<-names(X2)
    colnames(beta2.rnd.new)<-X2.names
    for(i in 1:m)
    {
      dtmp<-data.frame(score=log(score[i,1:nvec[i]]),X2[[i]])
      fit.tmp<-lm(score~.,data=dtmp)
      beta0.rnd.new[i]<-fit.tmp$coefficients[1]
      beta2.rnd.new[i,]<-fit.tmp$coefficients[-1]
    }
    beta0.new<-mean(beta0.rnd.new,na.rm=TRUE)
    sigma2.new<-mean((beta0.rnd.new-beta0.new)^2,na.rm=TRUE)
    beta2.new<-apply(beta2.rnd.new,2,mean,na.rm=TRUE)
    Omega.new<-cov(beta2.rnd.new)*(m-1)/m
    if(Omega.diag==TRUE)
    {
      Omega.tmp<-matrix(0,q2,q2)
      colnames(Omega.tmp)=rownames(Omega.tmp)<-X2.names
      diag(Omega.tmp)<-diag(Omega.new)
      Omega.new<-Omega.tmp
    }
    
    obj0<-obj.func(X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd.new,gamma=gamma.new,kappa=kappa.new,
                   beta0.rnd=beta0.rnd.new,beta2.rnd=beta2.rnd.new,beta0=beta0.new,beta2=beta2.new,sigma2=sigma2.new,Omega=Omega.new)
    #-----------------------------
    
    #-----------------------------
    if(trace)
    {
      beta0.rnd.trace<-cbind(beta0.rnd.new)
      beta0.trace<-cbind(beta0.new)
      sigma2.trace<-cbind(sigma2.new)
      beta2.rnd.trace<-cbind(beta2.rnd.new)
      beta2.trace<-cbind(beta2.new)
      Omega.trace<-cbind(Omega.new)
      
      gamma.rnd.trace<-cbind(gamma.rnd.new)
      gamma.trace<-cbind(gamma.new)
      kappa.trace<-c(kappa.new)
      
      obj.trace<-c(obj0)
    }
    #-----------------------------
    
    #-----------------------------
    s<-0
    diff<-100
    while(s<=max.itr&diff>tol)
    {
      s<-s+1
      
      # calculate score
      score<-compute.scores(Y.cov,gamma.rnd.new)

      # update random beta0
      beta0.rnd.upd<-rep(NA,m)
      for(i in 1:m)
      {
        fit.tmp<-beta0.rnd.new[i]+X2[[i]]%*%beta2.rnd.new[i,]
        
        pt1<-sum((1-score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2))+(beta0.rnd.new[i]-beta0.new)/sigma2.new
        pt2<-sum(score[i,1:nvec[i]]*exp(-fit.tmp)*Tmat[i,1:nvec[i]]/(2))+1/sigma2.new
        
        beta0.rnd.upd[i]<-beta0.rnd.new[i]-pt1/pt2
      }
      # update beta0
      beta0.upd<-mean(beta0.rnd.upd,na.rm=TRUE)
      # update beta0 variance
      sigma2.upd<-mean((beta0.rnd.upd-beta0.upd)^2,na.rm=TRUE)
      
      # update beta2
      beta2.rnd.upd<-matrix(NA,m,q2)
      rownames(beta2.rnd.upd)<-names(X2)
      colnames(beta2.rnd.upd)<-X2.names
      for(i in 1:m)
      {
        fit.tmp<-beta0.rnd.upd[i]+X2[[i]]%*%beta2.rnd.new[i,]
        
        pt1<-rep(0,q2)
        pt2<-matrix(0,q2,q2)
        for(j in 1:nvec[i])
        {
          pt1<-pt1+((1-score[i,j]*exp(-fit.tmp[j]))*Tmat[i,j]/(2))*X2[[i]][j,]
          
          pt2<-pt2+(score[i,j]*exp(-fit.tmp[j])*Tmat[i,j]/(2))*(X2[[i]][j,]%*%t(X2[[i]][j,]))
        }
        pt1<-pt1+c(ginv(Omega.new)%*%(beta2.rnd.new[i,]-beta2.new))
        pt2<-pt2+ginv(Omega.new)
        
        beta2.rnd.upd[i,]<-beta2.rnd.new[i,]-ginv(pt2)%*%pt1
      }
      beta2.upd<-apply(beta2.rnd.upd,2,mean,na.rm=TRUE)
      Omega.upd<-cov(beta2.rnd.upd)*(m-1)/m
      if(Omega.diag==TRUE)
      {
        Omega.tmp<-matrix(0,q2,q2)
        colnames(Omega.tmp)=rownames(Omega.tmp)<-X2.names
        diag(Omega.tmp)<-diag(Omega.upd)
        Omega.upd<-Omega.tmp
      }
      
      # update gamma's
      gamma.rnd.upd<-matrix(NA,m,p)
      rownames(gamma.rnd.upd)<-cluster.names
      gamma.rnd.upd.std<-gamma.rnd.upd
      for(i in 1:m)
      {
        Atmp<-matrix(0,p,p)
        fit.tmp<-beta0.rnd.upd[i]+X2[[i]]%*%beta2.rnd.upd[i,]
        for(j in 1:nvec[i])
        {
          Atmp<-Atmp+(Tmat[i,j]/2)*exp(-fit.tmp[j])*Y.cov[[i]][,,j]
        }
        
        gamma.est.tmp.out<-NULL
        try(gamma.est.tmp.out<-gamma.solve(A=Atmp,kappa=kappa.new,gamma=gamma.new,H=H.mat[,,i]),silent=TRUE)
        if(is.null(gamma.est.tmp.out)==FALSE)
        {
          gamma.est.tmp<-gamma.est.tmp.out$gamma
        }else
        {
          gamma.est.tmp<-c(gamma.solve.QP(A=Atmp,H=H.mat[,,i]))
        }
        if(gamma.est.tmp[which.max(abs(gamma.est.tmp))]<0)
        {
          gamma.est.tmp<--gamma.est.tmp
        }
        gamma.rnd.upd[i,]<-gamma.est.tmp
        gamma.rnd.upd.std[i,]<-gamma.est.tmp/sqrt(sum(gamma.est.tmp^2))
      }
      gamma.rnd.upd.avg<-apply(gamma.rnd.upd.std,2,mean,na.rm=TRUE)
      gamma.rnd.upd.avg.norm<-sqrt(sum(gamma.rnd.upd.avg^2))
      gamma.upd<-gamma.rnd.upd.avg/gamma.rnd.upd.avg.norm           # unit norm update
      kappa.upd<-gamma.rnd.upd.avg.norm*(p-gamma.rnd.upd.avg.norm^2)/(1-gamma.rnd.upd.avg.norm^2)
      
      # calculate converge criterion
      diff<-max(abs(c(beta0.upd-beta0.new,beta2.upd-beta2.new)))
      
      # calculate objective function
      obj.upd<-obj.func(X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd.upd,gamma=gamma.upd,kappa=kappa.upd,
                        beta0.rnd=beta0.rnd.upd,beta2.rnd=beta2.rnd.upd,beta0=beta0.upd,beta2=beta2.upd,sigma2=sigma2.upd,Omega=Omega.upd)
      
      # if(obj.upd<obj0&obj.upd>0)
      # {
      # update parameters
      beta0.rnd.new<-beta0.rnd.upd
      beta0.new<-beta0.upd
      sigma2.new<-sigma2.upd
      beta2.rnd.new<-beta2.rnd.upd
      beta2.new<-beta2.upd
      Omega.new<-Omega.upd
      
      gamma.rnd.new<-gamma.rnd.upd.std                # unit norm update
      # gamma.rnd.new<-gamma.rnd.upd                      # gamma update (without standardizing)
      gamma.new<-gamma.upd
      kappa.new<-kappa.upd
      
      obj0<-obj.upd
      
      if(trace)
      {
        beta0.rnd.trace<-cbind(beta0.rnd.trace,beta0.rnd.new)
        beta0.trace<-cbind(beta0.trace,beta0.new)
        sigma2.trace<-cbind(sigma2.trace,sigma2.new)
        beta2.rnd.trace<-cbind(beta2.rnd.trace,beta2.rnd.new)
        beta2.trace<-cbind(beta2.trace,beta2.new)
        Omega.trace<-cbind(Omega.trace,Omega.new)
        
        gamma.rnd.trace<-cbind(gamma.rnd.trace,gamma.rnd.new)
        gamma.trace<-cbind(gamma.trace,gamma.new)
        kappa.trace<-c(kappa.trace,kappa.new)
        
        obj.trace<-c(obj.trace,obj.upd)
      }
      
      # print(c(diff,obj.upd))
      # }else
      # {
      #   break
      # }
    }
    #-----------------------------
    
    #-----------------------------
    # standardize gamma estimate
    for(i in 1:m)
    {
      gamma.rnd.new[i,]<-gamma.rnd.new[i,]/sqrt(sum(gamma.rnd.new[i,]^2))
      if(gamma.rnd.new[i,which.max(abs(gamma.rnd.new[i,]))]<0)
      {
        gamma.rnd.new[i,]<--gamma.rnd.new[i,]
      }
    }
    gamma.rnd.new.avg<-apply(gamma.rnd.new,2,mean,na.rm=TRUE)
    gamma.rnd.new.avg.norm<-sqrt(sum(gamma.rnd.new.avg^2))
    gamma.new<-gamma.rnd.new.avg/gamma.rnd.new.avg.norm
    if(gamma.new[which.max(abs(gamma.new))]<0)
    {
      gamma.new<--gamma.new
    }
    kappa.new<-gamma.rnd.new.avg.norm*(p-gamma.rnd.new.avg.norm^2)/(1-gamma.rnd.new.avg.norm^2)
    
    # reestimate beta
    beta.out<-lcap.beta.cov.rnd(Y.cov=Y.cov,X2=X2,Tmat=Tmat,gamma.rnd=gamma.rnd.new,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE)
    beta.est<-beta.out$beta

    # calculate score
    score<-compute.scores(Y.cov,gamma.rnd.new)
    rownames(score)<-cluster.names
    colnames(score)<-paste0("unit",1:ncol(score))
    #-----------------------------

    if(trace)
    {
      colnames(beta2.trace)=colnames(beta0.rnd.trace)=colnames(beta0.trace)=colnames(sigma2.trace)=colnames(gamma.trace)<-paste0("iteration",0:(ncol(beta0.trace)-1))

      beta.trace<-rbind(beta0.trace,beta2.trace)
      rownames(beta.trace)[1]<-"Intercept"

      if(score.return)
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,
                 beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,beta2.rnd=beta.out$beta2.rnd,beta2.Omega=beta.out$beta2.Omega,
                 convergence=(s<max.itr),score=score,
                 gamma.trace=gamma.trace,gamma.rnd.trace=gamma.rnd.trace,kappa.trace=kappa.trace,
                 beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,beta2.rnd.trace=beta2.rnd.trace,beta2.Omega.trace=Omega.trace,obj=obj.trace)
      }else
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,
                 beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,beta2.rnd=beta.out$beta2.rnd,beta2.Omega=beta.out$beta2.Omega,
                 convergence=(s<max.itr),
                 gamma.trace=gamma.trace,gamma.rnd.trace=gamma.rnd.trace,kappa.trace=kappa.trace,
                 beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,beta2.rnd.trace=beta2.rnd.trace,beta2.Omega.trace=Omega.trace,obj=obj.trace)
      }
    }else
    {
      if(score.return)
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,
                 beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,beta2.rnd=beta.out$beta2.rnd,beta2.Omega=beta.out$beta2.Omega,
                 convergence=(s<max.itr),score=score)
      }else
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,
                 beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,beta2.rnd=beta.out$beta2.rnd,beta2.Omega=beta.out$beta2.Omega,
                 convergence=(s<max.itr))
      }
    }

    return(re)
  }
}
lcap.cov.D1.rnd.opt<-function(Y.cov,X2,Tmat,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,
                              max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,kappa0=NULL,seed=100)
{
  # Y.cov: a list of covariance matrix of Y
  # X2: a list of covariates corresponding to random-effect beta
  # Tmat: m by nmax matrix of sample size
  
  #-----------------------------
  m<-length(Y.cov)
  nvec<-sapply(Y.cov,dim)[3,]
  p<-sapply(Y.cov,dim)[1,1]
  Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
  
  if(is.null(names(Y.cov))==FALSE)
  {
    cluster.names<-names(Y.cov)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  q1<-0
  q2<-ncol(X2[[1]])
  if(is.null(names(X2))==TRUE)
  {
    names(X2)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X2[[i]]))==TRUE)
    {
      X2.names<-paste0("X2",1:q2)
      colnames(X2[[i]])<-X2.names
    }else
    {
      X2.names<-colnames(X2[[i]])
    }
  }
  q<-q2+1
  #-----------------------------
  
  if(method[1]=="CAP")
  {
    # set initial values
    if(is.null(gamma0.mat))
    {
      #--------------------------------
      set.seed(seed)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      #--------------------------------
    }
    if(is.null(ninitial))
    {
      ninitial<-min(ncol(gamma0.mat),10)
    }else
    {
      if(ninitial>ncol(gamma0.mat))
      {
        ninitial<-ncol(gamma0.mat)
      }
    }
    set.seed(seed)
    gamma0.mat<-matrix(gamma0.mat[,sort(sample(1:ncol(gamma0.mat),ninitial,replace=FALSE))],ncol=ninitial)
    
    re.tmp<-vector("list",ncol(gamma0.mat))
    obj<-rep(NA,ncol(gamma0.mat))
    for(kk in 1:ncol(gamma0.mat))
    {
      try(re.tmp[[kk]]<-lcap.cov.D1.rnd(Y.cov=Y.cov,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,
                                        max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0=gamma0.mat[,kk],kappa0=kappa0))
      
      if(is.null(re.tmp[[kk]])==FALSE)
      {
        if(re.tmp[[kk]]$convergence==TRUE)
        {
          obj[kk]<-obj.func(X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=re.tmp[[kk]]$gamma.rnd,gamma=re.tmp[[kk]]$gamma,kappa=re.tmp[[kk]]$kappa,
                            beta0.rnd=re.tmp[[kk]]$beta0.rnd,beta2.rnd=re.tmp[[kk]]$beta2.rnd,beta0=re.tmp[[kk]]$beta[1],beta2=re.tmp[[kk]]$beta[(q1+1+1):q],
                            sigma2=re.tmp[[kk]]$beta0.sigma2,Omega=re.tmp[[kk]]$beta2.Omega)
        }
      }
    }
    opt.idx<-which.min(obj)
    re<-re.tmp[[opt.idx]]
  }
  
  return(re)
}

# has both fixed and random effects
lcap.cov.D1.mix<-function(Y.cov,X1,X2,Tmat,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0=NULL,kappa0=NULL)
{
  # Y.cov: a list of covariance matrix of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # Tmat: m by nmax matrix of sample size
  
  #-----------------------------
  m<-length(Y.cov)
  nvec<-sapply(Y.cov,dim)[3,]
  p<-sapply(Y.cov,dim)[1,1]
  Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
  
  if(is.null(names(Y.cov))==FALSE)
  {
    cluster.names<-names(Y.cov)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  q1<-ncol(X1[[1]])
  if(is.null(names(X1))==TRUE)
  {
    names(X1)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X1[[i]]))==TRUE)
    {
      X1.names<-paste0("X1",1:q1)
      colnames(X1[[i]])<-X1.names
    }else
    {
      X1.names<-colnames(X1[[i]])
    }
  }
  q2<-ncol(X2[[1]])
  if(is.null(names(X2))==TRUE)
  {
    names(X2)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X2[[i]]))==TRUE)
    {
      X2.names<-paste0("X2",1:q2)
      colnames(X2[[i]])<-X2.names
    }else
    {
      X2.names<-colnames(X2[[i]])
    }
  }
  q<-q1+q2+1
  #-----------------------------
  
  #-----------------------------
  H.mat<-array(NA,c(p,p,m))
  dimnames(H.mat)[[3]]<-cluster.names
  for(i in 1:m)
  {
    if(H.type[1]=="CAvgCov"|H.type[1]=="AvgCov")
    {
      Htmp<-matrix(0,p,p)
      for(j in 1:nvec[i])
      {
        Htmp<-Htmp+Y.cov[[i]][,,j]*Tmat[i,j]
      }
      Htmp<-Htmp/sum(Tmat[i,1:nvec[i]])
      H.mat[,,i]<-Htmp
    }
    if(H.type[1]=="Identity")
    {
      H.mat[,,i]<-diag(rep(1,p))
    }
  }
  if(H.type[1]=="AvgCov")
  {
    Htmp<-matrix(0,p,p)
    for(i in 1:m)
    {
      Htmp<-Htmp+H.mat[,,i]*sum(Tmat[i,1:nvec[i]])/sum(Tmat,na.rm=TRUE)
    }
    for(i in 1:m)
    {
      H.mat[,,i]<-Htmp
    }
  }
  #-----------------------------
  
  if(method[1]=="CAP")
  {
    
    #-----------------------------
    # set initial value of gamma0
    if(is.null(gamma0))
    {
      set.seed(100)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      gamma0<-gamma0.mat[,sample(ncol(gamma0.mat),1)]
    }
    if(is.null(kappa0))
    {
      kappa0<-1
    }
    # generate random gamma
    gamma.new<-gamma0
    kappa.new<-kappa0
    gamma.rnd.new<-rvmf(m,mu=gamma.new,k=kappa.new) 
    rownames(gamma.rnd.new)<-cluster.names
    
    # calculate score
    score<-compute.scores(Y.cov,gamma.rnd.new)
    rownames(score)<-cluster.names
    colnames(score)<-paste0("unit",1:ncol(score))
    
    # Initial estimate: using a mixed effects model
    dtmp<-NULL
    for(i in 1:m)
    {
      dtmp<-rbind(dtmp,data.frame(ID=rep(i,nvec[i]),score=log(score[i,1:nvec[i]]),X1[[i]],X2[[i]]))
    }
    colnames(dtmp)<-c("ID","score",X1.names,X2.names)
    fit.tmp<-NULL
    try(eval(parse(text=paste0("fit.tmp<-lme(score~",paste(c(X1.names,X2.names),collapse="+"),",random=~1+",paste(X2.names,collapse="+"),"|ID,data=dtmp,control=lmeControl(opt='optim'))"))),silent=TRUE)
    if(is.null(fit.tmp)==FALSE)
    {
      coef.fix<-fit.tmp$coefficients$fixed
      beta1.new<-coef.fix[2:(q1+1)]
      beta0.rnd.new<-coef.fix[1]+fit.tmp$coefficients$random$ID[,1]
      beta0.new<-mean(beta0.rnd.new,na.rm=TRUE)
      sigma2.new<-mean((beta0.rnd.new-beta0.new)^2,na.rm=TRUE)
      beta2.rnd.new<-matrix(NA,m,q2)
      colnames(beta2.rnd.new)<-X2.names
      rownames(beta2.rnd.new)<-names(X2)
      for(kk in 1:q2)
      {
        beta2.rnd.new[,kk]<-coef.fix[q1+1+kk]+fit.tmp$coefficients$random$ID[,kk+1]
      }
      beta2.new<-apply(beta2.rnd.new,2,mean,na.rm=TRUE)
      Omega.new<-cov(beta2.rnd.new)*(m-1)/m
      if(Omega.diag==TRUE)
      {
        Omega.tmp<-matrix(0,q2,q2)
        colnames(Omega.tmp)=rownames(Omega.tmp)<-X2.names
        diag(Omega.tmp)<-diag(Omega.new)
        Omega.new<-Omega.tmp
      }
    }else
    {
      dtmp<-NULL
      for(i in 1:m)
      {
        dtmp<-rbind(dtmp,data.frame(ID=rep(i,nvec[i]),score=log(score[i,1:nvec[i]]),X1[[i]],X2[[i]]))
      }
      colnames(dtmp)<-c("ID","score",X1.names,X2.names)
      # fixed effects
      eval(parse(text=paste0("fit.fix<-lm(score~",paste(c(X1.names),collapse="+"),",data=dtmp)")))
      beta1.new<-fit.fix$coefficients[-1]
      
      beta0.rnd.new<-rnorm(m,mean=0,sd=1)
      beta0.new<-mean(beta0.rnd.new,na.rm=TRUE)
      sigma2.new<-mean((beta0.rnd.new-beta0.new)^2,na.rm=TRUE)
      
      beta2.rnd.new<-matrix(rnorm(m*q2,mean=0,sd=1),m,q2)
      colnames(beta2.rnd.new)<-X2.names
      rownames(beta2.rnd.new)<-names(X2)
      beta2.new<-apply(beta2.rnd.new,2,mean,na.rm=TRUE)
      Omega.new<-cov(beta2.rnd.new)*(m-1)/m
      if(Omega.diag==TRUE)
      {
        Omega.tmp<-matrix(0,q2,q2)
        colnames(Omega.tmp)=rownames(Omega.tmp)<-X2.names
        diag(Omega.tmp)<-diag(Omega.new)
        Omega.new<-Omega.tmp
      }
    }
    
    obj0<-obj.func(X1=X1,X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd.new,gamma=gamma.new,kappa=kappa.new,
                   beta0.rnd=beta0.rnd.new,beta1=beta1.new,beta2.rnd=beta2.rnd.new,beta0=beta0.new,beta2=beta2.new,sigma2=sigma2.new,Omega=Omega.new)
    #-----------------------------
    
    #-----------------------------
    if(trace)
    {
      beta1.trace<-cbind(beta1.new)
      beta0.rnd.trace<-cbind(beta0.rnd.new)
      beta0.trace<-cbind(beta0.new)
      sigma2.trace<-cbind(sigma2.new)
      beta2.rnd.trace<-cbind(beta2.rnd.new)
      beta2.trace<-cbind(beta2.new)
      Omega.trace<-cbind(Omega.new)
      
      gamma.rnd.trace<-cbind(gamma.rnd.new)
      gamma.trace<-cbind(gamma.new)
      kappa.trace<-c(kappa.new)
      
      obj.trace<-c(obj0)
    }
    #-----------------------------
    
    #-----------------------------
    s<-0
    diff<-100
    while(s<=max.itr&diff>tol)
    {
      s<-s+1
      
      # calculate score
      score<-compute.scores(Y.cov,gamma.rnd.new)
      
      # update random beta0
      beta0.rnd.upd<-rep(NA,m)
      for(i in 1:m)
      {
        fit.tmp<-beta0.rnd.new[i]+X1[[i]]%*%beta1.new+X2[[i]]%*%beta2.rnd.new[i,]
        
        pt1<-sum((1-score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2))+(beta0.rnd.new[i]-beta0.new)/sigma2.new
        pt2<-sum(score[i,1:nvec[i]]*exp(-fit.tmp)*Tmat[i,1:nvec[i]]/(2))+1/sigma2.new
        
        beta0.rnd.upd[i]<-beta0.rnd.new[i]-pt1/pt2
      }
      # update beta0
      beta0.upd<-mean(beta0.rnd.upd,na.rm=TRUE)
      # update beta0 variance
      sigma2.upd<-mean((beta0.rnd.upd-beta0.upd)^2,na.rm=TRUE)
      
      # update beta1
      pt1<-rep(0,q1)
      pt2<-matrix(0,q1,q1)
      for(i in 1:m)
      {
        fit.tmp<-beta0.rnd.upd[i]+X1[[i]]%*%beta1.new+X2[[i]]%*%beta2.rnd.new[i,]
        
        tmp1<-(1-score[i,1:nvec[i]]*exp(-fit.tmp))*Tmat[i,1:nvec[i]]/(2)
        pt1<-pt1+apply(t(apply(cbind(tmp1,X1[[i]]),1,function(x){return(x[1]*x[-1])})),2,sum,na.rm=TRUE)
        
        tmp2<-score[i,1:nvec[i]]*exp(-fit.tmp)*Tmat[i,1:nvec[i]]/(2)
        tmp3<-matrix(0,q1,q1)
        colnames(tmp3)=rownames(tmp3)<-X1.names
        for(j in 1:nvec[i])
        {
          tmp3<-tmp3+tmp2[j]*(X1[[i]][j,]%*%t(X1[[i]][j,]))
        }
        pt2<-pt2+tmp3
      }
      beta1.upd<-c(beta1.new-ginv(pt2)%*%pt1)
      
      # update beta2
      beta2.rnd.upd<-matrix(NA,m,q2)
      rownames(beta2.rnd.upd)<-names(X2)
      colnames(beta2.rnd.upd)<-X2.names
      for(i in 1:m)
      {
        fit.tmp<-beta0.rnd.upd[i]+X1[[i]]%*%beta1.upd+X2[[i]]%*%beta2.rnd.new[i,]
        
        pt1<-rep(0,q2)
        pt2<-matrix(0,q2,q2)
        for(j in 1:nvec[i])
        {
          pt1<-pt1+((1-score[i,j]*exp(-fit.tmp[j]))*Tmat[i,j]/(2))*X2[[i]][j,]
          
          pt2<-pt2+(score[i,j]*exp(-fit.tmp[j])*Tmat[i,j]/(2))*(X2[[i]][j,]%*%t(X2[[i]][j,]))
        }
        pt1<-pt1+c(ginv(Omega.new)%*%(beta2.rnd.new[i,]-beta2.new))
        pt2<-pt2+ginv(Omega.new)
        
        beta2.rnd.upd[i,]<-beta2.rnd.new[i,]-ginv(pt2)%*%pt1
      }
      beta2.upd<-apply(beta2.rnd.upd,2,mean,na.rm=TRUE)
      Omega.upd<-cov(beta2.rnd.upd)*(m-1)/m
      if(Omega.diag==TRUE)
      {
        Omega.tmp<-matrix(0,q2,q2)
        colnames(Omega.tmp)=rownames(Omega.tmp)<-X2.names
        diag(Omega.tmp)<-diag(Omega.upd)
        Omega.upd<-Omega.tmp
      }
      
      # update gamma's
      gamma.rnd.upd<-matrix(NA,m,p)
      rownames(gamma.rnd.upd)<-cluster.names
      gamma.rnd.upd.std<-gamma.rnd.upd
      for(i in 1:m)
      {
        Atmp<-matrix(0,p,p)
        fit.tmp<-beta0.rnd.upd[i]+X1[[i]]%*%beta1.upd+X2[[i]]%*%beta2.rnd.upd[i,]
        for(j in 1:nvec[i])
        {
          Atmp<-Atmp+(Tmat[i,j]/2)*exp(-fit.tmp[j])*Y.cov[[i]][,,j]
        }
        
        gamma.est.tmp.out<-NULL
        try(gamma.est.tmp.out<-gamma.solve(A=Atmp,kappa=kappa.new,gamma=gamma.new,H=H.mat[,,i]),silent=TRUE)
        if(is.null(gamma.est.tmp.out)==FALSE)
        {
          gamma.est.tmp<-gamma.est.tmp.out$gamma
        }else
        {
          gamma.est.tmp<-c(gamma.solve.QP(A=Atmp,H=H.mat[,,i]))
        }
        if(gamma.est.tmp[which.max(abs(gamma.est.tmp))]<0)
        {
          gamma.est.tmp<--gamma.est.tmp
        }
        gamma.rnd.upd[i,]<-gamma.est.tmp
        gamma.rnd.upd.std[i,]<-gamma.est.tmp/sqrt(sum(gamma.est.tmp^2))
      }
      gamma.rnd.upd.avg<-apply(gamma.rnd.upd.std,2,mean,na.rm=TRUE)
      gamma.rnd.upd.avg.norm<-sqrt(sum(gamma.rnd.upd.avg^2))
      gamma.upd<-gamma.rnd.upd.avg/gamma.rnd.upd.avg.norm           # unit norm update
      kappa.upd<-gamma.rnd.upd.avg.norm*(p-gamma.rnd.upd.avg.norm^2)/(1-gamma.rnd.upd.avg.norm^2)
      
      # calculate converge criterion
      diff<-max(abs(c(beta0.upd-beta0.new,beta1.upd-beta1.new,beta2.upd-beta2.new)))
      
      # calculate objective function
      obj.upd<-obj.func(X1=X1,X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=gamma.rnd.upd,gamma=gamma.upd,kappa=kappa.upd,
                        beta0.rnd=beta0.rnd.upd,beta1=beta1.upd,beta2.rnd=beta2.rnd.upd,beta0=beta0.upd,beta2=beta2.upd,sigma2=sigma2.upd,Omega=Omega.upd)
      
      # if(obj.upd<obj0&obj.upd>0)
      # {
        # update parameters
        beta0.rnd.new<-beta0.rnd.upd
        beta0.new<-beta0.upd
        sigma2.new<-sigma2.upd
        beta1.new<-beta1.upd
        beta2.rnd.new<-beta2.rnd.upd
        beta2.new<-beta2.upd
        Omega.new<-Omega.upd
        
        gamma.rnd.new<-gamma.rnd.upd.std                # unit norm update
        # gamma.rnd.new<-gamma.rnd.upd                      # gamma update (without standardizing)
        gamma.new<-gamma.upd
        kappa.new<-kappa.upd
        
        obj0<-obj.upd
        
        if(trace)
        {
          beta1.trace<-cbind(beta1.trace,beta1.new)
          beta0.rnd.trace<-cbind(beta0.rnd.trace,beta0.rnd.new)
          beta0.trace<-cbind(beta0.trace,beta0.new)
          sigma2.trace<-cbind(sigma2.trace,sigma2.new)
          beta2.rnd.trace<-cbind(beta2.rnd.trace,beta2.rnd.new)
          beta2.trace<-cbind(beta2.trace,beta2.new)
          Omega.trace<-cbind(Omega.trace,Omega.new)
          
          gamma.rnd.trace<-cbind(gamma.rnd.trace,gamma.rnd.new)
          gamma.trace<-cbind(gamma.trace,gamma.new)
          kappa.trace<-c(kappa.trace,kappa.new)
          
          obj.trace<-c(obj.trace,obj.upd)
        }
        
        # print(c(diff,obj.upd))
      # }else
      # {
      #   break
      # }
    }
    #-----------------------------
    
    #-----------------------------
    # standardize gamma estimate
    for(i in 1:m)
    {
      gamma.rnd.new[i,]<-gamma.rnd.new[i,]/sqrt(sum(gamma.rnd.new[i,]^2))
      if(gamma.rnd.new[i,which.max(abs(gamma.rnd.new[i,]))]<0)
      {
        gamma.rnd.new[i,]<--gamma.rnd.new[i,]
      }
    }
    gamma.rnd.new.avg<-apply(gamma.rnd.new,2,mean,na.rm=TRUE)
    gamma.rnd.new.avg.norm<-sqrt(sum(gamma.rnd.new.avg^2))
    gamma.new<-gamma.rnd.new.avg/gamma.rnd.new.avg.norm
    if(gamma.new[which.max(abs(gamma.new))]<0)
    {
      gamma.new<--gamma.new
    }
    kappa.new<-gamma.rnd.new.avg.norm*(p-gamma.rnd.new.avg.norm^2)/(1-gamma.rnd.new.avg.norm^2)
    
    # reestimate beta
    beta.out<-lcap.beta.cov.mix(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,gamma.rnd=gamma.rnd.new,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE)
    beta.est<-beta.out$beta
    
    # calculate score
    score<-compute.scores(Y.cov,gamma.rnd.new)
    rownames(score)<-cluster.names
    colnames(score)<-paste0("unit",1:ncol(score))
    #-----------------------------
    
    if(trace)
    {
      colnames(beta2.trace)=colnames(beta1.trace)=colnames(beta0.rnd.trace)=colnames(beta0.trace)=colnames(sigma2.trace)=colnames(gamma.trace)<-paste0("iteration",0:(ncol(beta0.trace)-1))
      
      beta.trace<-rbind(beta0.trace,beta1.trace,beta2.trace)
      rownames(beta.trace)[1]<-"Intercept"
      
      if(score.return)
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,
                 beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,beta2.rnd=beta.out$beta2.rnd,beta2.Omega=beta.out$beta2.Omega,
                 convergence=(s<max.itr),score=score,
                 gamma.trace=gamma.trace,gamma.rnd.trace=gamma.rnd.trace,kappa.trace=kappa.trace,
                 beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,beta2.rnd.trace=beta2.rnd.trace,beta2.Omega.trace=Omega.trace,obj=obj.trace)
      }else
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,
                 beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,beta2.rnd=beta.out$beta2.rnd,beta2.Omega=beta.out$beta2.Omega,
                 convergence=(s<max.itr),
                 gamma.trace=gamma.trace,gamma.rnd.trace=gamma.rnd.trace,kappa.trace=kappa.trace,
                 beta.trace=beta.trace,beta0.rnd.trace=beta0.rnd.trace,beta0.sigma2.trace=sigma2.trace,beta2.rnd.trace=beta2.rnd.trace,beta2.Omega.trace=Omega.trace,obj=obj.trace)
      }
    }else
    {
      if(score.return)
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,
                 beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,beta2.rnd=beta.out$beta2.rnd,beta2.Omega=beta.out$beta2.Omega,
                 convergence=(s<max.itr),score=score)
      }else
      {
        re<-list(gamma=gamma.new,kappa=kappa.new,gamma.rnd=gamma.rnd.new,
                 beta=beta.est,beta0.rnd=beta.out$beta0.rnd,beta0.sigma2=beta.out$beta0.sigma2,beta2.rnd=beta.out$beta2.rnd,beta2.Omega=beta.out$beta2.Omega,
                 convergence=(s<max.itr))
      }
    }
    
    return(re)
  }
}
lcap.cov.D1.mix.opt<-function(Y.cov,X1,X2,Tmat,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,
                              max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,kappa0=NULL,seed=100)
{
  # Y.cov: a list of covariance matrix of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # Tmat: m by nmax matrix of sample size
  
  #-----------------------------
  m<-length(Y.cov)
  nvec<-sapply(Y.cov,dim)[3,]
  p<-sapply(Y.cov,dim)[1,1]
  Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
  
  if(is.null(names(Y.cov))==FALSE)
  {
    cluster.names<-names(Y.cov)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  q1<-ncol(X1[[1]])
  if(is.null(names(X1))==TRUE)
  {
    names(X1)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X1[[i]]))==TRUE)
    {
      X1.names<-paste0("X1",1:q1)
      colnames(X1[[i]])<-X1.names
    }else
    {
      X1.names<-colnames(X1[[i]])
    }
  }
  q2<-ncol(X2[[1]])
  if(is.null(names(X2))==TRUE)
  {
    names(X2)<-cluster.names
  }
  for(i in 1:m)
  {
    if(is.null(colnames(X2[[i]]))==TRUE)
    {
      X2.names<-paste0("X2",1:q2)
      colnames(X2[[i]])<-X2.names
    }else
    {
      X2.names<-colnames(X2[[i]])
    }
  }
  q<-q1+q2+1
  #-----------------------------
  
  if(method[1]=="CAP")
  {
    # set initial values
    if(is.null(gamma0.mat))
    {
      #--------------------------------
      set.seed(seed)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      #--------------------------------
    }
    if(is.null(ninitial))
    {
      ninitial<-min(ncol(gamma0.mat),10)
    }else
    {
      if(ninitial>ncol(gamma0.mat))
      {
        ninitial<-ncol(gamma0.mat)
      }
    }
    set.seed(seed)
    gamma0.mat<-matrix(gamma0.mat[,sort(sample(1:ncol(gamma0.mat),ninitial,replace=FALSE))],ncol=ninitial)
    
    re.tmp<-vector("list",ncol(gamma0.mat))
    obj<-rep(NA,ncol(gamma0.mat))
    for(kk in 1:ncol(gamma0.mat))
    {
      try(re.tmp[[kk]]<-lcap.cov.D1.mix(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,
                                        max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0=gamma0.mat[,kk],kappa0=kappa0))
      
      if(is.null(re.tmp[[kk]])==FALSE)
      {
        if(re.tmp[[kk]]$convergence==TRUE)
        {
          obj[kk]<-obj.func(X1=X1,X2=X2,Tmat=Tmat,Sigma=Y.cov,gamma.rnd=re.tmp[[kk]]$gamma.rnd,gamma=re.tmp[[kk]]$gamma,kappa=re.tmp[[kk]]$kappa,
                            beta0.rnd=re.tmp[[kk]]$beta0.rnd,beta1=re.tmp[[kk]]$beta[2:(q1+1)],beta2.rnd=re.tmp[[kk]]$beta2.rnd,beta0=re.tmp[[kk]]$beta[1],beta2=re.tmp[[kk]]$beta[(q1+1+1):q],
                            sigma2=re.tmp[[kk]]$beta0.sigma2,Omega=re.tmp[[kk]]$beta2.Omega)
        }
      }
    }
    opt.idx<-which.min(obj)
    re<-re.tmp[[opt.idx]]
  }
  
  return(re)
}

# compile the functions
lcap.cov.D1<-function(Y.cov,X1=NULL,X2=NULL,Tmat,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0=NULL,kappa0=NULL)
{
  # Y.cov: a list of covariance matrix of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # Tmat: m by nmax matrix of sample size
  
  if(is.null(X1)==FALSE&is.null(X2)==FALSE)
  {
    re<-lcap.cov.D1.mix(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0=gamma0,kappa0=kappa0)
  }else
    if(is.null(X1)==FALSE)
    {
      re<-lcap.cov.D1.fix(Y.cov=Y.cov,X1=X1,Tmat=Tmat,method=method,H.type=H.type,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0=gamma0,kappa0=kappa0)
    }else
      if(is.null(X2)==FALSE)
      {
        re<-lcap.cov.D1.rnd(Y.cov=Y.cov,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0=gamma0,kappa0=kappa0)
      }
  
  return(re)
}
lcap.cov.D1.opt<-function(Y.cov,X1=NULL,X2=NULL,Tmat,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,
                          max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,kappa0=NULL,seed=100)
{
  # Y.cov: a list of covariance matrix of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # Tmat: m by nmax matrix of sample size
  
  if(is.null(X1)==FALSE&is.null(X2)==FALSE)
  {
    re<-lcap.cov.D1.mix.opt(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,
                            gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)
  }else
    if(is.null(X1)==FALSE)
    {
      re<-lcap.cov.D1.fix.opt(Y.cov=Y.cov,X1=X1,Tmat=Tmat,method=method,H.type=H.type,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,
                              gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)
    }else
      if(is.null(X2)==FALSE)
      {
        re<-lcap.cov.D1.rnd.opt(Y.cov=Y.cov,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,
                                gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)
      }
  
  return(re)
}

lcap.D1<-function(Y,X1=NULL,X2=NULL,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0=NULL,kappa0=NULL)
{
  # Y: a list of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  
  #-----------------------------
  m<-length(Y)
  nvec<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  if(is.null(names(Y))==FALSE)
  {
    cluster.names<-names(Y)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  Y.cov<-vector("list",length=m)
  names(Y.cov)<-cluster.names
  Tmat<-matrix(NA,m,max(nvec))
  rownames(Tmat)<-cluster.names
  for(i in 1:m)
  {
    Y.cov[[i]]<-array(NA,c(p,p,nvec[i]))
    for(j in 1:nvec[i])
    {
      Tmat[i,j]<-nrow(Y[[i]][[j]])
      Y.cov[[i]][,,j]<-cov(Y[[i]][[j]])*(Tmat[i,j]-1)/Tmat[i,j]
    }
  }
  
  re<-lcap.cov.D1(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0=gamma0,kappa0=kappa0)
  
  return(re)
}
lcap.D1.opt<-function(Y,X1=NULL,X2=NULL,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,
                      max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,kappa0=NULL,seed=100)
{
  # Y: a list of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  
  #-----------------------------
  m<-length(Y)
  nvec<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  if(is.null(names(Y))==FALSE)
  {
    cluster.names<-names(Y)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  Y.cov<-vector("list",length=m)
  names(Y.cov)<-cluster.names
  Tmat<-matrix(NA,m,max(nvec))
  rownames(Tmat)<-cluster.names
  for(i in 1:m)
  {
    Y.cov[[i]]<-array(NA,c(p,p,nvec[i]))
    for(j in 1:nvec[i])
    {
      Tmat[i,j]<-nrow(Y[[i]][[j]])
      Y.cov[[i]][,,j]<-cov(Y[[i]][[j]])*(Tmat[i,j]-1)/Tmat[i,j]
    }
  }
  
  re<-lcap.cov.D1.opt(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,
                      max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)
  
  return(re)
}
#################################################

#################################################
# Higher order components
lcap.cov.Dk<-function(Y.cov,X1=NULL,X2=NULL,Tmat,Gamma0.rnd=NULL,beta0.rnd=NULL,method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,
                      max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,kappa0=NULL,seed=100)
{
  # Y.cov: a list of covariance matrix of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # Tmat: m by nmax matrix of sample size
  # Gamma0.rnd: identified components of each cluster
  # beta0.rnd: estimated random effect of intercept
  
  if(is.null(Gamma0.rnd))
  {
    return(lcap.cov.D1.opt(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,
                           score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed))
  }else
  {
    #-----------------------------
    m<-length(Y.cov)
    nvec<-sapply(Y.cov,dim)[3,]
    p<-sapply(Y.cov,dim)[1,1]
    Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
    
    if(is.null(names(Y.cov))==FALSE)
    {
      cluster.names<-names(Y.cov)
    }else
    {
      cluster.names<-paste0("G",1:m)
    }
    #-----------------------------
    
    #-----------------------------
    p0<-dim(Gamma0.rnd)[2]
    # find the null space of Gamma0.rnd and estimate beta
    Theta.rnd<-array(NA,c(p,p-p0,m))
    Pi.rnd<-array(NA,c(p,p,m))
    dimnames(Theta.rnd)[[3]]=dimnames(Pi.rnd)[[3]]<-cluster.names
    for(i in 1:m)
    {
      Theta.rnd[,,i]<-Null(Gamma0.rnd[,,i])
      Pi.rnd[,,i]<-cbind(Theta.rnd[,,i],Gamma0.rnd[,,i])
    }
    if(is.null(beta0.rnd)==TRUE)
    {
      # beta estimate
      beta0.est<-matrix(NA,m,p0)
      rownames(beta0.est)<-cluster.names
      for(kk in 1:p0)
      {
        otmp<-lcap.beta.cov(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,gamma.rnd=t(Gamma0.rnd[,kk,]),Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE)
        beta0.est[,kk]<-otmp$beta0.rnd
      }
    }else
    {
      beta0.est<-beta0.rnd
      rownames(beta0.est)<-cluster.names
    }
    # generate new covariance outcomes
    Y.cov.new<-vector("list",length=m)
    names(Y.cov.new)<-cluster.names
    for(i in 1:m)
    {
      Y.cov.new[[i]]<-array(NA,c(p,p,nvec[i]))
      for(j in 1:nvec[i])
      {
        Delta.tmp<-t(Pi.rnd[,,i])%*%Y.cov[[i]][,,j]%*%Pi.rnd[,,i]
        Delta.new<-matrix(0,p,p)
        Delta.new[1:(p-p0),1:(p-p0)]<-Delta.tmp[1:(p-p0),1:(p-p0)]
        Delta.new[(p-p0+1):p,(p-p0+1):p]<-diag(exp(beta0.est[i,]),nrow=p0)
        Y.cov.new[[i]][,,j]<-Pi.rnd[,,i]%*%Delta.new%*%t(Pi.rnd[,,i])
      }
    }
    #-----------------------------
    
    if(method[1]=="CAP")
    {
      re<-lcap.cov.D1.opt(Y.cov=Y.cov.new,X1=X1,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,
                          score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)
      
      # orthogonality
      re$orthogonal.rnd<-matrix(NA,m,p0)
      rownames(re$orthogonal.rnd)<-cluster.names
      colnames(re$orthogonal.rnd)<-paste0("C",1:p0)
      for(i in 1:m)
      {
        re$orthogonal.rnd[i,]<-re$gamma.rnd[i,]%*%Gamma0.rnd[,,i]
      }
      
      return(re)
    }
  }
}

lcap.Dk<-function(Y,X1=NULL,X2=NULL,Gamma0.rnd=NULL,beta0.rnd=NULL,data.type=c("Y","Y.cov"),method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,
                  max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,kappa0=NULL,seed=100)
{
  # Y: a list of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # Gamma0.rnd: identified components of each cluster
  # beta0.rnd: estimated random effect of intercept
  
  #-----------------------------
  m<-length(Y)
  nvec<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  if(is.null(names(Y))==FALSE)
  {
    cluster.names<-names(Y)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  #-----------------------------
  
  if(data.type=="Y")
  {
    if(is.null(Gamma0.rnd))
    {
      return(lcap.D1.opt(Y=Y,X1=X1,X2=X2,method=method,H.type=H.type,Omega.diag=Omega.diag,
                         max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed))
    }else
    {
      #-----------------------------
      p0<-dim(Gamma0.rnd)[2]
      if(is.null(beta0.rnd)==TRUE)
      {
        # beta estimate
        beta0.est<-matrix(NA,m,p0)
        rownames(beta0.est)<-cluster.names
        for(kk in 1:p0)
        {
          otmp<-lcap.beta(Y=Y,X1=X1,X2=X2,gamma.rnd=t(Gamma0.rnd[,kk,]),Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE)
          beta0.est[,kk]<-otmp$beta0.rnd
        }
      }else
      {
        beta0.est<-beta0.rnd
        rownames(beta0.est)<-cluster.names
      }
      Ynew<-vector("list",length=m)
      names(Ynew)<-cluster.names
      for(i in 1:m)
      {
        Ynew[[i]]<-vector("list",length=nvec[i])
        names(Ynew[[i]])<-names(Y[[i]])
        for(j in 1:nvec[i])
        {
          Ytmp<-Y[[i]][[j]]-Y[[i]][[j]]%*%Gamma0.rnd[,,i]%*%t(Gamma0.rnd[,,i])
          Ytmp.svd<-svd(Ytmp)
          Ynew[[i]][[j]]<-Ytmp.svd$u%*%diag(c(Ytmp.svd$d[1:(p-p0)],sqrt(exp(beta0.est[i,])*nrow(Y[[i]][[j]]))))%*%t(Ytmp.svd$v)
        }
      }
      #-----------------------------
      
      if(method[1]=="CAP")
      {
        re<-lcap.D1.opt(Y=Ynew,X1=X1,X2=X2,method=method,H.type=H.type,Omega.diag=Omega.diag,
                        max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)
        
        # orthogonality
        re$orthogonal.rnd<-matrix(NA,m,p0)
        rownames(re$orthogonal.rnd)<-cluster.names
        colnames(re$orthogonal.rnd)<-paste0("C",1:p0)
        for(i in 1:m)
        {
          re$orthogonal.rnd[i,]<-re$gamma.rnd[i,]%*%Gamma0.rnd[,,i]
        }
        
        return(re)
      }
    }
  }
  if(data.type=="Y.cov")
  {
    Y.cov<-vector("list",length=m)
    names(Y.cov)<-cluster.names
    Tmat<-matrix(NA,m,max(nvec))
    rownames(Tmat)<-cluster.names
    for(i in 1:m)
    {
      Y.cov[[i]]<-array(NA,c(p,p,nvec[i]))
      for(j in 1:nvec[i])
      {
        Tmat[i,j]<-nrow(Y[[i]][[j]])
        Y.cov[[i]][,,j]<-cov(Y[[i]][[j]])*(Tmat[i,j]-1)/Tmat[i,j]
      }
    }
    
    re<-lcap.cov.Dk(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,Gamma0.rnd=Gamma0.rnd,method=method,H.type=H.type,Omega.diag=Omega.diag,
                    max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)
    
    return(re)
    #-----------------------------
  }
}
#################################################

#################################################
# deviation of diagonalization
DfD.cov<-function(Y.cov,Tmat,Gamma.rnd)
{
  # Y.cov: a list of covariance matrix of Y
  # Tmat: m by nmax matrix of sample size
  # Gamma.rnd: identified components of each cluster
  
  if(is.null(dim(Gamma.rnd))|dim(Gamma.rnd)[2]==1)
  {
    return("Dimension of Gamma is less than 2!")
  }else
  {
    #-----------------------------
    m<-length(Y.cov)
    nvec<-sapply(Y.cov,dim)[3,]
    p<-sapply(Y.cov,dim)[1,1]
    Tvec<-apply(Tmat,1,sum,na.rm=TRUE)
    
    if(is.null(names(Y.cov))==FALSE)
    {
      cluster.names<-names(Y.cov)
    }else
    {
      cluster.names<-paste0("G",1:m)
    }
    #-----------------------------
    
    #-----------------------------
    nD<-dim(Gamma.rnd)[2]
    ntotal<-sum(Tmat,na.rm=TRUE)
    dfd.mat<-array(NA,c(m,max(nvec),nD))
    dimnames(dfd.mat)[[1]]<-cluster.names
    dimnames(dfd.mat)[[3]]<-paste0("C",1:nD)
    for(i in 1:m)
    {
      for(j in 1:nvec[i])
      {
        dfd.mat[i,j,1]<-1
        for(kk in 2:nD)
        {
          gamma.tmp<-Gamma.rnd[,1:kk,i]
          mat.tmp<-t(gamma.tmp)%*%Y.cov[[i]][,,j]%*%gamma.tmp
          dfd.mat[i,j,kk]<-(det(diag(diag(mat.tmp)))/det(mat.tmp))^(Tmat[i,j]/ntotal)
        }
      }
    }
    pmean<-apply(dfd.mat,3,function(x){return(prod(x,na.rm=TRUE))})
    #-----------------------------
    
    re<-list(DfD.avg=pmean,DfD.unit=dfd.mat)
    
    return(re)
  }
}
DfD<-function(Y,Gamma.rnd)
{
  # Y: a list of Y
  # Gamma.rnd: identified components of each cluster
  
  #-----------------------------
  m<-length(Y)
  nvec<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  if(is.null(names(Y))==FALSE)
  {
    cluster.names<-names(Y)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  Y.cov<-vector("list",length=m)
  names(Y.cov)<-cluster.names
  Tmat<-matrix(NA,m,max(nvec))
  rownames(Tmat)<-cluster.names
  for(i in 1:m)
  {
    Y.cov[[i]]<-array(NA,c(p,p,nvec[i]))
    for(j in 1:nvec[i])
    {
      Tmat[i,j]<-nrow(Y[[i]][[j]])
      Y.cov[[i]][,,j]<-cov(Y[[i]][[j]])*(Tmat[i,j]-1)/Tmat[i,j]
    }
  }
  #-----------------------------
  
  re<-DfD.cov(Y.cov=Y.cov,Tmat=Tmat,Gamma.rnd=Gamma.rnd)
  
  return(re)
}
#################################################

#################################################
# lcap function: longitudinal cap regression
lcapReg<-function(Y,X1=NULL,X2=NULL,stop.crt=c("DfD","nD"),DfD.thred=2,nD=NULL,data.type=c("Y","Y.cov"),method=c("CAP"),H.type=c("CAvgCov","AvgCov","Identity"),Omega.diag=TRUE,
                  max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,kappa0=NULL,seed=100,verbose=TRUE)
{
  # Y: a list of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  
  #-----------------------------
  m<-length(Y)
  nvec<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  if(is.null(names(Y))==FALSE)
  {
    cluster.names<-names(Y)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }
  
  if(is.null(X1)==FALSE)
  {
    q1<-ncol(X1[[1]])
    if(is.null(names(X1))==TRUE)
    {
      names(X1)<-cluster.names
    }
    for(i in 1:m)
    {
      if(is.null(colnames(X1[[i]]))==TRUE)
      {
        X1.names<-paste0("X1",1:q1)
        colnames(X1[[i]])<-X1.names
      }else
      {
        X1.names<-colnames(X1[[i]])
      }
    }
  }
  if(is.null(X2)==FALSE)
  {
    q2<-ncol(X2[[1]])
    if(is.null(names(X2))==TRUE)
    {
      names(X2)<-cluster.names
    }
    for(i in 1:m)
    {
      if(is.null(colnames(X2[[i]]))==TRUE)
      {
        X2.names<-paste0("X2",1:q2)
        colnames(X2[[i]])<-X2.names
      }else
      {
        X2.names<-colnames(X2[[i]])
      }
    }
  }
  #-----------------------------
  
  #-----------------------------
  if(stop.crt[1]=="nD"&is.null(nD))
  {
    stop.crt<-"DfD"
  }
  #-----------------------------
  
  #-----------------------------
  if(method[1]=="CAP")
  {
    if(data.type=="Y")
    {
      #-----------------------------
      # First component
      tm1<-system.time(re1<-lcap.D1.opt(Y,X1=X1,X2=X2,method=method,H.type=H.type,Omega.diag=Omega.diag,
                                        max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed))
      
      Gamma.est<-matrix(re1$gamma,ncol=1)
      Gamma.rnd.est<-array(NA,c(p,1,m))
      dimnames(Gamma.rnd.est)[[3]]<-cluster.names
      Gamma.rnd.est[,1,]<-t(re1$gamma.rnd)
      kappa.est<-matrix(re1$kappa,ncol=1)
      rownames(kappa.est)<-"kappa"
      beta.est<-matrix(re1$beta,ncol=1)
      rownames(beta.est)<-names(re1$beta)
      beta0.rnd.est<-matrix(re1$beta0.rnd,ncol=1)
      rownames(beta0.rnd.est)<-cluster.names
      beta0.sigma2.est<-matrix(re1$beta0.sigma2,ncol=1)
      rownames(beta0.sigma2.est)<-"beta0.sigma2"
      if(is.null(X2)==FALSE)
      {
        beta2.rnd.est<-array(re1$beta2.rnd,c(m,q2,1))
        dimnames(beta2.rnd.est)[[1]]<-cluster.names
        dimnames(beta2.rnd.est)[[2]]<-X2.names
        
        beta2.Omega.est<-array(re1$beta2.Omega,c(q2,q2,1))
        dimnames(beta2.Omega.est)[[1]]=dimnames(beta2.Omega.est)[[2]]<-X2.names
      }
      
      cp.time<-matrix(as.numeric(tm1[1:3]),ncol=1)
      rownames(cp.time)<-c("user","system","elapsed")
      
      if(verbose)
      {
        print(paste0("Component ",ncol(Gamma.est)))
      }
      #-----------------------------
      
      #-----------------------------
      # Higher-order component
      if(stop.crt[1]=="DfD")
      {
        nD<-1
        if(score.return)
        {
          score<-array(NA,c(m,max(nvec),nD))
          dimnames(score)[[1]]<-cluster.names
          dimnames(score)[[3]]<-paste0("C",1:nD)
          
          score[,,1]<-re1$score
          score.tmp<-score
        }
        
        DfD.tmp<-1
        while(DfD.tmp<DfD.thred)
        {
          re.tmp<-NULL
          try(tm.tmp<-system.time(re.tmp<-lcap.Dk(Y,X1=X1,X2=X2,Gamma0.rnd=Gamma.rnd.est,beta0.rnd=beta0.rnd.est,data.type="Y",method=method,H.type=H.type,Omega.diag=Omega.diag,
                                                  max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)))
          
          if(is.null(re.tmp)==FALSE)
          {
            nD<-nD+1
            
            Gamma.rnd.new<-array(NA,c(p,nD,m))
            dimnames(Gamma.rnd.new)[[3]]<-cluster.names
            Gamma.rnd.new[,1:(nD-1),]<-Gamma.rnd.est
            Gamma.rnd.new[,nD,]<-t(re.tmp$gamma.rnd)
            
            DfD.out<-DfD(Y,Gamma.rnd.new)
            DfD.tmp<-DfD.out$DfD.avg[nD]
            
            if(DfD.tmp<DfD.thred)
            {
              Gamma.est<-cbind(Gamma.est,re.tmp$gamma)
              Gamma.rnd.est<-Gamma.rnd.new
              kappa.est<-cbind(kappa.est,re.tmp$kappa)
              beta.est<-cbind(beta.est,re.tmp$beta)
              beta0.rnd.est<-cbind(beta0.rnd.est,re.tmp$beta0.rnd)
              beta0.sigma2.est<-cbind(beta0.sigma2.est,re.tmp$beta0.sigma2)
              if(is.null(X2)==FALSE)
              {
                beta2.rnd.tmp<-array(NA,c(m,q2,nD))
                dimnames(beta2.rnd.tmp)[[1]]<-cluster.names
                dimnames(beta2.rnd.tmp)[[2]]<-X2.names
                dimnames(beta2.rnd.tmp)[[3]]<-paste0("C",1:nD)
                beta2.rnd.tmp[,,1:(nD-1)]<-beta2.rnd.est
                beta2.rnd.tmp[,,nD]<-re.tmp$beta2.rnd
                beta2.rnd.est<-beta2.rnd.tmp
                
                beta2.Omega.est.tmp<-beta2.Omega.est
                beta2.Omega.est<-array(NA,c(q2,q2,nD))
                dimnames(beta2.Omega.est)[[1]]=dimnames(beta2.Omega.est)[[2]]<-X2.names
                beta2.Omega.est[,,1:(nD-1)]<-beta2.Omega.est.tmp
                beta2.Omega.est[,,nD]<-re.tmp$beta2.Omega
              }
              
              cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
              
              if(verbose)
              {
                print(paste0("Component ",nD))
              }
              
              if(score.return)
              {
                score<-array(NA,c(m,max(nvec),nD))
                dimnames(score)[[1]]<-cluster.names
                
                score[,,1:(nD-1)]<-score.tmp
                score[,,nD]<-re.tmp$score
                
                score.tmp<-score
              }
            }
          }else
          {
            break
          }
        }
        
        nD<-ncol(Gamma.est)
        colnames(Gamma.est)=dimnames(Gamma.rnd.est)[[2]]=colnames(kappa.est)=colnames(beta.est)=colnames(beta0.rnd.est)=colnames(beta0.sigma2.est)<-paste0("C",1:nD)
        cp.time<-cbind(cp.time,apply(cp.time,1,sum))
        colnames(cp.time)<-c(paste0("C",1:nD),"Total")
        if(is.null(X2)==FALSE)
        {
          dimnames(beta2.rnd.est)[[3]]=dimnames(beta2.Omega.est)[[3]]<-paste0("C",1:nD)
        }
        if(is.null(score))
        {
          dimnames(score)[[3]]<-paste0("C",1:nD)
        }
        
        if(nD>1)
        {
          DfD.out<-DfD(Y,Gamma.rnd.est)
        }else
        {
          DfD.out<-list(DfD.avg=1)
        }
      }
      if(stop.crt[1]=="nD")
      {
        for(kk in 2:nD)
        {
          if(score.return)
          {
            score<-array(NA,c(m,max(nvec),nD))
            dimnames(score)[[1]]<-cluster.names
            
            score[,,1]<-re1$score
          }
          
          re.tmp<-NULL
          try(tm.tmp<-system.time(re.tmp<-lcap.Dk(Y,X1=X1,X2=X2,Gamma0.rnd=Gamma.rnd.est,beta0.rnd=beta0.rnd.est,data.type="Y",method=method,H.type=H.type,Omega.diag=Omega.diag,
                                                  max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)))
          
          if(is.null(re.tmp)==FALSE)
          {
            Gamma.est<-cbind(Gamma.est,re.tmp$gamma)
            Gamma.rnd.tmp<-Gamma.rnd.est
            Gamma.rnd.est<-array(NA,c(p,ncol(Gamma.est),m))
            dimnames(Gamma.rnd.est)[[3]]<-cluster.names
            Gamma.rnd.est[,1:(ncol(Gamma.est)-1),]<-Gamma.rnd.tmp
            Gamma.rnd.est[,ncol(Gamma.est),]<-t(re.tmp$gamma.rnd)
            kappa.est<-cbind(kappa.est,re.tmp$kappa)
            beta.est<-cbind(beta.est,re.tmp$beta)
            beta0.rnd.est<-cbind(beta0.rnd.est,re.tmp$beta0.rnd)
            beta0.sigma2.est<-cbind(beta0.sigma2.est,re.tmp$beta0.sigma2)
            if(is.null(X2)==FALSE)
            {
              beta2.rnd.tmp<-array(NA,c(m,q2,ncol(Gamma.est)))
              dimnames(beta2.rnd.tmp)[[1]]<-cluster.names
              dimnames(beta2.rnd.tmp)[[2]]<-X2.names
              dimnames(beta2.rnd.tmp)[[3]]<-paste0("C",1:ncol(Gamma.est))
              beta2.rnd.tmp[,,1:(ncol(Gamma.est)-1)]<-beta2.rnd.est
              beta2.rnd.tmp[,,ncol(Gamma.est)]<-re.tmp$beta2.rnd
              beta2.rnd.est<-beta2.rnd.tmp
              
              beta2.Omega.est.tmp<-beta2.Omega.est
              beta2.Omega.est<-array(NA,c(q2,q2,ncol(Gamma.est)))
              dimnames(beta2.Omega.est)[[1]]=dimnames(beta2.Omega.est)[[2]]<-X2.names
              beta2.Omega.est[,,1:(ncol(Gamma.est)-1)]<-beta2.Omega.est.tmp
              beta2.Omega.est[,,ncol(Gamma.est)]<-re.tmp$beta2.Omega
            }
            
            cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
            
            if(verbose)
            {
              print(paste0("Component ",ncol(Gamma.est)))
            }
            
            if(score.return)
            {
              score[,,kk]<-re.tmp$score
            }
          }else
          {
            break
          }
        }
        
        nD<-ncol(Gamma.est)
        colnames(Gamma.est)=dimnames(Gamma.rnd.est)[[2]]=colnames(kappa.est)=colnames(beta.est)=colnames(beta0.rnd.est)=colnames(beta0.sigma2.est)<-paste0("C",1:nD)
        cp.time<-cbind(cp.time,apply(cp.time,1,sum))
        colnames(cp.time)<-c(paste0("C",1:nD),"Total")
        if(is.null(X2)==FALSE)
        {
          dimnames(beta2.rnd.est)[[3]]=dimnames(beta2.Omega.est)[[3]]<-paste0("C",1:nD)
        }
        if(is.null(score))
        {
          dimnames(score)[[3]]<-paste0("C",1:nD)
        }
        
        if(nD>1)
        {
          DfD.out<-DfD(Y,Gamma.rnd.est)
        }else
        {
          DfD.out<-list(DfD.avg=1)
        }
      }
      #-----------------------------
    }
    if(data.type=="Y.cov")
    {
      Y.cov<-vector("list",length=m)
      names(Y.cov)<-cluster.names
      Tmat<-matrix(NA,m,max(nvec))
      rownames(Tmat)<-cluster.names
      for(i in 1:m)
      {
        Y.cov[[i]]<-array(NA,c(p,p,nvec[i]))
        for(j in 1:nvec[i])
        {
          Tmat[i,j]<-nrow(Y[[i]][[j]])
          Y.cov[[i]][,,j]<-cov(Y[[i]][[j]])*(Tmat[i,j]-1)/Tmat[i,j]
        }
      }
      
      #-----------------------------
      # First component
      tm1<-system.time(re1<-lcap.cov.D1.opt(Y.cov=Y.cov,X1=X1,X2=X2,Tmat=Tmat,method=method,H.type=H.type,Omega.diag=Omega.diag,
                                            max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed))
      
      Gamma.est<-matrix(re1$gamma,ncol=1)
      Gamma.rnd.est<-array(NA,c(p,1,m))
      dimnames(Gamma.rnd.est)[[3]]<-cluster.names
      Gamma.rnd.est[,1,]<-t(re1$gamma.rnd)
      kappa.est<-matrix(re1$kappa,ncol=1)
      rownames(kappa.est)<-"kappa"
      beta.est<-matrix(re1$beta,ncol=1)
      rownames(beta.est)<-names(re1$beta)
      beta0.rnd.est<-matrix(re1$beta0.rnd,ncol=1)
      rownames(beta0.rnd.est)<-cluster.names
      beta0.sigma2.est<-matrix(re1$beta0.sigma2,ncol=1)
      rownames(beta0.sigma2.est)<-"beta0.sigma2"
      if(is.null(X2)==FALSE)
      {
        beta2.rnd.est<-array(re1$beta2.rnd,c(m,q2,1))
        dimnames(beta2.rnd.est)[[1]]<-cluster.names
        dimnames(beta2.rnd.est)[[2]]<-X2.names
        
        beta2.Omega.est<-array(re1$beta2.Omega,c(q2,q2,1))
        dimnames(beta2.Omega.est)[[1]]=dimnames(beta2.Omega.est)[[2]]<-X2.names
      }
      
      cp.time<-matrix(as.numeric(tm1[1:3]),ncol=1)
      rownames(cp.time)<-c("user","system","elapsed")
      
      if(verbose)
      {
        print(paste0("Component ",ncol(Gamma.est)))
      }
      #-----------------------------
      
      #-----------------------------
      # Higher-order component
      if(stop.crt[1]=="DfD")
      {
        nD<-1
        if(score.return)
        {
          score<-array(NA,c(m,max(nvec),nD))
          dimnames(score)[[1]]<-cluster.names
          
          score[,,1]<-re1$score
          score.tmp<-score
        }
        
        DfD.tmp<-1
        while(DfD.tmp<DfD.thred)
        {
          re.tmp<-NULL
          try(tm.tmp<-system.time(re.tmp<-lcap.cov.Dk(Y=Y.cov,X1=X1,X2=X2,Tmat=Tmat,Gamma0.rnd=Gamma.rnd.est,beta0.rnd=beta0.rnd.est,method=method,H.type=H.type,Omega.diag=Omega.diag,
                                                      max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)))
          
          if(is.null(re.tmp)==FALSE)
          {
            nD<-nD+1
            
            Gamma.rnd.new<-array(NA,c(p,nD,m))
            dimnames(Gamma.rnd.new)[[3]]<-cluster.names
            Gamma.rnd.new[,1:(nD-1),]<-Gamma.rnd.est
            Gamma.rnd.new[,nD,]<-t(re.tmp$gamma.rnd)
            
            DfD.out<-DfD(Y,Gamma.rnd.new)
            DfD.tmp<-DfD.out$DfD.avg[nD]
            
            if(DfD.tmp<DfD.thred)
            {
              Gamma.est<-cbind(Gamma.est,re.tmp$gamma)
              Gamma.rnd.est<-Gamma.rnd.new
              kappa.est<-cbind(kappa.est,re.tmp$kappa)
              beta.est<-cbind(beta.est,re.tmp$beta)
              beta0.rnd.est<-cbind(beta0.rnd.est,re.tmp$beta0.rnd)
              beta0.sigma2.est<-cbind(beta0.sigma2.est,re.tmp$beta0.sigma2)
              if(is.null(X2)==FALSE)
              {
                beta2.rnd.tmp<-array(NA,c(m,q2,nD))
                dimnames(beta2.rnd.tmp)[[1]]<-cluster.names
                dimnames(beta2.rnd.tmp)[[2]]<-X2.names
                dimnames(beta2.rnd.tmp)[[3]]<-paste0("C",1:nD)
                beta2.rnd.tmp[,,1:(nD-1)]<-beta2.rnd.est
                beta2.rnd.tmp[,,nD]<-re.tmp$beta2.rnd
                beta2.rnd.est<-beta2.rnd.tmp
                
                beta2.Omega.est.tmp<-beta2.Omega.est
                beta2.Omega.est<-array(NA,c(q2,q2,nD))
                dimnames(beta2.Omega.est)[[1]]=dimnames(beta2.Omega.est)[[2]]<-X2.names
                beta2.Omega.est[,,1:(nD-1)]<-beta2.Omega.est.tmp
                beta2.Omega.est[,,nD]<-re.tmp$beta2.Omega
              }
              
              cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
              
              if(verbose)
              {
                print(paste0("Component ",nD))
              }
              
              if(score.return)
              {
                score<-array(NA,c(m,max(nvec),nD))
                dimnames(score)[[1]]<-cluster.names
                
                score[,,1:(nD-1)]<-score.tmp
                score[,,nD]<-re.tmp$score
                
                score.tmp<-score
              }
            }
          }else
          {
            break
          }
        }
        
        nD<-ncol(Gamma.est)
        colnames(Gamma.est)=dimnames(Gamma.rnd.est)[[2]]=colnames(kappa.est)=colnames(beta.est)=colnames(beta0.rnd.est)=colnames(beta0.sigma2.est)<-paste0("C",1:nD)
        cp.time<-cbind(cp.time,apply(cp.time,1,sum))
        colnames(cp.time)<-c(paste0("C",1:nD),"Total")
        if(is.null(X2)==FALSE)
        {
          dimnames(beta2.rnd.est)[[3]]=dimnames(beta2.Omega.est)[[3]]<-paste0("C",1:nD)
        }
        if(is.null(score))
        {
          dimnames(score)[[3]]<-paste0("C",1:nD)
        }
        
        if(nD>1)
        {
          DfD.out<-DfD(Y,Gamma.rnd.est)
        }else
        {
          DfD.out<-list(DfD.avg=1)
        }
      }
      if(stop.crt[1]=="nD")
      {
        for(kk in 2:nD)
        {
          if(score.return)
          {
            score<-array(NA,c(m,max(nvec),nD))
            dimnames(score)[[1]]<-cluster.names
            
            score[,,1]<-re1$score
          }
          
          re.tmp<-NULL
          try(tm.tmp<-system.time(re.tmp<-lcap.cov.Dk(Y=Y.cov,X1=X1,X2=X2,Tmat=Tmat,Gamma0.rnd=Gamma.rnd.est,beta0.rnd=beta0.rnd.est,method=method,H.type=H.type,Omega.diag=Omega.diag,
                                                      max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,kappa0=kappa0,seed=seed)))
          
          if(is.null(re.tmp)==FALSE)
          {
            Gamma.est<-cbind(Gamma.est,re.tmp$gamma)
            Gamma.rnd.tmp<-Gamma.rnd.est
            Gamma.rnd.est<-array(NA,c(p,ncol(Gamma.est),m))
            dimnames(Gamma.rnd.est)[[3]]<-cluster.names
            Gamma.rnd.est[,1:(ncol(Gamma.est)-1),]<-Gamma.rnd.tmp
            Gamma.rnd.est[,ncol(Gamma.est),]<-t(re.tmp$gamma.rnd)
            kappa.est<-cbind(kappa.est,re.tmp$kappa)
            beta.est<-cbind(beta.est,re.tmp$beta)
            beta0.rnd.est<-cbind(beta0.rnd.est,re.tmp$beta0.rnd)
            beta0.sigma2.est<-cbind(beta0.sigma2.est,re.tmp$beta0.sigma2)
            if(is.null(X2)==FALSE)
            {
              beta2.rnd.tmp<-array(NA,c(m,q2,ncol(Gamma.est)))
              dimnames(beta2.rnd.tmp)[[1]]<-cluster.names
              dimnames(beta2.rnd.tmp)[[2]]<-X2.names
              dimnames(beta2.rnd.tmp)[[3]]<-paste0("C",1:ncol(Gamma.est))
              beta2.rnd.tmp[,,1:(ncol(Gamma.est)-1)]<-beta2.rnd.est
              beta2.rnd.tmp[,,ncol(Gamma.est)]<-re.tmp$beta2.rnd
              beta2.rnd.est<-beta2.rnd.tmp
              
              beta2.Omega.est.tmp<-beta2.Omega.est
              beta2.Omega.est<-array(NA,c(q2,q2,ncol(Gamma.est)))
              dimnames(beta2.Omega.est)[[1]]=dimnames(beta2.Omega.est)[[2]]<-X2.names
              beta2.Omega.est[,,1:(ncol(Gamma.est)-1)]<-beta2.Omega.est.tmp
              beta2.Omega.est[,,ncol(Gamma.est)]<-re.tmp$beta2.Omega
            }
            
            cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
            
            if(verbose)
            {
              print(paste0("Component ",ncol(Gamma.est)))
            }
            
            if(score.return)
            {
              score[,,kk]<-re.tmp$score
            }
          }else
          {
            break
          }
        }
        
        nD<-ncol(Gamma.est)
        colnames(Gamma.est)=dimnames(Gamma.rnd.est)[[2]]=colnames(kappa.est)=colnames(beta.est)=colnames(beta0.rnd.est)=colnames(beta0.sigma2.est)<-paste0("C",1:nD)
        cp.time<-cbind(cp.time,apply(cp.time,1,sum))
        colnames(cp.time)<-c(paste0("C",1:nD),"Total")
        if(is.null(X2)==FALSE)
        {
          dimnames(beta2.rnd.est)[[3]]=dimnames(beta2.Omega.est)[[3]]<-paste0("C",1:nD)
        }
        if(is.null(score))
        {
          dimnames(score)[[3]]<-paste0("C",1:nD)
        }
        
        if(nD>1)
        {
          DfD.out<-DfD(Y,Gamma.rnd.est)
        }else
        {
          DfD.out<-list(DfD.avg=1)
        }
      }
      #-----------------------------
    }
    
    Gamma.rnd.otho<-array(NA,c(nD,nD,m))
    dimnames(Gamma.rnd.otho)[[1]]=dimnames(Gamma.rnd.otho)[[2]]<-paste0("C",1:nD)
    dimnames(Gamma.rnd.otho)[[3]]<-cluster.names
    for(i in 1:m)
    {
      Gamma.rnd.otho[,,i]<-t(Gamma.rnd.est[,,i])%*%Gamma.rnd.est[,,i]
    }
    Gamma.otho<-t(Gamma.est)%*%Gamma.est
    rownames(Gamma.otho)=colnames(Gamma.otho)<-paste0("C",1:nD)
    
    if(score.return)
    {
      if(is.null(X2)==FALSE)
      {
        re<-list(gamma=Gamma.est,kappa=kappa.est,gamma.rnd=Gamma.rnd.est,beta=beta.est,beta0.rnd=beta0.rnd.est,beta0.sigma2=beta0.sigma2.est,beta2.rnd=beta2.rnd.est,beta2.Omega=beta2.Omega.est,
                 DfD=DfD.out,gamma.othogonality=Gamma.otho,gamma.rnd.othogonality=Gamma.rnd.otho,score=score,time=cp.time)
      }else
      {
        re<-list(gamma=Gamma.est,kappa=kappa.est,gamma.rnd=Gamma.rnd.est,beta=beta.est,beta0.rnd=beta0.rnd.est,beta0.sigma2=beta0.sigma2.est,
                 DfD=DfD.out,gamma.othogonality=Gamma.otho,gamma.rnd.othogonality=Gamma.rnd.otho,score=score,time=cp.time)
      }
    }else
    {
      if(is.null(X2)==FALSE)
      {
        re<-list(gamma=Gamma.est,kappa=kappa.est,gamma.rnd=Gamma.rnd.est,beta=beta.est,beta0.rnd=beta0.rnd.est,beta0.sigma2=beta0.sigma2.est,beta2.rnd=beta2.rnd.est,beta2.Omega=beta2.Omega.est,
                 DfD=DfD.out,gamma.othogonality=Gamma.otho,gamma.rnd.othogonality=Gamma.rnd.otho,time=cp.time)
      }else
      {
        re<-list(gamma=Gamma.est,kappa=kappa.est,gamma.rnd=Gamma.rnd.est,beta=beta.est,beta0.rnd=beta0.rnd.est,beta0.sigma2=beta0.sigma2.est,
                 DfD=DfD.out,gamma.othogonality=Gamma.otho,gamma.rnd.othogonality=Gamma.rnd.otho,time=cp.time)
      }
    }
    
    return(re)
  }
  #-----------------------------
}
#################################################

#################################################
# Bootstrap inference on beta
lcap.beta.boot<-function(Y,X1=NULL,X2=NULL,gamma.rnd=NULL,boot=TRUE,boot.type=c("both","bycluster","withincluster"),sims=1000,boot.ci.type=c("boot.se","perc"),conf.level=0.95,boot.seed=100,verbose=TRUE,
                         Omega.diag=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE)
{
  # Y: a list of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # gamma.rnd: m by p matrix
  
  if(is.null(gamma.rnd))
  {
    stop("Error! Need gamma (random) values.")
  }else
  {
    if(boot)
    {
      #-----------------------------
      m<-length(Y)
      nvec<-sapply(Y,length)
      p<-ncol(Y[[1]][[1]])
      
      if(is.null(names(Y))==FALSE)
      {
        cluster.names<-names(Y)
      }else
      {
        cluster.names<-paste0("G",1:m)
      }
      
      if(is.null(X1)==FALSE)
      {
        q1<-ncol(X1[[1]])
        if(is.null(names(X1))==TRUE)
        {
          names(X1)<-cluster.names
        }
        for(i in 1:m)
        {
          if(is.null(colnames(X1[[i]]))==TRUE)
          {
            X1.names<-paste0("X1",1:q1)
            colnames(X1[[i]])<-X1.names
          }else
          {
            X1.names<-colnames(X1[[i]])
          }
        }
      }else
      {
        q1<-0
        X1.names<-NULL
      }
      if(is.null(X2)==FALSE)
      {
        q2<-ncol(X2[[1]])
        if(is.null(names(X2))==TRUE)
        {
          names(X2)<-cluster.names
        }
        for(i in 1:m)
        {
          if(is.null(colnames(X2[[i]]))==TRUE)
          {
            X2.names<-paste0("X2",1:q2)
            colnames(X2[[i]])<-X2.names
          }else
          {
            X2.names<-colnames(X2[[i]])
          }
        }
      }else
      {
        q2<-0
        X2.names<-NULL
      }
      q<-1+q1+q2
      X.names<-c("Intercept",X1.names,X2.names)
      #-----------------------------
      
      #-----------------------------
      beta.boot<-matrix(NA,q,sims)
      rownames(beta.boot)<-X.names
      beta0.rnd.boot<-matrix(NA,m,sims)
      rownames(beta0.rnd.boot)<-cluster.names
      beta0.sigma2.boot<-rep(NA,sims)
      if(is.null(X2)==FALSE)
      {
        beta2.rnd.boot<-array(NA,c(m,q2,sims))
        dimnames(beta2.rnd.boot)[[1]]<-cluster.names
        dimnames(beta2.rnd.boot)[[2]]<-X2.names
        beta2.Omega.boot<-array(NA,c(q2,q2,sims))
        dimnames(beta2.Omega.boot)[[1]]=dimnames(beta2.Omega.boot)[[2]]<-X2.names
      }
      for(b in 1:sims)
      {
        if(boot.type[1]=="both")
        {
          set.seed(boot.seed+b)
          idx.clst<-sample(1:m,m,replace=TRUE)
          idx.tmp<-vector("list",length=m)
          names(idx.tmp)<-cluster.names[idx.clst]
          for(i in 1:m)
          {
            idx.tmp[[i]]<-sample(1:nvec[idx.clst[i]],nvec[idx.clst[i]],replace=TRUE)
          }
          
          gamma.rnd.tmp<-gamma.rnd[idx.clst,]
          
          Ytmp<-vector("list",length=m)
          names(Ytmp)<-cluster.names[idx.clst]
          for(i in 1:m)
          {
            Ytmp[[i]]<-Y[[idx.clst[i]]][idx.tmp[[i]]]
          }
          if(is.null(X1)==FALSE)
          {
            X1.tmp<-vector("list",length=m)
            names(X1.tmp)<-cluster.names[idx.clst]
            for(i in 1:m)
            {
              X1.tmp[[i]]<-matrix(X1[[idx.clst[i]]][idx.tmp[[i]],],ncol=ncol(X1[[idx.clst[i]]]))
              colnames(X1.tmp[[i]])<-colnames(X1[[idx.clst[i]]])
              rownames(X1.tmp[[i]])<-rownames(X1[[idx.clst[i]]])[idx.tmp[[i]],]
            }
          }else
          {
            X1.tmp<-NULL
          }
          if(is.null(X2)==FALSE)
          {
            X2.tmp<-vector("list",length=m)
            names(X2.tmp)<-cluster.names[idx.clst]
            for(i in 1:m)
            {
              X2.tmp[[i]]<-matrix(X2[[idx.clst[i]]][idx.tmp[[i]],],ncol=ncol(X2[[idx.clst[i]]]))
              colnames(X2.tmp[[i]])<-colnames(X2[[idx.clst[i]]])
              rownames(X2.tmp[[i]])<-rownames(X2[[idx.clst[i]]])[idx.tmp[[i]],]
            }
          }else
          {
            X2.tmp<-NULL
          }
        }
        if(boot.type[1]=="bycluster")
        {
          set.seed(boot.seed+b)
          idx.clst<-sample(1:m,m,replace=TRUE)
          
          gamma.rnd.tmp<-gamma.rnd[idx.clst,]
          
          Ytmp<-Y[idx.clst]
          if(is.null(X1)==FALSE)
          {
            X1.tmp<-X1[idx.clst]
          }else
          {
            X1.tmp<-NULL
          }
          if(is.null(X2)==FALSE)
          {
            X2.tmp<-X2[idx.clst]
          }else
          {
            X2.tmp<-NULL
          }
        }
        if(boot.type[1]=="withincluster")
        {
          set.seed(boot.seed+b)
          idx.tmp<-vector("list",length=m)
          names(idx.tmp)<-cluster.names
          for(i in 1:m)
          {
            idx.tmp[[i]]<-sample(1:nvec[i],nvec[i],replace=TRUE)
          }
          
          gamma.rnd.tmp<-gamma.rnd
          
          Ytmp<-vector("list",length=m)
          names(Ytmp)<-cluster.names
          for(i in 1:m)
          {
            Ytmp[[i]]<-Y[[i]][idx.tmp[[i]]]
          }
          if(is.null(X1)==FALSE)
          {
            X1.tmp<-vector("list",length=m)
            names(X1.tmp)<-cluster.names
            for(i in 1:m)
            {
              X1.tmp[[i]]<-matrix(X1[[i]][idx.tmp[[i]],],ncol=ncol(X1[[i]]))
              colnames(X1.tmp[[i]])<-colnames(X1[[i]])
              rownames(X1.tmp[[i]])<-rownames(X1[[i]])[idx.tmp[[i]],]
            }
          }else
          {
            X1.tmp<-NULL
          }
          if(is.null(X2)==FALSE)
          {
            X2.tmp<-vector("list",length=m)
            names(X2.tmp)<-cluster.names
            for(i in 1:m)
            {
              X2.tmp[[i]]<-matrix(X2[[i]][idx.tmp[[i]],],ncol=ncol(X2[[i]]))
              colnames(X2.tmp[[i]])<-colnames(X2[[i]])
              rownames(X2.tmp[[i]])<-rownames(X2[[i]])[idx.tmp[[i]],]
            }
          }else
          {
            X2.tmp<-NULL
          }
        }
        
        re.tmp<-NULL
        try(re.tmp<-lcap.beta(Y=Ytmp,X1=X1.tmp,X2=X2.tmp,gamma.rnd=gamma.rnd.tmp,Omega.diag=TRUE,max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE))
        if(is.null(re.tmp)==FALSE)
        {
          if(re.tmp$convergence==TRUE)
          {
            beta.boot[,b]<-re.tmp$beta
            beta0.rnd.boot[,b]<-re.tmp$beta0.rnd
            beta0.sigma2.boot[b]<-re.tmp$beta0.sigma2
            if(is.null(X2)==FALSE)
            {
              beta2.rnd.boot[,,b]<-re.tmp$beta2.rnd
              beta2.Omega.boot[,,b]<-re.tmp$beta2.Omega
            }
          }
        }
        
        if(verbose)
        {
          print(paste0("Bootstrap sample ",b))
        }
      }
      
      for(j in 1:q)
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
      
      # summary
      beta.out<-cbind(Est=apply(beta.boot,1,mean,na.rm=TRUE),SE=apply(beta.boot,1,sd,na.rm=TRUE))
      beta.out<-cbind(beta.out,Stat=beta.out[,"Est"]/beta.out[,"SE"])
      beta.out<-cbind(beta.out,pvalue=(1-pnorm(abs(beta.out[,"Stat"])))*2)
      
      beta0.rnd.out<-cbind(Est=apply(beta0.rnd.boot,1,mean,na.rm=TRUE),SE=apply(beta0.rnd.boot,1,sd,na.rm=TRUE))
      beta0.rnd.out<-cbind(beta0.rnd.out,Stat=beta0.rnd.out[,"Est"]/beta0.rnd.out[,"SE"])
      beta0.rnd.out<-cbind(beta0.rnd.out,pvalue=(1-pnorm(abs(beta0.rnd.out[,"Stat"])))*2)
      
      beta0.sigma2.out<-cbind(Est=mean(beta0.sigma2.boot,na.rm=TRUE),SE=sd(beta0.sigma2.boot,na.rm=TRUE))
      beta0.sigma2.out<-cbind(beta0.sigma2.out,Stat=beta0.sigma2.out[,"Est"]/beta0.sigma2.out[,"SE"])
      beta0.sigma2.out<-cbind(beta0.sigma2.out,pvalue=(1-pnorm(abs(beta0.sigma2.out[,"Stat"])))*2)
      rownames(beta0.sigma2.out)<-"beta0.sigma2"
      
      if(is.null(X2)==FALSE)
      {
        beta2.rnd.out<-array(NA,c(q2,6,m))
        dimnames(beta2.rnd.out)[[1]]<-X2.names
        dimnames(beta2.rnd.out)[[2]]<-c("Est","SE","Stat","pvalue","LB","UB")
        dimnames(beta2.rnd.out)[[3]]<-cluster.names
        
        beta2.rnd.out[,1,]<-t(apply(beta2.rnd.boot,c(1,2),mean,na.rm=TRUE))
        beta2.rnd.out[,2,]<-t(apply(beta2.rnd.boot,c(1,2),sd,na.rm=TRUE))
        beta2.rnd.out[,3,]<-beta2.rnd.out[,1,]/beta2.rnd.out[,2,]
        beta2.rnd.out[,4,]<-(1-pnorm(abs(beta2.rnd.out[,3,])))*2
        
        beta2.Omega.out<-array(NA,c(q2,q2,6))
        dimnames(beta2.Omega.out)[[1]]=dimnames(beta2.Omega.out)[[2]]<-X2.names
        dimnames(beta2.Omega.out)[[3]]<-c("Est","SE","Stat","pvalue","LB","UB")
        beta2.Omega.out[,,1]<-apply(beta2.Omega.boot,c(1,2),mean,na.rm=TRUE)
        beta2.Omega.out[,,2]<-apply(beta2.Omega.boot,c(1,2),sd,na.rm=TRUE)
        beta2.Omega.out[,,3]<-beta2.Omega.out[,,1]/beta2.Omega.out[,,2]
        beta2.Omega.out[,,4]<-(1-pnorm(abs(beta2.Omega.out[,,3])))*2
      }
      
      if(boot.ci.type=="boot.se")
      {
        zv<-qnorm(1-(1-conf.level)/2)
        beta.out<-cbind(beta.out,LB=beta.out[,"Est"]-zv*beta.out[,"SE"],UB=beta.out[,"Est"]+zv*beta.out[,"SE"])
        beta0.rnd.out<-cbind(beta0.rnd.out,LB=beta0.rnd.out[,"Est"]-zv*beta0.rnd.out[,"SE"],UB=beta0.rnd.out[,"Est"]+zv*beta0.rnd.out[,"SE"])
        beta0.sigma2.out<-cbind(beta0.sigma2.out,LB=beta0.sigma2.out[,"Est"]-zv*beta0.sigma2.out[,"SE"],UB=beta0.sigma2.out[,"Est"]+zv*beta0.sigma2.out[,"SE"])
        
        if(is.null(X2)==FALSE)
        {
          beta2.rnd.out[,5,]<-beta2.rnd.out[,1,]-zv*beta2.rnd.out[,2,]
          beta2.rnd.out[,6,]<-beta2.rnd.out[,1,]+zv*beta2.rnd.out[,2,]
          
          beta2.Omega.out[,,5]<-beta2.Omega.out[,,1]-zv*beta2.Omega.out[,,2]
          beta2.Omega.out[,,6]<-beta2.Omega.out[,,1]+zv*beta2.Omega.out[,,2]
        }
      }
      if(boot.ci.type=="perc")
      {
        pbs<-c((1-conf.level)/2,1-(1-conf.level)/2)
        beta.out<-cbind(beta.out,LB=apply(beta.boot,1,quantile,probs=pbs[1],na.rm=TRUE),UB=apply(beta.boot,1,quantile,probs=pbs[2],na.rm=TRUE))
        beta0.rnd.out<-cbind(beta0.rnd.out,LB=apply(beta0.rnd.boot,1,quantile,probs=pbs[1],na.rm=TRUE),UB=apply(beta0.rnd.boot,1,quantile,probs=pbs[2],na.rm=TRUE))
        beta0.sigma2.out<-cbind(beta0.sigma2.out,LB=quantile(beta0.sigma2.boot,probs=pbs[1],na.rm=TRUE),UB=quantile(beta0.sigma2.boot,probs=pbs[2],na.rm=TRUE))
        
        if(is.null(X2)==FALSE)
        {
          beta2.rnd.out[,5,]<-t(apply(beta2.rnd.boot,c(1,2),quantile,probs=pbs[1],na.rm=TRUE))
          beta2.rnd.out[,6,]<-t(apply(beta2.rnd.boot,c(1,2),quantile,probs=pbs[2],na.rm=TRUE))
          
          beta2.Omega.out[,,5]<-t(apply(beta2.Omega.boot,c(1,2),quantile,probs=pbs[1],na.rm=TRUE))
          beta2.Omega.out[,,6]<-t(apply(beta2.Omega.boot,c(1,2),quantile,probs=pbs[2],na.rm=TRUE))
        }
      }
      #-----------------------------
      
      #-----------------------------
      if(is.null(X2)==FALSE)
      {
        re.inf<-list(beta=beta.out,beta0.rnd=beta0.rnd.out,beta0.sigma2=beta0.sigma2.out,beta2.rnd=beta2.rnd.out,beta2.Omega=beta2.Omega.out)
        re.boot<-list(beta=beta.boot,beta0.rnd=beta0.rnd.boot,beta0.sigma2=beta0.sigma2.boot,beta2.rnd=beta2.rnd.boot,beta2.Omega=beta2.Omega.boot)
      }else
      {
        re.inf<-list(beta=beta.out,beta0.rnd=beta0.rnd.out,beta0.sigma2=beta0.sigma2.out)
        re.boot<-list(beta=beta.boot,beta0.rnd=beta0.rnd.boot,beta0.sigma2=beta0.sigma2.boot)
      }
      
      re<-list(Inference=re.inf,boot=re.boot)
      
      return(re)
      #-----------------------------
    }
  }
}
#################################################

#################################################
# Asymptotic inference on beta (Theorem 1, Section 2.4)
#
# Asymptotic variances (known gamma):
#   (i)   sqrt(M_n)(beta1.hat - beta1*) -> N(0, 2*Q^{-1})
#         where Q = lim m^{-1} sum_i n_i^{-1} sum_j x_{1i(j)} x_{1i(j)}^T
#   (ii)  sqrt(m)(beta0.hat - beta0*)   -> N(0, sigma^{2*})
#         sqrt(m)(sigma2.hat - sigma^{2*}) -> N(0, 2*sigma^{4*})
#   (iii) sqrt(m)(beta2.hat - beta2*)   -> N(0, Omega*)
#         sqrt(m)*vec(Omega.hat - Omega*) -> N(0, 2*Omega* %x% Omega*)
#
lcap.beta.asym<-function(Y,X1=NULL,X2=NULL,gamma.rnd=NULL,conf.level=0.95,
                         Omega.diag=TRUE,max.itr=1000,tol=1e-4)
{
  # Y: a list of Y
  # X1: a list of covariates corresponding to fixed-effect beta
  # X2: a list of covariates corresponding to random-effect beta
  # gamma.rnd: m by p matrix (estimated cluster-specific gamma, fixed)

  if(is.null(gamma.rnd))
  {
    stop("Error! Need gamma (random) values.")
  }

  #-----------------------------
  m<-length(Y)
  nvec<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])

  if(is.null(names(Y))==FALSE)
  {
    cluster.names<-names(Y)
  }else
  {
    cluster.names<-paste0("G",1:m)
  }

  if(is.null(X1)==FALSE)
  {
    q1<-ncol(X1[[1]])
    if(is.null(names(X1))==TRUE)
    {
      names(X1)<-cluster.names
    }
    for(i in 1:m)
    {
      if(is.null(colnames(X1[[i]]))==TRUE)
      {
        X1.names<-paste0("X1",1:q1)
        colnames(X1[[i]])<-X1.names
      }else
      {
        X1.names<-colnames(X1[[i]])
      }
    }
  }else
  {
    q1<-0
    X1.names<-NULL
  }
  if(is.null(X2)==FALSE)
  {
    q2<-ncol(X2[[1]])
    if(is.null(names(X2))==TRUE)
    {
      names(X2)<-cluster.names
    }
    for(i in 1:m)
    {
      if(is.null(colnames(X2[[i]]))==TRUE)
      {
        X2.names<-paste0("X2",1:q2)
        colnames(X2[[i]])<-X2.names
      }else
      {
        X2.names<-colnames(X2[[i]])
      }
    }
  }else
  {
    q2<-0
    X2.names<-NULL
  }
  q<-1+q1+q2
  X.names<-c("Intercept",X1.names,X2.names)
  #-----------------------------

  #-----------------------------
  # Step 1: estimate regression parameters given fixed gamma.rnd
  beta.fit<-lcap.beta(Y=Y,X1=X1,X2=X2,gamma.rnd=gamma.rnd,Omega.diag=Omega.diag,max.itr=max.itr,tol=tol,score.return=FALSE,trace=FALSE)

  beta.est<-beta.fit$beta         # named vector: (Intercept, X1.names, X2.names)
  beta0.rnd<-beta.fit$beta0.rnd   # length m
  sigma2.est<-beta.fit$beta0.sigma2
  beta0.est<-beta.est[1]
  if(q1>0) beta1.est<-beta.est[2:(q1+1)]
  if(q2>0)
  {
    beta2.rnd<-beta.fit$beta2.rnd    # m by q2
    Omega.est<-beta.fit$beta2.Omega  # q2 by q2
    beta2.est<-beta.est[(q1+2):q]
  }
  #-----------------------------

  #-----------------------------
  # Step 2: compute M_n (total observations)
  Mn<-0
  for(i in 1:m)
  {
    for(j in 1:nvec[i])
    {
      Mn<-Mn+nrow(Y[[i]][[j]])
    }
  }
  #-----------------------------

  #-----------------------------
  # Step 3: asymptotic SE for beta0 and sigma2
  #   Var(beta0.hat) = sigma2 / m
  #   Var(sigma2.hat) = 2 * sigma2^2 / m
  se.beta0<-sqrt(sigma2.est/m)
  se.sigma2<-sqrt(2*sigma2.est^2/m)
  #-----------------------------

  #-----------------------------
  # Step 4: asymptotic SE for beta1 (if exists)
  #   Var(beta1.hat) = 2 * Q^{-1} / M_n
  #   Q = m^{-1} sum_i n_i^{-1} sum_j x_{1i(j)} x_{1i(j)}^T
  if(q1>0)
  {
    Q.mat<-matrix(0,q1,q1)
    for(i in 1:m)
    {
      Qtmp<-matrix(0,q1,q1)
      for(j in 1:nvec[i])
      {
        Qtmp<-Qtmp+X1[[i]][j,]%*%t(X1[[i]][j,])
      }
      Q.mat<-Q.mat+Qtmp/nvec[i]
    }
    Q.mat<-Q.mat/m

    var.beta1<-2*ginv(Q.mat)/Mn
    se.beta1<-sqrt(diag(var.beta1))
  }
  #-----------------------------

  #-----------------------------
  # Step 5: asymptotic SE for beta2 and Omega (if exists)
  #   Var(beta2.hat) = Omega / m
  #   Var(vec(Omega.hat)) = 2 * (Omega %x% Omega) / m
  if(q2>0)
  {
    var.beta2<-Omega.est/m
    se.beta2<-sqrt(diag(as.matrix(var.beta2)))

    var.vecOmega<-2*kronecker(Omega.est,Omega.est)/m
    se.Omega<-matrix(sqrt(diag(var.vecOmega)),q2,q2)
    rownames(se.Omega)=colnames(se.Omega)<-X2.names
  }
  #-----------------------------

  #-----------------------------
  # Step 6: assemble output tables (matching bootstrap output format)
  zv<-qnorm(1-(1-conf.level)/2)

  # beta = (Intercept, beta1, beta2) combined table
  beta.se<-rep(NA,q)
  beta.se[1]<-se.beta0
  if(q1>0) beta.se[2:(q1+1)]<-se.beta1
  if(q2>0) beta.se[(q1+2):q]<-se.beta2

  beta.out<-cbind(Est=beta.est,SE=beta.se)
  beta.out<-cbind(beta.out,Stat=beta.out[,"Est"]/beta.out[,"SE"])
  beta.out<-cbind(beta.out,pvalue=(1-pnorm(abs(beta.out[,"Stat"])))*2)
  beta.out<-cbind(beta.out,LB=beta.out[,"Est"]-zv*beta.out[,"SE"],UB=beta.out[,"Est"]+zv*beta.out[,"SE"])
  rownames(beta.out)<-X.names

  # sigma2
  beta0.sigma2.out<-cbind(Est=sigma2.est,SE=se.sigma2)
  beta0.sigma2.out<-cbind(beta0.sigma2.out,Stat=beta0.sigma2.out[,"Est"]/beta0.sigma2.out[,"SE"])
  beta0.sigma2.out<-cbind(beta0.sigma2.out,pvalue=(1-pnorm(abs(beta0.sigma2.out[,"Stat"])))*2)
  beta0.sigma2.out<-cbind(beta0.sigma2.out,LB=beta0.sigma2.out[,"Est"]-zv*beta0.sigma2.out[,"SE"],UB=beta0.sigma2.out[,"Est"]+zv*beta0.sigma2.out[,"SE"])
  rownames(beta0.sigma2.out)<-"beta0.sigma2"

  if(q2>0)
  {
    # Omega
    beta2.Omega.out<-array(NA,c(q2,q2,6))
    dimnames(beta2.Omega.out)[[1]]=dimnames(beta2.Omega.out)[[2]]<-X2.names
    dimnames(beta2.Omega.out)[[3]]<-c("Est","SE","Stat","pvalue","LB","UB")
    beta2.Omega.out[,,1]<-Omega.est
    beta2.Omega.out[,,2]<-se.Omega
    beta2.Omega.out[,,3]<-beta2.Omega.out[,,1]/beta2.Omega.out[,,2]
    beta2.Omega.out[,,4]<-(1-pnorm(abs(beta2.Omega.out[,,3])))*2
    beta2.Omega.out[,,5]<-beta2.Omega.out[,,1]-zv*beta2.Omega.out[,,2]
    beta2.Omega.out[,,6]<-beta2.Omega.out[,,1]+zv*beta2.Omega.out[,,2]

    re<-list(beta=beta.out,beta0.sigma2=beta0.sigma2.out,beta2.Omega=beta2.Omega.out)
  }else
  {
    re<-list(beta=beta.out,beta0.sigma2=beta0.sigma2.out)
  }

  return(re)
  #-----------------------------
}
#################################################




