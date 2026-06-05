#################################################
# CAP mediation
# add random error to mediator model
# hierarchical likelihood based estimator
#
# V4: RcppArmadillo-accelerated update of V3. Same interface and outputs; the hot
# loops (projected scores, the theta-update accumulation, eigen.solve, and the
# per-subject sample covariance) call a compiled kernel. The mixed-effects (lme)
# fit stays in R. See CAP-mediation/README.md.
#################################################



#################################################
# Load the compiled kernel (med_kernels.cpp -> shared object), building it once
# if necessary. Set option "med.v4.dir" or variable MED_V4_DIR to point at the
# folder containing med_kernels.cpp; defaults to "CAP-mediation/V4".


# helper: per-subject sample covariance cube (== cov(M[[i]]) for all i)
.med_cov <- function(M) .Call("med_cov_cpp", lapply(M, as.matrix))

#################################################
# objective function
obj.func<-function(X,Sigma,Y,nT,theta,alpha0,alpha0.rnd,alpha,beta,gamma0,gamma,tau2,sigma2)
{
  # X: (X, W) matrix
  # Sigma: covariance estimate of M
  # Y: n by 1 outcome vector
  # nT: n by 1 vector, # of observation of each subject
  # theta: projection vector
  # alpha0: fixed intercept of M model
  # alpha0.rnd: n by 1 vector, random effect of alpha0
  # alpha: coefficient of M model
  # beta: coefficient of M in Y model
  # gamma0: intercept of Y model
  # gamma: coefficient of Y model
  # tau2: M model error standard deviation
  # sigma2: Y model error standard deviation
  
  n<-length(nT)

  # t(theta)%*%Sigma%*%theta  (kernel)
  score<-as.numeric(.Call("med_score_cpp",Sigma,as.numeric(theta)))

  ll1<-sum(((alpha0.rnd+X%*%alpha)+score*exp(-alpha0.rnd-X%*%alpha))*nT)/2
  ll2<-sum((Y-gamma0-X%*%gamma-beta*log(score))^2/sigma2+log(sigma2))/2
  ll3<-sum((alpha0.rnd-alpha0)^2/tau2+log(tau2))/2
  
  ll<-ll1+ll2+ll3
  return(ll)
}
#################################################

#################################################
# given theta, estimate model coefficient parameters

# given an estimate of the covariance matrices
CAPMediation_coef_Mcov<-function(X,M.cov,Y,theta)
{
  # X: (X, W) matrix
  # M.cov: p by p by n array, covariance matrix of M
  # Y: n by 1 vector
  # theta: projection vector
  
  n<-length(Y)
  
  q<-ncol(X)-1
  if(is.null(colnames(X)))
  {
    if(q>0)
    {
      colnames(X)<-c("X",paste0("W",1:q)) 
    }else
    {
      colnames(X)<-c("X")
    }
  }
  
  score<-as.numeric(.Call("med_score_cpp",M.cov,as.numeric(theta)))
  
  # beta and gamma estimate
  Z<-cbind(rep(1,n),X,log(score))
  colnames(Z)<-c("Intercept",colnames(X),"M")
  mu.est<-c(ginv(t(Z)%*%Z)%*%t(Z)%*%Y)
  
  gamma0.est<-mu.est[1]
  gamma.est<-mu.est[2:(length(mu.est)-1)]
  beta.est<-mu.est[length(mu.est)]
  
  # sigma2 estimate
  sigma2.est<-mean((Y-Z%*%mu.est)^2,na.rm=TRUE)
  
  # estimate alpha0, alpha0.rnd, alpha using mixed effects model
  dtmp<-data.frame(ID=1:n,score=log(score),X)
  eval(parse(text=paste0("fit.tmp<-lme(score~",paste(colnames(X),collapse="+"),",random=~1|ID,data=dtmp,control=lmeControl(opt='optim'))")))
  alpha0.est<-fit.tmp$coefficients$fixed[1]
  alpha.est<-fit.tmp$coefficients$fixed[-1]
  alpha0.rnd.est<-c(fit.tmp$coefficients$random$ID+fit.tmp$coefficients$fixed[1])
  tau2.est<-mean((alpha0.rnd.est-alpha0.est)^2,na.rm=TRUE)

  IE.est<-alpha.est[1]*beta.est
  
  names(alpha0.est)<-NULL
  names(alpha.est)=names(gamma.est)<-colnames(X)
  re<-list(theta=theta,alpha=alpha.est,beta=beta.est,gamma=gamma.est,IE=IE.est,alpha0=alpha0.est,alpha0.rnd=alpha0.rnd.est,gamma0=gamma0.est,tau2=tau2.est,sigma2=sigma2.est)
  
  return(re)
}
# sample covariance of M
CAPMediation_coef<-function(X,M,Y,theta)
{
  # X: (1, X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # theta: projection vector
  
  M.cov<-.med_cov(M)
  
  re<-CAPMediation_coef_Mcov(X,M.cov,Y,theta)
  
  return(re)
}
#################################################

#################################################
# eigenvectors and eigenvalues of A with respect to H
# H positive definite and symmetric
eigen.solve<-function(A,H)
{
  # smallest generalized eigenvector of A w.r.t. H  (RcppArmadillo kernel)
  as.numeric(.Call("med_eigen_solve_cpp", as.matrix(A), as.matrix(H)))
}
#################################################

#################################################
# given an estimate of the covariance matrices
CAPMediation_D1_Mcov<-function(X,M.cov,Y,nT,H=NULL,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0=NULL)
{
  # X: (X, W) matrix
  # M.cov: p by p by n array, covariance matrix of M
  # Y: n by 1 vector
  # nT: n by 1 vector, # of observation of each subject
  # H: a p by p positive definite matrix for theta constraint
  
  n<-length(Y)
  
  q<-ncol(X)-1
  if(is.null(colnames(X)))
  {
    if(q>0)
    {
      colnames(X)<-c("X",paste0("W",1:q)) 
    }else
    {
      colnames(X)<-c("X")
    }
  }
  
  p<-dim(M.cov)[1]
  
  if(is.null(H))
  {
    H<-apply(M.cov,c(1,2),mean,na.rm=TRUE)
  }
  
  if(is.null(theta0))
  {
    theta0<-rep(1/sqrt(p),p)
  }
  
  # initial values
  otmp0<-CAPMediation_coef_Mcov(X,M.cov,Y,theta0)
  alpha.ini<-otmp0$alpha
  beta.ini<-otmp0$beta
  gamma.ini<-otmp0$gamma
  alpha0.ini<-otmp0$alpha0
  alpha0.rnd.ini<-otmp0$alpha0.rnd
  gamma0.ini<-otmp0$gamma0
  tau2.ini<-otmp0$tau2
  sigma2.ini<-otmp0$sigma2
  
  if(trace)
  {
    alpha.trace<-NULL
    beta.trace<-NULL
    gamma.trace<-NULL
    alpha0.trace<-NULL
    alpha0.rnd.trace<-NULL
    gamma0.trace<-NULL
    theta.trace<-NULL
    tau2.trace<-NULL
    sigma2.trace<-NULL
    
    obj<-NULL
  }
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    obj0<-obj.func(X,M.cov,Y,nT,theta0,alpha0.ini,alpha0.rnd.ini,alpha.ini,beta.ini,gamma0.ini,gamma.ini,tau2.ini,sigma2.ini)
    score<-as.numeric(.Call("med_score_cpp",M.cov,as.numeric(theta0)))
    
    # update theta
    U<-c(exp(-alpha0.rnd.ini-X%*%alpha.ini))
    V<-c(Y-gamma0.ini-X%*%gamma.ini)
    wA<-nT*U-(2*beta.ini*(V-beta.ini*log(score)))/(sigma2.ini*score)
    A<-.Call("med_accum_cpp",M.cov,as.numeric(wA))
    theta.new<-eigen.solve(A,H)
    
    otmp.new<-CAPMediation_coef_Mcov(X,M.cov,Y,theta.new)
    alpha.new<-otmp.new$alpha
    beta.new<-otmp.new$beta
    gamma.new<-otmp.new$gamma
    alpha0.new<-otmp.new$alpha0
    alpha0.rnd.new<-otmp.new$alpha0.rnd
    gamma0.new<-otmp.new$gamma0
    tau2.new<-otmp.new$tau2
    sigma2.new<-otmp.new$sigma2
    
    obj.new<-obj.func(X,M.cov,Y,nT,theta.new,alpha0.new,alpha0.rnd.new,alpha.new,beta.new,gamma0.new,gamma.new,tau2.new,sigma2.new)
    
    diff<-max(abs(c(alpha.new-alpha.ini,beta.new-beta.ini,gamma.new-gamma.ini))) 
    
    alpha.ini<-alpha.new
    beta.ini<-beta.new
    gamma.ini<-gamma.new
    alpha0.ini<-alpha0.new
    alpha0.rnd.ini<-alpha0.rnd.new
    gamma0.ini<-gamma0.new
    tau2.ini<-tau2.new
    sigma2.ini<-sigma2.new
    theta0<-theta.new
    
    if(trace)
    {
      alpha.trace<-cbind(alpha.trace,alpha.ini)
      beta.trace<-c(beta.trace,beta.ini)
      gamma.trace<-cbind(gamma.trace,gamma.ini)
      alpha0.trace<-cbind(alpha0.trace,alpha0.ini)
      alpha0.rnd.trace<-cbind(alpha0.rnd.trace,alpha0.rnd.ini)
      gamma0.trace<-cbind(gamma0.trace,gamma0.ini)
      tau2.trace<-c(tau2.trace,tau2.ini)
      sigma2.trace<-c(sigma2.trace,sigma2.ini)
      theta.trace<-cbind(theta.trace,theta0)
      
      obj<-c(obj,obj.new)
    }
    
    # print(c(diff,obj.new))
  }
  
  theta.est<-theta0/sqrt(sum(theta0^2))
  if(theta.est[which.max(abs(theta.est))]<0)
  {
    theta.est<--theta.est
  }
  otmp.est<-CAPMediation_coef_Mcov(X,M.cov,Y,theta.est)
  score<-as.numeric(.Call("med_score_cpp",M.cov,as.numeric(theta.est)))
  obj.est<-obj.func(X,M.cov,Y,nT,theta.est,alpha0=otmp.est$alpha0,alpha0.rnd=otmp.est$alpha0.rnd,alpha=otmp.est$alpha,beta=otmp.est$beta,gamma0=otmp.est$gamma0,gamma=otmp.est$gamma,
                    tau2=otmp.est$tau2,sigma2=otmp.est$sigma2)
  
  re<-otmp.est
  if(score.return)
  {
    re$score<-score
  }
  re$obj<-obj.est
  
  if(trace)
  {
    rownames(alpha.trace)=rownames(gamma.trace)<-colnames(X)
    
    re$theta.trace<-theta.trace
    re$alpha.trace<-alpha.trace
    re$beta.trace<-beta.trace
    re$gamma.trace<-gamma.trace
    re$alpha0.trace<-alpha0.trace
    re$alpha0.rnd.trace<-alpha0.rnd.trace
    re$gamma0.trace<-gamma0.trace
    re$tau2.trace<-tau2.trace
    re$sigma2.trace<-sigma2.trace
    re$obj.trace<-obj
  }
  
  return(re)
}
CAPMediation_D1<-function(X,M,Y,H=NULL,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0=NULL)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # H: a p by p positive definite matrix for theta constraint
  
  n<-length(Y)
  p<-ncol(M[[1]])
  
  nT<-sapply(M,nrow)
  M.cov<-.med_cov(M)
  
  re<-CAPMediation_D1_Mcov(X=X,M.cov=M.cov,Y=Y,nT=nT,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0=theta0)
  
  return(re)
}

# multiple initial values
CAPMediation_D1_opt_Mcov<-function(X,M.cov,Y,nT,H=NULL,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0.mat=NULL,ninitial=NULL,seed=100)
{
  # X: (X, W) matrix
  # M.cov: p by p by n array, covariance matrix of M
  # Y: n by 1 vector
  # nT: n by 1 vector, # of observation of each subject
  # H: a p by p positive definite matrix for theta constraint
  
  n<-nrow(X)
  
  q<-ncol(X)-1
  if(is.null(colnames(X)))
  {
    if(q>0)
    {
      colnames(X)<-c("X",paste0("W",1:q)) 
    }else
    {
      colnames(X)<-c("X")
    }
  }
  
  p<-dim(M.cov)[1]
  
  #----------------------------------
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
    if(is.null(ninitial))
    {
      ninitial<-min(p,10)
    }
  
  # theta initial
  if(is.null(theta0.mat))
  {
    set.seed(seed)
    theta.tmp<-matrix(rnorm((max(p,ninitial)+1+5)*p,mean=0,sd=1),nrow=p)
    theta0.mat<-apply(theta.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
  }
  set.seed(seed)
  theta0.mat<-matrix(theta0.mat[,sort(sample(1:ncol(theta0.mat),ninitial,replace=FALSE))],ncol=ninitial)
  #----------------------------------
  
  #----------------------------------
  # try different initial values with the lowest objective function
  re.tmp<-vector("list",ninitial)
  obj<-rep(NA,ninitial)
  for(kk in 1:ninitial)
  {
    try(re.tmp[[kk]]<-CAPMediation_D1_Mcov(X,M.cov,Y,nT=nT,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0=theta0.mat[,kk]))
    
    if(is.null(re.tmp[[kk]])==FALSE)
    {
      theta.unscale<-re.tmp[[kk]]$theta/sqrt(t(re.tmp[[kk]]$theta)%*%H%*%re.tmp[[kk]]$theta)[1,1]
      
      try(coef.tmp<-CAPMediation_coef_Mcov(X,M.cov,Y,theta=theta.unscale))
      try(obj[kk]<-obj.func(X,M.cov,Y,nT=nT,theta=theta.unscale,alpha0=coef.tmp$alpha0,alpha0.rnd=coef.tmp$alpha0.rnd,alpha=coef.tmp$alpha,beta=coef.tmp$beta,gamma0=coef.tmp$gamma0,gamma=coef.tmp$gamma,
                            tau2=coef.tmp$tau2,sigma2=coef.tmp$sigma2))
    }
  }
  opt.idx<-which.min(obj)
  re<-re.tmp[[opt.idx]]
  #----------------------------------
  
  return(re)
}
CAPMediation_D1_opt<-function(X,M,Y,H=NULL,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0.mat=NULL,ninitial=NULL,seed=100)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # H: a p by p positive definite matrix for theta constraint
  
  n<-length(Y)
  p<-ncol(M[[1]])
  
  nT<-sapply(M,nrow)
  M.cov<-.med_cov(M)
  
  re<-CAPMediation_D1_opt_Mcov(X=X,M.cov=M.cov,Y=Y,nT=nT,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed)
  
  return(re)
}
#################################################

#################################################
# higher order components
# option: 
# 1. remove identified components from M direct
# 2. remove identified components from M ~ X + W residuals

CAPMediation_Dk<-function(X,M,Y,H=NULL,Theta0=NULL,Y.remove=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0.mat=NULL,ninitial=NULL,seed=100)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # H: a p by p positive definite matrix for theta constraint
  # Theta0: p by k matrix, identified components
  
  n<-length(Y)
  p<-ncol(M[[1]])
  nT<-sapply(M,nrow)
  
  if(is.null(Theta0))
  {
    re<-CAPMediation_D1_opt(X=X,M=M,Y=Y,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed)
  }else
  {
    p0<-ncol(Theta0)
    
    M.cov<-.med_cov(M)
    score<-matrix(NA,n,p0)
    colnames(score)<-paste0("C",1:p0)
    for(kk in 1:p0)
    {
      score[,kk]<-as.numeric(.Call("med_score_cpp",M.cov,as.numeric(Theta0[,kk])))
    }
    
    # estimate alpha0
    alpha0.est<-rep(NA,p0)
    alpha.est<-matrix(NA,ncol(X),p0)
    rownames(alpha.est)<-colnames(X)
    colnames(alpha.est)<-paste0("C",1:p0)
    beta.est<-rep(NA,p0)
    alpha0.rnd.est<-matrix(NA,n,p0)
    colnames(alpha0.rnd.est)<-paste0("C",1:p0)
    for(kk in 1:p0)
    {
      otmp<-CAPMediation_coef(X,M,Y,theta=Theta0[,kk])
      alpha0.est[kk]<-otmp$alpha0
      alpha0.rnd.est[,kk]<-otmp$alpha0.rnd
      alpha.est[,kk]<-otmp$alpha
      beta.est[kk]<-otmp$beta
    }
    
    Mnew<-vector("list",length=n)
    names(Mnew)<-names(M)
    for(i in 1:n)
    {
      Mtmp<-M[[i]]-M[[i]]%*%Theta0%*%t(Theta0)
      Mtmp.svd<-svd(Mtmp)
      # dnew<-c(Mtmp.svd$d[1:(p-p0)],sqrt(exp(alpha0.rnd.est[i,])*nT[i]))
      dnew<-c(Mtmp.svd$d[1:(p-p0)],sqrt(exp(alpha0.rnd.est[i,]-alpha0.est)*nT[i]))
      Mnew[[i]]<-Mtmp.svd$u%*%diag(dnew)%*%t(Mtmp.svd$v)
    }
    
    if(Y.remove)
    {
      Ynew<-Y-log(score)%*%beta.est 
    }else
    {
      Ynew<-Y
    }
    
    re<-CAPMediation_D1_opt(X=X,M=Mnew,Y=Ynew,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed)
    re$orthogonal<-c(t(re$theta)%*%Theta0)
  }
  
  return(re)
}
#################################################

#################################################
# level of diagonalization
diag.level<-function(M,Theta)
{
  # M: a list of length n, M
  # Theta: p by k matrix, identified components
  
  n<-length(M)
  
  if(ncol(Theta)==1)
  {
    re<-list(avg.level=1,sub.level=rep(1,n))
  }else
  {
    p<-ncol(M[[1]])
    nT<-sapply(M,nrow)
    ps<-ncol(Theta)
    
    dl.sub<-matrix(NA,n,ps)
    colnames(dl.sub)<-paste0("C",1:ps)
    dl.sub[,1]<-1
    for(i in 1:n)
    {
      cov.tmp<-cov(M[[i]])
      
      for(j in 2:ps)
      {
        theta.tmp<-Theta[,1:j]
        mat.tmp<-t(theta.tmp)%*%cov.tmp%*%theta.tmp
        dl.sub[i,j]<-det(diag(diag(mat.tmp)))/det(mat.tmp)
      }
    }
    
    pmean<-apply(dl.sub,2,function(y){return(prod(apply(cbind(y,nT),1,function(x){return(x[1]^(x[2]/sum(nT)))})))})
    
    re<-list(avg.level=pmean,sub.level=dl.sub)
  }
  
  return(re)
}
#################################################

#################################################
# CAP mediation function: finding k components
CAPMediation<-function(X,M,Y,H=NULL,stop.crt=c("nD","DfD"),nD=NULL,DfD.thred=2,Y.remove=TRUE,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,theta0.mat=NULL,ninitial=NULL,seed=100,verbose=TRUE)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # H: a p by p positive definite matrix for theta constraint
  
  if(stop.crt[1]=="nD"&is.null(nD))
  {
    stop.crt<-"DfD"
  }
  
  n<-length(Y)
  p<-ncol(M[[1]])
  
  q<-ncol(X)-1
  if(is.null(colnames(X)))
  {
    if(q>0)
    {
      X.names<-c("X",paste0("W",1:q)) 
    }else
    {
      X.names<-c("X")
    }
    colnames(X)<-X.names
  }else
  {
    X.names<-colnames(X)
  }
  
  M.cov<-.med_cov(M)
  
  if(is.null(H))
  {
    H<-apply(M.cov,c(1,2),mean,na.rm=TRUE)
  }
  
  #--------------------------------------------
  # First direction
  tm1<-system.time(re1<-CAPMediation_D1_opt(X=X,M=M,Y=Y,H=H,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,theta0.mat=theta0.mat,ninitial=ninitial,seed=seed))
  
  Theta.est<-matrix(re1$theta,ncol=1)
  coef.est<-matrix(c(re1$alpha[1],re1$beta,re1$gamma[1],re1$IE,re1$gamma[1]),ncol=1)
  rownames(coef.est)<-c("alpha","beta","gamma","IE","DE")
  if(q>0)
  {
    ancoef.est<-matrix(c(re1$alpha0,re1$alpha[-1],re1$gamma0,re1$gamma[-1]),ncol=1)
    rownames(ancoef.est)<-paste0(rep(c("M","Y"),each=q+1),"-",rep(c("Intercept",X.names[-1]),2))
  }else
  {
    ancoef.est<-matrix(c(re1$alpha0,re1$gamma0),ncol=1)
    rownames(ancoef.est)<-paste0(c("M","Y"),"-","Intercept")
  }
  cp.time<-matrix(as.numeric(tm1[1:3]),ncol=1)
  rownames(cp.time)<-c("user","system","elapsed")
  
  if(score.return)
  {
    score<-matrix(re1$score,ncol=1)
  }
  
  if(verbose)
  {
    print(paste0("Component ",ncol(Theta.est)))
  }
  #--------------------------------------------
  
  if(stop.crt[1]=="nD")
  {
    if(nD>1)
    {
      for(j in 2:nD)
      {
        re.tmp<-NULL
        try(tm.tmp<-system.time(re.tmp<-CAPMediation_Dk(X=X,M=M,Y=Y,H=H,Theta0=Theta.est,Y.remove=Y.remove,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,
                                                        theta0.mat=theta0.mat,ninitial=ninitial,seed=seed)))
        
        if(is.null(re.tmp)==FALSE)
        {
          Theta.est<-cbind(Theta.est,re.tmp$theta)
          coef.est<-cbind(coef.est,c(re.tmp$alpha[1],re.tmp$beta,re.tmp$gamma[1],re.tmp$IE,re.tmp$gamma[1]))
          if(q>0)
          {
            ancoef.est<-cbind(ancoef.est,c(re.tmp$alpha0,re.tmp$alpha[-1],re.tmp$gamma0,re.tmp$gamma[-1]))
          }else
          {
            ancoef.est<-cbind(ancoef.est,c(re.tmp$alpha0,re.tmp$gamma0))
          }
          cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
          
          if(score.return)
          {
            score<-cbind(score,re.tmp$score)
          }
          
          if(verbose)
          {
            print(paste0("Component ",ncol(Theta.est)))
          }
        }else
        {
          break
        }
      }
    }
    
    colnames(Theta.est)=colnames(coef.est)=colnames(ancoef.est)=colnames(cp.time)<-paste0("C",1:ncol(Theta.est))
    
    cp.time<-cbind(cp.time,apply(cp.time,1,sum))
    colnames(cp.time)[ncol(cp.time)]<-"Total"
    
    if(score.return)
    {
      colnames(score)<-paste0("C",1:ncol(Theta.est))
    }
    
    DfD.out<-diag.level(M,Theta.est)
  }
  if(stop.crt[1]=="DfD")
  {
    nD<-1
    
    DfD.tmp<-1
    while(DfD.tmp<DfD.thred)
    {
      re.tmp<-NULL
      try(tm.tmp<-system.time(re.tmp<-CAPMediation_Dk(X=X,M=M,Y=Y,H=H,Theta0=Theta.est,Y.remove=Y.remove,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,
                                                      theta0.mat=theta0.mat,ninitial=ninitial,seed=seed)))
      
      if(is.null(re.tmp)==FALSE)
      {
        nD<-nD+1
        
        DfD.out<-diag.level(M,cbind(Theta.est,re.tmp$theta))
        DfD.tmp<-DfD.out$avg.level[nD]
        
        if(DfD.tmp<DfD.thred)
        {
          Theta.est<-cbind(Theta.est,re.tmp$theta)
          coef.est<-cbind(coef.est,c(re.tmp$alpha[1],re.tmp$beta,re.tmp$gamma[1],re.tmp$IE,re.tmp$gamma[1]))
          if(q>0)
          {
            ancoef.est<-cbind(ancoef.est,c(re.tmp$alpha0,re.tmp$alpha[-1],re.tmp$gamma0,re.tmp$gamma[-1]))
          }else
          {
            ancoef.est<-cbind(ancoef.est,c(re.tmp$alpha0,re.tmp$gamma0))
          }
          cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
          
          if(score.return)
          {
            score<-cbind(score,re.tmp$score)
          }
          
          if(verbose)
          {
            print(paste0("Component ",ncol(Theta.est)))
          }
        }
      }else
      {
        break
      }
    }
    
    colnames(Theta.est)=colnames(coef.est)=colnames(ancoef.est)=colnames(cp.time)<-paste0("C",1:ncol(Theta.est))
    
    cp.time<-cbind(cp.time,apply(cp.time,1,sum))
    colnames(cp.time)[ncol(cp.time)]<-"Total"
    
    if(score.return)
    {
      colnames(score)<-paste0("C",1:ncol(Theta.est))
    }
    
    DfD.out<-diag.level(M,Theta.est)
  }
  
  theta.orth<-t(Theta.est)%*%Theta.est
  
  if(score.return)
  {
    re<-list(theta=Theta.est,coef=coef.est,coef.other=ancoef.est,orthogonality=theta.orth,DfD=DfD.out,score=score,time=cp.time)
  }else
  {
    re<-list(theta=Theta.est,coef=coef.est,coef.other=ancoef.est,orthogonality=theta.orth,DfD=DfD.out,time=cp.time)
  }
  
  return(re)
}
#################################################

#################################################
# Bootstrap inference
CAPMediation_boot<-function(X,M,Y,theta=NULL,H=NULL,boot=TRUE,sims=1000,boot.ci.type=c("se","perc"),conf.level=0.95,seed.boot=100,verbose=TRUE)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # theta: components
  # H: a p by p positive definite matrix for theta constraint
  
  n<-length(Y)
  p<-ncol(M[[1]])
  
  nX<-ncol(X)
  q<-ncol(X)-1
  if(is.null(colnames(X)))
  {
    if(q>0)
    {
      X.names<-c("X",paste0("W",1:q)) 
    }else
    {
      X.names<-c("X")
    }
    colnames(X)<-X.names
  }else
  {
    X.names<-colnames(X)
  }
  
  M.cov<-.med_cov(M)
  
  if(is.null(H))
  {
    H<-apply(M.cov,c(1,2),mean,na.rm=TRUE)
  }
  
  if(boot)
  {
    if(is.null(theta))
    {
      stop("Error! Need theta value.")
    }else
    {
      coef.boot<-matrix(NA,5,sims)
      rownames(coef.boot)<-c("alpha","beta","gamma","IE","DE")
      if(q>0)
      {
        coef.other.boot<-matrix(NA,(q+1)*2,sims) 
        rownames(coef.other.boot)<-paste0(rep(c("M","Y"),each=q+1),"-",rep(c("Intercept",X.names[-1]),2))
      }else
      {
        coef.other.boot<-matrix(NA,2,sims)
        rownames(coef.other.boot)<-paste0(c("M","Y"),"-","Intercept")
      }
      
      for(b in 1:sims)
      {
        set.seed(seed.boot+b)
        idx.tmp<-sample(1:n,n,replace=TRUE)
        
        Ytmp<-Y[idx.tmp]
        Xtmp<-matrix(X[idx.tmp,],ncol=nX)
        colnames(Xtmp)<-X.names
        Mtmp<-M[idx.tmp]
        
        re.tmp<-NULL
        try(re.tmp<-CAPMediation_coef(X=Xtmp,M=Mtmp,Y=Ytmp,theta=theta))
        if(is.null(re.tmp)==FALSE)
        {
          coef.boot[,b]<-c(re.tmp$alpha[1],re.tmp$beta,re.tmp$gamma[1],re.tmp$IE,re.tmp$gamma[1])
          coef.other.boot[,b]<-c(re.tmp$alpha0,re.tmp$alpha[-1],re.tmp$gamma0,re.tmp$gamma[-1])
        }
        
        if(verbose)
        {
          print(paste0("Bootstrap sample ",b))
        }
      }
      
      coef.est<-matrix(NA,nrow(coef.boot),6)
      rownames(coef.est)<-rownames(coef.boot)
      colnames(coef.est)<-c("Estimate","SE","statistics","pvalue","LB","UB")
      coef.est[,1]<-apply(coef.boot,1,mean,na.rm=TRUE)
      coef.est[,2]<-apply(coef.boot,1,sd,na.rm=TRUE)
      coef.est[,3]<-coef.est[,1]/coef.est[,2]
      coef.est[,4]<-(1-pnorm(abs(coef.est[,3])))*2
      
      coef.other.est<-matrix(NA,nrow(coef.other.boot),6)
      rownames(coef.other.est)<-rownames(coef.other.boot)
      colnames(coef.other.est)<-c("Estimate","SE","statistics","pvalue","LB","UB")
      coef.other.est[,1]<-apply(coef.other.boot,1,mean,na.rm=TRUE)
      coef.other.est[,2]<-apply(coef.other.boot,1,sd,na.rm=TRUE)
      coef.other.est[,3]<-coef.other.est[,1]/coef.other.est[,2]
      coef.other.est[,4]<-(1-pnorm(abs(coef.other.est[,3])))*2
      
      if(boot.ci.type[1]=="se")
      {
        coef.est[,5]<-coef.est[,1]-qnorm(1-(1-conf.level)/2)*coef.est[,2]
        coef.est[,6]<-coef.est[,1]+qnorm(1-(1-conf.level)/2)*coef.est[,2]
        
        coef.other.est[,5]<-coef.other.est[,1]-qnorm(1-(1-conf.level)/2)*coef.other.est[,2]
        coef.other.est[,6]<-coef.other.est[,1]+qnorm(1-(1-conf.level)/2)*coef.other.est[,2]
      }
      if(boot.ci.type[1]=="perc")
      {
        coef.est[,5]<-apply(coef.boot,1,quantile,probs=(1-conf.level)/2)
        coef.est[,6]<-apply(coef.boot,1,quantile,probs=1-(1-conf.level)/2)
        
        coef.other.est[,5]<-apply(coef.other.boot,1,quantile,probs=(1-conf.level)/2)
        coef.other.est[,6]<-apply(coef.other.boot,1,quantile,probs=1-(1-conf.level)/2)
      }
    }
    
    re<-list(coef=coef.est,coef.other=coef.other.est,coef.boot=coef.boot,coef.other.boot=coef.other.boot)
    
    return(re)
  }else
  {
    stop("Error!")
  }
}
#################################################

#################################################
# CAP mediation refit with k components
CAPMediation_refit_Mcov<-function(X,M.cov,Y,Theta)
{
  # X: (X, W) matrix
  # M.cov: p by p by n array, covariance matrix of M
  # Y: n by 1 vector
  # Theta: p by k projection matrix
  
  k<-ncol(Theta)
  if(k==1)
  {
    re<-CAPMediation_coef_Mcov(X,M.cov,Y,theta=Theta)
  }else
  {
    n<-length(Y)
    
    q<-ncol(X)-1
    if(is.null(colnames(X)))
    {
      if(q>0)
      {
        colnames(X)<-c("X",paste0("W",1:q)) 
      }else
      {
        colnames(X)<-c("X")
      }
    }
    
    score<-matrix(NA,n,k)
    colnames(score)<-paste0("C",1:k)
    for(jj in 1:k)
    {
      score[,jj]<-as.numeric(.Call("med_score_cpp",M.cov,as.numeric(Theta[,jj])))
    }
    
    # beta and gamma estimate
    Z<-cbind(rep(1,n),X,log(score))
    colnames(Z)<-c("Intercept",colnames(X),paste0("M",1:k))
    mu.est<-c(ginv(t(Z)%*%Z)%*%t(Z)%*%Y)
    
    gamma0.est<-mu.est[1]
    gamma.est<-mu.est[2:(length(mu.est)-k)]
    beta.est<-mu.est[(length(mu.est)-k+1):length(mu.est)]
    
    # sigma2 estimate
    sigma2.est<-mean((Y-Z%*%mu.est)^2,na.rm=TRUE)
    
    # estimate alpha0, alpha0.rnd, alpha using mixed effects model
    alpha0.est<-rep(NA,k)
    alpha.est<-matrix(NA,ncol(X),k)
    rownames(alpha.est)<-colnames(X)
    colnames(alpha.est)<-paste0("M",1:k)
    alpha0.rnd.est<-matrix(NA,n,k)
    colnames(alpha0.rnd.est)<-paste0("M",1:k)
    tau2.est<-rep(NA,k)
    for(jj in 1:k)
    {
      dtmp<-data.frame(ID=1:n,score=log(score[,jj]),X)
      eval(parse(text=paste0("fit.tmp<-lme(score~",paste(colnames(X),collapse="+"),",random=~1|ID,data=dtmp,control=lmeControl(opt='optim'))")))
      
      alpha0.est[jj]<-fit.tmp$coefficients$fixed[1]
      alpha.est[,jj]<-fit.tmp$coefficients$fixed[-1]
      alpha0.rnd.est[,jj]<-c(fit.tmp$coefficients$random$ID+fit.tmp$coefficients$fixed[1])
      tau2.est[jj]<-mean((alpha0.rnd.est[,jj]-alpha0.est[jj])^2,na.rm=TRUE)
    }
    
    IE.est<-alpha.est[1,]*beta.est
    
    names(alpha0.est)=names(beta.est)=names(tau2.est)<-paste0("M",1:k)
    names(gamma.est)<-colnames(X)
    re<-list(theta=Theta,alpha=alpha.est,beta=beta.est,gamma=gamma.est,IE=IE.est,alpha0=alpha0.est,alpha0.rnd=alpha0.rnd.est,gamma0=gamma0.est,tau2=tau2.est,sigma2=sigma2.est)
    
    return(re)
  }
}
CAPMediation_refit<-function(X,M,Y,Theta)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # Theta: p by k projection matrix
  
  M.cov<-.med_cov(M)
  
  re<-CAPMediation_refit_Mcov(X,M.cov,Y,Theta)
  
  return(re)
}
#################################################

#################################################
# CAP mediation refit bootstrap inference
CAPMediation_refit_boot<-function(X,M,Y,Theta=NULL,H=NULL,boot=TRUE,sims=1000,boot.ci.type=c("se","perc"),conf.level=0.95,seed.boot=100,verbose=TRUE)
{
  # X: (X, W) matrix
  # M: a list of length n, M
  # Y: n by 1 vector
  # Theta: p by k projection matrix
  # H: a p by p positive definite matrix for theta constraint
  
  if(is.null(Theta)|ncol(Theta)==1)
  {
    re<-CAPMediation_boot(X,M,Y,theta=Theta,H=H,boot=boot,sims=sims,boot.ci.type=boot.ci.type,conf.level=conf.level,seed.boot=seed.boot,verbose=verbose)
  }else
  {
    k<-ncol(Theta)
    
    n<-length(Y)
    p<-ncol(M[[1]])
    
    nX<-ncol(X)
    q<-ncol(X)-1
    if(is.null(colnames(X)))
    {
      if(q>0)
      {
        X.names<-c("X",paste0("W",1:q)) 
      }else
      {
        X.names<-c("X")
      }
      colnames(X)<-X.names
    }else
    {
      X.names<-colnames(X)
    }
    
    M.cov<-.med_cov(M)
    
    if(is.null(H))
    {
      H<-apply(M.cov,c(1,2),mean,na.rm=TRUE)
    }
    
    if(boot)
    {
      coef.boot<-matrix(NA,3*k+2,sims)
      rownames(coef.boot)<-c(paste0("alpha_M",1:k),paste0("beta_M",1:k),"gamma",paste0("IE_M",1:k),"DE")
      if(q>0)
      {
        coef.other.boot<-matrix(NA,(q+1)*(k+1),sims) 
        rownames(coef.other.boot)<-paste0(rep(c(paste0("M",1:k),"Y"),each=q+1),"-",rep(c("Intercept",X.names[-1]),k+1))
      }else
      {
        coef.other.boot<-matrix(NA,k+1,sims)
        rownames(coef.other.boot)<-paste0(c(paste0("M",1:k),"Y"),"-","Intercept")
      }
      
      for(b in 1:sims)
      {
        set.seed(seed.boot+b)
        idx.tmp<-sample(1:n,n,replace=TRUE)
        
        Ytmp<-Y[idx.tmp]
        Xtmp<-matrix(X[idx.tmp,],ncol=nX)
        colnames(Xtmp)<-X.names
        Mtmp<-M[idx.tmp]
        
        re.tmp<-NULL
        try(re.tmp<-CAPMediation_refit(X=Xtmp,M=Mtmp,Y=Ytmp,Theta=Theta))
        if(is.null(re.tmp)==FALSE)
        {
          coef.boot[,b]<-c(re.tmp$alpha[1,],re.tmp$beta,re.tmp$gamma[1],re.tmp$IE,re.tmp$gamma[1])
          coef.other.boot[,b]<-c(c(rbind(re.tmp$alpha0,re.tmp$alpha[-1,])),re.tmp$gamma0,re.tmp$gamma[-1])
        }
        
        if(verbose)
        {
          print(paste0("Bootstrap sample ",b))
        }
      }
      
      coef.est<-matrix(NA,nrow(coef.boot),6)
      rownames(coef.est)<-rownames(coef.boot)
      colnames(coef.est)<-c("Estimate","SE","statistics","pvalue","LB","UB")
      coef.est[,1]<-apply(coef.boot,1,mean,na.rm=TRUE)
      coef.est[,2]<-apply(coef.boot,1,sd,na.rm=TRUE)
      coef.est[,3]<-coef.est[,1]/coef.est[,2]
      coef.est[,4]<-(1-pnorm(abs(coef.est[,3])))*2
      
      coef.other.est<-matrix(NA,nrow(coef.other.boot),6)
      rownames(coef.other.est)<-rownames(coef.other.boot)
      colnames(coef.other.est)<-c("Estimate","SE","statistics","pvalue","LB","UB")
      coef.other.est[,1]<-apply(coef.other.boot,1,mean,na.rm=TRUE)
      coef.other.est[,2]<-apply(coef.other.boot,1,sd,na.rm=TRUE)
      coef.other.est[,3]<-coef.other.est[,1]/coef.other.est[,2]
      coef.other.est[,4]<-(1-pnorm(abs(coef.other.est[,3])))*2
      
      if(boot.ci.type[1]=="se")
      {
        coef.est[,5]<-coef.est[,1]-qnorm(1-(1-conf.level)/2)*coef.est[,2]
        coef.est[,6]<-coef.est[,1]+qnorm(1-(1-conf.level)/2)*coef.est[,2]
        
        coef.other.est[,5]<-coef.other.est[,1]-qnorm(1-(1-conf.level)/2)*coef.other.est[,2]
        coef.other.est[,6]<-coef.other.est[,1]+qnorm(1-(1-conf.level)/2)*coef.other.est[,2]
      }
      if(boot.ci.type[1]=="perc")
      {
        coef.est[,5]<-apply(coef.boot,1,quantile,probs=(1-conf.level)/2)
        coef.est[,6]<-apply(coef.boot,1,quantile,probs=1-(1-conf.level)/2)
        
        coef.other.est[,5]<-apply(coef.other.boot,1,quantile,probs=(1-conf.level)/2)
        coef.other.est[,6]<-apply(coef.other.boot,1,quantile,probs=1-(1-conf.level)/2)
      }
      
      re<-list(coef=coef.est,coef.other=coef.other.est,coef.boot=coef.boot,coef.other.boot=coef.other.boot)
      
      return(re)
    }else
    {
      stop("Error!")
    }
  }
}
#################################################




