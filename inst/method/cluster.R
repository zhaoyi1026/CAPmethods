############################
# CAP parsimonious clustering (PCL)
#
# V4: RcppArmadillo-accelerated update of V3. Same interface and results; the hot
# loops (per-subject second moment Smat, projected scores, the gamma-update
# weighted accumulation, eigen.solve) call a compiled kernel. The multinomial
# membership fit (brglm2::brmultinom) stays in R. See CAP-clustering/README.md.
############################

# library("nnet")       # multinomial regression



############################
# Load the compiled kernel (cluster_kernels.cpp -> shared object), building it
# once if necessary. Set option "cluster.v4.dir" or variable CLUSTER_V4_DIR to
# point at the folder with cluster_kernels.cpp; defaults to "CAP-clustering/V4".
# per-subject second moment cube; projected scores; weighted accumulation
.cluster_smat  <- function(Y) .Call("cluster_smat_cpp", lapply(Y, as.matrix))
.cluster_score <- function(Smat, v) as.numeric(.Call("cluster_score_cpp", Smat, as.numeric(v)))
.cluster_accum <- function(Smat, w) .Call("cluster_accum_cpp", Smat, as.numeric(w))

# if(!require("pacman"))
# {
#   install.packages("pacman")
# }
# library("pacman")
# pacman::p_load(brglm2,MASS)

############################
# full likelihood function (with latent cluster indicator)

ll.full<-function(Y,X,W,eta.ind,gamma,beta,alpha)
{
  n<-nrow(X)
  p<-ncol(Y[[1]])
  q1<-ncol(X)
  q2<-ncol(W)
  Tvec<-sapply(Y,nrow)
  K<-ncol(eta.ind)
  
  Smat<-.cluster_smat(Y)$Smat
  pi.mat<-matrix(NA,n,K)
  for(i in 1:n)
  {
    otmp<-rep(NA,K)
    for(k in 1:K)
    {
      otmp[k]<-exp(t(W[i,])%*%alpha[,k])
    }
    pi.mat[i,]<-otmp/sum(otmp)
  }
  
  ll1=ll2<-0
  for(k in 1:K)
  {
    ll1<-ll1+sum(Tvec*eta.ind[,k]*log(pi.mat[,k]))
    
    ll2<-ll2+sum(Tvec*eta.ind[,k]*(-X%*%beta[,k]-exp(-X%*%beta[,k])*.cluster_score(Smat, gamma))/2)
  }
  
  ll<-ll1+ll2
  return(ll)
}
############################

############################
# given gamma
# estimate model coefficients and cluster
capPCL_coef<-function(Y,X,W,gamma,ncluster=2,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE)
{
  n<-nrow(X)
  p<-ncol(Y[[1]])
  q1<-ncol(X)
  q2<-ncol(W)
  Tvec<-sapply(Y,nrow)
  
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",0:(q1-1))
  }
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-paste0("W",0:(q2-1))
  }
  
  # Estimate covariance matrix
  Smat<-.cluster_smat(Y)$Smat
  Zmat<-matrix(NA,n,max(Tvec))
  for(i in 1:n)
  {
    Zmat[i,1:Tvec[i]]<-Y[[i]]%*%gamma
  }

  score<-.cluster_score(Smat, gamma)
  
  # set initial values (some random number)
  set.seed(100)
  beta0<-matrix(rnorm(q1*ncluster,mean=0,sd=1),q1,ncluster)
  rownames(beta0)<-colnames(X)
  alpha0<-matrix(rnorm(q2*ncluster,mean=0,sd=1),q2,ncluster)
  rownames(alpha0)<-colnames(W)
  colnames(beta0)=colnames(alpha0)<-paste0("cluster",1:ncluster)
  alpha0[,1]<-0         # first class coefficients zero
  
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
    
    #------------------------
    # E-step (updating cluster indicator eta)
    pi.mat<-matrix(NA,n,ncluster)
    phi.mat<-matrix(NA,n,ncluster)
    for(i in 1:n)
    {
      # cluster probability estimation: pi
      otmp<-rep(NA,ncluster)
      for(k in 1:ncluster)
      {
        otmp[k]<-exp(t(W[i,])%*%alpha0[,k])
      }
      pi.mat[i,]<-otmp/sum(otmp)
      
      # normal densities
      for(k in 1:ncluster)
      {
        phi.mat[i,k]<-exp(sum(log(dnorm(Zmat[i,1:Tvec[i]],mean=0,sd=sqrt(exp(t(X[i,])%*%beta0[,k]))))))
      }
    }
    
    eta.mat<-matrix(NA,n,ncluster)
    for(k in 1:ncluster)
    {
      eta.mat[,k]<-pi.mat[,k]*phi.mat[,k]/apply(pi.mat*phi.mat,1,sum)
      eta.mat[which(is.nan(eta.mat[,k])==TRUE),k]<-1/ncluster
    }
    #------------------------
    
    #------------------------
    # M-step 
    
    # update alpha: using multinomial log-linear regression
    cls.tmp<-apply(eta.mat,1,which.max)
    if(length(unique(cls.tmp))>1)
    {
      fit.tmp<-brmultinom(cls.tmp~0+W,weights=Tvec,ref=1)
      alpha.new<-matrix(0,q2,ncluster)
      rownames(alpha.new)<-colnames(W)
      colnames(alpha.new)<-paste0("cluster",1:ncluster)
      alpha.new[,2:ncluster]<-t(coef(fit.tmp))
      if(fit.tmp$converged==FALSE)
      {
        s<-max.itr+1
      }
    }else
    {
      alpha.new<-alpha0     # if only one cluster, do not update alpha
    }
    
    # update beta:
    beta.new<-matrix(NA,q1,ncluster)
    rownames(beta.new)<-colnames(X)
    colnames(beta.new)<-paste0("cluster",1:ncluster)
    for(k in 1:ncluster)
    {
      d1.tmp<-apply(t(apply(cbind(Tvec*eta.mat[,k]*(1-exp(-X%*%beta0[,k])*score)/2,X),1,function(x){return(x[1]*x[-1])})),2,sum)
      
      d2.tmp<-t(X)%*%diag(c(Tvec*eta.mat[,k]*exp(-X%*%beta0[,k])*score/2))%*%X
      
      beta.new[,k]<-beta0[,k]-ginv(d2.tmp)%*%d1.tmp
    }
    #------------------------
    
    if(trace)
    {
      alpha.trace<-cbind(alpha.trace,alpha.new)
      beta.trace<-cbind(beta.trace,beta.new)
      
      obj<-c(obj,ll.full(Y=Y,X=X,W=W,eta.ind=eta.mat,gamma=gamma,beta=beta.new,alpha=alpha.new))
    }
    
    diff<-max(c(c(abs(beta.new-beta0)),c(abs(alpha.new-alpha0))))
    
    beta0<-beta.new
    alpha0<-alpha.new
    
    # print(diff)
  }
  
  # an integer of clustering result eta
  eta.int<-t(apply(eta.mat,1,function(x){
    y<-rep(0,length(x))
    y[which.max(x)]<-1
    return(y)
  }))
  colnames(eta.mat)=colnames(eta.int)<-paste0("cluster",1:ncluster)
  # cluster indicator
  cl.ind<-apply(eta.int,1,function(x){return(which(x==1))})
  
  ll<-ll.full(Y=Y,X=X,W=W,eta.ind=eta.int,gamma=gamma,beta=beta.new,alpha=alpha.new)
  
  if(trace)
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma),beta=beta.new,alpha=alpha.new,class=cl.ind,eta=eta.int,eta.est=eta.mat,logLik=ll,nitr=s,convergence=(s<max.itr),score=score,
               beta.trace=beta.trace,alpha.trace=alpha.trace,logLik.trace=obj)
    }else
    {
      re<-list(gamma=c(gamma),beta=beta.new,alpha=alpha.new,class=cl.ind,eta=eta.int,eta.est=eta.mat,logLik=ll,nitr=s,convergence=(s<max.itr),
               beta.trace=beta.trace,alpha.trace=alpha.trace,logLik.trace=obj)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma),beta=beta.new,alpha=alpha.new,class=cl.ind,eta=eta.int,eta.est=eta.mat,logLik=ll,nitr=s,convergence=(s<max.itr),score=score)
    }else
    {
      re<-list(gamma=c(gamma),beta=beta.new,alpha=alpha.new,class=cl.ind,eta=eta.int,eta.est=eta.mat,logLik=ll,nitr=s,convergence=(s<max.itr))
    }
  }
  
  return(re)
}
############################

############################
# estimate both gamma and beta

# minimize t(x)%*%A%*%x such that t(x)%*%H%*%x=1
# eigenvectors and eigenvalues of A with respect to H
# H positive definite and symmetric
eigen.solve<-function(A,H)
{
  # smallest generalized eigenvector of A w.r.t. H  (RcppArmadillo kernel)
  as.numeric(.Call("cluster_eigen_solve_cpp", as.matrix(A), as.matrix(H)))
}

# finding the first component
capPCL_D1<-function(Y,X,W,H=NULL,ncluster=2,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,gamma0=NULL)
{
  n<-nrow(X)
  p<-ncol(Y[[1]])
  q1<-ncol(X)
  q2<-ncol(W)
  Tvec<-sapply(Y,nrow)
  
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",0:(q1-1))
  }
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-paste0("W",0:(q2-1))
  }
  
  # Estimate covariance matrix
  Smat<-.cluster_smat(Y)$Smat

  # H matrix
  if(is.null(H)==TRUE)
  {
    H<-.cluster_accum(Smat,Tvec)/sum(Tvec)
  }

  # set initial values (some random number)
  set.seed(100)
  beta0<-matrix(rnorm(q1*ncluster,mean=0,sd=1),q1,ncluster)
  rownames(beta0)<-colnames(X)
  alpha0<-matrix(rnorm(q2*ncluster,mean=0,sd=1),q2,ncluster)
  rownames(alpha0)<-colnames(W)
  colnames(beta0)=colnames(alpha0)<-paste0("cluster",1:ncluster)
  alpha0[,1]<-0         # first class coefficients zero
  if(is.null(gamma0))
  {
    H.eigen<-eigen(H)
    gamma0<-H.eigen$vectors[,round(runif(1,min=1,max=p))]
  }
  
  if(trace)
  {
    alpha.trace<-NULL
    beta.trace<-NULL
    gamma.trace<-NULL
    obj<-NULL
  }
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    #------------------------
    # update alpha, beta, and eta
    otmp<-capPCL_coef(Y=Y,X=X,W=W,gamma=gamma0,ncluster=ncluster,max.itr=max.itr,tol=tol,trace=FALSE,score.return=TRUE)
    if(otmp$convergence==FALSE)
    {
      s<-max.itr+1
    }
    
    eta.mat<-otmp$eta.est
    beta.new<-otmp$beta
    alpha.new<-otmp$alpha
    #------------------------
    
    #------------------------
    # update gamma:
    # per-subject weight = sum_k Tvec_i * eta_ik * exp(-x_i'beta_k) / 2
    wA<-(Tvec/2)*rowSums(eta.mat*exp(-X%*%beta.new))
    A<-.cluster_accum(Smat,wA)
    gamma.new<-eigen.solve(A,H)
    #------------------------
    
    if(trace)
    {
      alpha.trace<-cbind(alpha.trace,alpha.new)
      beta.trace<-cbind(beta.trace,beta.new)
      gamma.trace<-cbind(gamma.trace,gamma.new)
      
      obj<-c(obj,ll.full(Y=Y,X=X,W=W,eta.ind=eta.mat,gamma=gamma.new,beta=beta.new,alpha=alpha.new))
    }
    
    diff<-max(c(c(abs(beta.new-beta0)),c(abs(alpha.new-alpha0))))
    
    beta0<-beta.new
    alpha0<-alpha.new
    gamma0<-gamma.new
    
    # print(diff)
    # print(s)
  }
  
  # renormalize gamma vector
  gamma.new<-c(gamma.new)/sqrt(sum(gamma.new^2))
  if(gamma.new[which.max(abs(gamma.new))]<0)
  {
    gamma.new<--gamma.new
  }
  # estimate coefficients using the new gamma
  otmp<-capPCL_coef(Y=Y,X=X,W=W,gamma=gamma.new,ncluster=ncluster,max.itr=max.itr,tol=tol,trace=FALSE,score.return=TRUE)
  beta.new<-otmp$beta
  alpha.new<-otmp$alpha
  score<-.cluster_score(Smat, gamma.new)
  
  # an integer of clustering result eta
  eta.mat<-otmp$eta.est
  eta.int<-t(apply(eta.mat,1,function(x){
    y<-rep(0,length(x))
    y[which.max(x)]<-1
    return(y)
  }))
  colnames(eta.mat)=colnames(eta.int)<-paste0("cluster",1:ncluster)
  # cluster indicator
  cl.ind<-apply(eta.int,1,function(x){return(which(x==1))})
  
  ll<-otmp$logLik
  
  if(trace)
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma.new),beta=beta.new,alpha=alpha.new,class=cl.ind,eta=eta.int,eta.est=eta.mat,logLik=ll,nitr=s,convergence=(s<max.itr),score=score,
               beta.trace=beta.trace,alpha.trace=alpha.trace,logLik.trace=obj)
    }else
    {
      re<-list(gamma=c(gamma.new),beta=beta.new,alpha=alpha.new,class=cl.ind,eta=eta.int,eta.est=eta.mat,logLik=ll,nitr=s,convergence=(s<max.itr),
               beta.trace=beta.trace,alpha.trace=alpha.trace,logLik.trace=obj)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma.new),beta=beta.new,alpha=alpha.new,class=cl.ind,eta=eta.int,eta.est=eta.mat,logLik=ll,nitr=s,convergence=(s<max.itr),score=score)
    }else
    {
      re<-list(gamma=c(gamma.new),beta=beta.new,alpha=alpha.new,class=cl.ind,eta=eta.int,eta.est=eta.mat,logLik=ll,nitr=s,convergence=(s<max.itr))
    }
  }
  
  return(re)
}

# try several initial value of gamma and maximize the log-likelihood function
capPCL_D1_opt<-function(Y,X,W,H=NULL,ncluster=2,max.itr=1000,tol=1e-4,trace=TRUE,score.return=TRUE,gamma0.mat=NULL,ninitial=NULL,seed=100)
{
  n<-nrow(X)
  p<-ncol(Y[[1]])
  q1<-ncol(X)
  q2<-ncol(W)
  Tvec<-sapply(Y,nrow)
  
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",0:(q1-1))
  }
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-paste0("W",0:(q2-1))
  }
  
  if(is.null(ninitial))
  {
    ninitial<-min(p,10)
  }else
  {
    if(ninitial>p)
    {
      ninitial<-p
    }
  }
  
  if(is.null(gamma0.mat))
  {
    # Estimate covariance matrix
    Smat<-.cluster_smat(Y)$Smat

    # H matrix
    Savg<-.cluster_accum(Smat,Tvec)/sum(Tvec)

    Savg.eigen<-eigen(Savg)
    gamma0.mat<-Savg.eigen$vectors
  }
  set.seed(seed)
  gamma0.mat<-matrix(gamma0.mat[,sort(sample(1:ncol(gamma0.mat),ninitial,replace=FALSE))],ncol=ninitial)
  
  re.tmp<-vector("list",length=ncol(gamma0.mat))
  ll.vec<-matrix(NA,ncol(gamma0.mat),2)
  colnames(ll.vec)<-c("scaled","unsclaed")
  for(kk in 1:ncol(gamma0.mat))
  {
    try(re.tmp[[kk]]<-capPCL_D1(Y=Y,X=X,W=W,H=H,ncluster=ncluster,max.itr=max.itr,tol=tol,trace=TRUE,score.return=score.return,gamma0=gamma0.mat[,kk]))
    
    if(is.null(re.tmp[[kk]])==FALSE)
    {
      if(re.tmp[[kk]]$convergence==TRUE)
      {
        ll.vec[kk,1]<-re.tmp[[kk]]$logLik
        ll.vec[kk,2]<-re.tmp[[kk]]$logLik.trace[length(re.tmp[[kk]]$logLik.trace)]
      }
    }
    # print(kk)
  }
  opt.idx<-which.max(ll.vec[,2])
  re<-re.tmp[[opt.idx]]
  
  return(re)
}

# higher-order components
capPCL_Dk<-function(Y,X,W,Gamma0=NULL,H=NULL,ncluster=2,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,gamma0.mat=NULL,ninitial=NULL,seed=100)
{
  if(is.null(Gamma0))
  {
    re<-capPCL_D1_opt(Y=Y,X=X,W=W,H=H,ncluster=ncluster,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed)
  }else
  {
    n<-nrow(X)
    p<-ncol(Y[[1]])
    q1<-ncol(X)
    q2<-ncol(W)
    Tvec<-sapply(Y,nrow)
    
    if(is.null(colnames(X))==TRUE)
    {
      colnames(X)<-paste0("X",0:(q1-1))
    }
    if(is.null(colnames(W))==TRUE)
    {
      colnames(W)<-paste0("W",0:(q2-1))
    }
    
    # estimate model coefficients
    p0<-ncol(Gamma0)
    beta0<-array(NA,c(q1,ncluster,p0))
    dimnames(beta0)[[1]]<-colnames(X)
    alpha0<-array(NA,c(q2,ncluster,p0))
    dimnames(alpha0)[[1]]<-colnames(W)
    dimnames(beta0)[[2]]=dimnames(alpha0)[[2]]<-paste0("cluster",1:ncluster)
    dimnames(beta0)[[3]]=dimnames(alpha0)[[3]]<-paste0("C",1:p0)
    for(j in 1:p0)
    {
      otmp<-capPCL_coef(Y=Y,X=X,W=W,gamma=Gamma0[,j],ncluster=ncluster,max.itr=max.itr,tol=tol,trace=FALSE,score.return=TRUE)
      beta0[,,j]<-otmp$beta
      alpha0[,,j]<-otmp$alpha
    }
    # create new data
    Ytmp<-vector("list",length=n)
    for(i in 1:n)
    {
      Ynew<-Y[[i]]-Y[[i]]%*%Gamma0%*%t(Gamma0)
      Ynew.svd<-svd(Ynew)
      Dnew<-c(Ynew.svd$d[1:(p-p0)],sqrt(exp(beta0[1,1,])*nrow(Ynew)))
      Ytmp[[i]]<-Ynew.svd$u%*%diag(Dnew)%*%t(Ynew.svd$v)
    }
    
    if(p0>=2)
    {
      re<-capPCL_D1_opt(Y=Ytmp,X=X,W=W,H=H,ncluster=ncluster,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0.mat=gamma0.mat,
                        ninitial=max(min(10,p-p0-1),3),seed=seed)
    }else
    {
      re<-capPCL_D1_opt(Y=Ytmp,X=X,W=W,H=H,ncluster=ncluster,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0.mat=gamma0.mat,
                        ninitial=ninitial,seed=seed)
    }

    
    re$orthogonality<-t(re$gamma)%*%Gamma0
  }
  
  return(re)
}
############################

############################
# level of diagonalization
diag.level<-function(Y,Gamma)
{
  n<-length(Y)
  
  if(ncol(Gamma)==1)
  {
    re<-list(avg.level=1,sub.level=rep(1,n))
  }else
  {
    p<-ncol(Y[[1]])
    Tvec<-sapply(Y,nrow)
    r<-ncol(Gamma)
    
    Smat<-.cluster_smat(Y)$Smat
    dl.sub<-matrix(NA,n,r)
    colnames(dl.sub)<-paste0("C",1:r)
    dl.sub[,1]<-1
    for(i in 1:n)
    {
      cov.tmp<-Smat[,,i]

      for(j in 2:r)
      {
        gamma.tmp<-Gamma[,1:j]
        mat.tmp<-t(gamma.tmp)%*%cov.tmp%*%gamma.tmp
        dl.sub[i,j]<-det(diag(diag(mat.tmp)))/det(mat.tmp)
      }
    }
    
    pmean<-apply(dl.sub,2,function(y){return(prod(apply(cbind(y,Tvec),1,function(x){return(x[1]^(x[2]/sum(Tvec)))})))})
    
    re<-list(avg.level=pmean,sub.level=dl.sub)
  }
  
  return(re)
}
############################

############################
# finding first r directions
capPCL<-function(Y,X,W,ncluster=2,stop.crt=c("nD","DfD"),DfD.thred=2,nD=NULL,
                 max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE,gamma0.mat=NULL,ninitial=NULL,seed=100,verbose=TRUE)
{
  if(stop.crt[1]=="nD"&is.null(nD))
  {
    stop.crt<-"DfD"
  }
  
  n<-nrow(X)
  p<-ncol(Y[[1]])
  q1<-ncol(X)
  q2<-ncol(W)
  Tvec<-sapply(Y,nrow)
  
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",0:(q1-1))
  }
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-paste0("W",0:(q2-1))
  }
  
  #--------------------------------------------
  # First direction
  tm1<-system.time(re1<-capPCL_D1_opt(Y=Y,X=X,W=W,ncluster=ncluster,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed))
  
  Gamma.est<-matrix(re1$gamma,ncol=1)
  beta.est<-list()
  alpha.est<-list()
  beta.est[[1]]<-re1$beta
  alpha.est[[1]]<-re1$alpha
  class.est<-matrix(re1$class,ncol=1)
  eta.est<-list()
  eta.est[[1]]<-re1$eta.est
  logLik.est<-c(re1$logLik)
  nitr<-re1$nitr
  convergence<-re1$convergence
  
  cp.time<-matrix(as.numeric(tm1[1:3]),ncol=1)
  rownames(cp.time)<-c("user","system","elapsed")
  
  if(score.return)
  {
    score<-matrix(re1$score,ncol=1)
  }
  
  if(verbose)
  {
    print(paste0("Component ",ncol(Gamma.est)))
  }
  #--------------------------------------------
  
  #--------------------------------------------
  # nD stopping
  if(stop.crt[1]=="nD")
  {
    if(nD>1)
    {
      for(j in 2:nD)
      {
        re.tmp<-NULL
        try(tm.tmp<-system.time(re.tmp<-capPCL_Dk(Y=Y,X=X,W=W,Gamma0=Gamma.est,ncluster=ncluster,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,
                                                  gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed)))
        
        if(is.null(re.tmp)==FALSE)
        {
          Gamma.est<-cbind(Gamma.est,re.tmp$gamma)
          beta.est[[j]]<-re.tmp$beta
          alpha.est[[j]]<-re.tmp$alpha
          class.est<-cbind(class.est,re.tmp$class)
          eta.est[[j]]<-re.tmp$eta.est
          logLik.est<-c(logLik.est,re.tmp$logLik)
          nitr<-c(nitr,re.tmp$nitr)
          convergence<-c(convergence,re.tmp$convergence)
          
          cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
          
          if(score.return)
          {
            score<-cbind(score,re.tmp$score)
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
    
    colnames(Gamma.est)=names(beta.est)=names(alpha.est)=colnames(class.est)=names(eta.est)=names(logLik.est)=names(nitr)=names(convergence)=colnames(cp.time)<-paste0("C",1:ncol(Gamma.est))
    if(score.return)
    {
      colnames(score)<-paste0("C",1:ncol(Gamma.est))
    }
    DfD<-diag.level(Y,Gamma.est)
  }
  #--------------------------------------------
  
  #--------------------------------------------
  # DfD stopping
  if(stop.crt[1]=="DfD")
  {
    nD<-1
    
    DfD.tmp<-1
    while(DfD.tmp<DfD.thred)
    {
      re.tmp<-NULL
      try(tm.tmp<-system.time(re.tmp<-capPCL_Dk(Y=Y,X=X,W=W,Gamma0=Gamma.est,ncluster=ncluster,max.itr=max.itr,tol=tol,trace=trace,score.return=score.return,
                                                gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed)))
      
      if(is.null(re.tmp)==FALSE)
      {
        nD<-nD+1
        
        DfD<-diag.level(Y,cbind(Gamma.est,re.tmp$gamma))
        DfD.tmp<-DfD$avg.level[nD]
        
        if(is.nan(DfD.tmp)==FALSE)
        {
          if(DfD.tmp<DfD.thred)
          {
            Gamma.est<-cbind(Gamma.est,re.tmp$gamma)
            beta.est[[nD]]<-re.tmp$beta
            alpha.est[[nD]]<-re.tmp$alpha
            class.est<-cbind(class.est,re.tmp$class)
            eta.est[[nD]]<-re.tmp$eta.est
            logLik.est<-c(logLik.est,re.tmp$logLik)
            nitr<-c(nitr,re.tmp$nitr)
            convergence<-c(convergence,re.tmp$convergence)
            
            cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
            
            if(score.return)
            {
              score<-cbind(score,re.tmp$score)
            }
            
            if(verbose)
            {
              print(paste0("Component ",ncol(Gamma.est)))
            }
          }
        }else
        {
          DfD.tmp<-DfD.thred+100
        }
      }else
      {
        break
      }
      
      colnames(Gamma.est)=names(beta.est)=names(alpha.est)=colnames(class.est)=names(eta.est)=names(logLik.est)=names(nitr)=names(convergence)=colnames(cp.time)<-paste0("C",1:ncol(Gamma.est))
      if(score.return)
      {
        colnames(score)<-paste0("C",1:ncol(Gamma.est))
      }
      DfD<-diag.level(Y,Gamma.est)
    }
  }
  #--------------------------------------------
  
  cp.time.new<-cbind(cp.time,apply(cp.time,1,sum))
  rownames(cp.time.new)<-rownames(cp.time)
  colnames(cp.time.new)<-c(colnames(cp.time),"Total")
  
  if(score.return)
  {
    re<-list(gamma=Gamma.est,beta=beta.est,alpha=alpha.est,class=class.est,eta.est=eta.est,logLik=logLik.est,nitr=nitr,convergence=convergence,score=score,DfD=DfD)
    re$orthogonality<-t(Gamma.est)%*%Gamma.est
  }else
  {
    re<-list(gamma=Gamma.est,beta=beta.est,alpha=alpha.est,class=class.est,eta.est=eta.est,logLik=logLik.est,nitr=nitr,convergence=convergence,DfD=DfD)
    re$orthogonality<-t(Gamma.est)%*%Gamma.est
  }
  
  return(re)
}
############################

############################
# estimate of beta and alpha with given classification
capPCL_coef_class<-function(Y,X,W,gamma,eta,max.itr=1000,tol=1e-4,trace=FALSE,score.return=TRUE)
{
  n<-nrow(X)
  p<-ncol(Y[[1]])
  q1<-ncol(X)
  q2<-ncol(W)
  Tvec<-sapply(Y,nrow)
  
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",0:(q1-1))
  }
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-paste0("W",0:(q2-1))
  }
  
  # Estimate covariance matrix
  Smat<-.cluster_smat(Y)$Smat
  Zmat<-matrix(NA,n,max(Tvec))
  for(i in 1:n)
  {
    Zmat[i,1:Tvec[i]]<-Y[[i]]%*%gamma
  }

  score<-.cluster_score(Smat, gamma)
  
  ncluster<-ncol(eta)
  
  # estimate of alpha
  cls.tmp<-apply(eta,1,which.max)
  fit.tmp<-brmultinom(cls.tmp~0+W,weights=Tvec,ref=1)
  alpha.est<-matrix(0,q2,ncluster)
  rownames(alpha.est)<-colnames(W)
  colnames(alpha.est)<-paste0("cluster",1:ncluster)
  alpha.est[,2:ncluster]<-t(coef(fit.tmp))
  
  # estimate of beta
  # set initial values (some random number)
  set.seed(100)
  beta0<-matrix(rnorm(q1*ncluster,mean=0,sd=1),q1,ncluster)
  rownames(beta0)<-colnames(X)
  colnames(beta0)<-paste0("cluster",1:ncluster)
  
  if(trace)
  {
    beta.trace<-NULL
    obj<-NULL
  }
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update beta:
    beta.new<-matrix(NA,q1,ncluster)
    rownames(beta.new)<-colnames(X)
    colnames(beta.new)<-paste0("cluster",1:ncluster)
    for(k in 1:ncluster)
    {
      d1.tmp<-apply(t(apply(cbind(Tvec*eta[,k]*(1-exp(-X%*%beta0[,k])*score)/2,X),1,function(x){return(x[1]*x[-1])})),2,sum)
      
      d2.tmp<-t(X)%*%diag(c(Tvec*eta[,k]*exp(-X%*%beta0[,k])*score/2))%*%X
      
      beta.new[,k]<-beta0[,k]-ginv(d2.tmp)%*%d1.tmp
    }
    
    if(trace)
    {
      beta.trace<-cbind(beta.trace,beta.new)
      
      obj<-c(obj,ll.full(Y=Y,X=X,W=W,eta.ind=eta,gamma=gamma,beta=beta.new,alpha=alpha.est))
    }
    
    diff<-max(abs(beta.new-beta0))
    
    beta0<-beta.new
    
    # print(diff)
  }
  
  # an integer of clustering result eta
  eta.int<-t(apply(eta,1,function(x){
    y<-rep(0,length(x))
    y[which.max(x)]<-1
    return(y)
  }))
  colnames(eta.int)<-paste0("cluster",1:ncluster)
  # cluster indicator
  cl.ind<-apply(eta.int,1,function(x){return(which(x==1))})
  
  ll<-ll.full(Y=Y,X=X,W=W,eta.ind=eta.int,gamma=gamma,beta=beta.new,alpha=alpha.est)
  
  if(trace)
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma),beta=beta.new,alpha=alpha.est,class=cl.ind,eta=eta,logLik=ll,nitr=s,convergence=(s<max.itr),score=score,
               beta.trace=beta.trace,alpha.trace=alpha.trace,logLik.trace=obj)
    }else
    {
      re<-list(gamma=c(gamma),beta=beta.new,alpha=alpha.est,class=cl.ind,eta=eta,logLik=ll,nitr=s,convergence=(s<max.itr),
               beta.trace=beta.trace,alpha.trace=alpha.trace,logLik.trace=obj)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma),beta=beta.new,alpha=alpha.est,class=cl.ind,eta=eta,logLik=ll,nitr=s,convergence=(s<max.itr),score=score)
    }else
    {
      re<-list(gamma=c(gamma),beta=beta.new,alpha=alpha.est,class=cl.ind,eta=eta,logLik=ll,nitr=s,convergence=(s<max.itr))
    }
  }
  
  return(re)
}

# Inference of beta and alpha
# based on bootstrap
capPCL_coef_boot<-function(Y,X,W,gamma,eta,boot=TRUE,sims=1000,boot.ci.type=c("se","perc"),conf.level=0.95,max.itr=1000,tol=1e-4,verbose=TRUE)
{
  n<-nrow(X)
  p<-ncol(Y[[1]])
  q1<-ncol(X)
  q2<-ncol(W)
  Tvec<-sapply(Y,nrow)
  
  if(is.null(colnames(X))==TRUE)
  {
    colnames(X)<-paste0("X",0:(q1-1))
  }
  if(is.null(colnames(W))==TRUE)
  {
    colnames(W)<-paste0("W",0:(q2-1))
  }
  
  ncluster<-ncol(eta)
  
  if(boot)
  {
    beta.boot<-array(NA,c(q1,ncluster,sims))
    dimnames(beta.boot)[[1]]<-colnames(X)
    alpha.boot<-array(NA,c(q2,ncluster,sims))
    dimnames(alpha.boot)[[1]]<-colnames(W)
    dimnames(beta.boot)[[2]]=dimnames(alpha.boot)[[2]]<-paste0("cluster",1:ncluster)
    
    for(b in 1:sims)
    {
      set.seed(100+b)
      idx.tmp<-sample(1:n,n,replace=TRUE)
      
      Ytmp<-Y[idx.tmp]
      Xtmp<-X[idx.tmp,]
      Wtmp<-W[idx.tmp,]
      eta.tmp<-eta[idx.tmp,]
      
      re.tmp<-NULL
      try(re.tmp<-capPCL_coef_class(Y=Ytmp,X=Xtmp,W=Wtmp,gamma=gamma,eta=eta.tmp,max.itr=max.itr,tol=tol,trace=FALSE,score.return=FALSE))
      
      if(is.null(re.tmp)==FALSE)
      {
        if(re.tmp$convergence)
        {
          beta.boot[,,b]<-re.tmp$beta
          alpha.boot[,,b]<-re.tmp$alpha
        }
      }
      
      if(verbose)
      {
        print(paste0("Bootstrap sample ",b))
      }
    }
    
    # dealing with outliers and missing values (non-convergence sample)
    for(k in 1:ncluster)
    {
      for(j in 1:q1)
      {
        if(sum(is.na(beta.boot[j,k,]))>0)
        {
          itmp.nna<-which(is.na(beta.boot[j,k,])==FALSE)
          dis.cook<-cooks.distance(lm(beta.boot[j,k,]~1))
          cook.thred<-4/((length(itmp.nna)-2-2))
          itmp<-itmp.nna[which(dis.cook>cook.thred)]
          beta.boot[j,k,itmp]<-NA
          
          # print(length(itmp))
          
          while(length(itmp)>0)
          {
            itmp.nna<-which(is.na(beta.boot[j,k,])==FALSE)
            dis.cook<-cooks.distance(lm(beta.boot[j,k,]~1))
            cook.thred<-4/((length(itmp.nna)-2-2))
            itmp<-itmp.nna[which(dis.cook>cook.thred)]
            beta.boot[j,k,itmp]<-NA
            
            # print(length(itmp))
          }
        }
      }
      
      for(j in 1:q2)
      {
        if(sum(is.na(alpha.boot[j,k,]))>0)
        {
          itmp.nna<-which(is.na(alpha.boot[j,k,])==FALSE)
          dis.cook<-cooks.distance(lm(alpha.boot[j,k,]~1))
          cook.thred<-4/((length(itmp.nna)-2-2))
          itmp<-itmp.nna[which(dis.cook>cook.thred)]
          alpha.boot[j,k,itmp]<-NA
          
          # print(length(itmp))
          
          while(length(itmp)>0)
          {
            itmp.nna<-which(is.na(alpha.boot[j,k,])==FALSE)
            dis.cook<-cooks.distance(lm(alpha.boot[j,k,]~1))
            cook.thred<-4/((length(itmp.nna)-2-2))
            itmp<-itmp.nna[which(dis.cook>cook.thred)]
            alpha.boot[j,k,itmp]<-NA
            
            # print(length(itmp))
          }
        }
      }
    }
    
    beta.out<-array(NA,c(q1,6,ncluster))
    dimnames(beta.out)[[1]]<-colnames(X)
    alpha.out<-array(NA,c(q2,6,ncluster))
    dimnames(alpha.out)[[1]]<-colnames(W)
    dimnames(beta.out)[[2]]=dimnames(alpha.out)[[2]]<-c("Estimate","SE","statistics","pvalue","LB","UB")
    dimnames(beta.out)[[3]]=dimnames(alpha.out)[[3]]<-paste0("cluster",1:ncluster)
    
    beta.out[,1,]<-apply(beta.boot,c(1,2),mean,na.rm=TRUE)
    beta.out[,2,]<-apply(beta.boot,c(1,2),sd,na.rm=TRUE)
    beta.out[,3,]<-beta.out[,1,]/beta.out[,2,]
    beta.out[,4,]<-(1-pnorm(abs(beta.out[,3,])))*2
    
    alpha.out[,1,]<-apply(alpha.boot,c(1,2),mean,na.rm=TRUE)
    alpha.out[,2,]<-apply(alpha.boot,c(1,2),sd,na.rm=TRUE)
    alpha.out[,3,]<-alpha.out[,1,]/alpha.out[,2,]
    alpha.out[,3,1]<-0
    alpha.out[,4,]<-(1-pnorm(abs(alpha.out[,3,])))*2
    
    if(boot.ci.type[1]=="se")
    {
      zv<-qnorm(1-(1-conf.level)/2)
      
      beta.out[,5,]<-beta.out[,1,]-zv*beta.out[,2,]
      beta.out[,6,]<-beta.out[,1,]+zv*beta.out[,2,]
      
      alpha.out[,5,]<-alpha.out[,1,]-zv*alpha.out[,2,]
      alpha.out[,6,]<-alpha.out[,1,]+zv*alpha.out[,2,]
    }
    if(boot.ci.type[1]=="perc")
    {
      beta.out[,5,]<-apply(beta.boot,c(1,2),quantile,probs=(1-conf.level)/2,na.rm=TRUE)
      beta.out[,6,]<-apply(beta.boot,c(1,2),quantile,probs=1-(1-conf.level)/2,na.rm=TRUE)
      
      alpha.out[,5,]<-apply(alpha.boot,c(1,2),quantile,probs=(1-conf.level)/2,na.rm=TRUE)
      alpha.out[,6,]<-apply(alpha.boot,c(1,2),quantile,probs=1-(1-conf.level)/2,na.rm=TRUE)
    }
    
    re<-list(beta=beta.out,alpha=alpha.out,beta.boot=beta.boot,alpha.boot=alpha.boot)
    
    return(re)
  }else
  {
    stop("Error!")
  }
}
############################








