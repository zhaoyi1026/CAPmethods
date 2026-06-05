#################################################
# rvmf function in Rfast package
#################################################

rvmf_h<-function(n,ca,d1,x0,m,k,b)
{
  ta=u=z=tmp<-0
  w<-rep(NA,n)
  
  for(i in 1:n)
  {
    ta<--1000
    u<-1
    while(ta-ca<log(u))
    {
      z<-rbeta(1,m,m)
      u<-runif(1,0,1)
      tmp<-(1-(1+b)*z)/(1-(1-b)*z)
      ta<-k*tmp+d1*log(1-x0*tmp)
    }
    w[i]<-tmp
  }
  
  return(w)
}

rvmf<-function (n, mu, k) 
{
  rotation <- function(a, b) {
    p <- length(a)
    ab <- sum(a * b)
    ca <- a - b * ab
    ca <- ca/sqrt(sum(ca^2))
    A <- b %*% t(ca)
    A <- A - t(A)
    theta <- acos(ab)
    diag(p) + sin(theta) * A + (cos(theta) - 1) * (b %*% 
                                                     t(b) + ca %*% t(ca))
  }
  d <- length(mu)
  if (k > 0) {
    mu <- mu/sqrt(sum(mu^2))
    ini <- c(numeric(d - 1), 1)
    d1 <- d - 1
    
    # v1 <- Rfast::matrnorm(n, d1)
    # v <- v1/sqrt(Rfast::rowsums(v1^2))
    v1<-matrix(rnorm(n*d1),nrow=n)
    v<-t(apply(v1,1,function(x){return(x/sqrt(sum(x^2)))}))
    
    b <- (-2 * k + sqrt(4 * k^2 + d1^2))/d1
    x0 <- (1 - b)/(1 + b)
    m <- 0.5 * d1
    ca <- k * x0 + (d - 1) * log(1 - x0^2)
    
    # w <- .Call("Rfast_rvmf_h", PACKAGE = "Rfast", n, ca, 
    #            d1, x0, m, k, b)
    w<-rvmf_h(n,ca,d1,x0,m,k,b)
    
    S <- cbind(sqrt(1 - w^2) * v, w)
    if (isTRUE(all.equal(ini, mu, check.attributes = FALSE))) {
      x <- S
    }else 
      if (isTRUE(all.equal(-ini, mu, check.attributes = FALSE))) {
      x <- -S
      }else 
        {
          A <- rotation(ini, mu)
          x <- tcrossprod(S, A)
        }
  }
  else {
    # x1 <- Rfast::matrnorm(n, d)
    # x <- x1/sqrt(Rfast::rowsums(x1^2))
    x1<-matrix(rnorm(n*d),nrow=n)
    x<-t(apply(x1,1,function(x){return(x/sqrt(sum(x^2)))}))
  }
  colnames(x) <- names(mu)
  
  return(x)
}

