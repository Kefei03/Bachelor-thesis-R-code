
S <- function(sigma,delta,gamma,x,y,eta,alpha) {
  
  s=1/n*sum(delta*(sigma*y-gamma)*x)+eta*alpha
  
}


genXZ <- function(n,p,rho) {
  mu <- rep(0, p)
  Sigma <- diag(p)
  
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j] <- rho^abs(i-j)
    }}
  
  X <- matrix(rnorm(p*n), n)
  X <- scale(X, TRUE, FALSE)
  X <- X %*% svd(X, nu = 0)$v
  X <- scale(X, FALSE, TRUE)
  
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))){
    nm <- dn[[1L]]}
  dimnames(X) <- list(nm, NULL)
  if (n == 1){
    cX = drop(X)
  } else {cX = t(X)}
  
  cX <- scale(cX)
  colnames(cX) <- paste("X", 1:p,sep = ".")
  
  output <- cX
}

genY <- function(x,z,alpha,beta,K,prob,ssd,umax){
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  q <- dim(z)[2]
  
  coef1 <- matrix(alpha, nrow = K, ncol = p, byrow = TRUE)
  coef2 <- matrix(beta, nrow = K, ncol = p, byrow = TRUE)
  
  prob0 <- cumsum(prob)
  
  yobs <-c()
  c <- rep()
  dlt <- c()
  u <- c()
  tobs <- c()
  
  for(i in 1:n){
    epss <- rnorm(1, 0, 1)
    u1 <- runif(1)
    k = length(which(prob0<=u1)) + 1
    u[i] = k
    tobs[i] <- coef1[k,] %*% x[i,] + coef2[k,] %*% z[i,] + ssd[k] * epss
    
    c[i] <- log(runif(1, 0, umax))
    yobs[i] <- exp(min(tobs[i],c[i]))
    dlt[i] <- (tobs[i] < c[i])*1
  }
  
  output <- list(y=yobs,Delta=dlt,class_true=u)
}


L<-function(y,Delta,x,z,alpha,beta,mu,mysigma){
  
  logy <- log(y)
  
  b <- 0
  for (i in 1:n) {
    a <- 0
    for (k in 1:K) {
      temp <- logy[i] - x[i,k]*alpha[k] - z[i,k]*beta[k]
      f <- dnorm(temp, 0, sd = mysigma[k])
      S <- 1 - pnorm(temp, 0, sd = mysigma[k])
      
      a = a + mu[k] * f^Delta[i] * S^(1 - Delta[i])
    }
    b <- b + log(a+1e-20)
  }
  return(b/n)
}


GetFPTP<-function(theta,theta_hat){
  # to get TNR (True Negative Rate ) and TPR (True Positive Rate) 
  thea = abs(theta) > 0   # transform coefficients to binary values
  thea_hat = abs(theta_hat) > 1e-8  # convert estimated coefficients to binary values
  A = sum((!thea)*(!thea_hat))  # A: TN
  B = sum((!thea)*thea_hat)   # B: FP
  C = sum(thea*(!thea_hat))   # C: FN
  D = sum(thea*thea_hat)    # D: TP
  TPR = D/(D+C)    # TPR=TP/(TP+FN)  true positive rate (TPR) sensitivity
  FPR = B/(B+A)    # FPR=FP/(TN+FP)  false positive rate     
  result=list(TPR= TPR, FPR = FPR, TP=D, FP=B)
  return(result)
}


aft_corr<-function(x,z,y,Delta,K,lambda1,lambda2,class_old,ksi,tau,maxiter,iter_v){
  
  n<-dim(x)[1]
  p<-dim(x)[2]
  q<-dim(z)[2]
  
  # y <- y+0.01
  logy <- log(y)
  
  
  delta <- matrix(0,n,K)
  for (i in 1:n) {
    delta[i,class_old[i]] = 1
  }
  
  
  miu.old <- colSums(delta) / n
  miu.new <- miu.old
  
  
  #alpha/beta init (method: MCP)
  #sigma init (method: MLE)
  
  alpha.old <- matrix(0,p,K)
  beta.old <- matrix(0,q,K)
  sigma <- matrix(1,K,1)
  
  for (k in 1:K){
    model <- ncvsurv(cbind(x[class_old==k,],z[class_old==k,]),Surv(y[class_old==k], Delta[class_old==k]),
                     penalty = "MCP",lambda=0.08)

    temp <- model$beta

    alpha.old[,k]<-as.matrix(temp[1:(p)],p,1)
    beta.old[,k]<-as.matrix(temp[(p+1):(p+q)],q,1)
    # sigma[k] <- model$icoef[2]
  }

  alpha.old[is.na(alpha.old)] <- 0
  beta.old[is.na(beta.old)] <- 0

  
  c <- cor(x,z,method = 'pearson')
  c <- abs(c)*(abs(c)>0.15)
  
  
  sigma <- 1/sigma
  for (k in 1:K) {
    alpha.old[,k] <- alpha.old[,k]*sigma[k]
    beta.old[,k] <- beta.old[,k]*sigma[k]
  }
  
  
  gamma <- matrix(0,n,K)
  for (k in 1:K){
    gamma[,k] <- x%*%alpha.old[,k] + z%*%beta.old[,k]
  }
  
  
  alpha <- alpha.old
  beta <- beta.old
  
  diff_v = 1
  t = 1
  
  ################################  
  while ((diff_v>iter_v) && (t<maxiter)){
    
    alpha.old <- alpha
    beta.old <- beta
    miu.old <- miu.new
    
    
    ############# start E-step ####################
    
    
    #delta new
    ss <- matrix(1,n,K)
    tt <- matrix(1,n,K)
    pp <- matrix(1,n,K)
    for(k in 1:K){
      ss[,k] <- rep(0,n)
      for(l in 1:K){
        tt[,l] <- exp(-1/2*((sigma[l]*logy-gamma[,l])^2-(sigma[k]*logy-gamma[,k])^2))
        pp[,l] <- (1-pnorm(sigma[l]*logy-gamma[,l])) / (1-pnorm(sigma[k]*logy-gamma[,k])+1e-6)
        ss[,k] <- ss[,k] + miu.old[l]*(sigma[l]/sigma[k]) * tt[,l]^Delta * pp[,l]^(1-Delta)
      }
      if(miu.old[k]==0) {delta[,k] = 0} else{
        delta[,k] <- miu.old[k]/ss[,k]
      }
    }
    
    
    ############# end of E-step ####################
    
    #mu new
    miu.new <- colSums(delta) / n
    
    
    #time simulation
    v <- matrix(0,n,K)
    Aik <- matrix(0,n,K)
    omega <- matrix(0,n,K)
    for (k in 1:K) {
      omega[,k] <- sigma[k]*logy - gamma[,k]
      Aik[,k] <- dnorm(omega[,k]) / (1e-6 + 1-pnorm(omega[,k]))
      
      v[,k] <- Delta*logy + (1-Delta) * (x%*%alpha[,k] + z%*%beta[,k] + Aik[,k])/sigma[k]
    }
    
    for (k in 1:K) {
      omega[,k] <- sigma[k]*v[,k] - gamma[,k]
    }
    
    
    ############# start M-step ####################
    
    #sigma new
    a<-rep(1,2)
    b<-rep(1,2)
    c2<-rep(1,2)
    
    for (k in 1:K) {
      a[k]=sum(delta[,k] * v[,k]^2)
      b[k]=sum(delta[,k] * gamma[,k] * v[,k])
      c2[k]=sum(delta[,k])
      sigma[k]= (b[k]+sqrt(b[k]^2+4*a[k]*c2[k]))/(2*a[k])
    }
    
    
    #alpha new
    
    eta <- matrix(0,p,2)
    for (k in 1:K) {
      u <- matrix(0,p,K)
      for (j in 1:p) {
        eta[j,k] = sum(delta[,k]*x[,j]^2)/n
        s=S(sigma = sigma[k],delta[,k],gamma = gamma[,k],x[,j],y=v[,k],eta[j,k],alpha[j,k])
        
        u[j,k]=2/(tau)*exp(-alpha[j,k]^2/tau)*sum(c[j,which(beta[,k]!=0)])
        
        if(abs(s)>lambda1*ksi*(eta[j,k]+lambda2*u[j,k])){
          alpha[j,k]=(s)/(eta[j,k]+lambda2*u[j,k])
        } else {
          if (abs(s)>lambda1){
            alpha[j,k]=(s-sign(s)*lambda1)/(eta[j,k]+lambda2*u[j,k]-1/ksi)
          } else {
            alpha[j,k]=0
          }
        }
        
        alpha[abs(alpha)<1e-3]=0
        gamma[,k]=gamma[,k]+x[,j]*alpha[j,k]-x[,j]*alpha.old[j,k]
      }
      
      
      #beta new
      u <- matrix(0,q,K)
      for (l in 1:q) {
        eta[l,k] = sum(delta[,k]*z[,l]^2)/n
        s=S(sigma = sigma[k],delta[,k],gamma = gamma[,k],z[,l],y=v[,k],eta[l,k],beta[l,k])
        
        u[l,k]=2/(tau)*exp(-beta[l,k]^2/tau)*sum(c[which(alpha[,k]!=0),l])
        
        if(abs(s)>lambda1*ksi*(eta[l,k]+lambda2*u[l,k])){
          beta[l,k]=(s)/(eta[l,k]+lambda2*u[l,k])
        } else {
          if(abs(s)>lambda1){
            beta[l,k]=(s-sign(s)*lambda1)/(eta[l,k]+lambda2*u[l,k]-1/ksi)
          } else {
            beta[l,k]=0
          }
        }
        
        beta[abs(beta)<1e-3]=0
        
        gamma[,k]=gamma[,k]+z[,l]*beta[l,k]-z[,l]*beta.old[l,k]
      }
    }
    
    ############# end of M-step ####################
    
    maxmiu <- max(abs(miu.new-miu.old))
    maxalpha <- max(abs(alpha-alpha.old))
    maxbeta <- max(abs(beta-beta.old))
    
    diff_v <- max(maxmiu, maxalpha, maxbeta)
    t <- t+1
  }
  
  
  for (k in 1:K) {
    alpha[,k]<-alpha[,k]/sigma[k]
    beta[,k]<- beta[,k]/sigma[k]
  }
  
  alpha[abs(alpha)<1e-3] = 0
  beta[abs(beta)<1e-3] = 0
  
  sigma = 1/sigma
  
  output<-list(alpha=alpha,beta=beta,sigma=sigma,delta=delta,mu=miu.new)
}





matrix_vec <- function(P, cov.str, rho){
  tmp <- matrix(0, P, P)
  for (i in 1:P) {
    for (j in 1:P) {
      if (i == j) {
        tmp[i, j] <- 1
      } else if (cov.str == "ar") {
        tmp[i, j] <- rho ^ abs(i - j)
      } else if (cov.str == "flat") {
        tmp[i, j] <- rho
      }
    }
  }
  return(tmp)
}





