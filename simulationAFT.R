rm(list=ls())

library('ncvreg')
library('fmrs')
library('survival')

source("~/aft_functions.R")

for (case in 1:36){
  
  repli<-50
  case1<-matrix(0,repli*5+5,12)
  colnames(case1)<-c('tp','df','alpha.mse','alpha.tp','alpha.fp','alpha.tpr','alpha.fpr',
                     'beta.mse','beta.tp','beta.fp','beta.tpr','beta.fpr')
  repli.out<-matrix(0,5,13)
  colnames(repli.out)<-c('method',colnames(case1))
  #model specifications
  
  n <- 1000
  p <- 50
  q <- 50
  theta <- matrix(0,30,30)
  
  if (case%%6 == 1) {
    rho = 0.7
    for (i in 1:3){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)] = matrix_vec(10,cov.str = "ar",rho = 0.7)
    }
    
  } else if (case%%6 == 2) {
    ### same as case 1 except for that rho=0.5
    rho = 0.5
    for (i in 1:3){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)] = matrix_vec(10,cov.str = "ar",rho = 0.7)
    }
    
  } else if (case%%6 == 3) {
    ### same as case 1 except for that rho=0.3
    rho = 0.3
    for (i in 1:3){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)] = matrix_vec(10,cov.str = "ar",rho = 0.7)
    }
    
    
  } else if (case%%6 == 4){
    ### same as case 1 except for that each block of theta is AR(0.7).
    rho=0.7
    for (i in 1:3){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)] = matrix_vec(10,cov.str = "ar",rho = 0.3)
    }
    
  } else if (case%%6 == 5){
    rho=0.5
    for (i in 1:3){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)] = matrix_vec(10,cov.str = "ar",rho = 0.3)
    }
    
  } else if (case%%6 == 0){
    rho=0.3
    for (i in 1:3){
      theta[((i-1)*10+1):(i*10),((i-1)*10+1):(i*10)] = matrix_vec(10,cov.str = "ar",rho = 0.3)
    }
    
  }
  
  if (case<=6 || (case<=24 && case>18)){
    # two clusters have the same important markers but with different magnitudes
    cc<-11:40
    alpha.set <- cbind(c(rep(1.5,5),rep(0.5,5),rep(0,p-10)),c(rep(-0.5,10),rep(0,p-10)))
    beta.set<-cbind(c(rep(1.5,2),rep(0.5,8),rep(0,q-10)),c(rep(-0.5,10),rep(0,q-10)))
  } else if ((case<=12 && case>6) || (case<=30 && case>24)) {
    cc<-11:40
    # two clusters have the different important markers 
    alpha.set <- cbind(c(rep(1.5,5),rep(0.5,5),rep(0,p-10)),c(rep(0,3),rep(-1.5,10),rep(0,p-13)))
    beta.set<-cbind(c(rep(1.5,2),rep(0.5,8),rep(0,q-10)),c(rep(0,3),rep(-1.5,10),rep(0,q-13)))
  } else if ((case<=18 && case>12)|| (case<=36 && case>30)) {
    cc<-11:40
    # two clusters have some important markers with the same magnitudes
    alpha.set <- cbind(c(rep(1.5,5),rep(0.5,5),rep(0,p-10)),c(rep(-0.5,5),rep(0.5,2),rep(-0.5,3),rep(0,p-10)))
    beta.set<-cbind(c(rep(1.5,2),rep(0.5,8),rep(0,q-10)),c(rep(-0.5,5),rep(0.5,2),rep(-0.5,3),rep(0,q-10)))
  }
  
  if (case<=18){
    prob<-c(0.5,0.5)
  } else {
    prob<-c(0.4,0.6)
  }
  
  
  ssd <- c(0.5,0.5)
  prob <- c(0.5,0.5)
  K <- length(prob)
  
  
  ##################################
  for (re in 1:repli) {
    print(paste0('The ',re,'th replicate'))
    
    x <- genXZ(n,p,rho)
    z <- genXZ(n,q,rho)
    z[,cc]=x[,cc]%*%theta
    
    #generate data
    dat <- genY(x,z, alpha.set,beta.set,K,prob,ssd,umax = 100)
    y <- as.matrix(dat$y)
    Delta <- as.matrix(dat$Delta)
    sum(Delta)/n
    class_true <- dat$class_true
    
    maxiter <- 300
    ksi <- 6
    tau <- 1e-4
    iter_v <- 1e-4
    init_time <- 20
    
    method_result<-list()
    init_set<-list()
    result<-list()
    
    method <- c('proprosed','AFT_MCP','AFT-X','AFT-Z','MCP-XZ')
    
    lambda1_set <- c(seq(from=0.01,to=0.1,by=0.01))
    lambda2_set <-c(seq(from=0.001,to=0.01,by=0.001))
    l1=length(lambda1_set)
    l2=length(lambda2_set)
    
    
    for (tt in 1:5){
      
      if (tt==1){
        
        # init_class selection
        for (m in 1:init_time){
          # print(paste0('The ',m,'th initialization'))
          if (m==1){
            class_old<-sample(1:2,n,replace=T)
            init_set[[1]]<-class_old
          } else {
            aa=0
            while (aa<0.5*n){
              for (ii in 1:(m-1)){
                class_old<-sample(1:2,n,replace=T)
                aa=sum(abs(class_old-init_set[[ii]]))
                if (aa>0.5*n){
                  break
                }
              }
            }
          }
          
          init_set[[m]]<-class_old
        }
        
        lambda1=lambda1_set[ceiling(l1/2)]
        lambda2=lambda2_set[ceiling(l2/2)]
        bic_init<-matrix(1e+10,init_time,1)
        for (m in 1:init_time){
          print(paste0('The ',m,'th initialization'))
          class_old=init_set[[m]]
          temp<-aft_corr(x,z,y,Delta,K,lambda1,lambda2,class_old,ksi=ksi,tau=tau,maxiter=maxiter,iter_v=iter_v)
          delta=temp$delta
          alpha=temp$alpha
          beta=temp$beta
          mu=temp$mu
          sigma=temp$sigma
          lik<-L(y,Delta,x,z,alpha,beta,mu,sigma)
          bic_init[m] <- (-2)*lik+(sum(beta!=0)+sum(alpha!=0))*log(n)/n
        }
        id=which.min(bic_init)
        class_old=init_set[[id]]
        
        # lambda selection
        bic=matrix(1e+10,l1*l2,1)
        result_temp<-list()
        lll=1
        for (m1 in 1:l1){
          print(paste0('The ',m1,'th lambda'))
          lambda1=lambda1_set[m1]
          for (m2 in 1:l2){
            lambda2=lambda2_set[m2]
            temp<-aft_corr(x,z,y,Delta,K,lambda1,lambda2,class_old,ksi=ksi,tau=tau,maxiter=maxiter,iter_v=iter_v)
            result_temp[[lll]]<-temp
            delta=temp$delta
            alpha=temp$alpha
            beta=temp$beta
            mu=temp$mu
            sigma=temp$sigma
            lik<-L(y,Delta,x,z,alpha,beta,mu,sigma)
            bic[lll] <- (-2)*lik+(sum(beta!=0)+sum(alpha!=0))*log(n)/n
            lll=lll+1
          }
        }
        id=which.min(bic)
        result[[tt]]<-result_temp[[id]]
        
        delta = result[[tt]]$delta
        class_hat <- rep(0,n)
        for (i in 1:n) {
          if(delta[i,1]>delta[i,2]){class_hat[i]=1}else{class_hat[i]=2}
        }
        
        result[[tt]]$class_hat=class_hat
        
      } else if (tt==2) {
        
        res.mle <- fmrs.mle(y=y,x=cbind(x,z),delta=Delta,nComp=K,disFamily="lnorm")
        res.lam <- fmrs.tunsel(y=y,x=cbind(x,z),delta=Delta,
                               nComp=ncomp(res.mle), disFamily="lnorm",
                               initCoeff=c(coefficients(res.mle)),
                               initDispersion = dispersion(res.mle),
                               initmixProp = mixProp(res.mle),
                               penFamily = "adplasso")
        res.var <- fmrs.varsel(y=y,x=cbind(x,z),delta=Delta,
                               nComp = ncomp(res.mle), disFamily = "lnorm",
                               initCoeff=c(coefficients(res.mle)),
                               initDispersion = dispersion(res.mle),
                               initmixProp = mixProp(res.mle),
                               penFamily = "adplasso",
                               lambPen = slot(res.lam, "lambPen"))
        result[[tt]] <- list()
        
        coeffs <- res.var@coefficients
        coeffs[abs(coeffs)<1e-3] <- 0
        
        alpha <- coeffs[2:(p+1),]
        beta <- coeffs[(p+2):(p+q+1),]
        result[[tt]]$alpha <- alpha
        result[[tt]]$beta <- beta
        
        delta = res.var@weights
        class_hat <- rep(0,n)
        for (i in 1:n) {
          if(delta[i,1]>delta[i,2]){class_hat[i]=1}else{class_hat[i]=2}
        }
        
        result[[tt]]$class_hat=class_hat
        
      } else if (tt==3) {
        
        res.mle <- fmrs.mle(y=y,x=x,delta=Delta,nComp=K,disFamily="lnorm")
        res.lam <- fmrs.tunsel(y=y,x=x,delta=Delta,
                               nComp=ncomp(res.mle), disFamily="lnorm",
                               initCoeff=c(coefficients(res.mle)),
                               initDispersion = dispersion(res.mle),
                               initmixProp = mixProp(res.mle),
                               penFamily = "adplasso")
        res.var <- fmrs.varsel(y=y,x=x,delta=Delta,
                               nComp = ncomp(res.mle), disFamily = "lnorm",
                               initCoeff=c(coefficients(res.mle)),
                               initDispersion = dispersion(res.mle),
                               initmixProp = mixProp(res.mle),
                               penFamily = "adplasso",
                               lambPen = slot(res.lam, "lambPen"))
        result[[tt]] <- list()
        
        coeffs <- res.var@coefficients
        coeffs[abs(coeffs)<1e-3] <- 0
        
        alpha <- coeffs[2:(p+1),]
        beta <- matrix(0,q,K) #coeffs[(p+2):(p+q+1),]
        result[[tt]]$alpha <- alpha
        result[[tt]]$beta <- beta
        
        delta = res.var@weights
        class_hat <- rep(0,n)
        for (i in 1:n) {
          if(delta[i,1]>delta[i,2]){class_hat[i]=1}else{class_hat[i]=2}
        }
        
        result[[tt]]$class_hat=class_hat
        
      } else if (tt==4) {
        
        res.mle <- fmrs.mle(y=y,x=z,delta=Delta,nComp=K,disFamily="lnorm")
        res.lam <- fmrs.tunsel(y=y,x=z,delta=Delta,
                               nComp=ncomp(res.mle), disFamily="lnorm",
                               initCoeff=c(coefficients(res.mle)),
                               initDispersion = dispersion(res.mle),
                               initmixProp = mixProp(res.mle),
                               penFamily = "adplasso")
        res.var <- fmrs.varsel(y=y,x=z,delta=Delta,
                               nComp = ncomp(res.mle), disFamily = "lnorm",
                               initCoeff=c(coefficients(res.mle)),
                               initDispersion = dispersion(res.mle),
                               initmixProp = mixProp(res.mle),
                               penFamily = "adplasso",
                               lambPen = slot(res.lam, "lambPen"))
        result[[tt]] <- list()
        
        coeffs <- res.var@coefficients
        coeffs[abs(coeffs)<1e-3] <- 0
        
        alpha <- matrix(0,p,K)
        beta <- coeffs[2:(q+1),]
        result[[tt]]$alpha <- alpha
        result[[tt]]$beta <- beta
        
        delta = res.var@weights
        class_hat <- rep(0,n)
        for (i in 1:n) {
          if(delta[i,1]>delta[i,2]){class_hat[i]=1}else{class_hat[i]=2}
        }
        
        result[[tt]]$class_hat=class_hat

      } else if (tt==5) {
        
        X <- cbind(x,z)
        fit <- ncvsurv(X,Surv(y,Delta),penalty = "MCP", lambda = 0.05)
        
        #id <- which.min(fit3$lambda)
        temp1 = fit$beta
        
        
        alpha = matrix(0,p,K)
        beta = matrix(0,q,K)
        
        alpha[,1] = temp1[1:p]
        alpha[,2] = alpha[,1]
        
        beta[,1] = temp1[(p+1):(p+q)]
        beta[,2] = beta[,1]
        
        class_hat <- rep(1,n)
        
        result[[tt]]<-list(alpha=alpha,beta=beta,class_hat=class_hat)
        
      }
      
      alpha <- result[[tt]]$alpha
      beta <- result[[tt]]$beta
      class_hat <- result[[tt]]$class_hat
      
      ##True positive rate of delta
      indicator <- rep(0,n)
      for (i in 1:n) {
        if(class_hat[i]==dat$class_true[i]){indicator[i]=1}else{indicator[i]=0}
      }
      
      if (sum(indicator)/n < 1-sum(indicator)/n){
        alpha = alpha[,c(2,1)]
        beta = beta[,c(2,1)]
      }
      
      df<-sum(beta!=0)+sum(alpha!=0)
      tp<-max(sum(indicator)/n, 1-sum(indicator)/n)
      alpha.mse<-sqrt(sum((alpha-alpha.set)^2)) ####################
      beta.mse<-sqrt(sum((beta-beta.set)^2)) #############################
      alpha.tp<-GetFPTP(as.vector(alpha.set),as.vector(alpha))
      beta.tp<-GetFPTP(as.vector(beta.set),as.vector(beta))
      
      case1[(tt-1)*repli+re,]<-c(tp,df,alpha.mse,alpha.tp$TP,alpha.tp$FP,alpha.tp$TPR,alpha.tp$FPR,
                                 beta.mse,beta.tp$TP,beta.tp$FP,beta.tp$TPR,beta.tp$FPR)
      
    }
    print('---------------------------------')
    for (tt in 1:5){
      print(c(method[tt],round(case1[(tt-1)*repli+re,],2)),quote = FALSE)
      repli.out[tt,]<-c(method[tt],round(case1[(tt-1)*repli+re,],2))
    }
    # write.csv(repli.out,file =  paste0("repli",re,'case',case,".csv") )
    print('---------------------------------')
  }
  
  
  case1[5*repli+1,]=colMeans(case1[1:repli,])
  case1[5*repli+2,]=colMeans(case1[(repli+1):(2*repli),])
  case1[5*repli+3,]=colMeans(case1[(2*repli+1):(3*repli),])
  case1[5*repli+4,]=colMeans(case1[(3*repli+1):(4*repli),])
  case1[5*repli+5,]=colMeans(case1[(4*repli+1):(5*repli),])
  
  print('---------------------------------')
  print('FINAL RESULTS')
  for (tt in 1:5){
    print(c(method[tt],round(case1[5*repli+tt,],2)),quote = FALSE)
    
  }
  
  write.csv(case1,file=paste0('result_case',case,'.csv'))
}
