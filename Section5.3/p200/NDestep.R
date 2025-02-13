###
# x<-X.target
# y<-y.target
# betaw<-beta_coef
# esp1=1e-4
# T1=100
# T2=100
# t1=1
# t2=1
# lambda=0.1
  

  DeStep_feng<-function(x,y,init=NULL,betaw,lambda,esp1=1e-3,T1=100,T2=100){
    n<-nrow(x)
    K <- length(unique(y))
    p <- ncol(X)
    betaw<-betaw[-1,]
    init<-matrix(0, nrow = p, ncol = K-1)
    newdelta<-init
    ##
    convergence=0
    for (t1 in 1:T1) {
      
      Index<-newdelta
      
      for (k in 1:(K-1)) {
        for (t2 in 1:T2) {
          
          # eta<-exp(x%*%(betaw+newdelta[,-K]))
          # probs <- eta /(1+rowSums(eta))
          
          probs <- exp(x%*%(betaw+newdelta[,-K])) /(1+rowSums(exp(x%*%(betaw+newdelta[,-K]))))
          id=c(y==k)
          yk<-rep(0,n)
          yk[id]<-1
          
          index<-newdelta[,k]
          for (j in 1:p) {
            
            r<-yk-probs[,k]
            v<-probs[,k]
            #v<-probs[,k]*(1-probs[,k])
            #v<-probs[,k]*(1-probs[,k])
            SS<-sign((r%*%x[,j]+(v%*% x[,j]^2)*newdelta[,k][j]))*(max(abs((r%*%x[,j]+(v%*% x[,j]^2)*newdelta[,k][j]))/n-lambda,0))
            newdelta[,k][j]<-SS/( v%*% x[,j]^2 /n)
            
          } #end j
          
          if(max(abs(index-newdelta[,k]))<esp1){
            break
          }
          # if(sum((index-newdelta[,k])^2)<esp2){
          #   break
          # }
          
        } #end T2
        
      }   #end K
      
      if(max(abs(Index-newdelta))<esp1){
        convergence=1
        break
      }
      # if(sum((Index[,-K]-newdelta[,-K])^2)<esp1){
      #   convergence=1
      #   break
      # }
      
    }#end T1
    
    result<-list()
    result[[1]]<-newdelta
    result[[2]]<-convergence
    result[[3]]<-max(abs(Index-newdelta))
    result[[4]]<-c(t1,t2)
    return(result)
    
    
  }#end function
  
  
#####################################################
  
  nloglike<-function(X,Y,beta){
    n<-nrow(X)
    K <- length(unique(y))
    p <- ncol(X)
    L<-0
    for (k in 1:(K-1)) {
      id=c(Y==k)
      yk<-rep(0,n)
      yk[id]<-1
      L<-L-yk%*%(cbind(1,X)%*%beta[,k])  
    }
    L<-L+ sum( log( 1+rowSums( exp(cbind(1,X)%*%(beta)) ) ) )
    L<-L/n
    return(L)
  }
  
  
  
#####################################################
  cv.Destep<-function(x,y,init=NULL,betaw,lambda,kfold=3,lambdaseq=seq(0.01,0.3,length.out=10),esp1=1e-3,T1=100,T2=100){
    losses<-matrix(NA,length(lambdaseq),kfold)
    N<-nrow(x)
    for(lamb in 1:length(lambdaseq)){
      lam<-lambdaseq[lamb]
      loci<-createFolds(1:N,k=kfold)
      for(j in 1:kfold){
        train.x<-x[-as.vector(loci[[j]]),]
        train.y<-y[-as.vector(loci[[j]])]
        test.x<-x[as.vector(loci[[j]]),]
        test.y<-y[as.vector(loci[[j]])]
        
        deresult<-DeStep_feng(x=train.x,y=train.y,betaw = betaw,lambda = lam,T1=50,T2=50)
        delta<-deresult[[1]]
        
        losses[lamb,j]<-nloglike(X=test.x,Y=test.y,beta=delta+betaw[-1,])
        
      }
    }
    lam.min<-lambdaseq[which.min(apply(losses,1,mean))]
    return(list(losses,lam.min))
}
  
  
  

######################With intercept############
  
  
  DeStep_feng2<-function(x,y,init=NULL,betaw,lambda,esp1=1e-3,T1=100,T2=100){
    n<-nrow(x)
    K <- length(unique(y))
    p <- ncol(X)
    Betaw<-betaw
    
    init<-matrix(0, nrow = p+1, ncol = K)
    
    newdelta<-init
    ##
    convergence=0
    for (t1 in 1:T1) {
      
      Index<-newdelta
      
      for (k in 1:(K-1)) {
        for (t2 in 1:T2) {
          
          # eta<-exp(x%*%(betaw+newdelta[,-K]))
          # probs <- eta /(1+rowSums(eta))
          
          probs <- exp(cbind(1,x)%*%(Betaw+newdelta[,-K])) /(1+rowSums(exp(cbind(1,x)%*%(Betaw+newdelta[,-K]))))
          id=c(y==k)
          yk<-rep(0,n)
          yk[id]<-1
          
          index<-newdelta[,k]
          for (j in 1:(p+1)) {
            
            r<-yk-probs[,k]
            v<-probs[,k]
            #v<-probs[,k]*(1-probs[,k])
            #v<-probs[,k]*(1-probs[,k])
            if(j==1){
              SS<-sign( mean(r)+mean(v*newdelta[,k][j]) )*( max(abs(mean(r)+mean(v*newdelta[,k][j]))/n -lambda,0))
              newdelta[,k][j]<-SS/mean(v)
            }else{
            SS<-sign((r%*%x[,j-1]+(v%*% x[,j-1]^2)*newdelta[,k][j]))*(max(abs((r%*%x[,j-1]+(v%*% x[,j-1]^2)*newdelta[,k][j]))/n-lambda,0))
            newdelta[,k][j]<-SS/( v%*% x[,j-1]^2 /n)
            }
            
          } #end j
          
          if(max(abs(index-newdelta[,k]))<esp1){
            break
          }
          # if(sum((index-newdelta[,k])^2)<esp2){
          #   break
          # }
          
        } #end T2
        
      }   #end K
      
      if(max(abs(Index-newdelta))<esp1){
        convergence=1
        break
      }
      # if(sum((Index[,-K]-newdelta[,-K])^2)<esp1){
      #   convergence=1
      #   break
      # }
      
    }#end T1
    
    result<-list()
    result[[1]]<-newdelta[,-K]
    result[[2]]<-convergence
    result[[3]]<-max(abs(Index-newdelta))
    result[[4]]<-c(t1,t2)
    return(result)
    
    
  }#end function
  
  ####################################################
  nloglike2<-function(X,Y,beta){
    n<-nrow(X)
    K <- length(unique(y))
    p <- ncol(X)
    L<-0
    for (k in 1:(K-1)) {
      id=c(Y==k)
      yk<-rep(0,n)
      yk[id]<-1
      L<-L-yk%*%(cbind(1,X)%*%beta[,k])  
    }
    L<-L+ sum( log( 1+rowSums( exp(cbind(1,X)%*%(beta)) ) ) )
    L<-L/n
    return(L)
  }
  

  #####################################################
  cv.Destep2<-function(x,y,init=NULL,betaw,lambda,kfold=3,lambdaseq=seq(0.01,0.3,length.out=10),esp1=1e-3,T1=100,T2=100){
    losses<-matrix(NA,length(lambdaseq),kfold)
    N<-nrow(x)
    for(lamb in 1:length(lambdaseq)){
      lam<-lambdaseq[lamb]
      loci<-createFolds(1:N,k=kfold)
      for(j in 1:kfold){
        train.x<-x[-as.vector(loci[[j]]),]
        train.y<-y[-as.vector(loci[[j]])]
        test.x<-x[as.vector(loci[[j]]),]
        test.y<-y[as.vector(loci[[j]])]
        
        deresult<-DeStep_feng2(x=train.x,y=train.y,betaw = betaw,lambda = lam,T1=50,T2=50)
        delta<-deresult[[1]]
        
        losses[lamb,j]<-nloglike2(X=test.x,Y=test.y,beta=delta+betaw)
        
      }
    }
    lam.min<-lambdaseq[which.min(apply(losses,1,mean))]
    return(list(losses,lam.min))
  }
  
  
  
  
  ########################################################
  
  
  
  
  
  
  
  
  
  
  
######################################Detection###########################
  
  
  # X0=X.target
  # Y0=y.target
  # XK=X.source
  # YK=y.source
  # M=M
  # c=0.1
  # fold=3
  
  sdet_feng = function(X0,Y0,XK,YK,M,c,fold=3){
    L0 = rep(0,fold)
    LM = matrix(0,M,fold)
    ntarget=nrow(X0)
    for(r in 1:fold){
      index = ( (r-1)*ntarget/fold+1):(r*ntarget/fold)
      fit.tar <- cv.pemultinom(X0[-index,],Y0[-index],ncores = 1)
      beta_tar <- fit.tar$beta.min
      
      L0[r] = nloglike(X0[index,],Y0[index],beta_tar)
      
      for(m in 1:M){
        fit_m =  cv.pemultinom(rbind(X0[-index,],XK[[m]]),c(Y0[-index],YK[,m]),ncores = 1)
        beta_m = fit_m$beta.min
        LM[m,r] = nloglike(X0[index,],Y0[index],beta_m)
      }
    }
    
    Loss_target = mean(L0)
    Loss_source = apply(LM,1,mean)
    sigma = sqrt(sum((L0-Loss_target)^2)/2)
    threshold = Loss_target + c*max(sigma,0.01)
    id = which(Loss_source <= threshold)
    
    rr<-list(id,Loss_target,Loss_source,sigma)
    return(rr)
  }
  
  
 
  
  
  
######  For  backup
  # DeStep<-function(x,y,init=NULL,betaw,lambda,esp1=1e-4,T1=100,T2=100){
  #   n<-nrow(x)
  #   K <- length(unique(y))
  #   p <- ncol(X)
  #   betaw<-betaw[-1,]
  #   init<-matrix(0, nrow = p, ncol = K)
  #   newdelta<-init
  #   ##
  #   convergence=0
  #   for (t1 in 1:T1) {
  #     
  #     Index<-newdelta
  #     
  #     for (k in 1:(K-1)) {
  #       for (t2 in 1:T2) {
  #         
  #         # eta<-exp(x%*%(betaw+newdelta[,-K]))
  #         # probs <- eta /(1+rowSums(eta))
  #         
  #         probs <- exp(x%*%(betaw+newdelta[,-K])) /(1+rowSums(exp(x%*%(betaw+newdelta[,-K]))))
  #         id=c(y==k)
  #         yk<-rep(0,n)
  #         yk[id]<-1
  #         
  #         index<-newdelta[,k]
  #         for (j in 1:p) {
  #           
  #           r<-yk-probs[,k]
  #           #v<-probs[,k]*((1+rowSums(exp(x%*%(betaw+newdelta[,-K]))))-probs[,k])
  #           v<-probs[,k]
  #           #v<-probs[,k]*(1-probs[,k])
  #           newdelta[,k][j]<-newdelta[,k][j]+(r%*%x[,j]+lambda*sign(newdelta[,k][j])*n)/( v%*% x[,j]^2 )
  #           
  #         } #end j
  #         
  #         if(max(abs(index-newdelta[,k]))<esp1){
  #           break
  #         }
  #         # if(sum((index-newdelta[,k])^2)<esp2){
  #         #   break
  #         # }
  #         
  #       } #end T2
  #       
  #     }   #end K
  #     
  #     if(max(abs(Index[,-K]-newdelta[,-K]))<esp1){
  #       convergence=1
  #       break
  #     }
  #     # if(sum((Index[,-K]-newdelta[,-K])^2)<esp1){
  #     #   convergence=1
  #     #   break
  #     # }
  #     
  #   }#end T1
  #   
  #   result<-list()
  #   result[[1]]<-newdelta
  #   result[[2]]<-convergence
  #   result[[3]]<-max(abs(Index[,-K]-newdelta[,-K]))
  #   result[[4]]<-c(t1,t2)
  #   return(result)
  #   
  #   
  # }#end function
  
  # DeStep2<-function(x,y,init=NULL,betaw,lambda,esp1=1e-4,T1=100){
  #   n<-nrow(x)
  #   K <- length(unique(y))
  #   p <- ncol(X)
  #   betaw<-betaw[-1,]
  #   init<-matrix(0, nrow = p, ncol = K)
  #   newdelta<-init
  #   ##
  #   convergence=0
  #   for (t1 in 1:T1) {
  #     
  #     Index<-newdelta
  #     
  #     for (k in 1:(K-1)) {
  #       
  #         eta<-exp(x%*%(betaw+newdelta[,-K]))
  #         probs <- eta /(1+rowSums(eta))
  #         id=c(y==k)
  #         yk<-rep(0,n)
  #         yk[id]<-1
  #         
  #         for (j in 1:p) {
  #           
  #           r<-yk-probs[,k]
  #           v<-probs[,k]*((1+rowSums(exp(x%*%(betaw+newdelta[,-K]))))-probs[,k])
  #           newdelta[,k][j]<-newdelta[,k][j]+(r%*%x[,j]+lambda*sign(newdelta[,k][j]))/( v%*% x[,j]^2 )
  #           
  #         } #end j
  #         
  #       
  #     }   #end K
  #     
  #     if(max(abs(Index[,-K]-newdelta[,-K]))<esp1){
  #       convergence=1
  #       break
  #     }
  #     # if(sum((Index[,-K]-newdelta[,-K])^2)<esp1){
  #     #   convergence=1
  #     #   break
  #     # }
  #     
  #   }#end T1
  #   
  #   result<-list()
  #   result[[1]]<-newdelta
  #   result[[2]]<-convergence
  #   result[[3]]<-max(abs(Index[,-K]-newdelta[,-K]))
  #   result[[4]]<-c(t1)
  #   return(result)
  #   
  #   
  # }#end function