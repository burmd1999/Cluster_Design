#All functions used in Extrapolation simulation

#####Calculate MSE's and bias, variance for designs
MSEextrapolation <- function(n,N,X,x0,sigma2,k2,k3,model){
  if (model == 1){
    f <- function(x,k2,k3){1 + x + k2*x^2 + k3*x^3}
  }else if (model == 2){
    f <- function(x,k2,k3){1 + x + x^2 + k2*x^3 + k3*x^4}
  }
  
  #Vector for N extrapolations at x0
  yhat <- vector(mode = "double",length = N)
  if (model == 1){
    betahat <- matrix(NA,ncol = N,nrow = 2)
  }else if (model == 2){
    betahat <- matrix(NA,ncol = N,nrow = 3)
  }
  
  for (j in 1:N){
    #Sample data
    y <- vector(mode="double",length = n)
    for (i in 1:n){
      y[i] <- f(X[i,2],k2,k3) + rnorm(1,0,sd = sigma2)
    }
    #Estimate parameters
    Beta <- solve(t(X)%*%X)%*%t(X)%*%y
    #Estimate extrapolation
    if (model == 1){
      yhat[j] <- Beta[1] + Beta[2]*x0
    }else if (model == 2){
      yhat[j] <- Beta[1] + Beta[2]*x0 + Beta[3]*x0^2
    }
    betahat[,j] <- Beta 
  }
  
  #Extrapolation Bias, Variance and MSE
  if (model == 1){
    #Asymptotic formulas
    bias <- c(1,x0)%*%solve(t(X)%*%X)%*%t(X)%*%f(X[,2],k2,k3) - f(x0,k2,k3)
    variance <- c(1,x0)%*%solve(t(X)%*%X)%*%c(1,x0)*sigma2
    # #Simulated formulas
    # bias <- mean(yhat) - f(x0,k2,k3)
    # variance <- var(yhat)
  }else if (model == 2){
    #Asymptotic formulas
    bias <- c(1,x0,x0^2)%*%solve(t(X)%*%X)%*%t(X)%*%f(X[,2],k2,k3) - f(x0,k2,k3)
    variance <- c(1,x0,x0^2)%*%solve(t(X)%*%X)%*%c(1,x0,x0^2)*sigma2
    # #Simulated formulas
    # bias <- mean(yhat) - f(x0,k2,k3)
    # variance <- var(yhat)
  }
  MSE <- variance + bias^2 
  
  betahat <- rowMeans(betahat)
  
  return(c(bias,variance,MSE,betahat))
}


#####Creates design matrix according to HL design
HLdesign <- function(n,x0,model){
  if (model == 1){
    x02 <- 2*x0
    w <- c((x02-1)/(2*x02),(x02+1)/(2*x02))
    xi <- c(rep(-0.5,n*w[1]),rep(0.5,n*w[2]))
  }else if (model == 2){
    x02 <- 2*x0
    w <- c(x02*(x02-1)/(2*(2*x02^2 - 1)),
           (x02-1)*(x02+1)/(2*x02^2 - 1),
           x02*(x02+1)/(2*(2*x02^2 - 1)))
    xi <- c(rep(-0.5,n*w[1]),rep(0,n*w[2]),rep(0.5,n*w[3]))
  }
  if (model == 1){
    X <- cbind(rep(1,n),xi)
  }else if (model == 2){
    X <- cbind(rep(1,n),xi,xi^2)
  }
  
  return(X)
}

#####Creates design matrix according to WX design
WXdesign <- function(n,x0,m,model){
  if (model == 1){
    if (x0 == 0.75){
      a1 <- 14.06
      a2 <- 3.18
      a3 <- 11.23
      a4 <- 1
      a5 <- -16.52
      a6 <- 0.0032
    }else if (x0 == 2.5){
      a1 <- 148.64
      a2 <- 9.73
      a3 <- 80.38
      a4 <- 1
      a5 <- -525.93
      a6 <- 0.000122
    }
    f <- function(x){
      pmax(((a1*x*2 + a2)*(a3*x*2 + a4) + a5),0)/((a3*x*2 + a4)^2 + a6*(a3*x*2 + a4)^4)
    }
    totint <- integrate(f,-0.5,0.5)$value
    f <- function(x){
      (pmax(((a1*x*2 + a2)*(a3*x*2 + a4) + a5),0)/((a3*x*2 + a4)^2 + a6*(a3*x*2 + a4)^4))/totint
    }
  }else if (model == 2){
    if (x0 == 0.75){
      b0 <- 1.23
      b1 <- -0.858
      b2 <- -4.10
      a0 <- 1
      a1 <- -0.0710
      a2 <- -2.17
      c <- -0.396
      d <- 0.627
    }else if (x0 == 2.5){
      b0 <- 1.69
      b1 <- -0.255
      b2 <- -5.01
      a0 <- 1
      a1 <- -0.0079
      a2 <- -2.11
      c <- -0.482
      d <- 0.949
    }
    f <- function(x){
      pmax((a0 + a1*x*2 + a2*x^2*4)*(b0 + b1*x*2 + b2*x^2*4) + c,0)/((a0 + a1*x*2 + a2*x^2*4)^2 + d*(a0 + a1*x*2 + a2*x^2*4)^4)
    }
    totint <- integrate(f,-0.5,0.5)$value
    f <- function(x){
      pmax((a0 + a1*x*2 + a2*x^2*4)*(b0 + b1*x*2 + b2*x^2*4) + c,0)/((a0 + a1*x*2 + a2*x^2*4)^2 + d*(a0 + a1*x*2 + a2*x^2*4)^4)/totint
    }
  }

  #integrate until Huber's equation is satisfied
  h <- 0.0001
  xj <- vector(mode="double",length = m)
  
  for(j in 1:m){
    i <- -0.5
    val <- 0
    while(val < (j - 0.5)/m){
      val <- integrate(f,-0.5,i)$value
      i <- i + h
    }
    xj[j] <- i
  }
  if (model == 1){
    X <- cbind(rep(1,n),xj)
  }else if (model == 2){
    X <- cbind(rep(1,n),xj,xj^2)
  }
  
  return(X)
}

#####Clustered Design for extrapolation
Cldesign <- function(n,x0,p,model){
  if (model == 1){
    x02 <- 2*x0
    w <- c((x02-1)/(2*x02),(x02+1)/(2*x02))
    xi <- c(seq(-0.5,-0.5+p/2,by = (p/2)/(n*w[1] - 1)),
            seq(0.5-p/2,0.5, by = (p/2)/(n*w[2] - 1)))
  }else if (model == 2){
    x02 <- 2*x0
    w <- c(x02*(x02-1)/(2*(2*x02^2 - 1)),
           (x02-1)*(x02+1)/(2*x02^2 - 1),
           x02*(x02+1)/(2*(2*x02^2 - 1)))
    xi <- c(seq(-0.5,-0.5+p/3,by = (p/3)/(n*w[1] - 1)),
            seq(-p/6,p/6,by = (p/3)/(n*w[2] - 1)),
            seq(0.5-p/3,0.5,by = (p/3)/(n*w[3] - 1)))
  }
  if (model == 1){
    X <- cbind(rep(1,n),xi)
  }else if (model == 2){
    X <- cbind(rep(1,n),xi,xi^2)
  }
  
  
  return(X)
}

#####Find optimal p for extrapolation
Optp <- function(n,N,x0,sigma2vector,k2vector,k3vector,model){
  prange <- seq(0.05,0.95,by=0.01)
  opt <- vector(mode="double",length = length(sigma2vector))
  #Matrix for storing estimated MSE for each p for each sigma^2
  pMSE <- matrix(0,nrow = length(prange), ncol = length(sigma2vector))
  colnames(pMSE) <- c("0.01","0.05","0.1","0.2","0.5")
  
  for (joe in 1:10){
    print(joe)
    for (j in 1:length(prange)){
      p <- prange[j]
      X <- Cldesign(n,x0,p,model)
      
      #Estimate det(MSE) for each sigma^2 for this p
      for (i in 1:length(sigma2vector)){
        pMSE[j,i] <- pMSE[j,i] + MSEextrapolation(n,N,X,x0,sigma2vector[i],k2vector[i],k3vector[i],model = model)[[3]]
        #pMSE[j,i] <- pMSE[j,i] + MSEextrapolation(n,N,X,x0,sigma2vector[i],k2vector[i],k3vector[i],model = model)[[1]]
      }
    }
  }
  pMSE <- pMSE/10
  
  #Plots
  for(i in 1:length(sigma2vector)){
    #Code for printing plots
    if (model == 1){
      plot(prange,pMSE[,i],xlab = "p",ylab = "MSE(yhat(x0))", pch = 19,col = "gray45",
           main = paste("Model ",model,".  sigma^2 =",sigma2vector[i],"\nk2 =",k2vector[i],
                        ", k3 =",k3vector[i],sep = ""))
    }else if (model == 2){
      plot(prange,pMSE[,i],xlab = "p",ylab = "MSE(yhat(x0))", pch = 19,col = "gray45",
           main = paste("Model ",model,".  sigma^2 =",sigma2vector[i],"\nk3 =",k2vector[i],
                        ", k4 =",k3vector[i],sep = ""))
    }
    
    fhat <- ksmooth(prange,pMSE[,i],kernel = "normal", bandwidth = 0.15,
                    range.x = range(prange),n.points = 100, x.points = prange)
    
    #Code for saving plots
    lines(fhat)
    dev.copy(png,paste("model",model,"s",sigma2vector[i],"x0_",x0,"EP.png",sep=""))
    dev.off()
    
    opt[i] <- fhat$x[which.min(fhat$y)]
  }
  return(opt)
}

#####Results tables function 
resultstable <- function(X,n,N,x0,kresults,Clust,model){
  sigma2vector <- kresults[,1]
  k2vector <- kresults[,4+model]
  k3vector <- kresults[,5+model]
  pvector <- kresults[,6+model]
  
  results <- matrix(NA,nrow = length(sigma2vector),ncol = 6)
  if (model == 1){
    colnames(results) <- c("Sigma^2","k2","k3","EBias","EVariance","EMSE")
  }else if (model == 2){
    colnames(results) <- c("Sigma^2","k3","k4","EBias","EVariance","EMSE")
  }
  results[,1] <- sigma2vector
  results[,2] <- k2vector
  results[,3] <- k3vector
  for (i in 1:length(sigma2vector)){
    if (Clust == T){
      p <- pvector[i]
      X <- Cldesign(n,x0,p,model)
    }
    design <- MSEextrapolation(n,N,X,x0,sigma2 = sigma2vector[i],k2 = k2vector[i],k3 = k3vector[i],model = model)
    results[i,4] <- design[[1]] #EBias
    results[i,5] <- design[[2]] #EVariance
    results[i,6] <- design[[3]] #EMSE
  }
  
  if (Clust == T){
    results <- cbind(results,pvector)
    colnames(results)[7] <- "p"
  }
  return(results)
}


#####The following results were obtained in Maple
kL2 <- function(x0,model,sigma2vector,n,maxk){
  values <- matrix(NA,nrow = 5,ncol = 3+model)
  
  if (model == 1){
    if (x0 == 0.75){     #Second extrapolation bias condition
      values[,1] <- c(1.0198,1.0443,1.0626,1.0885,1.1399) #beta0
      values[,2] <- c(0.9645,0.9206,0.8877,0.8412,0.7489) #beta1
      values[,3] <- c(0.2375,0.5310,0.7510,1.0621,1.6793) #k2
      values[,4] <- c(-0.2367,-0.5293,-0.7485,-1.0586,-1.6738) #k3
    }else if (x0 == 2.5){
      values[,1] <- c(1.0021,1.0040,1.0066,1.0093,1.0140) #beta0
      values[,2] <- c(0.9982,0.9963,0.9943,0.9924,0.9874) #beta1
      values[,3] <- c(0.025,0.048,0.079,0.111,0.168) #k2
      values[,4] <- c(-0.012,-0.025,-0.038,-0.051,-0.084) #k3
    }
    
    
  }else if (model == 2){
    if (x0 == 0.75){     #Second extrapolation bias condition
      values[,1] <- c(0.9901,0.9778,0.9686,0.9556,0.9298) #beta0
      values[,2] <- c(0.8674,0.7034,0.5806,0.4069,0.0622) #beta1
      values[,3] <- c(1.3972,1.8883,2.2562,2.7765,3.8089) #beta2
      values[,4] <- c(-0.8842,-1.9771,-2.7961,-3.9543,-6.2522) #k3
      values[,5] <- c(1.8538,4.1452,5.8622,8.2904,13.1083) #k4
    }else if (x0 == 2.5){
      values[,1] <- c(1.0000,1.0000,1.0001,1.0001,1.0002) #beta0
      values[,2] <- c(1.0005,1.0033,1.0053,1.0080,1.0122) #beta1
      values[,3] <- c(0.9996,0.9981,0.9964,0.9951,0.9916) #beta2
      values[,4] <- c(0.003,0.022,0.035,0.053,0.081) #k3
      values[,5] <- c(-0.002,-0.009,-0.017,-0.023,-0.039) #k4
    }
    
  }
  
  return(values)
}               