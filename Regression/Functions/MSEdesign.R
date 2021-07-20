#Inputs
#model = 1: Linear
#model = 2: quadratic
#model = 3: MLR

#model = 1 or 2:
#X: nx2 matrix of vectors (1,x)
#model = 3:
#X: nx3 matrix of vectors (1,x1,x2)

#model = 4,5,6 are the same as 1,2,3 except with 3 missing true parameters in regression model

#n: sample size
#N: Number of iterations (for estimating Bhat)
#sigma2: Variance of ei's
#model 1: k is true value of B2
#model 2: k is true value of B3
#model 3: k is true value of B12

#Outputs
#output: list containing biases, variances and determinants of interest

MSEdesign <- function(X,n,N,sigma2,k,model){
  B0 <- 1
  B1 <- 1
  
  #####LINEAR MODEL
  if (model == 1|model == 4){
    BetaN <- matrix(0,N,2)
    
    for (j in 1:N){
      #Simulate residuals and y values
      y <- vector(mode="double",length = n)
      if (model == 1){
        for (i in 1:n){
          y[i] <- B0 + B1*X[i,2] + k*X[i,2]^2 + rnorm(1, mean = 0, sd = sqrt(sigma2))
        }
      }else if (model == 4){
        for (i in 1:n){
          y[i] <- B0 + B1*X[i,2] + k*(X[i,2]^2 + X[i,2]^3 + X[i,2]^4) + rnorm(1, mean = 0, sd = sqrt(sigma2))
        }
      }

      #Estimate Beta for each simulated dataset
      BetaN[j,] <- solve(t(X)%*%X)%*%t(X)%*%y
    }
    
    #Compute MSE estimate
    BiasB <- colMeans(BetaN) - c(1,1)
    VarB <- c(var(BetaN[,1]),var(BetaN[,2]))
  }
  
  #####QUADRATIC AND MLR MODEL
  if (model == 2|model == 3|model == 5|model == 6){
    B2 <- 1
    BetaN <- matrix(0,N,3)
    
    if (model == 2|model == 5){
      X <- cbind(X,X[,2]^2)
    }
    
    for (j in 1:N){
      #Simulate residuals and y values
      y <- vector(mode="double",length = n)
      if (model == 2|model == 3){
        for (i in 1:n){
          #Note that this works for both quadratic model and MLR model
          y[i] <- B0 + B1*X[i,2] + B2*X[i,3] + k*X[i,2]*X[i,3] + rnorm(1, mean = 0, sd = sqrt(sigma2))
        }
      }else if (model == 5){
        for (i in 1:n){
          y[i] <- 1 + X[i,2] + X[i,2]^2 + k*(X[i,2]^3 + X[i,2]^4 + X[i,2]^5) + rnorm(1,mean=0,sd=sigma2)
        }
      }else if (model == 6){
        for (i in 1:n){
          y[i] <- 1 + X[i,2] + X[i,3] + k*(X[i,2]*X[i,3] + X[i,2]^2 + X[i,3]^2) + rnorm(1,mean=0,sd=sigma2)
        }
      }
      
      #Estimate Beta for each simulated dataset
      BetaN[j,] <- solve(t(X)%*%X)%*%t(X)%*%y
    }
    
    #Compute bias and variances of beta
    BiasB <- colMeans(BetaN) - c(1,1,1)
    VarB <- c(var(BetaN[,1]),var(BetaN[,2]),var(BetaN[,3]))
  }
  
  
  ####Compute Cov and MSE and determinants and output list of values
  CovB <- sigma2*solve(t(X)%*%X)
  #CovB <- cov(BetaN)
  detCovB <- det(CovB)
  MSEBhat <- CovB + BiasB%*%t(BiasB)
  detMSEBhat <- det(MSEBhat)
  
  
  output <- list(detMSEBhat,detCovB,BiasB,VarB)
  return(output)
}