#Finds the max value of k for which Ho: B2 = 0 is accepted for all sigmas inputted

#Inputs
#model = 1: Linear
#model = 2: quadratic
#model = 3: MLR

#model = 1 or 2:
#X: nx2 matrix of vectors (1,x)
#model = 3:
#X: nx3 matrix of vectors (1,x1,x2)

#n: sample size
#N: Number of hypothesis tests
#sigma2vector: vector of sigma^2's (variance) to find k for
#model 1: maxBplus is max true value of B2 to test
#model 2: maxBplus is max true value of B3 to test
#model 3: maxBplus is max true value of B12 to test
#step: amount to decrease from maxBplus by until k is found


#Output
#kresults: value of k for each sigma in sigmavector

kaccept <- function(X,sigma2vector,n,N,maxBplus,step,model){
  B0 <- 1
  B1 <- 1
  B2 <- 1
  
  if (model == 1){
    X <- cbind(X,X[,2]^2)
  }else if (model == 2){
    X <- cbind(X,X[,2]^2,X[,2]^3)
  }else if (model == 3){
    X <- cbind(X,X[,2]*X[,3])
  }
  
  #Create results vector
  kresults <- vector(mode="double",length = length(sigma2vector))
  
  for (joe in 1:100){
    print(joe)
    #Count variable for inputting results
    sigmacount <- 0
    #Start looping through sigma^2 values
    for (sigma2 in sigma2vector){
      sigmacount <- sigmacount + 1
      
      #Loop through k values starting large and working way down until H0 is accepted
      for (k in seq(maxBplus,0,by=-step)){
        #k is the true value of B2
        Bplus <- k
        
        #Responses for each N trial
        accept <- vector(mode="double",length=N)
        #Loop through N hypothesis tests of this comb. of X,k and sigma
        for (j in 1:N){
          #Create y (response) vector
          y <- vector(mode="double",length = n)
          
          #Generate simulated data and values needed for hypothesis test based on model chosen
          if (model == 1){
            for (i in 1:n){
              y[i] <- B0 + B1*X[i,2] + Bplus*X[i,3] + rnorm(1, mean = 0, sd = sqrt(sigma2))
            }
            Bhat <- solve(t(X)%*%X)%*%t(X)%*%y
            CovB <- sigma2*solve(t(X)%*%X)
            VarBplushat <- CovB[3,3]
            Bplushat <- Bhat[3]
          }else if (model == 2|model == 3){
            for (i in 1:n){
              y[i] <- B0 + B1*X[i,2] + B2*X[i,3] + Bplus*X[i,4] + rnorm(1, mean = 0, sd = sqrt(sigma2))
            }
            Bhat <- solve(t(X)%*%X)%*%t(X)%*%y
            CovB <- sigma2*solve(t(X)%*%X)
            VarBplushat <- CovB[4,4]
            Bplushat <- Bhat[4]
          }
          
          #If Ho is rejected, indicate with a 1
          if (abs(Bplushat) <= 1.96*sqrt(VarBplushat)){
            accept[j] <- 1
          }
        }
        
        #Check if Ho was rejected less than 95% of times
        #If so, we have found the max k value and exit loop
        if (sum(accept) > 0.95*N){
          kresults[sigmacount] <- kresults[sigmacount] + k
          break
        }
      }
    }
  }
  kresults <- kresults/100
  
  return(kresults)
}