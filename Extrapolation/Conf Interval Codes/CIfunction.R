
#Computes all measurements of interest for extrapolation
#Including: relative bias, ECI, ACI, PI, interval lengths, simulated SE, average ASD

CIfunction <- function(n,N,X,x0,sigma2,k2,k3,model){
  if (model == 1){
    f <- function(x,k2,k3){1 + x + k2*x^2 + k3*x^3}
    z <- c(1,x0)
  }else if (model == 2){
    f <- function(x,k2,k3){1 + x + x^2 + k2*x^3 + k3*x^4}
    z <- c(1,x0,x0^2)
  }
  
  #Vector for N extrapolations at x0
  yhat <- vector(mode = "double", length = N)
  betahat <- matrix(NA,nrow = model + 1,ncol = N)
  sigma2hat <- vector(mode = "double", length = N)
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
    
    sigma2hat[j] <- sum((y - yhat[j])^2)/(n - (model + 1))
  }
  
  output <- matrix(NA,nrow = 1,ncol = 11)
  colnames(output) <- c("rel bias (%)","ECI CP","ACI CP","PI CP","ECI avg length","ACI avg length","PI avg length",
                        "True y(x_0)","ybar(x_0)","sim SE","ASD")
  
  y_x0 <- f(x0,k2,k3)
  ###Relative Bias
  output[1] <- (mean(yhat) - y_x0)/y_x0*100
  
  ###Coverage Percentage and CI lengths
  ###Empirical CI (ECI)
  ECIs <- matrix(NA,nrow = 2,ncol = N)
  ECIs[1,] <- yhat + qt(0.05/2, df = n-(model+1),lower.tail = F)*sd(yhat)
  ECIs[2,] <- yhat - qt(0.05/2, df = n-(model+1),lower.tail = F)*sd(yhat)
  count <- 0
  for (j in 1:N){
    if (y_x0 >= ECIs[2,j] && y_x0 <= ECIs[1,j]){
      count <- count + 1
    }
  }
  output[2] <- count/N*100
  output[5] <- mean(ECIs[1,] - ECIs[2,])
  
  ###Asymptotic CI (ACI)
  ACIs <- matrix(NA,nrow = 2,ncol = N)
  ACIs[1,] <- yhat + qt(0.05/2, df = n-(model+1),lower.tail = F)*sqrt(sigma2hat)*c(sqrt(z%*%solve(t(X)%*%X)%*%z))
  ACIs[2,] <- yhat - qt(0.05/2, df = n-(model+1),lower.tail = F)*sqrt(sigma2hat)*c(sqrt(z%*%solve(t(X)%*%X)%*%z))
  count <- 0
  for (j in 1:N){
    if (y_x0 >= ACIs[2,j] && y_x0 <= ACIs[1,j]){
      count <- count + 1
    }
  }
  output[3] <- count/N*100
  output[6] <- mean(ACIs[1,] - ACIs[2,])
  
  ###Prediction Interval (PI)
  PIs <- matrix(NA,nrow = 2,ncol = N)
  PIs[1,] <- yhat + qt(0.05/2, df = n-(model+1),lower.tail = F)*sqrt(sigma2hat)*c(sqrt(z%*%solve(t(X)%*%X)%*%z + 1))
  PIs[2,] <- yhat - qt(0.05/2, df = n-(model+1),lower.tail = F)*sqrt(sigma2hat)*c(sqrt(z%*%solve(t(X)%*%X)%*%z + 1))
  count <- 0
  for (j in 1:N){
    if (y_x0 >= PIs[2,j] && y_x0 <= PIs[1,j]){
      count <- count + 1
    }
  }
  output[4] <- count/N*100
  output[7] <- mean(PIs[1,] - PIs[2,])
  
  output[8] <- y_x0     #True y(x0)
  output[9] <- mean(yhat)  #ybar(x0)
  output[10] <- sd(yhat)   #simulated SE
  output[11] <- mean(sqrt(sigma2hat))*c(sqrt(z%*%solve(t(X)%*%X)%*%z))   #average Asymptotic Stand. Dev. (ASD)
  
  return(output)
}