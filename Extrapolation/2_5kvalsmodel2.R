rm(list=ls())
cat('\f')

sigma2vector <- c(0.01,0.05,0.1,0.2,0.5)
x0 <- 2.5
n <- 49

kresults <- matrix(NA,nrow = 9,ncol = 5)
colnames(kresults) <- c("0.01","0.05","0.1","0.2","0.5")
rownames(kresults) <- c("beta0","beta1","beta2","k3","k4","er1","er2","er1+er2","sigma2/n")

kmin <- 0.001
kmax <- 3
kby <- kmin


######er1 only
for (j in 1:length(sigma2vector)){
  sigma <- sqrt(sigma2vector[j])
  print(sigma2vector[j])
  
  #table for error results for each k2,k3
  table <- matrix(NA,nrow = length(seq(-kmax,kmax,by = kby)),ncol = 5)
  colnames(table) <- c("k3","k4","er1=eta2-intS","er2=eta2-intExtrap","er1+er2")
  table[,2] <- sort(seq(-kmax,kmax,by = kby))
  
  for (i in 1:nrow(table)){
    #Ignore cases with imaginary solutions
    if (is.na( (1.227736317*10^(-10)*(-1.762757143*10^10*table[i,2]*n + sqrt(-4.956216559*10^18*table[i,2]^2*n^2 + 8.145071429*10^17*n*sigma^2)))/n) == T){
      table[i,1] <- NA
      table[i,3] <- NA
      table[i,4] <- NA
      table[i,5] <- NA
    }else{
      #k3 values
      table[i,1] <- (1.227736317*10^(-10)*(-1.762757143*10^10*table[i,2]*n + sqrt(-4.956216559*10^18*table[i,2]^2*n^2 + 8.145071429*10^17*n*sigma^2)))/n
      #define function given k3,k4
      beta2 <- 1 + 3*table[i,2]/14
      beta0 <- 13/12 + table[i,2]/80 - beta2/12
      beta1 <- 1 + 3*table[i,1]/20


      f <- function(x){(1 + x + x^2 + table[i,1]*x^3 + table[i,2]*x^4 - beta0 - beta1*x - beta2*x^2)^2}
      
      #compute errors and sum of errors
      table[i,3] <- sigma^2/n - integrate(f,-0.5,0.5)$value
      table[i,4] <- sigma^2/n - integrate(f,0.5,x0)$value
      table[i,5] <- table[i,3] + table[i,4]
    }
  }
  table2 <- table[is.na(table[,1])==F,]
  table2 <- table2[table2[,3]>=0,]
  if (nrow(table2) == 0){
    print(paste("sigma2 = ",sigma^2," - no valid solutions"))
  }else{
    kresults[4,j] <- table2[which.min(table2[,3]),1] #k3
    kresults[5,j] <- table2[which.min(table2[,3]),2] #k4
    kresults[3,j] <- 1 + 3*kresults[5,j]/14 #beta2
    kresults[1,j] <- 13/12 + kresults[5,j]/80 - kresults[3,j]/12 #beta0
    kresults[2,j] <- 1 + 3*kresults[4,j]/20 #beta1
    
    #Errors
    f <- function(x){(1 + x + x^2 + kresults[4,j]*x^3 + kresults[5,j]*x^4 - kresults[1,j] - kresults[2,j]*x - kresults[3,j]*x^2)^2}
    kresults[6,j] <- sigma^2/n - integrate(f,-0.5,0.5)$value
    kresults[7,j] <- sigma^2/n - integrate(f,0.5,x0)$value
    kresults[8,j] <- kresults[6,j] + kresults[7,j]
    
    kresults[9,j] <- sigma^2/n
  }
}


######er2 only
for (j in 1:length(sigma2vector)){
  sigma <- sqrt(sigma2vector[j])
  print(sigma2vector[j])
  
  #table for error results for each k2,k3
  table <- matrix(NA,nrow = length(seq(-kmax,kmax,by = kby)),ncol = 5)
  colnames(table) <- c("k3","k4","er1=eta2-intS","er2=eta2-intExtrap","er1+er2")
  table[,2] <- sort(seq(-kmax,kmax,by = kby))
  
  for (i in 1:nrow(table)){
    #Ignore cases with imaginary solutions
    if (is.na( 2*sqrt(-7*n*(table[i,2]^2*n - 44100*sigma^2))/(21*n)) == T){
      table[i,1] <- NA
      table[i,3] <- NA
      table[i,4] <- NA
      table[i,5] <- NA
    }else{
      #k3 values
      table[i,1] <- 2*sqrt(-7*n*(table[i,2]^2*n - 44100*sigma^2))/(21*n)
      #define function given k3,k4
      beta2 <- 1 + 3*table[i,2]/14
      beta0 <- 13/12 + table[i,2]/80 - beta2/12
      beta1 <- 1 + 3*table[i,1]/20
      
      
      f <- function(x){(1 + x + x^2 + table[i,1]*x^3 + table[i,2]*x^4 - beta0 - beta1*x - beta2*x^2)^2}
      
      #compute errors and sum of errors
      table[i,3] <- sigma^2/n - integrate(f,-0.5,0.5)$value
      table[i,4] <- sigma^2/n - integrate(f,0.5,x0)$value
      table[i,5] <- table[i,3] + table[i,4]
    }
  }
  table2 <- table[is.na(table[,1])==F,]
  table2 <- table2[table2[,4]>=0,]
  if (nrow(table2) == 0){
    print(paste("sigma2 = ",sigma^2," - no valid solutions"))
  }else{
    kresults[4,j] <- table2[which.min(table2[,4]),1] #k3
    kresults[5,j] <- table2[which.min(table2[,4]),2] #k4
    kresults[3,j] <- 1 + 3*kresults[5,j]/14 #beta2
    kresults[1,j] <- 13/12 + kresults[5,j]/80 - kresults[3,j]/12 #beta0
    kresults[2,j] <- 1 + 3*kresults[4,j]/20 #beta1
    
    #Errors
    f <- function(x){(1 + x + x^2 + kresults[4,j]*x^3 + kresults[5,j]*x^4 - kresults[1,j] - kresults[2,j]*x - kresults[3,j]*x^2)^2}
    kresults[6,j] <- sigma^2/n - integrate(f,-0.5,0.5)$value
    kresults[7,j] <- sigma^2/n - integrate(f,0.5,x0)$value
    kresults[8,j] <- kresults[6,j] + kresults[7,j]
    
    kresults[9,j] <- sigma^2/n
  }
}

######er1 + er2
for (j in 1:length(sigma2vector)){
  ktot <- length(seq(kmin,kmax,by = kby))
  sigma <- sqrt(sigma2vector[j])
  print(sigma2vector[j])
  
  #table for error results for each k2,k3
  table <- matrix(NA,nrow = ktot^2,ncol = 5)
  colnames(table) <- c("k3","k4","er1=eta2-intS","er2=eta2-intExtrap","er1+er2")
  table[,1:2] <- cbind(sort(rep(seq(kmin,kmax,by = kby),ktot)),rep(seq(-kmin,-kmax,by = -kby),ktot))
  
  for (i in 1:nrow(table)){
    
    #define function given k3,k4
    beta2 <- 1 + 3*table[i,2]/14
    beta0 <- 13/12 + table[i,2]/80 - beta2/12
    beta1 <- 1 + 3*table[i,1]/20
    f <- function(x){(1 + x + x^2 + table[i,1]*x^3 + table[i,2]*x^4 - beta0 - beta1*x - beta2*x^2)^2}
      
    #compute errors and sum of errors
    table[i,3] <- sigma^2/n - integrate(f,-0.5,0.5)$value
    table[i,4] <- sigma^2/n - integrate(f,0.5,x0)$value
    table[i,5] <- table[i,3] + table[i,4]
  }
  
  #remove negative errors
  table2 <- table[table[,3]>=0,]
  table2 <- table2[table2[,4]>=0,]

  kresults[4,j] <- table2[which.min(table2[,5]),1] #k3
  kresults[5,j] <- table2[which.min(table2[,5]),2] #k4
  kresults[3,j] <- 1 + 3*kresults[5,j]/14 #beta2
  kresults[1,j] <- 13/12 + kresults[5,j]/80 - kresults[3,j]/12 #beta0
  kresults[2,j] <- 1 + 3*kresults[4,j]/20 #beta1
  
  #Errors
  f <- function(x){(1 + x + x^2 + kresults[4,j]*x^3 + kresults[5,j]*x^4 - kresults[1,j] - kresults[2,j]*x - kresults[3,j]*x^2)^2}
  kresults[6,j] <- sigma^2/n - integrate(f,-0.5,0.5)$value
  kresults[7,j] <- sigma^2/n - integrate(f,0.5,x0)$value
  kresults[8,j] <- kresults[6,j] + kresults[7,j]
  
  kresults[9,j] <- sigma^2/n
}





kresults <- round(kresults,digits = 6)
