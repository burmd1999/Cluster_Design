#Find optimal p values of clustered design for models:

#Inputs:
#model = 1: Linear
#model = 2: quadratic
#model = 3: MLR        ***n MUST be 36 for model = 3

#n: sample size
#N: Number of iterations (for estimating Bhat) to be submitted into MSEdesign function
#sigma2vector: vector of sigma^2's (variance) to find p for
#kvector: vector of k's found using kaccept.R function
#ktype: string for the saved-graph file name


Optp <- function(n,N,sigma2vector,kvector,model,ktype){
  
  prange <- seq(0.05,0.95,by=0.01)
  opt <- vector(mode="double",length = length(sigma2vector))
  #Matrix for storing estimated MSE for each p for each sigma^2
  pMSE <- matrix(0,nrow = length(prange), ncol = length(sigma2vector))
  colnames(pMSE) <- c("0.01","0.05","0.1","0.2","0.5")
  for (joe in 1:10){
    print(joe)
    for (j in 1:length(prange)){
      p <- prange[j]
      print(p)
      
      #Set up clustered design matrix
      if (model == 1|model == 4){
        X <- cbind(rep(1,n),c(seq(-0.5,-0.5 + p/2,by=(p/2)/(n/2-1)),seq(0.5 - p/2,0.5,by=(p/2)/(n/2-1))))
      }else if (model == 2|model == 5){
        X <- cbind(rep(1,n),c(seq(-0.5,-0.5+p/3,by = p/3/(n/3-1)),
                              seq(0-p/6,0+p/6,by = p/3/(n/3-1)),
                              seq(0.5-p/3,0.5,by = p/3/(n/3-1))))
      }else if (model == 3|model == 6){
        X <- cbind(rep(1,n),
                   c(rep(0.5,6),rep(0.5-0.25*sqrt(p),6),rep(0.5-0.5*sqrt(p),6),
                     rep(-0.5+0.5*sqrt(p),6),rep(-0.5+0.25*sqrt(p),6),rep(-0.5,6)),
                   rep(c(0.5,0.5-0.25*sqrt(p),0.5-0.5*sqrt(p),-0.5+0.5*sqrt(p),
                         -0.5+0.25*sqrt(p),-0.5),6)
        )
      }
      #Estimate det(MSE) for each sigma^2 for this p
      for (i in 1:length(sigma2vector)){
        pMSE[j,i] <- pMSE[j,i] + MSEdesign(X,n,N,sigma2vector[i],kvector[i],model = model)[[1]]
      }
    }
  }
  pMSE <- pMSE/10
  
  #Plots
  for(i in 1:length(sigma2vector)){
    plot(prange,pMSE[,i],xlab = "p",ylab = "det(MSE(betahat))", pch = 19,col = "gray45",
         main = paste("Model ",model,".  sigma^2 =",sigma2vector[i],", k =",kvector[i]))

    fhat <- ksmooth(prange,pMSE[,i],kernel = "normal", bandwidth = 0.15,
                    range.x = range(prange),n.points = 100, x.points = prange)
    lines(fhat)
    dev.copy(png,paste("model",model,"s",sigma2vector[i],ktype,".png",sep=""))
    dev.off()
    
    opt[i] <- fhat$x[which.min(fhat$y)]
  }

  return(opt)
}
