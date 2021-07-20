#Makes results table for designs

#Inputs:
#Cldesign: boolean stating whether its for the cluster design (T) or not (F)
#If Cldesign is false pvector can be anything
#If Cldesign is true then X can be anything
#Other variables are the same as other functions

resultstable <- function(X,n,N,sigma2vector,kvector,Cldesign,pvector,model){
  #####For linear model since table has less columns than model 2 or 3
  if (model == 1|model == 4){
    results <- matrix(NA,nrow = length(sigma2vector),ncol = 8)
    colnames(results) <- c("Sigma^2","k","BiasBo","BiasB1","VarB0","VarB1",
                           "det(Cov)","det(MSE)")
    results[,1] <- sigma2vector
    results[,2] <- kvector
    for (i in 1:length(sigma2vector)){
      if (Cldesign == T){
        p <- pvector[i]
        X <- cbind(rep(1,n),c(seq(-0.5,-0.5 + p/2,by=(p/2)/(n/2-1)),seq(0.5 - p/2,0.5,by=(p/2)/(n/2-1))))
      }
      design <- MSEdesign(X,n,N,sigma2 = sigma2vector[i],k = kvector[i],model = model)
      results[i,3] <- design[[3]][[1]] #BiasB0
      results[i,4] <- design[[3]][[2]] #BiasB1
      results[i,5] <- design[[4]][[1]] #VarB0
      results[i,6] <- design[[4]][[2]] #VarB1
      results[i,7] <- design[[2]] #det(cov)
      results[i,8] <- design[[1]] #det(MSE)
    }
    if (Cldesign == T){
      results <- cbind(results,pvector)
      colnames(results)[9] <- "p"
    }
    ##### For quadratic and MLR models
  }else if (model == 2|model == 3|model == 5|model == 6){
    results <- matrix(0,nrow = length(sigma2vector),ncol = 10)
    colnames(results) <- c("Sigma^2","k","BiasB0","BiasB1","BiasB2","VarB0","VarB1","VarB2",
                            "det(Cov)","det(MSE)")
    results[,1] <- sigma2vector
    results[,2] <- kvector
    for (i in 1:length(sigma2vector)){
      if (Cldesign == T){
        p <- pvector[i]
        if (model == 2|model == 5){
          X <- cbind(rep(1,n),c(seq(-0.5,-0.5+p/3,by = p/3/(n/3-1)),
                                seq(0-p/6,0+p/6,by = p/3/(n/3-1)),
                                seq(0.5-p/3,0.5,by = p/3/(n/3-1))))
        }else if (model == 3|model == 6){
          X <- cbind(rep(1,n),
                     c(rep(0.5,6),rep(0.5-0.25*sqrt(p),6),rep(0.5-0.5*sqrt(p),6),
                       rep(-0.5+0.5*sqrt(p),6),rep(-0.5+0.25*sqrt(p),6),rep(-0.5,6)),
                     rep(c(0.5,0.5-0.25*sqrt(p),0.5-0.5*sqrt(p),-0.5+0.5*sqrt(p),
                           -0.5+0.25*sqrt(p),-0.5),6)  )
        }
      }
      design <- MSEdesign(X,n,N,sigma2 = sigma2vector[i],k = kvector[i],model = model)
      results[i,3] <- design[[3]][[1]] #BiasB0
      results[i,4] <- design[[3]][[2]] #BiasB1
      results[i,5] <- design[[3]][[3]] #BiasB2
      results[i,6] <- design[[4]][[1]] #VarB0
      results[i,7] <- design[[4]][[2]] #VarB1
      results[i,8] <- design[[4]][[3]] #VarB2
      results[i,9] <- design[[2]] #det(cov)
      results[i,10] <- design[[1]] #det(MSE)
    }
    if (Cldesign == T){
      results <- cbind(results,pvector)
      colnames(results)[11] <- "p"
    }
  }
  
  return(results)
}