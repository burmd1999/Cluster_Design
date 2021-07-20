#Function kL2

#Function used to determine k values for model 4,5,6 (linear, quadratic, MLR where
#there are 3 missing true parameters)


kL2 <- function(n,sigma2vector,maxk,model){
  
  kresults <- vector(mode="double",length = length(sigma2vector))
  
  for (i in 1:length(sigma2vector)){
    sigma2 <- sigma2vector[i]
    print(sigma2)
    #eta2 - boundary condition
    eta2 <- sigma2/n
    
    if (model == 1|model == 4){
      for (k in seq(maxk,0,by=-0.00001)){
        
        Betahat <- vector(mode="double",length = 2)
        Betahat[1] <- 1 + 2*k*(0.5^3/3 + 0.5^5/5)
        Betahat[2] <- 1 + 3*k*0.5^2/5
        
        #integrate until integral is less than eta2
        integral <- integrate(function(x){
          (1 + x + k*(x^2 + x^3 + x^4) - Betahat[1] - Betahat[2]*x)^2
        },lower = -0.5,upper = 0.5)[[1]]
        
        val <- (1 + 0.75 + k*(0.75^2 + 0.75^3 + 0.75^4) - Betahat[1] - Betahat[2]*0.75)^2
        
        if (integral < eta2 && val < eta2){
          break
        }
      }
    }else if (model == 2|model == 5){
      for (k in seq(maxk,0,by=-0.00001)){
        
        Betahat <- vector(mode = "double",length = 3)
        Betahat[1] <- -0.5^2/3*(1+k*0.5^2*(1/5 - 3/7)/(1/3 - 3/5)) + 1 + 0.5^2/3 + k*0.5^4/5
        Betahat[2] <- 1 + 3*k*(0.5^2/5 + 0.5^4/7)
        Betahat[3] <- 1 + k*0.5^2*(1/5 - 3/7)/(1/3 - 3/5)
        
        #integrate until integral is less than eta2
        integral <- integrate(function(x){
          (1 + x + x^2 + k*(x^3 + x^4 + x^5) - Betahat[1] - Betahat[2]*x - Betahat[3]*x^2)^2
        },lower = -0.5,upper = 0.5)[[1]]
        
        if (integral < eta2){
          break
        }
      }
    }else if (model == 3|model == 6){
      for (k in seq(maxk,0,by=-0.00001)){
        
        Betahat <- vector(mode= "double",length = 3)
        Betahat[1] <- 1 + 2*k*0.5^2/3
        Betahat[2] <- 1
        Betahat[3] <- 1
        
        #integrate until integral is less than eta2
        integral <- integrate(function(y) {
          sapply(y,function(y){
            integrate(function(x){
              sapply(x,function(x){
                (1 + x + y + k*(x*y + x^2 + y^2) - Betahat[1] - Betahat[2]*x - Betahat[3]*y)^2
              })
            },-0.5,0.5)$value
          })
        },-0.5,0.5)$value
        
        if (integral < eta2){
          break
        }
      }
    }
    
    kresults[i] <- round(k,4)
  }
  
  return(kresults)
}