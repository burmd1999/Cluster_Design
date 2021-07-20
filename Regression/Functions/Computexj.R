#Find xj for Huber's implemented design

#Inputs:
#model = 1: Linear
#model = 2: quadratic
#model = 3: MLR

#m: number of distinct design points

#Output:
#design points

Computexj <- function(m,model){
  
  if (model == 4){
    model <- 1
  }else if (model == 5){
    model <- 2
  }else if (model == 6){
    model <- 3
  }
  
  
  #####LINEAR AND QUADRATIC INTEGRATION
  if (model == 1|model == 2){
    if (model == 1){
      f <- function(x){
        5.12*x^2 + 0.573
      }
    }else if (model == 2){
      f <- function(x){
        35.0934*((x^2-0.1487^2)*(x^2-0.1489^2)+0.0192)
      }
    }
    
    h <- 0.0001
    xj <- vector(mode="double",length = m)
    
    for(j in 1:m){
      i <- -0.5
      x <- 0
      while(x < (j - 0.5)/m){
        x <- integrate(f,-0.5,i)$value
        i <- i + h
      }
      xj[j] <- i
    }
  #####MLR DOUBLE INTEGRATION  
  }else if (model == 3){
    h <- 0.0001
    xj <- matrix(NA,nrow = m,ncol = 2)
    
    for(j in 1:(m/4)){
      i <- -0.5
      x <- 0
      while(x < (j - 0.5)/m){
        x <- integrate(function(y) {
          sapply(y,function(y){
            integrate(function(x){
              sapply(x,function(x){181.02*pmax((x^2-0.2333^2+0.044)+(y^2-0.2333^2+0.044),0)} )
            },-0.5,y)$value
          })
        },-0.5,i)$value
        i <- i + h
      }
      xj[m - j + 1,] <- c(i,i)
    }
    
    #Solve for remaining design points using symmetry
    if (m == 4){
      xj[1,] <- c(-i,-i)
      xj[2,] <- c(-i,i)
      xj[3,] <- c(i,-i)
    }else if (m > 4){
      xj[1:(m/4),] <- apply(-xj[(3*m/4+1):m,],2,rev)
      xj[(m/4+1):(2*m/4),] <- cbind(xj[1:(m/4),1],rev(xj[(3*m/4+1):m,2]))
      xj[(2*m/4+1):(3*m/4),] <- cbind(xj[(3*m/4+1):m,1],rev(xj[1:(m/4),2]))
    }

  }
  
  return(xj)
}
