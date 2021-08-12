#Comparison of Clustered Designs
rm(list = ls())
cat('\f')

library(xtable)

setwd('C:/Users/dburm/Desktop/School/MATH 4F90/')
source('EPfunctions.R')

N <- 5000
sigma2vector <- c(0.01,0.05,0.1,0.2,0.5)




####MODEL 1
model <- 1
x0 <- 2.5
n <- 30
results <- matrix(NA,nrow = length(sigma2vector),ncol = 6+model)
colnames(results) <- c("sigma2","sigma2/n","beta0","beta1","k2","k3","p")
results[,1] <- sigma2vector
results[,2] <- sigma2vector/n
results[,3:6] <- kL2(x0,model)
results[,7] <- rep(0.05,5)


p <- 0.05
x02 <- 2*x0
w <- c((x02-1)/(2*x02),(x02+1)/(2*x02))


#p = 1/2, 1/2, n = weights
xi <- c(seq(-0.5,-0.5+p/2, by = (p/2)/(n*w[1] - 1)),
        seq(0.5-p/2,0.5, by = (p/2)/(n*w[2] - 1)))
X <- cbind(rep(1,n),xi)
resultsCl1 <- resultstable(X,n,N,x0,kresults = results,Clust = F,model = 1)

#p = w1, w2, n = 1/2,1/2
xi <- c(seq(-0.5,-0.5+p*w[1],by = (p*w[1])/(n/2 - 1)),
        seq(0.5-p*w[2],0.5, by = (p*w[2])/(n/2 - 1)))
X <- cbind(rep(1,n),xi)
resultsCl2 <- resultstable(X,n,N,x0,kresults = results,Clust = F,model = 1)


#p = w1, w2, n = w1, w2
xi <- c(seq(-0.5,-0.5+p*w[1],by = (p*w[1])/(n*w[1] - 1)),
        seq(0.5-p*w[2],0.5, by = (p*w[2])/(n*w[2] - 1)))
X <- cbind(rep(1,n),xi)
resultsCl3 <- resultstable(X,n,N,x0,kresults = results,Clust = F,model = 1)

resultsCl1
resultsCl2
resultsCl3

table <- cbind(sigma2vector,resultsCl1[,6],resultsCl2[,6],resultsCl3[,6])
colnames(table) <- c("$\\sigma^2$","Method 1","Method 2","Method 3")
table1 <- xtable(table,caption = paste("Comparison of SMSE's for different methods of 
                                       constructing clustered design for extrapolation. Model = ",model,", $x_0$ = ",x0,sep = ""),
                 align = c("|c","|c","|c","|c","|c|"),digits = c(0,2,6,6,6))
print(table1, file = paste("methodcomp_model",model,"x0_",x0,".tex",sep=""),append = T,table.placement = "h",
      caption.placement = "bottom",hline.after = seq(from=-1,to=nrow(table1),by=1),include.rownames = F,sanitize.text.function=function(x){x})






####MODEL 2
model <- 2
x0 <- 0.75
n <- 28*3

# x0 <- 2.5
# n <- 49


results <- matrix(NA,nrow = length(sigma2vector),ncol = 8)
colnames(results) <- c("sigma2","sigma2/n","beta0","beta1","beta2","k3","k4","p")
results[,1] <- sigma2vector
results[,2] <- sigma2vector/n
results[,3:7] <- kL2(x0,model)
if (x0 == 0.75){
    results[,8] <- c(0.95,0.95,0.87,0.05,0.05)   
}else if (x0 == 2.5){
    results[,8] <- c(0.80,0.05,0.05,0.05,0.05)
}

x02 <- 2*x0
w <- c(x02*(x02-1)/(2*(2*x02^2 - 1)),
       (x02-1)*(x02+1)/(2*x02^2 - 1),
       x02*(x02+1)/(2*(2*x02^2 - 1)))

for (i in 1:5){
    p <- results[i,8]
    
    #p = 1/3, 1/3, 1/3, n = weights
    xi <- c(seq(-0.5,-0.5+p/3, by = (p/3)/(n*w[1] - 1)),
            seq(-p/6,p/6, by = (p/3)/(n*w[2] - 1)),
            seq(0.5-p/3,0.5, by = (p/3)/(n*w[3] - 1)))
    X <- cbind(rep(1,n),xi,xi^2)
    resultsCl1[i,] <- resultstable(X,n,N,x0,kresults = results,Clust = F,model = 2)[i,]
    
    
    #p = w1,w2,w3,  n = 1/3,1/3,1/3
    xi <- c(seq(-0.5,-0.5+p*w[1],by = (p*w[1])/(n/3 - 1)),
            seq(-p*w[2]/2,p*w[2]/2, by = (p*w[2])/(n/3 - 1)),
            seq(0.5-p*w[3],0.5, by = (p*w[3])/(n/3 - 1)))
    X <- cbind(rep(1,n),xi,xi^2)
    resultsCl2[i,] <- resultstable(X,n,N,x0,kresults = results,Clust = F,model = 2)[i,]
    
    #p = w1,w2,w3    n = w1,w2,w3
    xi <- c(seq(-0.5,-0.5+p*w[1],by = (p*w[1])/(n*w[1] - 1)),
            seq(-p*w[2]/2,p*w[2]/2,by = (p*w[2])/(n*w[2] - 1)),
            seq(0.5-p*w[3],0.5, by = (p*w[3])/(n*w[3] - 1)))
    X <- cbind(rep(1,n),xi,xi^2)
    resultsCl3[i,] <- resultstable(X,n,N,x0,kresults = results,Clust = F,model = 2)[i,]
}

resultsCl1
resultsCl2
resultsCl3
