#Model 4 - Linear with complex contamination function

cat('\f')
rm(list=ls())
#Retrieve functions
setwd('C:/Users/dburm/Desktop/School/MATH 4F90/')
source('MSEdesign.R')
source('Computexj.R')
source('kaccept.R')
source('Optp.R')
source('resultstable.R')
source('kL2.R')
library(xtable)


#####Determine k and p
n <- 30 
N <- 50000
model <- 4
#Create vector of sigma's to be tested
sigma2vector <- c(0.01,0.05,0.1,0.2,0.5)

results <- matrix(0,nrow = length(sigma2vector),ncol = 4)
colnames(results) <- c("Sigma^2","12Sigma^2/n","k","p")
results[,1] <- sigma2vector
results[,2] <- 12*results[,1]/n
results[,3] <- kL2(n,sigma2vector,maxk = 3,model = model)
results[,4] <- Optp(n,N=5000,sigma2vector,kvector = results[,3],model = model,ktype="")


#k-C
#Classical D-Opt Design
X <- cbind(rep(1,n),c(rep(-0.5,n/2),rep(0.5,n/2)))
resultsCDD.1 <- resultstable(X,n,N,sigma2vector,kvector=results[,3],Cldesign = F,pvector = F,model = model) 

#Huber's Implemented Design
X <- cbind(rep(1,n),c(rep(-0.396,n/3),rep(0,n/3),rep(0.396,n/3)))
resultsH3.1 <- resultstable(X,n,N,sigma2vector,kvector=results[,3],Cldesign = F,pvector = F,model = model)

#Uniform Design
X <- cbind(rep(1,n),seq(-0.5,0.5,by=1/(n-1)))
resultsU.1 <- resultstable(X,n,N,sigma2vector,kvector=results[,3],Cldesign = F,pvector = F,model = model)

#Clustered Design for model 1
resultsCl.1 <- resultstable(X=F,n,N,sigma2vector,kvector=results[,3],Cldesign = T,pvector = results[,4],model = model)

#Huber's Implemented Design
X <- cbind(rep(1,n),Computexj(m = 30,model = 1))
resultsH30.2 <- resultstable(X,n,N,sigma2vector,kvector = results[,3],Cldesign = F,pvector = F,model = model)


#Tables of interest
BiasVarTable.1 <- matrix(NA,nrow = length(sigma2vector)*6,ncol = 5)
colnames(BiasVarTable.1) <- c("BiasB0","BiasB1","VarB0","VarB1","det(MSE)")
for (i in 1:length(sigma2vector)){
  BiasVarTable.1[6*(i-1)+1, 1:4] <- c(results[i,c(1,3,4)],N)
  BiasVarTable.1[6*(i-1)+2, 1:5] <- c(resultsCDD.1[i,c(3,4,5,6,8)])
  BiasVarTable.1[6*(i-1)+3, 1:5] <- c(resultsH3.1[i,c(3,4,5,6,8)])
  BiasVarTable.1[6*(i-1)+4, 1:5] <- c(resultsH30.2[i,c(3,4,5,6,8)])
  BiasVarTable.1[6*(i-1)+5, 1:5] <- c(resultsCl.1[i,c(3,4,5,6,8)])
  BiasVarTable.1[6*i, 1:5] <- c(resultsU.1[i,c(3,4,5,6,8)])
}

EffTable.1 <- matrix(NA,nrow = length(sigma2vector),ncol = 7)
colnames(EffTable.1) <- c("sigma^2","k","p","CDD/CL","H3/CL","H30/CL","U/CL")
EffTable.1[,1] <- sigma2vector
EffTable.1[,2] <- results[,3]
EffTable.1[,3] <- results[,4]
EffTable.1[,4] <- resultsCDD.1[,8]/resultsCl.1[,8]
EffTable.1[,5] <- resultsH3.1[,8]/resultsCl.1[,8]
EffTable.1[,6] <- resultsH30.2[,8]/resultsCl.1[,8]
EffTable.1[,7] <- resultsU.1[,8]/resultsCl.1[,8]

#Save results as RDS
model4results <- list(BiasVarTable = BiasVarTable.1,
                      EffTable = EffTable.1,
                      results = results,
                      resultsCDD = resultsCDD.1,
                      resultsCl = resultsCl.1,
                      resultsH3 = resultsH3.1,
                      resultsH30 = resultsH30.2,
                      resultsU = resultsU.1)
saveRDS(model4results,file = paste("model",model,"results.rds",sep=""))







######Format and export tables for latex
customround <- function(numbers){ #Function for rounding numbers for latex tables
  numbers[abs(numbers) >= 0.0001] <- round(numbers[abs(numbers) >= 0.0001],digits = 5)
  numbers[abs(numbers) < 0.0001] <- format(signif(numbers[abs(numbers) < 0.0001],digits = 2),scientific = T)
  
  return(numbers)
}
#Table 1 - k and p values 
colnames(results) <- c("$\\sigma^2$","$12\\frac{\\sigma^2}{n}$","$k$","$p$")
table1 <- xtable(results,caption = paste("Estimated values of $k$ and $p$ when fitting Model ",model," with $n$ =",n,sep=""),
                 align = c("|c","|c","|c","|c","|c|"),digits = c(0,2,4,4,2))
print(table1,file=paste("model",model,"table1.tex",sep=""),append=T,table.placement="h",caption.placement = "bottom",
      hline.after = seq(from=-1,to=nrow(table1),by=1),include.rownames = F,sanitize.text.function=function(x){x})


#Table 2a - Bias, Variance, det table
for (i in 1:length(sigma2vector)){
  BiasVarTable.1[6*(i-1)+2, 1:5] <- customround(c(resultsCDD.1[i,c(3,4,5,6,8)]))
  BiasVarTable.1[6*(i-1)+3, 1:5] <- customround(c(resultsH3.1[i,c(3,4,5,6,8)]))
  BiasVarTable.1[6*(i-1)+4, 1:5] <- customround(c(resultsH30.2[i,c(3,4,5,6,8)]))
  BiasVarTable.1[6*(i-1)+5, 1:5] <- customround(c(resultsCl.1[i,c(3,4,5,6,8)]))
  BiasVarTable.1[6*i, 1:5] <- customround(c(resultsU.1[i,c(3,4,5,6,8)]))
}
BiasVarTable.1 <- cbind(rep(c(NA,"$\\xi_{CDD}^{(2)}$","$\\xi_{H}^{(3)}$","$\\xi_{H}^{(30)}$","$\\xi_{Cl}$","$\\xi_{U}$"),length(sigma2vector)),BiasVarTable.1)
colnames(BiasVarTable.1) <- c(" ","$Bias(\\hat{\\beta_0})$","$Bias(\\hat{\\beta_1})$","$Var(\\hat{\\beta_0})$",
                              "$Var(\\hat{\\beta_1})$","det(MSE($\\boldsymbol{\\hat{\\beta}}$))")
addtorow <- list()
for (i in 1:length(sigma2vector)){
  addtorow$pos[i] <- list(5*(i-1))
  addtorow$command[i] <- paste("& \\multicolumn{4}{l}{\\underline{$\\sigma^2 = ",BiasVarTable.1[(6*i-5),2],"$,
                               $k = ",BiasVarTable.1[(6*i-5),3],"$,
                               $p = ",BiasVarTable.1[(6*i-5),4],"$,
                               $n = ",n,"$}} & \\\\",sep = "")
}
BiasVarTable.1 <- BiasVarTable.1[-seq(1,length(sigma2vector)*6,by=6),]
table1 <- xtable(BiasVarTable.1,caption = "Table of Bias, Variance and D-optimality criteria.",align=c("|c","|c","c","c","c","c","c|"),
                 digits = c(0,0,5,5,5,5,5))
print(table1,file=paste("model",model,"table2a.tex",sep=""),append=T,table.placement="h",caption.placement="bottom",
      hline.after = seq(from=0,to=nrow(table1),by=5),include.rownames = F,sanitize.text.function=function(x){x},
      add.to.row = addtorow)

#Table 2b - efficiency comparison for kC
colnames(EffTable.1) <- c("$\\sigma^2$","$k$","$p$","Eff($\\xi_{Cl}$,$\\xi_{CDD}^{(2)}$)",
                          "Eff($\\xi_{Cl}$,$\\xi_{H}^{(3)}$)","Eff($\\xi_{Cl}$,$\\xi_{H}^{30}$)",
                          "Eff($\\xi_{Cl}$,$\\xi_{U}$)")
table1 <- xtable(EffTable.1,caption = "Efficiency Comparison of det(MSE($\\boldsymbol{\\hat{\\beta}}$)) between designs.",
                 align=c("|c","|c","|c","|c","|c","|c","|c","|c|"),digits = c(0,2,4,2,4,4,4,4))
print(table1,file=paste("model",model,"table2b.tex",sep=""),append=T,table.placement="h",caption.placemenet="bottom",
      hline.after = seq(from=-1,to=nrow(table1),by=1),include.rownames = F,sanitize.text.function=function(x){x})
