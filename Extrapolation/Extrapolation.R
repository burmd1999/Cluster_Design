#***For every function except kL2.R, for model 2 k3 and k4 still use variable names k2 and k3 respectively
#(Or k3vector and k4vector have variable names k2vector and k3vector respectively) for simplicity of code




#Determine k and optimal p
results <- matrix(NA,nrow = length(sigma2vector),ncol = 6+model)
if (model == 1){
  colnames(results) <- c("sigma2","sigma2/n","beta0","beta1","k2","k3","p")
}else if (model == 2){
  colnames(results) <- c("sigma2","sigma2/n","beta0","beta1","beta2","k3","k4","p")
}
results[,1] <- sigma2vector
results[,2] <- sigma2vector/n
if (model == 1){
  results[,3:6] <- kL2(x0,model)
}else if (model == 2){
  results[,3:7] <- kL2(x0,model)
}
results[,6+model] <- Optp(n,N=100,x0,sigma2vector,k2vector = results[,4+model],k3vector = results[,5+model],model = model)
#results[,6+model] <- c(0.05,0.05,0.05,0.05,0.05)


#Classical (HL) design results table
X <- HLdesign(n,x0,model = model)
resultsHL <- resultstable(X,n,N,x0,kresults = results,Clust=F,model = model)

#WX design (Huber's implemented)
X <- WXdesign(n,x0,m = n,model = model)
resultsWX <- resultstable(X,n,N,x0,kresults = results,Clust=F,model = model)

#Clustered design results table
resultsCl <- resultstable(X=F,n,N,x0,kresults = results,Clust=T,model = model)

#Uniform design results table
X <- cbind(rep(1,n),seq(-0.5,0.5,by=1/(n-1)))
if (model == 2){
  X <- cbind(X,X[,2]^2)
}
resultsU <- resultstable(X,n,N,x0,kresults = results,Clust=F,model = model)

#analysis table
comparisontable <- matrix(NA,nrow = length(sigma2vector)*5,ncol = 4)
colnames(comparisontable) <- c("Bias","Var","MSE","RE")
for (i in 1:length(sigma2vector)){
  comparisontable[5*(i-1)+1, 1:4] <- c(results[i,c(1,4+model,5+model,6+model)]) #sigma,k2,k3,p
  comparisontable[5*(i-1)+2, 1:4] <- c(resultsCl[i,c(4,5,6)],NA)
  comparisontable[5*(i-1)+3, 1:4] <- c(resultsHL[i,c(4,5,6)],resultsHL[i,6]/resultsCl[i,6])
  comparisontable[5*(i-1)+4, 1:4] <- c(resultsWX[i,c(4,5,6)],resultsWX[i,6]/resultsCl[i,6])
  comparisontable[5*(i-1)+5, 1:4] <- c(resultsU[i,c(4,5,6)],resultsU[i,6]/resultsCl[i,6])
}




#Save rds file
EPresults <- list(results = results,
                        resultsHL = resultsHL,
                        resultsCl = resultsCl,
                        resultsWX = resultsWX,
                        resultsU = resultsU,
                        comparisontable = comparisontable)
saveRDS(EPresults,file = paste("model",model,"EPx0_",x0,"results.rds",sep=""))






#Format and export table for latex
customround <- function(numbers){ #Function for rounding numbers for latex tables
  numbers[abs(numbers) >= 0.0001] <- round(numbers[abs(numbers) >= 0.0001],digits = 5)
  numbers[abs(numbers) < 0.0001] <- format(signif(numbers[abs(numbers) < 0.0001],digits = 2),scientific = T)

  return(numbers)
}
#Table 1 - k and p values
if(model == 1){
  colnames(results) <- c("$\\sigma^2$","$\\frac{\\sigma^2}{n}$","$\\beta_0^*$","$\\beta_1^*$","$k_2$","$k_3$","$p$")
  table1 <- xtable(results,caption = paste("Values of $k_2$, $k_3$ and $p$ when fitting Model ",model,
                                           " for extrapolation at $x_0 = ",x0,"$, with $n =",n,"$.",sep=""),
                   align = c("|c","|c","|c","|c","|c","|c","|c","|c|"),digits = c(0,2,4,4,4,4,4,2))
}else if (model == 2){
  colnames(results) <- c("$\\sigma^2$","$\\frac{\\sigma^2}{n}$","$\\beta_0^*$","$\\beta_1^*$","$\\beta_2^*$","$k_3$","$k_4$","$p$")
  table1 <- xtable(results,caption = paste("Values of $k_3$, $k_4$ and $p$ when fitting Model ",model,
                                           " for extrapolation at $x_0 = ",x0,"$, with $n =",n,"$.",sep=""),
                   align = c("|c","|c","|c","|c","|c","|c","|c","|c","|c|"),digits = c(0,2,4,4,4,4,4,4,2))
}

print(table1,file=paste("model",model,"EPx0_",x0,"table1.tex",sep=""),append=T,table.placement="h",caption.placement = "bottom",
      hline.after = seq(from=-1,to=nrow(table1),by=1),include.rownames = F,sanitize.text.function=function(x){x})


#Comparison table
for (i in 1:length(sigma2vector)){
  comparisontable[5*(i-1)+2, 1:3] <- customround(resultsCl[i,c(4,5,6)])
  comparisontable[5*(i-1)+3, 1:4] <- customround(c(resultsHL[i,c(4,5,6)],resultsHL[i,6]/resultsCl[i,6]))
  comparisontable[5*(i-1)+4, 1:4] <- customround(c(resultsWX[i,c(4,5,6)],resultsWX[i,6]/resultsCl[i,6]))
  comparisontable[5*(i-1)+5, 1:4] <- customround(c(resultsU[i,c(4,5,6)],resultsU[i,6]/resultsCl[i,6]))
}
comparisontable <- cbind(rep(c(NA,"$\\xi_{Cl}$","$\\xi_{HL}$","$\\xi_{WX}$","$\\xi_{U}$")),comparisontable)
colnames(comparisontable) <- c(" ","$Bias(\\hat{y}(x_0))$","$Var(\\hat{y}(x_0))$","$MSE(\\hat{y}(x_0))$","Eff($\\xi_{Cl}$,$\\xi$)")
addtorow <- list()
for (i in 1:length(sigma2vector)){
  addtorow$pos[i] <- list(4*(i-1))
  addtorow$command[i] <- paste("& \\multicolumn{3}{l}{\\underline{$\\sigma^2 = ",comparisontable[(5*(i-1)+1),2],"$,
                               $k_",model+1," = ",comparisontable[(5*(i-1)+1),3],"$,
                               $k_",model+2," = ",comparisontable[(5*(i-1)+1),4],"$,
                               $p = ",comparisontable[(5*(i-1)+1),5],"$}} & \\\\",sep = "")
}
comparisontable <- comparisontable[-seq(1,length(sigma2vector)*5,by=5),]
table1 <- xtable(comparisontable,caption = paste("Table of Bias, Variance and D-optimality criteria for Model ",model,
                                                 " extrapolation with $x_0 = ",x0,"$, and $n = ",n,"$.",sep=""),
                 align=c("|c","|c","|c"," c"," c"," c|"), digits = c(0,0,5,5,5,4))
print(table1,file=paste("model",model,"EPx0_",x0,"table2.tex",sep=""),append=T,table.placement="h",caption.placement="bottom",
      hline.after = seq(from=0,to=nrow(table1),by=4),include.rownames = F,sanitize.text.function=function(x){x},
      add.to.row = addtorow)

