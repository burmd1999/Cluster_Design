#Code for Relative Bias, Coverage Percentage and avg length of each conf. interval
#Computes all above values of interest and exports results into a latex table


comptable <- matrix(NA,nrow = 5*length(sigma2vector),ncol = 7)
colnames(comptable) <- c("rel bias (%)","ECI CP","ACI CP","PI CP","ECI avg length","ACI avg length","PI avg length")
rownames(comptable) <- rep(c("","HL","WX","Cl","U"),length(sigma2vector))

comptable2 <- matrix(NA,nrow = 5*length(sigma2vector),ncol = 4)
colnames(comptable2) <- c("True y(x0)","ybar(x0)","SE","ASD")
rownames(comptable2) <- rep(c("","HL","WX","Cl","U"),length(sigma2vector))

for (i in 1:length(sigma2vector)){
  comptable[5*(i-1)+1,1:4] <- c(results[i,c(1,4+model,5+model,6+model)])
  comptable2[5*(i-1)+1,1:4] <- c(results[i,c(1,4+model,5+model,6+model)])
}

#HL design
X <- HLdesign(n,x0,model = model)
for (i in 1:length(sigma2vector)){
  values <- CIfunction(n,N,X,x0,sigma2vector[i],results[i,4+model],results[i,5+model],model = model)
  comptable[5*(i-1)+2,] <- values[1:7]
  comptable2[5*(i-1)+2,] <- values[8:11]
}

#WX design
X <- WXdesign(n,x0,m = n,model = model)
for (i in 1:length(sigma2vector)){
  values <- CIfunction(n,N,X,x0,sigma2vector[i],results[i,4+model],results[i,5+model],model = model)
  comptable[5*(i-1)+3,] <- values[1:7]
  comptable2[5*(i-1)+3,] <- values[8:11]
}

#Clustered Design
for (i in 1:length(sigma2vector)){
  X <- Cldesign(n,x0,results[i,6+model],model=model)
  values <- CIfunction(n,N,X,x0,sigma2vector[i],results[i,4+model],results[i,5+model],model = model)
  comptable[5*(i-1)+4,] <- values[1:7]
  comptable2[5*(i-1)+4,] <- values[8:11]
}

#Uniform Design
X <- cbind(rep(1,n),seq(-0.5,0.5,by=1/(n-1)))
if (model == 2){
  X <- cbind(X,X[,2]^2)
}
for (i in 1:length(sigma2vector)){
  values <- CIfunction(n,N,X,x0,sigma2vector[i],results[i,4+model],results[i,5+model],model = model) 
  comptable[5*i,] <- values[1:7]
  comptable2[5*i,] <- values[8:11]
}


#Save results as RDS
saveitems <- list(comptable,comptable2)
saveRDS(saveitems,file = paste("model",model,"EPx0_",x0,"CI.rds",sep=""))



#Format and export table for latex
customround <- function(numbers,digits){ #Function for rounding numbers for latex tables
  numbers[abs(numbers) >= 0.0001] <- round(numbers[abs(numbers) >= 0.0001],digits = digits)
  numbers[abs(numbers) < 0.0001] <- format(signif(numbers[abs(numbers) < 0.0001],digits = 2),scientific = T)
  
  return(numbers)
}
comparisontable <- matrix(NA,nrow = 5*length(sigma2vector),ncol = 7)
for(i in 1:length(sigma2vector)){
  comparisontable[5*(i-1)+1,1:4] <- comptable[5*(i-1)+1,1:4]
  comparisontable[5*(i-1)+2,] <- customround(comptable[5*(i-1)+2,],digits=4)
  comparisontable[5*(i-1)+3,] <- customround(comptable[5*(i-1)+3,],digits=4)
  comparisontable[5*(i-1)+4,] <- customround(comptable[5*(i-1)+4,],digits=4)
  comparisontable[5*(i-1)+5,] <- customround(comptable[5*(i-1)+5,],digits=4)
}

comparisontable <- cbind(rep(c(NA,"$\\xi_{Cl}$","$\\xi_{HL}$","$\\xi_{WX}$","$\\xi_{U}$")),comparisontable)
colnames(comparisontable) <- c(" ","RB (\\%)","ECI CP","ACI CP","PI CP","ECI avg.l","ACI avg.l","PI avg.l")
addtorow <- list()
for (i in 1:length(sigma2vector)){
  addtorow$pos[i] <- list(4*(i-1))
  addtorow$command[i] <- paste("& \\multicolumn{6}{l}{\\underline{$\\sigma^2 = ",comparisontable[(5*(i-1)+1),2],"$,
                               $k_",model+1," = ",comparisontable[(5*(i-1)+1),3],"$,
                               $k_",model+2," = ",comparisontable[(5*(i-1)+1),4],"$,
                               $p = ",comparisontable[(5*(i-1)+1),5],"$}} & \\\\",sep = "")
}
comparisontable <- comparisontable[-seq(1,length(sigma2vector)*5,by=5),]
table1 <- xtable(comparisontable,caption = paste("Table comparing values of interest between designs for Model ",model,
                                                 " extrapolation with $x_0 = ",x0,"$, and $n = ",n,"$.  Abbreviations from left to right stand for:
                                                 Relative Bias, ECI Coverage Percentage, ACI Coverage Percentage, PI Coverage Percentage,
                                                 ECI average CI length, ACI average CI length, PI average CI length.
                                                 ",sep=""),align=c("|c","|c","|c"," c"," c"," c"," c"," c"," c|"), digits = c(0,0,2,2,2,2,4,4,4))
print(table1,file=paste("model",model,"EPx0_",x0,"CIcomp.tex",sep=""),append=T,table.placement="h",caption.placement="bottom",
      hline.after = seq(from=0,to=nrow(table1),by=4),include.rownames = F,sanitize.text.function=function(x){x},
      add.to.row = addtorow)




comparisontable2 <- matrix(NA,nrow = 5*length(sigma2vector),ncol = 4)
for(i in 1:length(sigma2vector)){
  comparisontable2[5*(i-1)+1,1:4] <- comptable2[5*(i-1)+1,1:4]
  comparisontable2[5*(i-1)+2,] <- customround(comptable2[5*(i-1)+2,],digits=4)
  comparisontable2[5*(i-1)+3,] <- customround(comptable2[5*(i-1)+3,],digits=4)
  comparisontable2[5*(i-1)+4,] <- customround(comptable2[5*(i-1)+4,],digits=4)
  comparisontable2[5*(i-1)+5,] <- customround(comptable2[5*(i-1)+5,],digits=4)
}

comparisontable2 <- cbind(rep(c(NA,"$\\xi_{Cl}$","$\\xi_{HL}$","$\\xi_{WX}$","$\\xi_{U}$")),comparisontable2)
colnames(comparisontable2) <- c(" ","$y(x_0)$","$\\overline{\\hat{y}(x_0)}$","sim SE","average ASD")
addtorow <- list()
for (i in 1:length(sigma2vector)){
  addtorow$pos[i] <- list(4*(i-1))
  addtorow$command[i] <- paste("& \\multicolumn{3}{l}{\\underline{$\\sigma^2 = ",comparisontable2[(5*(i-1)+1),2],"$}} & \\\\",sep = "")
}
comparisontable2 <- comparisontable2[-seq(1,length(sigma2vector)*5,by=5),]
table1 <- xtable(comparisontable2,caption = paste("Table comparing values of interest between designs for Model ",model,
                                                 " extrapolation with $x_0 = ",x0,"$, and $n = ",n,"$.  Lists from left to right, the true y value,
                                                 the average estimate of y at $x_0$, the simulated standard error of the estimate of y at $x_0$, and average asymptotic standard deviation.
                                                 ",sep=""),align=c("|c","|c","|c"," c"," c"," c|"), digits = c(0,0,4,4,4,4))
print(table1,file=paste("model",model,"EPx0_",x0,"CIcomp2.tex",sep=""),append=T,table.placement="h",caption.placement="bottom",
      hline.after = seq(from=0,to=nrow(table1),by=4),include.rownames = F,sanitize.text.function=function(x){x},
      add.to.row = addtorow)
