#Source code to run Extrapolation.R for the different models and x0's

#Values of n for each scenario:
#Model 1: x0 = 0.75 -> n = 30,   x0 = 2.5 -> n = 30
#Model 2: x0 = 0.75 -> n = 28,   x0 = 2.5 -> n = 49

library(xtable)
setwd('C:/Users/dburm/Desktop/School/MATH 4F90/')
source('EPfunctions.R')

n <- 30
N <- 50000
sigma2vector <- c(0.01,0.05,0.1,0.2,0.5)

#MODEL 1
model <- 1
x0 <- 0.75
source('Extrapolation.R',echo=T)

x0 <- 2.5
source('Extrapolation.R',echo=T)

#MODEL 2
model <- 2
n <- 28
x0 <- 0.75
source('Extrapolation.R',echo=T)

n <- 49
x0 <- 2.5
source('Extrapolation.R',echo=T)
