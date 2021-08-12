#Source code for running CIcomparison.R for all models and x0's

library(xtable)
setwd('C:/Users/dburm/Desktop/School/MATH 4F90/')
source('EPfunctions.R')
source('CIfunction.R')

N <- 50000
sigma2vector <- c(0.01,0.05,0.1,0.2,0.5)


model <- 1
n <- 30
x0 <- 0.75
rdsfile <- readRDS('rds files/model1EPx0_0.75results.rds')
results <- rdsfile$results
source('CIcomparison.R',echo=T)


x0 <- 2.5
rdsfile <- readRDS('rds files/model1EPx0_2.5results.rds')
results <- rdsfile$results
source('CIcomparison.R',echo=T)

model <- 2
n <- 28
x0 <- 0.75
rdsfile <- readRDS('rds files/model2EPx0_0.75results.rds')
results <- rdsfile$results
source('CIcomparison.R',echo=T)

n <- 49
x0 <- 2.5
rdsfile <- readRDS('rds files/model2EPx0_2.5results.rds')
results <- rdsfile$results
source('CIcomparison.R',echo=T)