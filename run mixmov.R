# rm(list=ls(all=TRUE))
library('MCMCpack')
library('Rcpp')
set.seed(2)

setwd('U:\\GIT_models\\mixture_movem')
source('mixmov_function.R')
source('mixmov_gibbs.R')
sourceCpp('aux1.cpp')
dat=read.csv('fake data.csv',as.is=T)

#prior
gamma1=0.1
alpha=0.1

#prepare for gibbs
ngibbs=10000
nburn=ngibbs/2
nmaxclust=10
model1=mixture_movement(dat=dat,gamma1=gamma1,alpha=alpha,
                        ngibbs=ngibbs,nmaxclust=nmaxclust,
                        nburn=nburn)


