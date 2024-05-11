
# -*- coding: utf-8 -*-

rm(list=ls())

r = 1

p = 50 # no. of mediators + exposure
q = 4 # no. of layers of X
pE = 0.15
n = 60

library(igraph)
library(RColorBrewer)
library(Matrix)
library(bama)
library(parallel)
source("/root/Project/BANSMA/code/BANSMA.R")
source("/root/Project/BANSMA/code/BANSMA_simulation.R")

path_result <- "/root/Project/BANSMA/simulation/p50q4n60/result"

load(paste(path_result,"/ER_Setting_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
Y=dat$Y
load(paste(path_result,"/ER_BANSMA_p",p,"q",q,"n",n,"r",r,".RData",sep=""))

##########
# BANSMA-probit outcome
set.seed(123)
lambda = 5
delta = 2
burnin.S = 20000
inf.S = 1000

t= q+1
y = Y[,chlist[[t]]]
X = Y[,palist[[t]]]

z = ifelse(y>median(y),1,0)
fitProbit.out = fit.Outcome(y=z,X=X,lambda=lambda,delta=delta,burnin.S=burnin.S,inf.S=inf.S,model="probit")
# fit.out = fit.Outcome(y=y,X=X,lambda=lambda,delta=delta,burnin.S=burnin.S,inf.S=inf.S,model="normal")

bansmaProbit.out = getIndEffect(palist=palist,chlist=chlist,fit.med=fit.med,fit.out=fitProbit.out)
# bansma.out = getIndEffect(palist=palist,chlist=chlist,fit.med=fit.med,fit.out=fit.out)

save(fit.med,fitProbit.out,bansmaProbit.out,file=paste(path_result,"/ER_BANSMAProbit_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
# save(fit.med,fit.out,fitProbit.out,bansma.out,bansmaProbit.out,file=paste(path_result,"/ER_BANSMAProbit_p",p,"q",q,"n",n,"r",r,".RData",sep=""))

