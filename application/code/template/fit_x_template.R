
# -*- coding: utf-8 -*-

rm(list=ls())

# 1 - 12
pp = 1

library(igraph)
library(Matrix)
library(BANS)
library(parallel)
source("/root/Project/BANSMA/code/BANSMA.R")

path_data <- "/root/Project/BANSMA/application/data"
path_result <- "/root/Project/BANSMA/application/result"

CNA = read.csv(file=paste0(path_data,"/CNA_BO.csv"),as.is=T)
mRNA = read.csv(file=paste0(path_data,"/mRNA_BO.csv"),as.is=T)
RPPA = read.csv(file=paste0(path_data,"/RPPA_BO.csv"),as.is=T)

load(paste0(path_data,"/StringV10_ppi_pathway.RData"))
padat = ppi$pathwaydat
padat[,2] = gsub("-","",padat[,2])
pname = unique(padat[,1])[pp]
genes = unlist(strsplit(padat[padat[,1] == pname,3],split=", "))

RPPA = RPPA[,colnames(RPPA)%in%genes]
mRNA = mRNA[,colnames(mRNA)%in%genes]
CNA = CNA[,colnames(CNA)%in%genes]

X = cbind(CNA,mRNA,RPPA)
colnames(X) = 
    c(
        paste("CNA_",colnames(CNA),sep=""),
        paste("mRNA_",colnames(mRNA),sep=""),
        paste("RPPA_",colnames(RPPA),sep="")
    )

X = X[rowSums(is.na(X))==0,]

p1 = ncol(CNA)
p2 = ncol(mRNA)
p3 = ncol(RPPA)
addr = cumsum(c(p1,p2,p3))

q = 3

chlist = vector("list",q)
chlist[[1]] = 1:addr[1]
chlist[[2]] = (addr[1]+1):addr[2]
chlist[[3]] = (addr[2]+1):addr[3]

palist = vector("list",q)
palist[[2]] = 1:addr[1]
palist[[3]] = 1:addr[2]

##########
# BANSMA
lambda = 5
delta = 2
burnin.S = 20000
inf.S = 10000

# Mediator Model Fitting
set.seed(1234, "L'Ecuyer")
fit = parallel::mclapply(1:q,function(t) fit.Mediator(v.ch=chlist[[t]],v.pa=palist[[t]],Y=X,lambda=lambda,delta=delta,burnin.S=burnin.S,inf.S=inf.S,eta.prob=0.5,gamma.prob=0.5), mc.cores = 3)

save(X,fit,chlist,palist,p1,p2,p3,inf.S,pname,file=paste(path_result,"/fit_x_pp",pp,".RData",sep=""))
