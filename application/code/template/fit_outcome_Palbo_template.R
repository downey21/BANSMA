
# -*- coding: utf-8 -*-

rm(list=ls())

# 1 - 12
pp = 1

library(igraph)
library(Matrix)
library(BANS)
source("/root/Project/BANSMA/code/BANSMA.R")

path_data <- "/root/Project/BANSMA/application/data"
path_result <- "/root/Project/BANSMA/application/result"

CNA = read.csv(file=paste0(path_data,"/CNA_BO.csv"),as.is=T)
mRNA = read.csv(file=paste0(path_data,"/mRNA_BO.csv"),as.is=T)
RPPA = read.csv(file=paste0(path_data,"/RPPA_BO.csv"),as.is=T)
Drug = read.csv(file=paste0(path_data,"/GDSC2_BO.csv"),as.is=T)
Drug = Drug[,c("X","Palbociclib")]
Drug = Drug[rowSums(is.na(Drug))==0,]

rnames= Drug[,1]
CNA = CNA[match(rnames,CNA[,1]),]
mRNA = mRNA[match(rnames,mRNA[,1]),]
RPPA = RPPA[match(rnames,RPPA[,1]),]

load(paste0(path_data,"/StringV10_ppi_pathway.RData"))
padat = ppi$pathwaydat
padat[,2] = gsub("-","",padat[,2])
pname = unique(padat[,1])[pp]
genes = unlist(strsplit(padat[padat[,1] == pname,3],split=", "))

RPPA = RPPA[,colnames(RPPA)%in%genes]
mRNA = mRNA[,colnames(mRNA)%in%genes]
CNA = CNA[,colnames(CNA)%in%genes]

X = cbind(CNA,mRNA,RPPA,Drug[,-1])
colnames(X) = 
    c(
        paste("CNA_",colnames(CNA),sep=""),
        paste("mRNA_",colnames(mRNA),sep=""),
        paste("RPPA_",colnames(RPPA),sep=""),
        paste("Drug_",colnames(Drug)[-1],sep="")
    )

X = X[rowSums(is.na(X))==0,]

p1 = ncol(CNA)
p2 = ncol(mRNA)
p3 = ncol(RPPA)
p4 = ncol(Drug)-1
addr = cumsum(c(p1,p2,p3,p4))

q = 4

chlist = vector("list",q)
chlist[[1]] = 1:addr[1]
chlist[[2]] = (addr[1]+1):addr[2]
chlist[[3]] = (addr[2]+1):addr[3]
chlist[[4]] = (addr[3]+1):addr[4]

palist = vector("list",q)
palist[[2]] = 1:addr[1]
palist[[3]] = 1:addr[2]
palist[[4]] = 1:addr[3]

##########
# BANSMA
lambda = 5
delta = 2
burnin.S = 20000
inf.S = 10000

# Mediator Model Fitting
t = 4
set.seed(1234)
X <- as.matrix(X)
fit = fit.Outcome(y=X[,chlist[[t]]],X=X[,palist[[t]]],lambda=lambda,delta=delta,burnin.S=burnin.S,inf.S=inf.S,gamma.prob=0.5)

save(X,fit,chlist,palist,p1,p2,p3,p4,inf.S,pname,file=paste(path_result,"/fit_outcome_Palbo_pp",pp,".RData",sep=""))
