
# -*- coding: utf-8 -*-

rm(list=ls())

r = 1

p = 50 # no. of mediators + exposure
q = 4 # no. of layers of X
pE = 0.15
n = 200

library(igraph)
library(RColorBrewer)
library(Matrix)
library(bama)
library(mediation)
library(parallel)
source("/root/Project/BANSMA/code/BANSMA.R")
source("/root/Project/BANSMA/code/BANSMA_simulation.R")

path_result <- "/root/Project/BANSMA/simulation/p50q4n200/result"

##########
# Graph generation
set.seed(77109)
sim  = simul.Chaingraph(p,q,pE,scale.free=T)

simAug = sim
A = matrix(0,p+1,p+1)
A[1:p,1:p] = sim$A
dd = rowSums(A)
w = sample(size=p*pE,x=which(dd>=2))
A[p+1,w] = 1

simAug$A = A
simAug$ch.addr[[q+1]] = p+1
ov = do.call(c,simAug$ch.addr)
layerinfo = as.character(do.call(c,sapply(1:(q+1),function(k)rep(k,length(simAug$ch.addr[[k]])))))

##########
# Draw graph
adjA = t(simAug$A)
und.w = which(adjA + t(adjA)==2,arr.ind=T)
und.w = unique(t(apply(und.w,1,sort)))
dir.w = which(adjA==1 & t(adjA)==0,arr.ind=T)
E = rbind(und.w,dir.w)

ov = do.call(c,simAug$ch.addr)
layerinfo = as.character(do.call(c,sapply(1:(q+1),function(k)rep(k,length(simAug$ch.addr[[k]])))))
w = which(!ov %in% unique(c(E)))
if (length(w)>0) {
    ov = ov[-w]
    layerinfo=layerinfo[-w]
}

E = data.frame(E,type=c(rep("undirected",nrow(und.w)),rep("directed",nrow(dir.w))),weight=c(rep(900,nrow(und.w)),rep(1000,nrow(dir.w))))
w= apply(E[,1:2],1,function(x)(p+1)%in%x)
E[w,"weight"] = E[w,"weight"] + 100
G = graph_from_data_frame(E)
G = set.vertex.attribute(G,"layer",as.character(ov),layerinfo)

E(G)$arrow.size= 0.5
E(G)$width=ifelse(E(G)$type=="undirected",6,2)
E(G)$color=ifelse(E(G)$type=="undirected","grey","black")

pdf("/root/Project/BANSMA/simulation/p50q4n200/figure/BA_p50q4.pdf")
plot(G,layout=layout_with_kk(G,weight=E(G)$weight),vertex.color=V(G)$layer,vertex.size=10,edge.width=E(G)$width,edge.color=E(G)$color)
dev.off()

##########
# Data generation 
set.seed(101 * (r+77177))
dat = rmvnorm.Chaingraph(A=simAug$A,ch.addr=simAug$ch.addr,n=n)
Y = dat$Y

chlist = simAug$ch.addr

palist = vector("list",q+1)
for (i in 1:q+1) {
    for (j in 1:(i-1)) {
        palist[[i]] = c(chlist[[j]],palist[[i]])
    }
}

layerinfo = do.call(c,sapply(1:(q+1),function(k)rep(k,length(simAug$ch.addr[[k]]))))
ov = do.call(c,chlist)

##########
# True
v.y = p+1
K = dat$K[-v.y,-v.y]
B = dat$B[-v.y,-v.y]
theta = dat$B[v.y,-v.y]

if (q-1 == 1) {
    Aarray = c(sapply(1:(q-1),function(a) chlist[[a]],simplify=T))
}else {
    Aarray = do.call(c,sapply(1:(q-1),function(a) chlist[[a]],simplify=T))
}
true.out = numeric(0)
for (a in Aarray) {
    a.layer= layerinfo[which(ov==a)]
    if ((a.layer+1)==q) {
        m= c(sapply((a.layer+1):q,function(a) chlist[[a]],simplify=T))
    }else {
        m = do.call(c,sapply((a.layer+1):q,function(a) chlist[[a]],simplify=T))
    }
    IED = theta[m] * B[m,a]
    IEC = sapply(m,function(j) theta[j]*sum(K[j,setdiff(m,j)]* B[setdiff(m,j),a]))
    IE = IED - IEC
    true.out = rbind(true.out,data.frame(rep(a,length(m)),m,IED=IED,IEC=IEC,IE=IE))
}
colnames(true.out) = c("Exposure","Mediator","IED","IEC","IE")

save(sim,dat,chlist,palist,true.out,file=paste(path_result,"/BA_Setting_p",p,"q",q,"n",n,"r",r,".RData",sep=""))

##########
# BANSMA
lambda = 5
delta = 2
burnin.S = 5000
inf.S = 1000

# Mediator Model Fitting
set.seed(123, "L'Ecuyer")
fit.med = parallel::mclapply(2:q,function(t) ch.chaingraph(v.ch=chlist[[t]],v.pa=palist[[t]],Y=Y,lambda=lambda,delta=delta,burnin.S=burnin.S,inf.S=inf.S), mc.cores = 3)

# Outcome Model Fitting
t = q+1
y = Y[,chlist[[t]]]
X = Y[,palist[[t]]]
fit.out= fit.Outcome(y=y,X=X,lambda=lambda,delta=delta,burnin.S=burnin.S,inf.S=inf.S)

# Indirect effects
bansma.out = getIndEffect(palist,chlist,fit.med,fit.out)

save(fit.med,fit.out,bansma.out,file=paste(path_result,"/BA_BANSMA_p",p,"q",q,"n",n,"r",r,".RData",sep=""))

##########
# BAMA
set.seed(123)
y = Y[,chlist[[q+1]]]
if (q-1 == 1) {
    Aarray = c(sapply(1:(q-1),function(a) chlist[[a]],simplify=T))
}else {
    Aarray = do.call(c,sapply(1:(q-1),function(a) chlist[[a]],simplify=T))
}

bama.out = numeric(0)
for (a in Aarray) {
    a.layer= layerinfo[which(ov==a)]
    A = Y[,a]
    if (a.layer!=1) {
        if (a.layer==1) {
        covs = c(sapply(1:a.layer,function(a) chlist[[a]],simplify=T))
        }else {
            covs= do.call(c,sapply(1:a.layer,function(a) chlist[[a]],simplify=T))
        }
        C1 = Y[,setdiff(covs,a),drop=F]
    } else {
        C1<- matrix(1, nrow(Y), 1)
    }
    if ((a.layer+1)==q) {
        m= c(sapply((a.layer+1):q,function(a) chlist[[a]],simplify=T))
    } else {
        m = do.call(c,sapply((a.layer+1):q,function(a) chlist[[a]],simplify=T))
    }

    M = Y[,m]
    beta.m  <- rep(0, ncol(M))
    alpha.a <- rep(0, ncol(M))
    C2 = C1
    out <- bama(Y = y, A = A, M = M, C1 = C1, C2 = C2, method = "BSLMM",
                burnin = 100, ndraws = 110, weights = NULL, inits = NULL, 
                control = list(k = 2, lm0 = 1e-04, lm1 = 1, lma1 = 1, l = 1))
    out <- summary(out)
    bama.out= rbind(bama.out,data.frame(rep(a,length(m)),m,out$estimate,out$ci.lower,out$ci.upper,out$pip))
}
colnames(bama.out) =c("Exposure","Mediator","Estimate","CI_lower","CI_upper","PIP")

save(bama.out,file=paste(path_result,"/BA_BAMA_p",p,"q",q,"n",n,"r",r,".RData",sep=""))

##########
# Mediation
set.seed(123)
y = Y[,chlist[[q+1]]]
if (q-1 == 1) {
    Aarray = c(sapply(1:(q-1),function(a) chlist[[a]],simplify=T))
} else {
    Aarray = do.call(c,sapply(1:(q-1),function(a) chlist[[a]],simplify=T))
}

med.single.out = numeric(0)
for (a in Aarray) {
    a.layer = layerinfo[which(ov==a)]
    A = Y[,a]
    if (a.layer!=1) {
        if (a.layer==1) {
            covs = c(sapply(1:a.layer,function(a) chlist[[a]],simplify=T))
        } else {
            covs= do.call(c,sapply(1:a.layer,function(a) chlist[[a]],simplify=T))
        }
        cc=setdiff(covs,a)
        C1 = Y[,cc,drop=F]
        colnames(C1) = paste("C",cc,sep="")
    } else {
        C1 <- matrix(1, nrow(Y), 1)
    }
    if ((a.layer+1)==q) {
        m = c(sapply((a.layer+1):q,function(a) chlist[[a]],simplify=T))
    } else {
        m = do.call(c,sapply((a.layer+1):q,function(a) chlist[[a]],simplify=T))
    }

    M = Y[,m]
    colnames(M) = paste("M",m,sep="")
    if (a.layer!=1) {
        dat = data.frame(y=y,A=A,M,C1)
    } else {
        dat = data.frame(y=y,A=A,M)
    }
    out.layer = data.frame(matrix(NA,nrow=length(m),ncol=5))
    out.layer[,1] = m
    colnames(out.layer) = c("Mediator","estimate","ci.lower","ci.upper","p-value")
    k = 0
        for (j in m) {
            k = k+1
            if (a.layer!=1) {
                nn = paste("C",cc,collapse="+",sep="")
                form.med = as.formula(paste("M",j,"~A+",nn,sep=""))
                form.out = as.formula(paste("y~A+M",j,"+",nn,sep=""))
            } else {
                form.med = as.formula(paste("M",j,"~A",sep=""))
                form.out = as.formula(paste("y~A+M",j,sep=""))
            }
            fit.med.single = lm(form.med,data=dat)
            fit.out.single = lm(form.out,data=dat)
            out <- mediate(fit.med.single, fit.out.single, treat = "A", mediator = paste("M",j,sep=""), robustSE=TRUE, sims = 100)
            out.layer[k,2] = out$d0
            out.layer[k,3:4] = out$d0.ci
            out.layer[k,5] = out$d0.p
        }
    med.single.out = rbind(med.single.out,data.frame(rep(a,length(m)),out.layer))
}
colnames(med.single.out) =c("Exposure","Mediator","Estimate","CI_lower","CI_upper","p-value")

save(med.single.out,file=paste(path_result,"/BA_Mediate_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
