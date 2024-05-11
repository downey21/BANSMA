
# -*- coding: utf-8 -*-

library(Matrix)
library(mvtnorm)

scaledMat <- function(x){
    newx=x/sqrt(diag(x) %*% t(diag(x)))
    return(newx)
}

add.alpha <- function(col, alpha=1){
    if(missing(col))
    stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2,
    function(x)
    rgb(x[1], x[2], x[3], alpha=alpha))
}

simul.Chaingraph <- function(p,q,pE,scale.free=F) {
    # ----------------------------------------
    # generate ER model for the skeleton of the chain graph
    # p : number of vertices
    # q : number of chain components
    # pE : connection probability of ER model
    # ----------------------------------------
    
    stopifnot(pE>0,pE<1,q>1)
    addr = cbind(1:p,1:p%%q+1)
    ch.addr = sapply(1:q,function(x)addr[addr[,2]==x,1],simplify=F)
    ov = do.call(c,ch.addr)
    
    ### edges in each chain
    A = matrix(0,p,p)
    for (i in 1:q) {
        qp = length(ch.addr[[i]])
        if (!scale.free) {
            Atmp = matrix(0,qp,qp)
            w = which(lower.tri(Atmp))
            Atmp[w] = rbinom(n=length(w),size=1,prob=pE)
            A[ch.addr[[i]],ch.addr[[i]]] = Atmp + t(Atmp)
        }else {
            g <- barabasi.game(qp, power = 1.2, m = NULL, out.dist = NULL, out.seq = NULL,
            out.pref = FALSE, zero.appeal = 1, directed = FALSE,
            algorithm ="psumtree", start.graph = NULL)
            Atmp = as_adjacency_matrix(g)
            A[ch.addr[[i]],ch.addr[[i]]] = as.matrix(Atmp)
        }
    }
    
    # edges between chains #
    for (i in 1:(q-1)){
        qp = length(ch.addr[[i+1]])
        qpp = length(ch.addr[[i]])
        tmp = matrix(rbinom(n=qp*qpp,size=1,prob=pE/2) ,ncol=qpp,nrow=qp)
        A[ch.addr[[i+1]],ch.addr[[i]]] = tmp
    }
    # Replace undirected edges to directed edges
    return(list(A=A,ch.addr=ch.addr))
}

simul.MedModel <- function(p,q,pE) {
# ----------------------------------------
# generate ER model for the skeleton of the chain graph
# p : number of vertices for mediators
# q : number of chain components for mediators
# pE : connection probability of ER model across all exposure, mediators and outcome
# ----------------------------------------
	
	stopifnot(pE>0,pE<1,q>=1)
	addr = rbind(c(1,1),cbind(1:p,1:p%%q+1)+1)
    q = q+1
    
	ch.addr = sapply(1:q,function(x)addr[addr[,2]==x,1],simplify=F)
	ov = do.call(c,ch.addr)
### edges in each chain
	A = matrix(0,p+1,p+1)
	for (i in 1:q) {
		qp = length(ch.addr[[i]])
		Atmp = matrix(0,qp,qp)
		w = which(lower.tri(Atmp))
		Atmp[w] = rbinom(n=length(w),size=1,prob=pE) 
		A[ch.addr[[i]],ch.addr[[i]]] = Atmp + t(Atmp)
	}
    
# edges between chains #
	for (i in 1:(q-1)){
		qp = length(ch.addr[[i+1]])
		qpp = length(ch.addr[[i]])
		tmp = matrix(rbinom(n=qp*qpp,size=1,prob=pE) ,ncol=qpp,nrow=qp)
		A[ch.addr[[i+1]],ch.addr[[i]]] = tmp
	}
# Replace undirected edges to directed edges
	return(list(A=A,ch.addr=ch.addr))
}

rmvnorm.Chaingraph <- function(A,ch.addr,n) {
# ----------------------------------------
# generate data given a chain graph
# A : adjacency matrix 
# ch.addr : a list for chain components
# n : sample size
# ----------------------------------------
	p = ncol(A)
	ov = do.call(c,ch.addr)
	stopifnot(sum(1:p%in%ov)==p)
	
	und.edges = which(A!=0 & t(A)!=0,arr.ind=T)
	und.edges = und.edges[!duplicated(t(apply(und.edges,1,sort))),]
	dir.edges = which(A!=0 & t(A)==0,arr.ind=T)
	
# Set undirected part
	K = matrix(0,p,p)
	K[und.edges] = sample(c(-1,1),nrow(und.edges),replace=T)*runif(n=nrow(und.edges),min=1,max=1.5)
	K = K + t(K)
	diag(K) = colSums(abs(K)) + 0.1
	R = scaledMat(K)
# Set
    B = matrix(0,p,p)
	B[dir.edges] = sample(c(-1,1),nrow(dir.edges),replace=T)*runif(n=nrow(dir.edges),min=1,max=1.5)
# Generate data 
	Ip = diag(p)
	Omega = t(Ip-B)%*%K%*%(Ip-B)
	Sigma = solve(Omega)
	Y = rmvnorm(n,mean=rep(0,p),sigma=Sigma)
	return(list(K=K,B=B,Y=Y,R=R))
}
