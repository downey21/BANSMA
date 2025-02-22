
# -*- coding: utf-8 -*-

# install.packages("rms", repos = "https://cran.yu.ac.kr/")
# install.packages("RcppArmadillo", repos = "https://cran.yu.ac.kr/")
# install.packages("RcppEigen", repos = "https://cran.yu.ac.kr/")
# install.packages(
#     '/root/Project/BANSMA/BANS_package/BANS_1.0.tar.gz',
#     repos = NULL,
#     type = 'source'
# )

library(BANS)
library(Matrix)

fit.Outcome <- function(y,X,gamma.prob=0.1,lambda,delta,burnin.S,inf.S,model="normal") {
    # Modified from ch.chaingraph function in the BANS package to allow one outcome
    # This function fits regression model for a chain component
    # Input
    # - y :nx1 vector for response (continuous or ordinal)
    # - X : nx(p+1) covariate matrix for exposure and mediators (ncol should be >2 with a exposure and a mediator)
    # - model: select "normal" or "probit"
    stopifnot(is.matrix(X),ncol(X)>=2)
    stopifnot(model%in% c("normal","probit"))
    n = nrow(X)
    S = burnin.S + inf.S
    
    if (model=="normal") {
        dat.C = scale(y)
    }else {
        zlevels = unique(y)
        K.z.levels = length(zlevels)
        z = as.numeric(factor(y,labels=1:K.z.levels))   ### Replace y to z
        gf = c(-Inf,qnorm(1:(K.z.levels-1)/K.z.levels),Inf) ### g_0,....g_K Initialize g parameters
    }
    pmat = 0.2
    
    dat.P = scale(X)
    pP = ncol(dat.P)
    Gamma.list = B.list = matrix(NA,nrow=pP,ncol=inf.S)
    CC =rep(1/lambda,pP)# hyper parameter for b for nonzero gamma
    qmat = matrix(rep(0.1,ncol(dat.P)),ncol=pP)
    
    kappa.list = rep(0,nrow=inf.S)
    ###########################################
    ####### Initialize all regression parameters #########
    ###########################################
    # Gamma (indicators for b)#
    Gamma = matrix(rbinom(pP,size=1,prob=gamma.prob),ncol=pP)
    # B (pxpP b) #
    B = matrix(qmat*Gamma,ncol=pP)
    # kappa (px1 vector) #
    kappa = 1
    
    
    s = 0
    while (s<S) {
        s = s+1
        CCinv = diag((1/CC) * kappa)
        tempB = B
        tempGamma = Gamma
        tempqmat = qmat
        
        if (model=="probit") {
            ###########################################
            ####### Sample Y|Z,B,kappa  (Hoff and Niu p.212)#########
            ###########################################
            ey = dat.P %*% c(tempB) #### nx1 expected response for the latent
            a = gf[z] # nx1 lower bound
            b = gf[z+1] # n
            u = sapply(1:n,function(i) max(min(runif(1,pnorm((a[i]-ey[i])/kappa),pnorm((b[i]-ey[i])/kappa)),0.999),0.0001))
            y = ey + kappa * qnorm(u) #### Sampled y
            # print(cbind(ey,a,b,u,y,z))
            ###########################################
            ####### Sample g_k|y,z, other g  (Hoff and Niu p.213) and Albert and Chib, 1993 #########
            ###########################################
            for (k in 2:K.z.levels) {
                a = max(max(y[z==(k-1)]),gf[k-1])
                b = max(min(min(y[z==k]),gf[k+1]),a+0.0001)
                gk = runif(1,a,b)
                #stopifnot(sum(is.nan(gk))==0)
                gf[k] =gk
            }
            #  print(a)
            # print(b)
            #print(gf)
            
            dat.C = y
        }
        
        up = updateDirected(y=dat.C,X=dat.P,lambda=lambda,delta=delta,CCinv=CCinv,B=tempB,Gamma=tempGamma,kappa=kappa,qmat=tempqmat,v.no.ue =0,v.aa=0,no.tau=1)
        if (up$is.move){
            Gamma[1,] = up$Gamma
            kappa = up$kappa
            B[1,] = up$B
        }
        ### Store values
        if (s>burnin.S) {
            ss = s - burnin.S
            Gamma.list[,ss] = Gamma; B.list[,ss]=B
            kappa.list[ss] = kappa
        }
    }#while (s<S)
    return(list(Gamma=Gamma.list,B=B.list,kappa=kappa.list))
}

fit.Mediator <- function(v.ch,v.pa,Y,eta.prob=0.1,gamma.prob=0.1,lambda,delta,burnin.S,inf.S) {
    # Modified from ch.chaingraph function in the BANS package
    # This function fits regression model for a chain component
    # Input
    # - v.ch : indices for the target chain component
    # - v.pa : indices for parents set
    # - Y : nxp data matrix
    stopifnot(length(v.ch)>2)
    S = burnin.S + inf.S
    dat.C = scale(Y[,v.ch,drop=F])
    n = nrow(dat.C)
    p = ncol(dat.C)
    pmat = matrix(0.2,ncol(dat.C),ncol(dat.C))
    diag(pmat) = 0
    
    if (is.null(v.pa)) {
        Gamma.list = NULL
        B.list=NULL
    }else{
        dat.P = scale(Y[,v.pa,drop=F])
        pP = ncol(dat.P)
        Gamma.list = B.list = array(0,dim=c(p,pP,inf.S))
        CC = matrix(1/lambda,ncol=pP,nrow=p)# hyper parameter for b for nonzero gamma
        qmat = matrix(0.1,ncol(dat.C),ncol(dat.P))
    }
    eta.list =A.list= array(0,dim=c(p,p,inf.S))
    kappa.list = matrix(0,nrow=inf.S,ncol=p)
    ###########################################
    ####### Initialize all parameters #########
    ###########################################
    # eta (indicators for alpha) #
    w.upper = which(upper.tri(diag(p)))
    eta = matrix(0,p,p)
    eta[w.upper] = rbinom(length(w.upper),size=1,prob=eta.prob)
    eta = eta + tCpp(eta)
    diag(eta) = 0
    # Gamma (indicators for b)#
    if (!is.null(v.pa)) {Gamma = matrix(rbinom(p*pP,size=1,prob=gamma.prob),nrow=p,ncol=pP)}
    # A (pxp alpha)$
    A = pmat*eta
    # B (pxpP b) #
    if (!is.null(v.pa)) {B = qmat*Gamma}
    # kappa (px1 vector) #
    kappa = rep(1,p)
    
    s = 0
    
    while (s<S) {
        s = s+1
        if (s%%100==0) cat("no. of samples=",s,"\n")
        for (v in sample(1:p)) { # in random order
            # Update eta, A, kappa
            if (!is.null(v.pa)) {
                tempDat = dat.C - tcrossprodCpp(dat.P,B)
                no.de = rowSums(B!=0)
                bCb = sapply(1:p,function(v)sum(B[v,]^2/CC[v,]))
            } else{tempDat=dat.C
                no.de=bCb = rep(0,p)
            }
            up = updateUndirected(v=v,dat=tempDat,lambda=lambda,delta=delta,Alpha=A,eta=eta,kappa=kappa,pmat=pmat,no.de=no.de,bCb=bCb,no.tau=p)
            if (up$is.move) {
                eta=up$eta
                kappa = up$kappa
                A =up$Alpha
            }
            if (!is.null(v.pa)) {
                # Update Gamma, B, kappa
                v.ne = which(A[v,]!=0)
                v.ne.l = length(v.ne)
                vv.ne = c(v,v.ne)
                l.cl = length(vv.ne)
                if (v.ne.l>0) {
                    alpha = A[v,v.ne,drop=F]
                    y = c(dat.C[,v] - prod(dat.C[,v.ne,drop=F],alpha))
                    X = sapply(alpha,function(k) -k*dat.P,simplify=F)
                    X = cbind(dat.P,do.call(cbind,X))
                }else{
                    y = dat.C[,v]
                    X = dat.P
                    alpha=0
                }
                CCinv = (1/CC[vv.ne,,drop=F]) * kappa[vv.ne]
                ##### Modified for BANSMA to cover one parent.
                if (!all(dim(CCinv)==c(1,1))){
                    if (ncol(CCinv)==1) {
                        CCinv = diag(c(CCinv))
                    } else{
                        CCinv = as.matrix(bdiag(sapply(1:l.cl,function(x)diag(CCinv[x,]),simplify=F)))
                    }
                }
                ###############################################
                tempB = B[vv.ne,,drop=F]
                tempGamma = Gamma[vv.ne,,drop=F]
                tempqmat = qmat[vv.ne,,drop=F]
                up = updateDirected(y=y,X=X,lambda=lambda,delta=delta,CCinv=CCinv,B=tempB,Gamma=tempGamma,kappa=kappa[v],qmat=tempqmat,v.no.ue =v.ne.l,v.aa=sum(alpha^2),no.tau=p)
                if (up$is.move){
                    Gamma[vv.ne,] = up$Gamma
                    kappa[v] = up$kappa
                    B[vv.ne,] = up$B
                }
            } #(!is.null(v.pa))
        }# for (v in sample(1:p))
        ### Store values
        if (s>burnin.S) {
            ss = s - burnin.S
            if (!is.null(v.pa)) {Gamma.list[,,ss] = Gamma; B.list[,,ss]=B}
            eta.list[,,ss] = eta
            A.list[,,ss] = A
            kappa.list[ss,] = kappa
        }
    }#while (s<S)
    return(list(Gamma=Gamma.list,eta=eta.list,A=A.list,B=B.list,kappa=kappa.list))
}

getIndEffect <- function(palist,chlist,fit.med,fit.out,ci=c(0.025,0.975)) {
    #### - palist : parent list (used for fitting with the variable order)
    #### - chlist : multilayer nodes
    #### - fit.med: q-1 mediator model. (fit.med$Gamma includes no. mediator x (no. exposure + no. covariate)) (from second layer, does not allow GGM in the first layer)
    #### - fit.out: outcome model
    #### - ci: length 2 vector that includes upper and lower credible interval to calculate
    
    ov = do.call(c,chlist)
    layerinfo = do.call(c,sapply(1:length(chlist),function(k)rep(k,length(chlist[[k]]))))
    
    q = length(fit.med)+1   ### because the first layer model is not included
    p = length(ov) - 1
    
    inf.S = ncol(fit.out$B)
    
    
    if (q-1 == 1) {
        A.array = c(sapply(1:(q-1),function(a) chlist[[a]],simplify=T))
    }else {
        Aarray = do.call(c,sapply(1:(q-1),function(a) chlist[[a]],simplify=F))
    }
    IEDout = IECout=IEout = numeric(0)
    for (a in Aarray) {
        a.layer= layerinfo[which(ov==a)]
        m.layer= (a.layer+1):q
        IED = IEC = IE = numeric(0)
        for (m in m.layer) {
            out = v.m.getIndEffect(v=a,m,palist,chlist,fit.med[[m-1]],fit.out,ci=c(0.025,0.975))
            IED = rbind(IED,data.frame(rep(a,length(chlist[[m]])),out$v.med,out$IED))
            IEC = rbind(IEC,data.frame(rep(a,length(chlist[[m]])),out$v.med,out$IEC))
            IE = rbind(IE,data.frame(rep(a,length(chlist[[m]])),out$v.med,out$IE))
        }
        IEDout = rbind(IEDout,IED)
        IECout = rbind(IECout,IEC)
        IEout = rbind(IEout,IE)
    }
    colnames(IEDout) = c("Exposure","Mediator","IED","CI_lower","CI_upper","PIP")
    colnames(IECout) = c("Exposure","Mediator","IEC","CI_lower","CI_upper","PIP")
    colnames(IEout) = c("Exposure","Mediator","IE","CI_lower","CI_upper","PIP")
    return(list(IED=IEDout,IEC=IECout,IE=IEout))
}

v.m.getIndEffect<- function(v,m,palist,chlist,fit.med,fit.out,ci=c(0.025,0.975)){
    # Input
    #### - v : index for exposure
    #### - m : position of layer of the whole mlGGM for fit.med
    #### - palist : parent list (used for fitting with the variable order)
    #### - chlist : multilayer nodes
    #### - fit.med: mediator model. (fit.med$Gamma includes no. mediator x (no. exposure + no. covariate))
    #### - fit.out: outcome model
    #### - ci: length 2 vector that includes upper and lower credible interval to calculate
   
   
    v.med = chlist[[m]] # nodes in the mediator layer
    p.med = palist[[m]] # parents of the mediator layer m
    p.out = palist[[length(palist)]] # covariates (parents) of outcome model
    pos.v.med = which(p.med==v) ## position of v in parent set of mediator model (v is in the parent set)
    
    
    stopifnot()
    d.med = dim(fit.med$Gamma)
    d.out = dim(fit.out$Gamma)
    stopifnot(d.med[3]==d.out[2])
    
    inf.S = d.med[3]
    ######### IED ############
    outIED = outIEC = outIE = matrix(NA,nrow=length(v.med),ncol=4)
    colnames(outIED) = c("IED","CI_lower","CI_upper","PIP")
    colnames(outIEC) = c("IEC","CI_lower","CI_upper","PIP")
    colnames(outIE) = c("IE","CI_lower","CI_upper","PIP")
    k = 0
    for (j in v.med) {
        k = k+1
        pos.j.med = which(v.med==j) ## position of mediator j in mediator model
        pos.j.out = which(p.out==j) ## position of mediator j in outcome model
        IED = sapply(1:inf.S,function(s) fit.out$B[pos.j.out,s] * fit.med$B[pos.j.med,pos.v.med,s] )  #IED
        IEC = sapply(1:inf.S,function(s) fit.out$B[pos.j.out,s] * (fit.med$A[pos.j.med,-pos.j.med,s]%*%fit.med$B[-pos.j.med,pos.v.med,s]))
        IE = IED - IEC
        
        outIED[k,1] = mean(IED)
        outIED[k,2:3] = quantile(IED,probs=ci)
        outIED[k,4] = mean(IED!=0)
        
        outIEC[k,1] = mean(IEC)
        outIEC[k,2:3] = quantile(IEC,probs=ci)
        outIEC[k,4] = mean(IEC!=0)
        
        outIE[k,1] = mean(IE)
        outIE[k,2:3] = quantile(IE,probs=ci)
        outIE[k,4] = mean(IE!=0)
    }
    return(list(IED=outIED,IEC=outIEC,IE=outIE,v.med=v.med))
}