
# -*- coding: utf-8 -*-

rm(list=ls())

path_result <- "/root/Project/BANSMA/application/result"
path_result_analysis <- "/root/Project/BANSMA/application/result_analysis"

pparray = (1:12)

pnodeattr = numeric(0)
pedgeattr = numeric(0)

for (pp in pparray){
    print(pp)
    
    load(paste(path_result,"/fit_x_pp",pp,".RData",sep=""))
    p = ncol(X)
    q = length(chlist)
    
    B = matrix(0,p,p)
    cB = numeric(0)
    for (j in 1:q) {
        tmp = apply(fit[[j]]$eta,c(1,2),mean)
        B[chlist[[j]],chlist[[j]]] = tmp
        cB = c(cB,tmp[lower.tri(tmp)])
        if (j!=1) {
            tmp = c(apply(fit[[j]]$Gamma,c(1,2),mean))
            B[as.matrix(expand.grid(chlist[[j]],palist[[j]]))] = tmp
            cB = c(cB,tmp)
        }
    }
    o = order(cB,decreasing=T)
    ii = max(which(cumsum(1-cB[o]) * (1/1:length(cB)) < 0.1))  # FDR at 0.1
    cut = cB[o][ii]
    
    nodes = colnames(X)
    matB = matrix(0,p,p)
    matB[B>=cut] = 1
    und.w = which(matB + t(matB)==2,arr.ind=T)
    dir.w = which(matB==1 & t(matB)==0,arr.ind=T)
    
    if (length(und.w) > 0) {
        und.w = unique(t(apply(und.w,1,sort)))
        undattr = cbind(nodes[und.w[,1]],nodes[und.w[,2]],"und",round(B[und.w],digit=2))
    } else {
        undattr = NULL
    }
    
    if (length(dir.w) > 0) {
        dirattr = cbind(nodes[dir.w[,2]],nodes[dir.w[,1]],"dir",round(B[dir.w],digit=2))
    } else {
        dirattr = NULL
    }
    
    edgeattr = rbind(undattr,dirattr)
    
    onodes = nodes
    nodes = gsub("CNA_","CNA:",nodes)
    nodes = gsub("mRNA_","mRNA:",nodes)
    nodes = gsub("RPPA_","RPPA:",nodes)
    nodeattr = cbind(onodes,matrix(unlist(strsplit(nodes,split=":")),byrow=T,ncol=2))
    
    # Palbo
    load(paste(path_result,"/fit_outcome_Palbo_pp",pp,".RData",sep=""))
    cB = rowMeans(fit$Gamma)
    o = order(cB,decreasing=T)
    ii = max(which(cumsum(1-cB[o]) * (1/1:length(cB)) < 0.1))  # FDR at 0.1
    cut = cB[o][ii]
    w = which(cB>=cut)
    
    if (length(w)>0) {
        dname = colnames(X)[chlist[[4]]]
        ddir = cbind(colnames(X)[w],rep(dname,length(w)),rep("dir",length(w)),cB[w])
        dnod = c(dname,unlist(strsplit(dname,"_")))
        edgeattr = rbind(edgeattr,ddir)
        nodeattr = rbind(nodeattr,dnod)
    }
    
    # Fulv
    load(paste(path_result,"/fit_outcome_Fulv_pp",pp,".RData",sep=""))
    cB = rowMeans(fit$Gamma)
    o = order(cB,decreasing=T)
    ii = max(which(cumsum(1-cB[o]) * (1/1:length(cB)) < 0.1))  ## FDR at 0.1
    cut = cB[o][ii]
    w = which(cB>=cut)
    
    if (length(w)>0) {
        dname = colnames(X)[chlist[[4]]]
        ddir = cbind(colnames(X)[w],rep(dname,length(w)),rep("dir",length(w)),cB[w])
        dnod = c(dname,unlist(strsplit(dname,"_")))
        edgeattr = rbind(edgeattr,ddir)
        nodeattr = rbind(nodeattr,dnod)
    }
    
    # Tamo
    load(paste(path_result,"/fit_outcome_Tamo_pp",pp,".RData",sep=""))
    cB = rowMeans(fit$Gamma)
    o = order(cB,decreasing=T)
    ii = max(which(cumsum(1-cB[o]) * (1/1:length(cB)) < 0.1))  ## FDR at 0.1
    cut = cB[o][ii]
    w = which(cB>=cut)
    
    if (length(w)>0) {
        dname = colnames(X)[chlist[[4]]]
        ddir = cbind(colnames(X)[w],rep(dname,length(w)),rep("dir",length(w)),cB[w])
        dnod = c(dname,unlist(strsplit(dname,"_")))
        edgeattr = rbind(edgeattr,ddir)
        nodeattr = rbind(nodeattr,dnod)
    }
    
    pedgeattr = rbind(pedgeattr,edgeattr)
    pnodeattr = rbind(pnodeattr,nodeattr)
    
    write.table(edgeattr,file=paste(path_result_analysis,"/whole_edgeattrSingleDrug_pp",pp,".txt",sep=""),col.names=F,row.names=F,quote=F)
    write.table(nodeattr,file=paste(path_result_analysis,"/whole_nodeattrSingleDrug_pp",pp,".txt",sep=""),col.names=F,row.names=F,quote=F)
    
}

write.table(pedgeattr,file=paste(path_result_analysis,"/whole_edgeattrSingleDrug_pp1to12",".txt",sep=""),col.names=F,row.names=F,quote=F)
write.table(pnodeattr,file=paste(path_result_analysis,"/whole_nodeattrSingleDrug_pp1to12",".txt",sep=""),col.names=F,row.names=F,quote=F)

