
# -*- coding: utf-8 -*-

rm(list=ls())

library(ROSE)

cols = c("darkorange","chartreuse3","darkviolet","skyblue2")

path_result <- "/root/Project/BANSMA/simulation/p50q4n200/result"
path_figure <- "/root/Project/BANSMA/simulation/p50q4n200/figure"

p = 50 # no. of mediators + exposure
q = 4 # no. of layers of X
n = 200

##########
# ER Probit plot

pbansmaIED = pbansmaIEC = pbansmaIE = bansmaIED = bansmaIEC = bansmaIE = bamaIED = medIED = truemat = numeric(0)

for (r in 1:25) {
    print(r)

    load(paste(path_result,"/ER_Setting_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
    truemat = rbind(truemat,true.out)

    load(paste(path_result,"/ER_BANSMA_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
    bansmaIE = rbind(bansmaIE,bansma.out$IE)
    bansmaIED = rbind(bansmaIED,bansma.out$IED)
    bansmaIEC = rbind(bansmaIEC,bansma.out$IEC)
    
    load(paste(path_result,"/ER_BANSMAProbit_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
    pbansmaIE = rbind(pbansmaIE,bansmaProbit.out$IE)
    pbansmaIED = rbind(pbansmaIED,bansmaProbit.out$IED)
    pbansmaIEC = rbind(pbansmaIEC,bansmaProbit.out$IEC)
    

    load(paste(path_result,"/ER_BAMA_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
    bamaIED = rbind(bamaIED,bama.out)

    load(paste(path_result,"/ER_Mediate_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
    medIED = rbind(medIED,med.single.out)
}

fit1.IED = ROSE::roc.curve(response=as.numeric(truemat[,"IED"]!=0),predicted=bansmaIED[,"PIP"],plotit=FALSE)
fit2.IED = ROSE::roc.curve(response=as.numeric(truemat[,"IED"]!=0),predicted=bamaIED[,"PIP"],plotit=FALSE)
fit3.IED = ROSE::roc.curve(response=as.numeric(truemat[,"IED"]!=0),predicted=-log10(medIED[,"p-value"]),plotit=FALSE)
fit4.IED = ROSE::roc.curve(response=as.numeric(truemat[,"IED"]!=0),predicted=pbansmaIED[,"PIP"],plotit=FALSE)

fit1.IE = ROSE::roc.curve(response=as.numeric(truemat[,"IE"]!=0),predicted=bansmaIE[,"PIP"],plotit=FALSE)
fit2.IE = ROSE::roc.curve(response=as.numeric(truemat[,"IE"]!=0),predicted=bamaIED[,"PIP"],plotit=FALSE)
fit3.IE = ROSE::roc.curve(response=as.numeric(truemat[,"IE"]!=0),predicted=-log10(medIED[,"p-value"]),plotit=FALSE)
fit4.IE = ROSE::roc.curve(response=as.numeric(truemat[,"IE"]!=0),predicted=pbansmaIE[,"PIP"],plotit=FALSE)

fit1.IEC = ROSE::roc.curve(response=as.numeric(truemat[,"IEC"]!=0),predicted=bansmaIEC[,"PIP"],plotit=FALSE)
fit2.IEC = ROSE::roc.curve(response=as.numeric(truemat[,"IEC"]!=0),predicted=bamaIED[,"PIP"],plotit=FALSE)
fit3.IEC = ROSE::roc.curve(response=as.numeric(truemat[,"IEC"]!=0),predicted=-log10(medIED[,"p-value"]),plotit=FALSE)
fit4.IEC = ROSE::roc.curve(response=as.numeric(truemat[,"IEC"]!=0),predicted=pbansmaIEC[,"PIP"],plotit=FALSE)

pdf(paste(path_figure,"/ER_Probit_Setting_p",p,"q",q,"n",n,".pdf",sep=""),width=12,height=4)

par(mfrow=c(1,3),mar=c(4.5,4.5,2,2),ps=20)
plot(as.numeric(fit1.IED$false.positive.rate),as.numeric(fit1.IED$true.positive.rate),col=cols[1],lwd=2,type="l",xlab="Specificity",ylab="Sensitivity")
lines(as.numeric(fit2.IED$false.positive.rate),as.numeric(fit2.IED$true.positive.rate), col=cols[2],lwd=2,lty=1)
lines(as.numeric(fit3.IED$false.positive.rate),as.numeric(fit3.IED$true.positive.rate), col=cols[3],lwd=2,lty=1)
lines(as.numeric(fit4.IED$false.positive.rate),as.numeric(fit4.IED$true.positive.rate), col=cols[4],lwd=2,lty=1)
legend("bottomright",legend=paste(c("BANSMA-IED","BAMA","SMA","pBANSMA-IED")," (AUC=",round(c(fit1.IED$auc,fit2.IED$auc,fit3.IED$auc,fit4.IED$auc),digit=2),")",sep=""),col=cols,lwd=2,lty=1,bty="n",y.intersp=1.2)
title("IED")
abline(0,1,col="grey")

plot(as.numeric(fit1.IEC$false.positive.rate),as.numeric(fit1.IEC$true.positive.rate),col=cols[1],lwd=2,type="l",xlab="Specificity",ylab="Sensitivity")
lines(as.numeric(fit2.IEC$false.positive.rate),as.numeric(fit2.IEC$true.positive.rate), col=cols[2],lwd=2,lty=1)
lines(as.numeric(fit3.IEC$false.positive.rate),as.numeric(fit3.IEC$true.positive.rate), col=cols[3],lwd=2,lty=1)
lines(as.numeric(fit4.IEC$false.positive.rate),as.numeric(fit4.IEC$true.positive.rate), col=cols[4],lwd=2,lty=1)
legend("bottomright",legend=paste(c("BANSMA-IEC","BAMA","SMA","pBANSMA-IEC")," (AUC=",round(c(fit1.IEC$auc,fit2.IEC$auc,fit3.IEC$auc,fit4.IEC$auc),digit=2),")",sep=""),col=cols,lwd=2,lty=1,bty="n",y.intersp=1.5)
title("IEC")
abline(0,1,col="grey")

plot(as.numeric(fit1.IE$false.positive.rate),as.numeric(fit1.IE$true.positive.rate),col=cols[1],lwd=2,type="l",xlab="Specificity",ylab="Sensitivity")
lines(as.numeric(fit2.IE$false.positive.rate),as.numeric(fit2.IE$true.positive.rate), col=cols[2],lwd=2,lty=1)
lines(as.numeric(fit3.IE$false.positive.rate),as.numeric(fit3.IE$true.positive.rate), col=cols[3],lwd=2,lty=1)
lines(as.numeric(fit4.IE$false.positive.rate),as.numeric(fit4.IE$true.positive.rate), col=cols[4],lwd=2,lty=1)
legend("bottomright",legend=paste(c("BANSMA-IE","BAMA","SMA","pBANSMA-IE")," (AUC=",round(c(fit1.IE$auc,fit2.IE$auc,fit3.IE$auc,fit4.IE$auc),digit=2),")",sep=""),col=cols,lwd=2,lty=1,bty="n",y.intersp=1.5)
title("IE")
abline(0,1,col="grey")

dev.off()

##########
# ER plot

cols = c("darkorange","chartreuse3","darkviolet")

bansmaIED = bansmaIEC = bansmaIE =bamaIED = medIED = truemat = numeric(0)

for (r in 1:25) {
    print(r)

    load(paste(path_result,"/ER_Setting_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
    truemat = rbind(truemat,true.out)

    load(paste(path_result,"/ER_BANSMA_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
    bansmaIE = rbind(bansmaIE,bansma.out$IE)
    bansmaIED = rbind(bansmaIED,bansma.out$IED)
    bansmaIEC = rbind(bansmaIEC,bansma.out$IEC)
    

    load(paste(path_result,"/ER_BAMA_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
    bamaIED = rbind(bamaIED,bama.out)

    load(paste(path_result,"/ER_Mediate_p",p,"q",q,"n",n,"r",r,".RData",sep=""))
    medIED = rbind(medIED,med.single.out)
}

fit1.IED = ROSE::roc.curve(response=as.numeric(truemat[,"IED"]!=0),predicted=bansmaIED[,"PIP"],plotit=FALSE)
fit2.IED = ROSE::roc.curve(response=as.numeric(truemat[,"IED"]!=0),predicted=bamaIED[,"PIP"],plotit=FALSE)
fit3.IED = ROSE::roc.curve(response=as.numeric(truemat[,"IED"]!=0),predicted=-log10(medIED[,"p-value"]),plotit=FALSE)

fit1.IE = ROSE::roc.curve(response=as.numeric(truemat[,"IE"]!=0),predicted=bansmaIE[,"PIP"],plotit=FALSE)
fit2.IE = ROSE::roc.curve(response=as.numeric(truemat[,"IE"]!=0),predicted=bamaIED[,"PIP"],plotit=FALSE)
fit3.IE = ROSE::roc.curve(response=as.numeric(truemat[,"IE"]!=0),predicted=-log10(medIED[,"p-value"]),plotit=FALSE)

fit1.IEC = ROSE::roc.curve(response=as.numeric(truemat[,"IEC"]!=0),predicted=bansmaIEC[,"PIP"],plotit=FALSE)
fit2.IEC = ROSE::roc.curve(response=as.numeric(truemat[,"IEC"]!=0),predicted=bamaIED[,"PIP"],plotit=FALSE)
fit3.IEC = ROSE::roc.curve(response=as.numeric(truemat[,"IEC"]!=0),predicted=-log10(medIED[,"p-value"]),plotit=FALSE)

pdf(paste(path_figure,"/ER_Setting_p",p,"q",q,"n",n,".pdf",sep=""),width=12,height=4)

par(mfrow=c(1,3),mar=c(4.5,4.5,2,2),ps=20)
plot(as.numeric(fit1.IED$false.positive.rate),as.numeric(fit1.IED$true.positive.rate),col=cols[1],lwd=2,type="l",xlab="Specificity",ylab="Sensitivity")
lines(as.numeric(fit2.IED$false.positive.rate),as.numeric(fit2.IED$true.positive.rate), col=cols[2],lwd=2,lty=1)
lines(as.numeric(fit3.IED$false.positive.rate),as.numeric(fit3.IED$true.positive.rate), col=cols[3],lwd=2,lty=1)
legend("bottomright",legend=paste(c("BANSMA-IED","BAMA","SMA")," (AUC=",round(c(fit1.IED$auc,fit2.IED$auc,fit3.IED$auc),digit=2),")",sep=""),col=cols,lwd=2,lty=1,bty="n",y.intersp=1.5)
title("IED")
abline(0,1,col="grey")

plot(as.numeric(fit1.IEC$false.positive.rate),as.numeric(fit1.IEC$true.positive.rate),col=cols[1],lwd=2,type="l",xlab="Specificity",ylab="Sensitivity")
lines(as.numeric(fit2.IEC$false.positive.rate),as.numeric(fit2.IEC$true.positive.rate), col=cols[2],lwd=2,lty=1)
lines(as.numeric(fit3.IEC$false.positive.rate),as.numeric(fit3.IEC$true.positive.rate), col=cols[3],lwd=2,lty=1)
legend("bottomright",legend=paste(c("BANSMA-IEC","BAMA","SMA")," (AUC=",round(c(fit1.IEC$auc,fit2.IEC$auc,fit3.IEC$auc),digit=2),")",sep=""),col=cols,lwd=2,lty=1,bty="n",y.intersp=1.5)
title("IEC")
abline(0,1,col="grey")

plot(as.numeric(fit1.IE$false.positive.rate),as.numeric(fit1.IE$true.positive.rate),col=cols[1],lwd=2,type="l",xlab="Specificity",ylab="Sensitivity")
lines(as.numeric(fit2.IE$false.positive.rate),as.numeric(fit2.IE$true.positive.rate), col=cols[2],lwd=2,lty=1)
lines(as.numeric(fit3.IE$false.positive.rate),as.numeric(fit3.IE$true.positive.rate), col=cols[3],lwd=2,lty=1)
legend("bottomright",legend=paste(c("BANSMA-IE","BAMA","SMA")," (AUC=",round(c(fit1.IE$auc,fit2.IE$auc,fit3.IEC$auc),digit=2),")",sep=""),col=cols,lwd=2,lty=1,bty="n",y.intersp=1.5)
title("IE")
abline(0,1,col="grey")

dev.off()
