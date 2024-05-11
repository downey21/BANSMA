
# -*- coding: utf-8 -*-

rm(list=ls())

source("/root/Project/BANSMA/code/BANSMA.R")

path_data <- "/root/Project/BANSMA/application/data"
path_result <- "/root/Project/BANSMA/application/result"
path_result_analysis <- "/root/Project/BANSMA/application/result_analysis"

load(paste0(path_data,"/StringV10_ppi_pathway.RData"))
padat = ppi$pathwaydat
padat[,2] = gsub("-","",padat[,2])
ppname = unique(padat[,1])

pparray = (1:12)

mat= NULL
for (pp in pparray){
    print(pp)
    
    load(paste(path_result,"/fit_x_pp",pp,".RData",sep=""))
    medfit = list()
    medfit[[1]] = fit[[2]]
    medfit[[2]] = fit[[3]]
    
    # Palbo
    load(paste(path_result,"/fit_outcome_Palbo_pp",pp,".RData",sep=""))
    nodenames=colnames(X)
    outfit=fit
    bansma.out = getIndEffect(palist=palist,chlist=chlist,fit.med=medfit,fit.out=outfit)
    IE = as.matrix(bansma.out$IE)
    IED = as.matrix(bansma.out$IED)
    IEC = as.matrix(bansma.out$IEC)
    
    w =  which(IED[,6]>0.5|IEC[,6]>0.5|IE[,6]>0.5)
    if (length(w)>0) {
        dmat = matrix(NA,nrow=length(w),ncol=10)
        dmat[,1] =rep(ppname[pp],length(w))
        dmat[,2] =rep("Palbociclib",length(w))
        dmat[,3:4] =  nodenames[IE[w,1:2]]
        
        dmat[,5] = paste(round(IED[w,"IED"],digit=2)," (",round(IED[w,"CI_lower"],digit=2),",",round(IED[w,"CI_upper"],digit=2),")",sep="")
        dmat[,6] = round(IED[w,"PIP"],digit=2)
        dmat[,7] = paste(round(IEC[w,"IEC"],digit=2)," (",round(IEC[w,"CI_lower"],digit=2),",",round(IEC[w,"CI_upper"],digit=2),")",sep="")
        dmat[,8] = round(IEC[w,"PIP"],digit=2)
        dmat[,9] = paste(round(IE[w,"IE"],digit=2)," (",round(IE[w,"CI_lower"],digit=2),",",round(IE[w,"CI_upper"],digit=2),")",sep="")
        dmat[,10] =round(IE[w,"PIP"],digit=2)
        
        mat = rbind(mat,dmat)
    }
    
    # TAMO
    load(paste(path_result,"/fit_outcome_Tamo_pp",pp,".RData",sep=""))
    outfit=fit
    bansma.out = getIndEffect(palist=palist,chlist=chlist,fit.med=medfit,fit.out=outfit)
    IE = as.matrix(bansma.out$IE)
    IED = as.matrix(bansma.out$IED)
    IEC = as.matrix(bansma.out$IEC)
    
    w =  which(IED[,6]>0.5|IEC[,6]>0.5|IE[,6]>0.5)
    if (length(w)>0) {
        dmat = matrix(NA,nrow=length(w),ncol=10)
        dmat[,1] =rep(ppname[pp],length(w))
        dmat[,2] =rep("Tamoxifen",length(w))
        dmat[,3:4] = nodenames[IE[w,1:2]]
        
        dmat[,5] = paste(round(IED[w,"IED"],digit=2)," (",round(IED[w,"CI_lower"],digit=2),",",round(IED[w,"CI_upper"],digit=2),")",sep="")
        dmat[,6] = round(IED[w,"PIP"],digit=2)
        dmat[,7] = paste(round(IEC[w,"IEC"],digit=2)," (",round(IEC[w,"CI_lower"],digit=2),",",round(IEC[w,"CI_upper"],digit=2),")",sep="")
        dmat[,8] = round(IEC[w,"PIP"],digit=2)
        dmat[,9] = paste(round(IE[w,"IE"],digit=2)," (",round(IE[w,"CI_lower"],digit=2),",",round(IE[w,"CI_upper"],digit=2),")",sep="")
        dmat[,10] = round(IE[w,"PIP"],digit=2)
        
        mat = rbind(mat,dmat)
    }
    
    # Fulv
    load(paste(path_result,"/fit_outcome_Fulv_pp",pp,".RData",sep=""))
    outfit=fit
    bansma.out = getIndEffect(palist=palist,chlist=chlist,fit.med=medfit,fit.out=outfit)
    IE = as.matrix(bansma.out$IE)
    IED = as.matrix(bansma.out$IED)
    IEC = as.matrix(bansma.out$IEC)
    
    w =  which(IED[,6]>0.5|IEC[,6]>0.5|IE[,6]>0.5)
    if (length(w)>0) {
        dmat = matrix(NA,nrow=length(w),ncol=10)
        dmat[,1] = rep(ppname[pp],length(w))
        dmat[,2] = rep("Fulvestrant",length(w))
        dmat[,3:4] = nodenames[IE[w,1:2]]
        
        dmat[,5] = paste(round(IED[w,"IED"],digit=2)," (",round(IED[w,"CI_lower"],digit=2),",",round(IED[w,"CI_upper"],digit=2),")",sep="")
        dmat[,6] = round(IED[w,"PIP"],digit=2)
        dmat[,7] = paste(round(IEC[w,"IEC"],digit=2)," (",round(IEC[w,"CI_lower"],digit=2),",",round(IEC[w,"CI_upper"],digit=2),")",sep="")
        dmat[,8] = round(IEC[w,"PIP"],digit=2)
        dmat[,9] = paste(round(IE[w,"IE"],digit=2)," (",round(IE[w,"CI_lower"],digit=2),",",round(IE[w,"CI_upper"],digit=2),")",sep="")
        dmat[,10] = round(IE[w,"PIP"],digit=2)
        
        mat = rbind(mat,dmat)
    }
    
}

colnames(mat) = c("Pathway","Drug","Exposure","Mediator","IED (95% CI)","PIP","IEC (95% CI)","PIP","IE (95% CI)","PIP")
write.csv(mat,file=paste0(path_result_analysis,"/IE_results.csv"))

##########
# Make network label
pdf(paste0(path_result_analysis,"/net_lable.pdf"),width=10,height=10)

plot(0,0.5,col="white",xlim=c(0,1),ylim=c(0,1),pch=16,cex=10)
points(0.1,0.7,col=rgb(red=204,green=204,blue=204,maxColorValue=255),pch=16,cex=10)
text(0.4,0.7,"Drug",cex=3)
points(0.1,0.5,col=rgb(red=102,green=204,blue=0,maxColorValue=255),pch=16,cex=10)
text(0.4,0.5,"Protein",cex=3)
points(0.1,0.3,col=rgb(red=102,green=153,blue=255,maxColorValue=255),pch=16,cex=10)
text(0.4,0.3,"mRNA",cex=3)
points(0.1,0.1,col=rgb(red=255,green=153,blue=102,maxColorValue=255),pch=16,cex=10)
text(0.4,0.1,"Copy number",cex=3)

dev.off()

##########
# Make table
mat = read.csv(paste0(path_result_analysis,"/IE_results.csv"))
w = order(mat[,"Drug"])
mat = mat[w,]
remat = mat[,c("Pathway","Drug","Exposure","Mediator","PIP","PIP.1","PIP.2")]
rremat = data.frame(remat[,2],"&",remat[,1],"&",remat[,3],"&",remat[,4],"&",remat[,5],"&",remat[,6],"&",remat[,7],"\\")
write.table(rremat,file=paste0(path_result_analysis,"/IE_out.txt"),col.names=F,row.names=F,quote=F)
