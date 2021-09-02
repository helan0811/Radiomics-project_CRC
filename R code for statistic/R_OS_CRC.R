#####   Writen by Lan He
##### Email: helan0811@126.com
##### Department of Radiology, Guangdong Gernal Hospital, Guangdong Academy of Medical Sciences, 106 Zhongshan Er road, Guangzhou 510080, China. 


library(glmnet)
library(survival)

rm(list=ls(all=TRUE))
dat=read.csv("dat_surgery_GGH.csv",header=TRUE)
dat2=read.csv("dat_surgery_YCH.csv",header=TRUE)

datt=dat[1:161,]
datv=dat[162:409,]

################ (4*4) for lasso cox
set.seed(1)
y=Surv(datt$OS_time,datt$OS_status) 
x=as.matrix(datt[113:157])
fit_cv=cv.glmnet(x,y,family="cox",nfolds=10)
plot(fit_cv,cex.lab=1.25,cex.axis=1.25,cex.main=1.25)
abline(v=log(fit_cv$lambda.1se),col="darkgray",lwd=2)
abline(v=log(fit_cv$lambda.min),col="darkgray",lwd=2)
axis=(3,cex.lab=1.25)


model.final=fit_cv$glmnet.fit
plot(model.final,xvar=c("lambda"),label=TRUE,cex.axis=1.25,cex.lab=1.25,cex.main=1.25)
abline(v=log(fit_cv$lambda.min),col="darkgray",lwd=2)


#beta value
c=round(coef(model.final,s=fit_cv$lambda.min),8)
cc=summary(c)
#selected variables
LassoM.coef=coef(model.final,s=fit_cv$lambda.min)
VV=names(LassoM.coef[as.vector(LassoM.coef[,1]!=0),])
VV

## building signature
attach(datt)
datt$PI=get(VV[1])*cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
datt$PI=datt$PI+ve
}

attach(datv)
datv$PI=get(VV[1])*cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
datv$PI=datv$PI+ve
}

attach(dat2)
dat2$PI=get(VV[1])*cc[3][1,1]
for(m in 2:length(VV))
{
a=get(VV[m])
b=cc[3][m,1]
ve=a*b
dat2$PI=dat2$PI+ve
}


#####association with OS
library(survivalROC)
cutoff=36
datt$OS_status=ifelse(datt$OS_status==1,1,0)
roc=survivalROC(Stime=datt$OS_time,status=datt$OS_status,marker=datt$PI,predict.time=cutoff,method="KM")
cut.op=roc$cut.values[which.max(roc$TP-roc$FP)]

plot(roc$FP,roc$TP,type="l",lwd=2,xlim=c(0,1),ylim=c(0,1),col='firebrick3',,xlab="1-Specificity",ylab="Sensitivity",main="Time-dependt ROC at 3 years in the training cohort")
abline(0,1,lty=3,lwd=2,col='dodgerblue4')
legend(0.65,0.05,"Best cut-off = -1.235",box.lty=0,cex=1)


######  KM curve for training cohort
datt$group_5=ifelse(datt$PI<=cut.op,1,2)
s<-survdiff(Surv(datt$OS_time,datt$OS_status)~group_5,data=datt,rho=1)
p.val<-1-pchisq(s$chisq,length(s$n)-1)
p.val
HR<-(s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
HR
up95<-exp(log(HR)+qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
up95
low95<-exp(log(HR)-qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
low95

km<-survfit(Surv(datt$OS_time,datt$OS_status)~group_5,data=datt,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Overall Survival",xlab="Time Since Surgery (months)",main="Training cohort",col=c('deepskyblue3','goldenrod1'),cex.lab=1.25,cex.main=1.25,bty="n",yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1)
axis(1,at=c(0,12,24,36,48,60,72,84,96,108,120),labels=c(0,12,24,36,48,60,72,84,96,108,120),cex.axis=1)
legend(80,1.05,c("Low EcoRad-score","High EcoRad-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','goldenrod1'),box.lty=0)

#### C-index for signature
library(Hmisc)
r=rcorrcens(Surv(datt$OS_time,datt$OS_status==1,type="right")~(datt$PI))
1-r

 ######  KM curve for internal testing cohort
datv$group_5=ifelse(datv$PI<=cut.op,1,2)
s<-survdiff(Surv(datv$OS_time,datv$OS_status)~group_5,data=datv,rho=1)
p.val<-1-pchisq(s$chisq,length(s$n)-1)
p.val
HR<-(s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
HR
up95<-exp(log(HR)+qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
up95
low95<-exp(log(HR)-qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
low95

km<-survfit(Surv(datv$OS_time,datv$OS_status)~group_5,data=datv,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Overall Survival",xlab="Time Since Surgery (months)",main="Internal validation cohort",col=c('deepskyblue3','goldenrod1'),cex.lab=1.25,cex.main=1.25,bty="n",yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1)
axis(1,at=c(0,12,24,36,48,60,72,84),labels=c(0,12,24,36,48,60,72,84),cex.axis=1)
legend(54,1.05,c("Low EcoRad-score","High EcoRad-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','goldenrod1'),box.lty=0)

##### C-index for internal testing cohort
library(Hmisc)
r=rcorrcens(Surv(datv$OS_time,datv$OS_status==1,type="right")~(datv$PI))
1-r

 ######  KM curve for extrernal testing cohort
dat2$group_5=ifelse(dat2$PI<=cut.op,1,2)
s<-survdiff(Surv(dat2$OS_time,dat2$OS_status)~group_5,data=dat2,rho=1)
p.val<-1-pchisq(s$chisq,length(s$n)-1)
p.val
HR<-(s$obs[2]/s$exp[2])/(s$obs[1]/s$exp[1])
HR
up95<-exp(log(HR)+qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
up95
low95<-exp(log(HR)-qnorm(0.975)*sqrt(1/s$exp[2]+1/s$exp[1]))
low95

km<-survfit(Surv(dat2$OS_time,dat2$OS_status)~group_5,data=dat2,conf.int=TRUE)
plot(km,lwd=c(4,4),lty=c(1,1),ylab="Overall Survival",xlab="Time Since Surgery (months)",main="External validation cohort",col=c('deepskyblue3','goldenrod1'),cex.lab=1.25,cex.main=1.25,,bty="n",yaxt="n",xaxt="n",conf.int=F)
axis(2,at=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=c("0.0","0.2","0.4","0.6","0.8","1.0"),cex.axis=1)
axis(1,at=c(0,12,24,36,48,60,72,84),labels=c(0,12,24,36,48,60,72,84),cex.axis=1)
legend(50,1.05,c("Low EcoRad-score","High EcoRad-score"),cex=c(1,1),lty=1,col=c('deepskyblue3','goldenrod1'),box.lty=0)

##### C-index for external testing cohort
library(Hmisc)
r=rcorrcens(Surv(dat2$OS_time,dat2$OS_status==1,type="right")~(dat2$PI))
1-r

######### Cox for model building
library(rms)
dd<-datadist(datt)
attach(datt)
datt$OS_status=ifelse(datt$OS_status==1,1,0)
f=Surv(datt$OS_time,datt$OS_status==1,type="right")

options(datadist='dd')
coxm<-cph(f ~ tumor_location+age+sex+stage+T_stage+N_stage+PI,x=T,y=T,surv=TRUE,time.inc=36)
coxm
 
 ###after variables selecting
coxm<-cph(f ~ T_stage+N_stage+PI,x=T,y=T,surv=TRUE,time.inc=36)
coxm

c=coef(coxm)
datt$score=datt$T_stage*c[1]+datt$N_stage*c[2]+datt$PI*c[3]
datv$score=datv$T_stage*c[1]+datv$N_stage*c[2]+datv$PI*c[3]
dat2$score=dat2$T_stage*c[1]+dat2$N_stage*c[2]+dat2$PI*c[3]

####### C-index for model
library(Hmisc)
r=rcorrcens(Surv(datt$OS_time,datt$OS_status==1,type="right")~(datt$score))
1-r

library(Hmisc)
r=rcorrcens(Surv(datv$OS_time,datv$OS_status==1,type="right")~(datv$score))
1-r

library(Hmisc)
r=rcorrcens(Surv(dat2$OS_time,dat2$OS_status==1,type="right")~(dat2$score))
1-r

##### nomogram
dd<-datadist(datt)
attach(datt)
datt$OS_status=ifelse(datt$OS_status==1,1,0)
f=Surv(datt$OS_time,datt$OS_status==1,type="right")

coxm<-cph(f ~ T_stage+N_stage+PI,x=T,y=T,surv=TRUE,time.inc=36)
coxm

options(datadist="dd")
surv<-Survival(coxm)
surv1<-function(x)surv(36,lp=x)
coxm1<- Newlabels(coxm, c(PI='EcoRad signature')) 
coxm2<- Newlabels(coxm1, c(N_stage='N stage'))
coxm3<- Newlabels(coxm2, c(T_stage='T stage'))
non<-nomogram(coxm3,fun=list(surv1),lp=F,funlabel=c("3-year Overall Survival"),maxscale=10,fun.at=c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))
plot(non,lwd=5,xfrac=.25,label.every=2,cex.axis=1,cex.var=1,lmgp=.2,vnames=c("PI","Rad score"),varname.label=F,force.label=F)

##### nomogram
dd<-datadist(datt)
attach(datt)
datt$OS_status=ifelse(datt$OS_status==1,1,0)
f=Surv(datt$OS_time,datt$OS_status==1,type="right")

coxm1<-cph(f ~ T_stage+N_stage,x=T,y=T,surv=TRUE,time.inc=36)
coxm1

options(datadist="dd")
surv<-Survival(coxm1)
surv1<-function(x)surv(36,lp=x)
coxm3<- Newlabels(coxm1, c(N_stage='N stage'))
coxm4<- Newlabels(coxm3, c(T_stage='T stage'))
non<-nomogram(coxm4,fun=list(surv1),lp=F,funlabel=c("Overall Survival"),maxscale=10,fun.at=c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0))
plot(non,lwd=5,xfrac=.25,label.every=2,cex.axis=1,cex.var=1,lmgp=.2,vnames=c("PI","Rad score"),varname.label=F,force.label=F)

c1=coef(coxm1)
datt$score1=datt$T_stage*c1[1]+datt$N_stage*c1[2]
datv$score1=datv$T_stage*c1[1]+datv$N_stage*c1[2]
dat2$score1=dat2$T_stage*c1[1]+dat2$N_stage*c1[2]


library(Hmisc)
r=rcorrcens(Surv(datt$OS_time,datt$OS_status==1,type="right")~(datt$score1))
1-r

library(Hmisc)
r=rcorrcens(Surv(datv$OS_time,datv$OS_status==1,type="right")~(datv$score1))
1-r

library(Hmisc)
r=rcorrcens(Surv(dat2$OS_time,dat2$OS_status==1,type="right")~(dat2$score1))
1-r

#calibrate for nomogram
dd<-datadist(datt)
attach(datt)
datt$OS_status=ifelse(datt$OS_status==1,1,0)
f=Surv(datt$OS_time,datt$OS_status==1,type="right")
f_1=Surv(datt$OS_time_5,datt$OS_status_5==1,type="right")

coxm<-cph(f ~ T_stage+N_stage+PI,x=T,y=T,surv=TRUE,time.inc=36)
coxm
coxm_1<-cph(f_1 ~ T_stage+N_stage+PI,x=T,y=T,surv=TRUE,time.inc=36)



cal<-calibrate(coxm,u=36,cmethod='KM',method='boot',m=38,B=1000,riskdist=F)
cal_1<-calibrate(coxm_1,u=60,cmethod='KM',method='boot',m=38,B=1000,riskdist=F)

plot(cal,4,lwd=3,lty=1,xlim=c(0,1),ylim=c(0,1.0),bty="n",yaxt="n",xaxt="n",subtitle=FALSE, conf.int=TRUE,errbar.col='deepskyblue3',xlab="Radiomics model predict survival probability",ylab="Observed fraction survival probability",col='deepskyblue3')
plot(cal_1,4,lwd=3,lty=1,xlim=c(0,1),ylim=c(0,1.0),bty="n",yaxt="n",xaxt="n",subtitle=FALSE,conf.int=TRUE,errbar.col='goldenrod1',xlab="Nomogram predict survival probability",ylab="Observed fraction survival probability",col='goldenrod1',add=TRUE,main='Training cohort')
lines(cal[,'mean.predicted'],cal[,'KM'],type="b",lwd=3,col='deepskyblue3',pch=17)
lines(cal_1[,'mean.predicted'],cal_1[,'KM'],type="b",lwd=3,col='goldenrod1',pch=15)
abline(0,1,lty=3,lwd=2)
axis(1,seq(0,1,0.2),seq(0,1,0.2))
axis(2,seq(0,1,0.2),seq(0,1,0.2))
legend(0,1.0,c("3-year Survival","5-year Survival"),pch=c(17,15),col=c('deepskyblue3','goldenrod1'),box.lty=0)

dd<-datadist(datv)
attach(datv)
datv$OS_status=ifelse(datv$OS_status==1,1,0)
f=Surv(datv$OS_time,datv$OS_status==1,type="right")
f_1=Surv(datv$OS_time_5,datv$OS_status_5==1,type="right")

coxm<-cph(f ~ score,x=T,y=T,surv=TRUE,time.inc=36)
coxm
coxm_1<-cph(f_1 ~ score,x=T,y=T,surv=TRUE,time.inc=36)



cal<-calibrate(coxm,u=36,cmethod='KM',method='boot',m=62,B=1000,riskdist=F)
cal_1<-calibrate(coxm_1,u=60,cmethod='KM',method='boot',m=62,B=1000,riskdist=F)

plot(cal,4,lwd=3,lty=1,xlim=c(0,1),ylim=c(0,1.0),bty="n",yaxt="n",xaxt="n",subtitle=FALSE, conf.int=TRUE,errbar.col='deepskyblue3',xlab="Radiomics model predict survival probability",ylab="Observed fraction survival probability",col='deepskyblue3')
plot(cal_1,4,lwd=3,lty=1,xlim=c(0,1),ylim=c(0,1.0),bty="n",yaxt="n",xaxt="n",subtitle=FALSE,conf.int=TRUE,errbar.col='goldenrod1',xlab="Radiomics model predict survival probability",ylab="Observed fraction survival probability",col='goldenrod1',add=TRUE,main='Training cohort')
lines(cal[,'mean.predicted'],cal[,'KM'],type="b",lwd=3,col='deepskyblue3',pch=17)
lines(cal_1[,'mean.predicted'],cal_1[,'KM'],type="b",lwd=3,col='goldenrod1',pch=15)
abline(0,1,lty=3,lwd=2)
axis(1,seq(0,1,0.2),seq(0,1,0.2))
axis(2,seq(0,1,0.2),seq(0,1,0.2))
legend(0,1.0,c("3-year Survival","5-year Survival"),pch=c(17,15),col=c('deepskyblue3','goldenrod1'),box.lty=0)


dd<-datadist(dat2)
attach(dat2)
dat2$OS_status=ifelse(dat2$OS_status==1,1,0)
f=Surv(dat2$OS_time,dat2$OS_status==1,type="right")
f_1=Surv(dat2$OS_time_5,dat2$OS_status_5==1,type="right")

coxm<-cph(f ~ score,x=T,y=T,surv=TRUE,time.inc=36)
coxm
coxm_1<-cph(f_1 ~ score,x=T,y=T,surv=TRUE,time.inc=36)



cal<-calibrate(coxm,u=36,cmethod='KM',method='boot',m=25,B=1000,riskdist=F)
cal_1<-calibrate(coxm_1,u=60,cmethod='KM',method='boot',m=25,B=1000,riskdist=F)

plot(cal,4,lwd=3,lty=1,xlim=c(0,1),ylim=c(0,1.0),bty="n",yaxt="n",xaxt="n",subtitle=FALSE, conf.int=TRUE,errbar.col='deepskyblue3',xlab="Radiomics model predict survival probability",ylab="Observed fraction survival probability",col='deepskyblue3')
plot(cal_1,4,lwd=3,lty=1,xlim=c(0,1),ylim=c(0,1.0),bty="n",yaxt="n",xaxt="n",subtitle=FALSE,conf.int=TRUE,errbar.col='goldenrod1',xlab="Radiomics model predict survival probability",ylab="Observed fraction survival probability",col='goldenrod1',add=TRUE,main='Training cohort')
lines(cal[,'mean.predicted'],cal[,'KM'],type="b",lwd=3,col='deepskyblue3',pch=17)
lines(cal_1[,'mean.predicted'],cal_1[,'KM'],type="b",lwd=3,col='goldenrod1',pch=15)
abline(0,1,lty=3,lwd=2)
axis(1,seq(0,1,0.2),seq(0,1,0.2))
axis(2,seq(0,1,0.2),seq(0,1,0.2))
legend(0,1.0,c("3-year Survival","5-year Survival"),pch=c(17,15),col=c('deepskyblue3','goldenrod1'),box.lty=0)



########decision curve analysis (DCA) ##################
source("stdca.R")
Srv=Surv(datt$OS_time,datt$OS_status)
coxm=coxph(Srv ~T_stage+N_stage+PI,data=datt)
coxm1=coxph(Srv ~T_stage+N_stage, data=datt)
datt$pr_failure18=c(1-(summary(survfit(coxm,newdata=datt),times=36)$surv))
datt$pr_failure18_1=c(1-(summary(survfit(coxm1,newdata=datt),times=36)$surv))
#Run the decision curve analysis (with a smoother)
attach(datt)
c1=stdca(data=datt, outcome="OS_status",ttoutcome="OS_time",timepoint=36, 
predictors="pr_failure18", xstop=0.9, smooth=TRUE)
c2=stdca(data=datt, outcome="OS_status",ttoutcome="OS_time",timepoint=36, 
predictors="pr_failure18_1", xstop=0.9, smooth=TRUE)
#Plotting the curves
plot(c1$net.benefit$threshold, c1$net.benefit$none, type = "l", lwd=1, xlim=c(0,0.7),
ylim=c(-0.1,0.25),xlab = "Threshold Probability",main="Decision Curve",ylab = "Net Benefit")
lines(c1$net.benefit$threshold, c1$net.benefit$all, type="l", col="MediumAquamarine", 
lwd=1)
# lines(c2$net.benefit$threshold, c2$net.benefit$all, type="l", col=10, lwd=3, lty=2)
lines(c1$net.benefit$threshold, c1$net.benefit$pr_failure18, type="l", col='firebrick',
lwd=1)
lines(c2$net.benefit$threshold, c2$net.benefit$pr_failure18_1, type="l", col ='deepskyblue3', lty=1,lwd=1)
legend("topright", cex=0.8, legend=c("None", "Intervention All", "Radiomics model", "Reference model"),
col=c(17, "MediumAquamarine",'firebrick','deepskyblue3'), lwd=c(1, 1, 1, 1, 1), lty=c(1, 1, 1, 1, 1),box.lty=0)



############## iAUC ##################################
library(MASS)
library(risksetROC)
library(survival)
surv.prob=unique(survfit(Surv(datt$OS_time,datt$OS_status)~1)$surv)
fit0=coxph(Surv(datt$OS_time,datt$OS_status)~datt$score)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(datt$OS_time[datt$OS_status==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,datt$OS_time,datt$OS_status,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC1=AUC
iAUC1=iAUC
utimes1=utimes

surv.prob=unique(survfit(Surv(datt$OS_time,datt$OS_status)~1)$surv)
fit0=coxph(Surv(datt$OS_time,datt$OS_status)~datt$score1)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(datt$OS_time[datt$OS_status==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,datt$OS_time,datt$OS_status,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC2=AUC
iAUC2=iAUC
utimes2=utimes


plot(utimes1,AUC1,type="b",xlim=c(0,60),ylim=c(0.5,1.0),lty=2,bty="n",xaxt="n",yaxt="n",xlab='Times(Months)',ylab='Time-dependent AUC for Overall survival',col='deepskyblue3',main="Training cohort",pch=19,cex=1,lwd=1.5)
lines(utimes2,AUC2,type="b",ylim=c(0.4,1.0),lty=2,col='goldenrod1',pch=19,lwd=1.5)
axis(1,seq(0,60,12),seq(0,60,12))
axis(2,seq(0.5,1.0,0.1),seq(0.5,1.0,0.1))
legend(40,0.6,c("Radiomics model","Reference model"),cex=c(1,1),pch=19,col=c('deepskyblue3','goldenrod1'),box.lty=0)
abline(h=0.5,type="l",ylim=c(0.4,1.0),lty=2,col='firebrick',lwd=1.5)


library(MASS)
library(risksetROC)
library(survival)
surv.prob=unique(survfit(Surv(datv$OS_time,datv$OS_status)~1)$surv)
fit0=coxph(Surv(datv$OS_time,datv$OS_status)~datv$score)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(datv$OS_time[datv$OS_status==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,datv$OS_time,datv$OS_status,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC1=AUC
iAUC1=iAUC
utimes1=utimes

surv.prob=unique(survfit(Surv(datv$OS_time,datv$OS_status)~1)$surv)
fit0=coxph(Surv(datv$OS_time,datv$OS_status)~datv$score1)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(datv$OS_time[datv$OS_status==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,datv$OS_time,datv$OS_status,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC2=AUC
iAUC2=iAUC
utimes2=utimes

plot(utimes1,AUC1,type="b",xlim=c(0,60),ylim=c(0.5,1.0),lty=2,bty="n",xaxt="n",yaxt="n",xlab='Times(Months)',ylab='Time-dependent AUC for Overall survival',col='deepskyblue3',main="Internal validation cohort",pch=19,cex=1,lwd=1.5)
lines(utimes2,AUC2,type="b",ylim=c(0.4,1.0),lty=2,col='goldenrod1',pch=19,lwd=1.5)
axis(1,seq(0,60,12),seq(0,60,12))
axis(2,seq(0.5,1.0,0.1),seq(0.5,1.0,0.1))
legend(40,0.6,c("Radiomics model","Reference model"),cex=c(1,1),pch=19,col=c('deepskyblue3','goldenrod1'),box.lty=0)
abline(h=0.5,type="l",ylim=c(0.4,1.0),lty=2,col='firebrick',lwd=1.5)

library(MASS)
library(risksetROC)
library(survival)
surv.prob=unique(survfit(Surv(dat2$OS_time,dat2$OS_status)~1)$surv)
fit0=coxph(Surv(dat2$OS_time,dat2$OS_status)~dat2$score)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat2$OS_time[dat2$OS_status==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat2$OS_time,dat2$OS_status,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC1=AUC
iAUC1=iAUC
utimes1=utimes

surv.prob=unique(survfit(Surv(dat2$OS_time,dat2$OS_status)~1)$surv)
fit0=coxph(Surv(dat2$OS_time,dat2$OS_status)~dat2$score1)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat2$OS_time[dat2$OS_status==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat2$OS_time,dat2$OS_status,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC2=AUC
iAUC2=iAUC
utimes2=utimes

plot(utimes1,AUC1,type="b",xlim=c(0,60),ylim=c(0.5,1.0),lty=2,bty="n",xaxt="n",yaxt="n",xlab='Times(Months)',ylab='Time-dependent AUC for Overall survival',col='deepskyblue3',main="External validation cohort",pch=19,cex=1,lwd=1.5)
lines(utimes2,AUC2,type="b",ylim=c(0.4,1.0),lty=2,col='goldenrod1',pch=19,lwd=1.5)
axis(1,seq(0,60,12),seq(0,60,12))
axis(2,seq(0.5,1.0,0.1),seq(0.5,1.0,0.1))
legend(40,0.6,c("Radiomics model","Reference model"),cex=c(1,1),pch=19,col=c('deepskyblue3','goldenrod1'),box.lty=0)
abline(h=0.5,type="l",ylim=c(0.4,1.0),lty=2,col='firebrick',lwd=1.5)


#################################NRI & IDI #####################################
#####datt
library(survIDINRI)
library(survC1)
attach(datt)
datt$OS_status=ifelse(datt$OS_status==1,1,0)
datt$OS_time=datt$OS_time*30
D=subset(datt,select=c("OS_time","OS_status","score"))
D=D[!is.na(apply(D,1,mean)),]
D1=subset(datt,select=c("OS_time","OS_status","score1"))
D1=D1[!is.na(apply(D1,1,mean)),]
t0=365*3
indata1=D
indata0=D1
covs1<-as.matrix(indata1[,c(-1,-2)])
covs0<-as.matrix(indata0[,c(-1,-2)])
x=IDI.INF(D[,1:2],covs0,covs1,t0,npert=300)
IDI.INF.OUT(x)
IDI.INF.GRAPH(x)
#CI
Delta=Inf.Cval.Delta(D[,1:2], covs0, covs1, t0, itr=200)
round(Delta, digits=3)

#####datv
library(survIDINRI)
library(survC1)
attach(datv)
datv$OS_status=ifelse(datv$OS_status==1,1,0)
datv$OS_time=datv$OS_time*30
D=subset(datv,select=c("OS_time","OS_status","score"))
D=D[!is.na(apply(D,1,mean)),]
D1=subset(datv,select=c("OS_time","OS_status","score1"))
D1=D1[!is.na(apply(D1,1,mean)),]
t0=365*3
indata1=D
indata0=D1
covs1<-as.matrix(indata1[,c(-1,-2)])
covs0<-as.matrix(indata0[,c(-1,-2)])
x=IDI.INF(D[,1:2],covs0,covs1,t0,npert=300)
IDI.INF.OUT(x)
IDI.INF.GRAPH(x)
#CI
Delta=Inf.Cval.Delta(D[,1:2], covs0, covs1, t0, itr=200)
round(Delta, digits=3)

#########dat2
library(survIDINRI)
library(survC1)
attach(dat2)
dat2$OS_status=ifelse(dat2$OS_status==1,1,0)
dat2$OS_time=dat2$OS_time*30
D=subset(dat2,select=c("OS_time","OS_status","score"))
D=D[!is.na(apply(D,1,mean)),]
D1=subset(dat2,select=c("OS_time","OS_status","score1"))
D1=D1[!is.na(apply(D1,1,mean)),]
t0=365*3
indata1=D
indata0=D1
covs1<-as.matrix(indata1[,c(-1,-2)])
covs0<-as.matrix(indata0[,c(-1,-2)])
x=IDI.INF(D[,1:2],covs0,covs1,t0,npert=300)
IDI.INF.OUT(x)
IDI.INF.GRAPH(x)
#CI
Delta=Inf.Cval.Delta(D[,1:2], covs0, covs1, t0, itr=200)
round(Delta, digits=3)

############# iAUC for signature
library(MASS)
library(risksetROC)
library(survival)
surv.prob=unique(survfit(Surv(datt$OS_time,datt$OS_status)~1)$surv)
fit0=coxph(Surv(datt$OS_time,datt$OS_status)~datt$PI)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(datt$OS_time[datt$OS_status==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,datt$OS_time,datt$OS_status,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC1=AUC
iAUC1=iAUC
utimes1=utimes

library(MASS)
library(risksetROC)
library(survival)
surv.prob=unique(survfit(Surv(datv$OS_time,datv$OS_status)~1)$surv)
fit0=coxph(Surv(datv$OS_time,datv$OS_status)~datv$PI)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(datv$OS_time[datv$OS_status==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,datv$OS_time,datv$OS_status,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC1=AUC
iAUC1=iAUC
utimes1=utimes

library(MASS)
library(risksetROC)
library(survival)
surv.prob=unique(survfit(Surv(dat2$OS_time,dat2$OS_status)~1)$surv)
fit0=coxph(Surv(dat2$OS_time,dat2$OS_status)~dat2$PI)
eta=fit0$linear.predictor
model.score=eta
utimes=unique(dat2$OS_time[dat2$OS_status==1])
utimes=utimes[order(utimes)]
AUC=rep(NA,length(utimes))
for (j in 1:length(utimes))
{
	out=CoxWeights(eta,dat2$OS_time,dat2$OS_status,utimes[j])
	AUC[j]=out$AUC
}
iAUC=IntegrateAUC(AUC,utimes,surv.prob,tmax=36)
iAUC
AUC1=AUC
iAUC1=iAUC
utimes1=utimes

#############  AUC at 3 years
library(survival)
library(timeROC)
ROC.bili=timeROC(T=datt$OS_time,delta=datt$OS_status,marker=datt$PI,cause=1,weighting="marginal",times=c(12,24,36,48,60),ROC=T,iid=T)
ROC.bili
confint(ROC.bili,level=0.95)

library(survival)
library(timeROC)
ROC.bili=timeROC(T=datv$OS_time,delta=datv$OS_status,marker=datv$PI,cause=1,weighting="marginal",times=c(12,24,36,48,60),ROC=T,iid=T)
ROC.bili
confint(ROC.bili,level=0.95)

library(survival)
library(timeROC)
ROC.bili=timeROC(T=dat2$OS_time,delta=dat2$OS_status,marker=dat2$PI,cause=1,weighting="marginal",times=c(12,24,36,48,60),ROC=T,iid=T)
ROC.bili
confint(ROC.bili,level=0.95)

library(survival)
library(timeROC)
ROC.bili=timeROC(T=datt$OS_time,delta=datt$OS_status,marker=datt$score,cause=1,weighting="marginal",times=c(12,24,36,48,60),ROC=T,iid=T)
ROC.bili
confint(ROC.bili,level=0.95)

library(survival)
library(timeROC)
ROC.bili=timeROC(T=datv$OS_time,delta=datv$OS_status,marker=datv$score,cause=1,weighting="marginal",times=c(12,24,36,48,60),ROC=T,iid=T)
ROC.bili
confint(ROC.bili,level=0.95)

library(survival)
library(timeROC)
ROC.bili=timeROC(T=dat2$OS_time,delta=dat2$OS_status,marker=dat2$score,cause=1,weighting="marginal",times=c(12,24,36,48,60),ROC=T,iid=T)
ROC.bili
confint(ROC.bili,level=0.95)

library(survival)
library(timeROC)
ROC.bili=timeROC(T=datt$OS_time,delta=datt$OS_status,marker=datt$score1,cause=1,weighting="marginal",times=c(12,24,36,48,60),ROC=T,iid=T)
ROC.bili
confint(ROC.bili,level=0.95)

library(survival)
library(timeROC)
ROC.bili=timeROC(T=datv$OS_time,delta=datv$OS_status,marker=datv$score1,cause=1,weighting="marginal",times=c(12,24,36,48,60),ROC=T,iid=T)
ROC.bili
confint(ROC.bili,level=0.95)

library(survival)
library(timeROC)
ROC.bili=timeROC(T=dat2$OS_time,delta=dat2$OS_status,marker=dat2$score1,cause=1,weighting="marginal",times=c(12,24,36,48,60),ROC=T,iid=T)
ROC.bili
confint(ROC.bili,level=0.95)


