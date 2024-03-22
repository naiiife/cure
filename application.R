library(smcure)

## Read data
dat = read.csv('leukemiaPKU.csv')
dat$TRANSPLANT = 1 - dat$TRANSPLANT
X = c('CR','ALL')
V = 'MRD'
W = NULL
TRT = 'TRANSPLANT'
Time = 'LFST'
status = 'LFS'

newx = rbind(c(1,1,1),c(1,0,1),c(0,1,1),c(0,0,1),
             c(1,1,0),c(1,0,0),c(0,1,0),c(0,0,0))
colnames(newx) = c('CR','ALL','MRD')


## Draw Suu(t), without conf int

fit = surv_uu(Time,status,X,V,W=NULL,TRT,dat,rho=1,newdata=NULL)
t1 = fit$time1
t0 = fit$time0
plot(t1[t1<3000],fit$Suu1[t1<3000],type='s',lwd=1.5,col='brown',xlab='Time',
     ylab='Survival Probability',main='Overall Population')
points(t0[t0<3000],fit$Suu0[t0<3000],type='s',lwd=1.5,col='darkcyan')
#points(t1[t1<3000],fit$Suc1[t1<3000],type='s',lwd=1.5,col='brown',lty=6)
#points(t0[t0<3000],fit$Scu0[t0<3000],type='s',lwd=1.5,col='darkcyan',lty=6)
#legend('topright',cex=0.8,col=rep(c('brown','darkcyan'),each=2),lty=c(1,6,1,6),
#       lwd=rep(1.5,4),legend=c('UU treated','UC treated','UU control','CU control'))

fit = surv_pi(Time,status,X,V,W=NULL,TRT,dat,rho=1,newdata=NULL)
t1 = fit$time1
t0 = fit$time0
points(t1[t1<3000],fit$Suu1[t1<3000],type='s',lwd=1.5,lty=4,col='brown')
points(t0[t0<3000],fit$Suu0[t0<3000],type='s',lwd=1.5,lty=4,col='darkcyan')

fit0 = smcure(Surv(LFST,LFS)~CR+ALL,~CR+ALL+MRD,data=dat[dat[,TRT]==0,],model='ph',Var=FALSE)
pre0 = predictsmcure(fit0, newX=dat[,X], newZ=dat[,c(X,V)], model='ph')
t0 = c(0,sort(fit0$Time))
s0 = c(1,sort(fit0$s, decreasing=TRUE))
fit1 = smcure(Surv(LFST,LFS)~CR+ALL,~CR+ALL+MRD,data=dat[dat[,TRT]==1,],model='ph',Var=FALSE)
pre1 = predictsmcure(fit0, newX=dat[,X], newZ=dat[,c(X,V)], model='ph')
t1 = c(0,sort(fit1$Time))
s1 = c(1,sort(fit1$s, decreasing=TRUE))
points(t1[t1<3000],s1[t1<3000],type='s',lwd=1.5,col='brown',lty=5)
points(t0[t0<3000],s0[t0<3000],type='s',lwd=1.5,col='darkcyan',lty=5)
legend('topright',cex=0.9,lwd=rep(1.5,6),col=rep(c('brown','darkcyan'),each=3),
       lty=rep(c(1,4,5),2),
       legend=c('MSDT (Prop-SV)','MSDT (Prop-PI)','MSDT (AMC)',
                'Haplo-SCT (Prop-SV)','Haplo-SCT (Prop-PI)','Haplo-SCT (AMC)'))


## Draw Suu(t), with conf int

# Substitutional variable
set.seed(2024)
fit = surv_uu_ci(Time,status,X,V,W=NULL,TRT,dat,rho=1,method='uu',nboot=200,newdata=NULL)
t1 = fit$time1
t0 = fit$time0
plot(t1[t1<3000],fit$Suu1[t1<3000],type='s',lwd=1.5,col='brown',xlab='Time',
     ylab='Survival Probability',main='Overall Population')
points(t0[t0<3000],fit$Suu0[t0<3000],type='s',lwd=1.5,col='darkcyan')
points(t1[t1<3000],fit$ci1_u[t1<3000],type='s',lty=2,col='brown')
points(t1[t1<3000],fit$ci1_l[t1<3000],type='s',lty=2,col='brown')
points(t0[t0<3000],fit$ci0_u[t0<3000],type='s',lty=2,col='darkcyan')
points(t0[t0<3000],fit$ci0_l[t0<3000],type='s',lty=2,col='darkcyan')
legend('topright',cex=0.9,lwd=c(1.5,1.5),col=c('brown','darkcyan'),
       legend=c('MSDT','Haplo-SCT'))
cat(paste0('Prop UU: ',round(fit$Puu,3),' (',round(fit$se.uu[1],3),', ',round(fit$se.uu[2],3),')\n'))
cat(paste0('Prop UC: ',round(fit$Puc,3),' (',round(fit$se.uc[1],3),', ',round(fit$se.uc[2],3),')\n'))
cat(paste0('Prop CU: ',round(fit$Pcu,3),' (',round(fit$se.cu[1],3),', ',round(fit$se.cu[2],3),')\n'))
cat(paste0('Prop CC: ',round(fit$Pcc,3),' (',round(fit$se.cc[1],3),', ',round(fit$se.cc[2],3),')\n'))
cat(paste0('TDU: ',round(fit$TDU,1),' (',round(fit$se.tdu[1],1),', ',round(fit$se.tdu[2],2),')\n'))

# Principal ignorability
set.seed(2024)
fit.adj = surv_uu_ci(Time,status,X=c(X,V),V=NULL,W=NULL,TRT,dat,rho=1,method='pi',nboot=200,newdata=NULL)
t1 = fit.adj$time1
t0 = fit.adj$time0
plot(t1[t1<3000],fit.adj$Suu1[t1<3000],type='s',lwd=1.5,col='brown',xlab='Time',
     ylab='Survival Probability',main='Overall Population')
points(t0[t0<3000],fit.adj$Suu0[t0<3000],type='s',lwd=1.5,col='darkcyan')
points(t1[t1<3000],fit.adj$ci1_u[t1<3000],type='s',lty=2,col='brown')
points(t1[t1<3000],fit.adj$ci1_l[t1<3000],type='s',lty=2,col='brown')
points(t0[t0<3000],fit.adj$ci0_u[t0<3000],type='s',lty=2,col='darkcyan')
points(t0[t0<3000],fit.adj$ci0_l[t0<3000],type='s',lty=2,col='darkcyan')
legend('topright',cex=0.9,lwd=c(1.5,1.5),col=c('brown','darkcyan'),
       legend=c('MSDT','Haplo-SCT'))
cat(paste0('Prop UU: ',round(fit.adj$Puu,3),' (',round(fit.adj$se.uu[1],3),', ',round(fit.adj$se.uu[2],3),')\n'))
cat(paste0('Prop UC: ',round(fit.adj$Puc,3),' (',round(fit.adj$se.uc[1],3),', ',round(fit.adj$se.uc[2],3),')\n'))
cat(paste0('Prop CU: ',round(fit.adj$Pcu,3),' (',round(fit.adj$se.cu[1],3),', ',round(fit.adj$se.cu[2],3),')\n'))
cat(paste0('Prop CC: ',round(fit.adj$Pcc,3),' (',round(fit.adj$se.cc[1],3),', ',round(fit.adj$se.cc[2],3),')\n'))
cat(paste0('TDU: ',round(fit.adj$TDU,2),' (',round(fit.adj$se.tdu[1],2),', ',round(fit.adj$se.tdu[2],2),')\n'))

# Conventional method (AMC)

fit0 = smcure(Surv(LFST,LFS)~CR+ALL+MRD,~CR+ALL+MRD,data=dat[dat[,TRT]==0,],model='ph',Var=FALSE)
pre0 = predictsmcure(fit0, newX=dat[,c(X,V)], newZ=dat[,c(X,V)], model='ph')
t0 = c(0,sort(fit0$Time))
s0 = c(1,sort(fit0$s, decreasing=TRUE))
fit1 = smcure(Surv(LFST,LFS)~CR+ALL+MRD,~CR+ALL+MRD,data=dat[dat[,TRT]==1,],model='ph',Var=FALSE)
pre1 = predictsmcure(fit0, newX=dat[,c(X,V)], newZ=dat[,c(X,V)], model='ph')
t1 = c(0,sort(fit1$Time))
s1 = c(1,sort(fit1$s, decreasing=TRUE))
points(t1[t1<3000],s1[t1<3000],type='s',lwd=1.5,col='brown',lty=5)
points(t0[t0<3000],s0[t0<3000],type='s',lwd=1.5,col='darkcyan',lty=5)


## Sensitivity analysis

# strata proportion and tdu versus rho
rseq = seq(0,1,length=31)
Puu = Puc = Pcu = Pcc = TDU.sv = TDU.pi = NULL
for (r in rseq){
  fit = surv_uu(Time,status,X,V,W=NULL,TRT,dat,rho=r,newdata=NULL)
  Puu = append(Puu, fit$Puu)
  Puc = append(Puc, fit$Puc)
  Pcu = append(Pcu, fit$Pcu)
  Pcc = append(Pcc, fit$Pcc)
  TDU.sv = append(TDU.sv, fit$TDU)
  fit = surv_pi(Time,status,X,V,W=NULL,TRT,dat,rho=r,newdata=NULL)
  TDU.pi = append(TDU.pi, fit$TDU)
  cat(r,'\n')
}
plot(rseq,Puu,type='l',lwd=2,xlab=expression(rho),ylab='Probability',ylim=c(0,.7),
     main='Strata Proportions',col='darkred')
points(rseq,Puc,type='l',lwd=2,col='darkblue')
points(rseq,Pcu,type='l',lwd=2,col='darkgreen')
points(rseq,Pcc,type='l',lwd=2,col='darkorange')
legend('topleft',cex=0.8,lwd=rep(2,4),
       col=c('darkred','darkblue','darkgreen','darkorange'),
       legend=c('UU','UC','CU','CC'))
plot(rseq,TDU.sv,type='l',lwd=2,xlab=expression(rho),ylab='TDU',ylim=c(-150,250),
     main='Mean Survival Time Difference',col='black')
points(rseq,TDU.pi,type='l',lwd=2,lty=4,col='black')
legend('topleft',cex=0.8,lwd=rep(2,2),lty=c(1,4),
       col=c('black','black'),legend=c('Prop-SV','Prop-PI'))
abline(h=0,lty=2)

# survival curve versus rho
plot(NULL,NULL,type='s',lwd=1.5,col='brown',ylim=c(0,1),xlim=c(0,3000),xlab='Time',
     ylab='Survival Probability',main='Sensitivity Analysis')
for (r in c(1,.5,0)){
  l = 3 - 2*r
  fit = surv_uu(Time,status,X,V,W=NULL,TRT,dat,rho=r,newdata=NULL)
  t1 = fit$time1
  t0 = fit$time0
  points(t1[t1<3000],fit$Suu1[t1<3000],type='s',lwd=1.5,lty=l,col='brown')
  points(t0[t0<3000],fit$Suu0[t0<3000],type='s',lwd=1.5,lty=l,col='darkcyan')
}
legend('topright',cex=0.8,legend=c('MSDT, rho = 1', 'MSDT, rho = 0.5', 'MSDT, rho = 0',
                                   'Haplo-SCT, rho = 1', 'Haplo-SCT, rho = 0.5', 'Haplo-SCT, rho = 0'), 
       lwd=rep(1.5,6), lty=c(1,2,3,1,2,3), col=rep(c('brown','darkcyan'),each=3))
