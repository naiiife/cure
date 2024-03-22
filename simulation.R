generate <- function(n,rho=1,homo=TRUE){
  X = rbinom(n,1,0.5)
  V = rbinom(n,1,0.5)
  W = rbinom(n,1,0.5)
  A = rbinom(n,1,expit(-1+X/2+V/2-W/2))
  pu1 = expit(0+X/4-V/2)  # uncured prob
  pu0 = expit(-1/2-X/3+V)   # uncured prob
  uu = round(rho*((pu1+pu0)-abs(pu1-pu0))/2 + (1-rho)*pu1*pu0, 4)
  uc = round(pu1 - uu, 4)
  cu = round(pu0 - uu, 4)
  cc = round(1 - uu - uc - cu, 4)
  strata = rep(0,n)
  prob = cbind(uu,uc,cu,cc)
  for (i in 1:n){
    gen = rmultinom(n=1, size=1, prob=prob[i,])
    strata[i] = rownames(gen)[gen>0]
  }
  b1uu = b1uc = c(0,-1,1)
  b0uu = b0cu = c(1,1,0)
  if (!homo) {
    b1uc = c(.5,-.5,.5)
    b0cu = c(.5,.5,.5)
  }
  # Positive X leads to higher uncured rate, but longer survival time if uncured
  ind.uu = which(strata=='uu')
  ind.uc = which(strata=='uc')
  ind.cu = which(strata=='cu')
  ind.cc = which(strata=='cc')
  T1uu = rweibull(length(ind.uu),2,exp(-cbind(1,X,W)[ind.uu,]%*%b1uu/2)*5)
  T0uu = rweibull(length(ind.uu),2,exp(-cbind(1,X,W)[ind.uu,]%*%b0uu/2)*5)
  T1uc = rweibull(length(ind.uc),2,exp(-cbind(1,X,W)[ind.uc,]%*%b1uc/2)*5)
  T0cu = rweibull(length(ind.cu),2,exp(-cbind(1,X,W)[ind.cu,]%*%b0cu/2)*5)
  T1 = T0 = rep(99,n)
  T1[ind.uu] = T1uu
  T1[ind.uc] = T1uc
  T0[ind.uu] = T0uu
  T0[ind.cu] = T0cu
  C1 = runif(n,6,30)
  C0 = runif(n,5,30)
  T = A*T1 + (1-A)*T0
  C = A*C1 + (1-A)*C0
  status = (C>=T)
  Time = status*T + (1-status)*C
  return(data.frame(A,X,W,V,strata,T1,T0,Time,status))
}

n = 2000
x = (0:20)/2
bias1.prop = bias1.naive = bias1.pi = bias0.prop = bias0.naive = bias0.pi = NULL
set.seed(999)
for (k in 1:100){
  cat(k,'\n')
  dat0 = generate(n,rho=1,homo=FALSE)
  x1 = c(0,sort(dat0$T1[dat0$strata=='uu']),99)
  x0 = c(0,sort(dat0$T0[dat0$strata=='uu']),99)
  s1 = c(seq(1,0,length=1+sum(dat0$strata=='uu')),0)
  s0 = c(seq(1,0,length=1+sum(dat0$strata=='uu')),0)
  ty1 = matchy(s1,x1,x)
  ty0 = matchy(s0,x0,x)
  fit = try(surv_uu('Time','status','X','V','W','A',data=dat0,rho=1))
  if ('try-error' %in% class(fit)) next
  bias1.prop = rbind(bias1.prop,matchy(fit$Suu1,fit$time1,x) - ty1)
  bias0.prop = rbind(bias0.prop,matchy(fit$Suu0,fit$time0,x) - ty0)
  #fit.mis = try(surv_uu('Time','status','X','V','W','A',data=dat0,rho=0))
  #if ('try-error' %in% class(fit.mis)) next
  #bias1.mis = rbind(bias1.mis,matchy(fit.mis$Suu1,fit.mis$time1,x) - ty1)
  #bias0.mis = rbind(bias0.mis,matchy(fit.mis$Suu0,fit.mis$time0,x) - ty0)
  fit.pi = try(surv_pi('Time','status','X','V','W','A',data=dat0,rho=1))
  if ('try-error' %in% class(fit.pi)) next
  bias1.pi = rbind(bias1.pi,matchy(fit.pi$Suu1,fit.pi$time1,x) - ty1)
  bias0.pi = rbind(bias0.pi,matchy(fit.pi$Suu0,fit.pi$time0,x) - ty0)
  fit.naive1 = smcure(Surv(Time,status)~X+V+W,~X+V,data=dat0[dat0$A==1,],model='ph',
                    Var=FALSE)
  y1 = cbind(fit.naive1$Time,fit.naive1$s)
  y1 = rbind(c(0,1),y1[order(y1[,1]),])
  fit.naive0 = smcure(Surv(Time,status)~X+V+W,~X+V,data=dat0[dat0$A==0,],model='ph',
                    Var=FALSE)
  y0 = cbind(fit.naive0$Time,fit.naive0$s)
  y0 = rbind(c(0,1),y0[order(y0[,1]),])
  bias1.naive = rbind(bias1.naive,matchy(y1[,2],y1[,1],x) - ty1)
  bias0.naive = rbind(bias0.naive,matchy(y0[,2],y0[,1],x) - ty0)
}
bias1.prop = na.omit(bias1.prop)
bias0.prop = na.omit(bias0.prop)
bias1.pi = na.omit(bias1.pi)
bias0.pi = na.omit(bias0.pi)
bias1.naive = na.omit(bias1.naive)
bias0.naive = na.omit(bias0.naive)
meanbias1.prop = colMeans(bias1.prop)
meanbias1.pi = colMeans(bias1.pi)
meanbias1.naive = colMeans(bias1.naive)
rmse1.prop = sqrt(colMeans(bias1.prop^2))
rmse1.pi = sqrt(colMeans(bias1.pi^2))
rmse1.naive = sqrt(colMeans(bias1.naive^2))
meanbias0.prop = colMeans(bias0.prop)
meanbias0.pi = colMeans(bias0.pi)
meanbias0.naive = colMeans(bias0.naive)
rmse0.prop = sqrt(colMeans(bias0.prop^2))
rmse0.pi = sqrt(colMeans(bias0.pi^2))
rmse0.naive = sqrt(colMeans(bias0.naive^2))

plot(x,meanbias1.prop,type='l',ylim=c(-.2,.3),col='brown',lty=1,lwd=1.5,
     ylab='',xlab='Time',main='Mean Bias')
points(x,meanbias1.pi,type='l',col='brown',lty=4,lwd=1.5)
points(x,meanbias1.naive,type='l',col='brown',lty=5,lwd=1.5)
points(x,meanbias0.prop,type='l',ylim=c(-.1,.1),col='darkcyan',lty=1,lwd=1.5)
points(x,meanbias0.pi,type='l',col='darkcyan',lty=4,lwd=1.5)
points(x,meanbias0.naive,type='l',col='darkcyan',lty=5,lwd=1.5)
abline(h=0,lty=2)
mtext('Not principal ignorable')
legend('topright',cex=0.8,lty=rep(c(1,4,5),2),col=rep(c('brown','darkcyan'),each=3),
       lwd=rep(1.5,6),
       legend=c('Prop-SV (Treated)','Prop-PI (Treated)', 'AMC (Treated)',
                'Prop-SV (Controlled)','Prop-PI (Controlled)','AMC (Controlled)'))
plot(x,rmse1.prop,type='l',ylim=c(0,.25),col='brown',lty=1,lwd=1.5,ylab='',
     xlab='Time',main='Root Mean Squared Error')
points(x,rmse1.pi,type='l',col='brown',lty=4,lwd=1.5)
points(x,rmse1.naive,type='l',col='brown',lty=5,lwd=1.5)
points(x,rmse0.prop,type='l',ylim=c(-.1,.1),col='darkcyan',lty=1,lwd=1.5)
points(x,rmse0.pi,type='l',col='darkcyan',lty=4,lwd=1.5)
points(x,rmse0.naive,type='l',col='darkcyan',lty=5,lwd=1.5)
abline(h=0,lty=2)
mtext('Not principal ignorable')
legend('topright',cex=0.8,lty=rep(c(1,4,5),2),col=rep(c('brown','darkcyan'),each=3),
       lwd=rep(1.5,6),
       legend=c('Prop-SV (Treated)','Prop-PI (Treated)', 'AMC (Treated)',
                'Prop-SV (Controlled)','Prop-PI (Controlled)','AMC (Controlled)'))

#plot(x1,s1,type='l',col='brown',lwd=2,xlim=c(0,10),ylab='Survival Probability',
#     xlab='Time',main='True Survival function in UU')
#points(x0,s0,type='l',lwd=2,col='darkcyan')
#legend('topright',cex=0.8,col=c('brown','darkcyan'),
#       lwd=c(2,2),legend=c('True treated','True controlled'))
#points(fit$time1,fit$Suu1,type='s',lwd=2,col='brown')
#points(fit$time0,fit$Suu0,type='s',lwd=2,col='darkcyan')
#points(y1[,1],y1[,2],type='s',lwd=1.5,lty=5,col='brown')
#points(y0[,1],y0[,2],type='s',lwd=1.5,lty=5,col='darkcyan')
#legend('topright',cex=0.7,col=rep(c('brown','darkcyan'),each=3),
#       lwd=c(1,2,1.5,1,2,1.5),lty=c(1,1,5,1,1,5),
#       legend=c('True treated','Proposed','AMC','True controlled','Proposed','AMC'))
#x1 = c(0,sort(dat0$T1[dat0$strata=='uu'|dat0$strata=='uc']),99)
#x0 = c(0,sort(dat0$T0[dat0$strata=='uu'|dat0$strata=='cu']),99)
#s1 = c(seq(1,0,length=1+sum(dat0$strata=='uu'|dat0$strata=='uc')),0)
#s0 = c(seq(1,0,length=1+sum(dat0$strata=='uu'|dat0$strata=='cu')),0)
#plot(x1,s1,type='l',col='brown',lwd=2,xlim=c(0,10),ylab='Survival Probability',
#     xlab='Time',main='True Survival function')
#points(x0,s0,type='l',lwd=2,col='darkcyan')
#points(fit$time1,fit$Suu1,type='s',lwd=2,col='brown')
#points(fit$time0,fit$Suu0,type='s',lwd=2,col='darkcyan')
#legend('topright',cex=0.8,col=c('brown','darkcyan'),
#       lwd=c(2,2),legend=c('True treated','True controlled'))


## Sensitivity analysis on TDU
n = 2000
x = (0:20)/2
set.seed(999)
rseq = seq(0,1,length=11)
bias.prop = bias.pi = bias.naive = NULL
for (r in rseq){
  bias.r = bias.a = bias.naive.r = 0
  for (k in 1:100){
    cat(k,'\n')
    dat0 = generate(n,rho=r,homo=FALSE)
    x1 = c(0,sort(dat0$T1[dat0$strata=='uu']),99)
    x0 = c(0,sort(dat0$T0[dat0$strata=='uu']),99)
    s1 = c(seq(1,0,length=1+sum(dat0$strata=='uu')),0)
    s0 = c(seq(1,0,length=1+sum(dat0$strata=='uu')),0)
    #ty1 = matchy(s1,x1,x)
    #ty0 = matchy(s0,x0,x)
    TDU = sum(diff(c(0,x1))*s1) - sum(diff(c(0,x0))*s0)
    fit.r = try(surv_uu('Time','status','X','V','W','A',data=dat0,rho=1))
    if ('try-error' %in% class(fit.r)) next
    fit.a = try(surv_pi('Time','status',c('X','V'),NULL,'W','A',data=dat0,rho=1))
    if ('try-error' %in% class(fit.a)) next
    TDUr = sum(diff(c(0,fit.r$time1))*fit.r$Suu1)-sum(diff(c(0,fit.r$time0))*fit.r$Suu0)
    bias.r = bias.r + TDUr - TDU
    TDUa = sum(diff(c(0,fit.a$time1))*fit.a$Suu1)-sum(diff(c(0,fit.a$time0))*fit.a$Suu0)
    bias.a = bias.a + TDUa - TDU
    fit.naive1 = smcure(Surv(Time,status)~X+V+W,~X+V,data=dat0[dat0$A==1,],model='ph',
                      Var=FALSE)
    TDU1amc = sum(diff(c(0,sort(fit.naive1$Time)))*sort(fit.naive1$s,decreasing=TRUE))
    fit.naive0 = smcure(Surv(Time,status)~X+V+W,~X+V,data=dat0[dat0$A==0,],model='ph',
                      Var=FALSE)
    TDU0amc = sum(diff(c(0,sort(fit.naive0$Time)))*sort(fit.naive0$s,decreasing=TRUE))
    bias.naive.r = bias.naive.r + TDU1amc - TDU0amc - TDU
  }
  bias.prop = append(bias.prop,bias.r/100)
  bias.pi = append(bias.pi,bias.a/100)
  bias.naive = append(bias.naive,bias.naive.r/100)
}
plot(rseq, bias.prop, type='l', lwd=1.5, ylab='Bias of TDU', xlab=expression(rho),
     ylim=c(-2,2), col='black', main='Sensitivity Analysis')
points(rseq,bias.pi, type='l', col='black', lwd=1.5, lty=4)
points(rseq,bias.naive, type='l', col='black', lwd=1.5, lty=5)
mtext('Not principal ignorable')
abline(h=0,lty=2)
legend('topright',lty=c(1,4,5),lwd=c(1.5,1.5,1.5),col=rep('black',3),
       legend=c('Prop-SV','Prop-PI','AMC'))

