library(smcure)

expit <- function(x) exp(x)/(1+exp(x))

matchy <- function(yvec,xvec,newx){
  ivec = sapply(newx, function(x) max(which(xvec<=x)))
  if (is.vector(yvec)) {
    return(yvec[ivec])
  } else {
    return(yvec[ivec,])
  }
}

## Function to estimate Suu(t)

# binary V
surv_uu <- function(Time,status,X,V,W=NULL,TRT,data,rho=1,newdata=NULL){
  timeform = formula(paste0('Surv(',Time,',',status,')~',paste0(c(X,W,V),collapse='+')))
  cureform = formula(paste0('~',paste0(c(X,V),collapse='+')))
  treated = data[,TRT] == 1
  n = nrow(data)
  Xm = as.matrix(data[,X])
  Vm = as.vector(data[,V])
  # Fit mixture cure model
  cure1 = smcure(timeform, cureform, data=data[treated,], model='ph', Var=FALSE)
  names(cure1$b) = cure1$bnm
  cure0 = smcure(timeform, cureform, data=data[!treated,], model='ph', Var=FALSE)
  names(cure0$b) = cure0$bnm
  # Estimate survival probability in UU
  if (is.null(newdata)){
    Pv = mean(data[,V])
    pi1.v1 = as.numeric(expit(cure1$b[1] + Xm%*%cure1$b[X] + cure1$b[V]))
    pi1.v0 = as.numeric(expit(cure1$b[1] + Xm%*%cure1$b[X]))
    pi1 = cbind(pi1.v1,pi1.v0)
    pi0.v1 = as.numeric(expit(cure0$b[1] + Xm%*%cure0$b[X] + cure0$b[V]))
    pi0.v0 = as.numeric(expit(cure0$b[1] + Xm%*%cure0$b[X]))
    pi0 = cbind(pi0.v1,pi0.v0)
    Pu1 = as.numeric(expit(cure1$b[1] + Xm%*%cure1$b[X] + Vm*cure1$b[V]))
    Pu0 = as.numeric(expit(cure0$b[1] + Xm%*%cure0$b[X] + Vm*cure0$b[V]))
    P.uu = rho*((Pu1+Pu0)-abs(Pu1-Pu0))/2 + (1-rho)*Pu1*Pu0
    P.uc = Pu1 - P.uu
    P.cu = Pu0 - P.uu
    P.cc = 1 - P.uu - P.uc - P.cu
    pi.uu = rho*((pi1+pi0)-abs(pi1-pi0))/2 + (1-rho)*pi1*pi0
    pi.uc = pi1 - pi.uu
    pi.cu = pi0 - pi.uu
    pi.cc = 1 - pi.uu - pi.uc - pi.cu
    Su1.v1 = predictsmcure(cure1, newX=data.frame(data[,c(X,W)],V=1), newZ=data.frame(data[,X],V=1), model='ph')
    Su1.v0 = predictsmcure(cure1, newX=data.frame(data[,c(X,W)],V=0), newZ=data.frame(data[,X],V=0), model='ph')
    Su0.v1 = predictsmcure(cure0, newX=data.frame(data[,c(X,W)],V=1), newZ=data.frame(data[,X],V=1), model='ph')
    Su0.v0 = predictsmcure(cure0, newX=data.frame(data[,c(X,W)],V=0), newZ=data.frame(data[,X],V=0), model='ph')
    time1 = c(0,Su1.v1$prediction[,1+n])
    time0 = c(0,Su0.v1$prediction[,1+n])
    od1 = order(time1)
    od0 = order(time0)
    time1 = time1[od1]
    time0 = time0[od0]
    Su1.v1 = (t(Su1.v1$prediction[,1:n]) - 1) / Su1.v1$newuncureprob[1,] + 1
    Su1.v0 = (t(Su1.v0$prediction[,1:n]) - 1) / Su1.v0$newuncureprob[1,] + 1
    Su0.v1 = (t(Su0.v1$prediction[,1:n]) - 1) / Su0.v1$newuncureprob[1,] + 1
    Su0.v0 = (t(Su0.v0$prediction[,1:n]) - 1) / Su0.v0$newuncureprob[1,] + 1
    Suu1 = round((pi.uc[,2]*pi1[,1]*Su1.v1 - pi.uc[,1]*pi1[,2]*Su1.v0) / (pi.uc[,2]*pi1[,1] - pi.uc[,1]*pi1[,2]), 8)
    Suu0 = round((pi.cu[,2]*pi0[,1]*Su0.v1 - pi.cu[,1]*pi0[,2]*Su0.v0) / (pi.cu[,2]*pi0[,1] - pi.cu[,1]*pi0[,2]), 8)
    Suu1[is.na(Suu1)] = 1
    Suu0[is.na(Suu0)] = 1
    rowSumspi.uu = pi.uu[,1]*Pv+pi.uu[,2]*(1-Pv)
    rowSumspi.uc = pi.uc[,1]*Pv+pi.uc[,2]*(1-Pv)
    rowSumspi.cu = pi.cu[,1]*Pv+pi.cu[,2]*(1-Pv)
    eps = max(abs(pi.uc[,2]*pi1[,1]-pi.uc[,1]*pi1[,2])[apply(Suu1,1,max)>1|apply(Suu1,1,min)<0])
    eps = max(eps, 0.001)
    ucnull = abs(pi.uc[,2]*pi1[,1]-pi.uc[,1]*pi1[,2])<=eps
    Suu1[ucnull,] = (Su1.v1*pi.uu[,1]*Pv/rowSumspi.uu)[ucnull,]+(Su1.v0*pi.uu[,2]*(1-Pv)/rowSumspi.uu)[ucnull,]
    eps = max(abs(pi.cu[,2]*pi0[,1]-pi.cu[,1]*pi0[,2])[apply(Suu0,1,max)>1|apply(Suu0,1,min)<0])
    eps = max(eps, 0.001)
    cunull = abs(pi.cu[,2]*pi0[,1]-pi.cu[,1]*pi0[,2])<=eps
    Suu0[cunull,] = (Su0.v1*pi.uu[,1]*Pv/rowSumspi.uu)[cunull,]+(Su0.v0*pi.uu[,2]*(1-Pv)/rowSumspi.uu)[cunull,]
    
    Suc1 = round((pi.uu[,2]*pi1[,1]*Su1.v1 - pi.uu[,1]*pi1[,2]*Su1.v0) / (pi.uu[,2]*pi1[,1] - pi.uu[,1]*pi1[,2]), 8)#
    Scu0 = round((pi.uu[,2]*pi0[,1]*Su0.v1 - pi.uu[,1]*pi0[,2]*Su0.v0) / (pi.uu[,2]*pi0[,1] - pi.uu[,1]*pi0[,2]), 8)#
    Suc1[is.na(Suc1)] = 1#
    Scu0[is.na(Scu0)] = 1#
    eps = max(abs(pi.uu[,2]*pi1[,1]-pi.uu[,1]*pi1[,2])[apply(Suc1,1,max)>1|apply(Suc1,1,min)<0])
    eps = max(eps, 0.001)
    uunull = abs(pi.uu[,2]*pi1[,1]-pi.uu[,1]*pi1[,2])<=eps
    Suc1[uunull,] = (Su1.v1*pi.uc[,1]*Pv/rowSumspi.uc)[uunull,]+(Su1.v0*pi.uc[,2]*(1-Pv)/rowSumspi.uc)[uunull,]
    eps = max(abs(pi.uu[,2]*pi0[,1]-pi.uu[,1]*pi0[,2])[apply(Scu0,1,max)>1|apply(Scu0,1,min)<0])
    eps = max(eps, 0.001)
    uunull = abs(pi.uu[,2]*pi0[,1]-pi.uu[,1]*pi0[,2])<=eps
    Scu0[uunull,] = (Su0.v1*pi.cu[,1]*Pv/rowSumspi.cu)[uunull,]+(Su0.v0*pi.cu[,2]*(1-Pv)/rowSumspi.cu)[uunull,]
    
    Puu = mean(P.uu)
    Puc = mean(P.uc)
    Pcu = mean(P.cu)
    Pcc = mean(P.cc)
    Suu1 = c(1, colMeans(Suu1*P.uu)/Puu)
    Suu0 = c(1, colMeans(Suu0*P.uu)/Puu)
    Suu1 = Suu1[od1]
    Suu0 = Suu0[od0]
    
    Suc1 = c(1, colMeans(Suc1*P.uc)/Puc)
    Scu0 = c(1, colMeans(Scu0*P.cu)/Pcu)
    Suc1 = Suc1[od1]
    Scu0 = Scu0[od0]
  } else {
    Pv = mean(data[colSums(t(data[,X])-newdata[X])==0 & colSums(t(data[,W])-newdata[W])==0, V])
    Su1.v1 = predictsmcure(cure1, newX=data.frame(t(newdata[c(X,W)]),V=1), newZ=data.frame(t(newdata[X]),V=1), model='ph')
    Su1.v0 = predictsmcure(cure1, newX=data.frame(t(newdata[c(X,W)]),V=0), newZ=data.frame(t(newdata[X]),V=0), model='ph')
    Su0.v1 = predictsmcure(cure0, newX=data.frame(t(newdata[c(X,W)]),V=1), newZ=data.frame(t(newdata[X]),V=1), model='ph')
    Su0.v0 = predictsmcure(cure0, newX=data.frame(t(newdata[c(X,W)]),V=0), newZ=data.frame(t(newdata[X]),V=0), model='ph')
    time1 = c(0,Su1.v1$prediction[,2])
    time0 = c(0,Su0.v1$prediction[,2])
    od1 = order(time1)
    od0 = order(time0)
    time1 = time1[od1]
    time0 = time0[od0]
    Pu1.v1 = Su1.v1$newuncureprob[1,]
    Pu1.v0 = Su1.v0$newuncureprob[1,]
    Pu0.v1 = Su0.v1$newuncureprob[1,]
    Pu0.v0 = Su0.v0$newuncureprob[1,]
    Su1.v1 = (Su1.v1$prediction[,1] - 1) / Pu1.v1 + 1
    Su1.v0 = (Su1.v0$prediction[,1] - 1) / Pu1.v0 + 1
    Su0.v1 = (Su0.v1$prediction[,1] - 1) / Pu0.v1 + 1
    Su0.v0 = (Su0.v0$prediction[,1] - 1) / Pu0.v0 + 1
    pi1 = c(Pu1.v1,Pu1.v0)
    pi0 = c(Pu0.v1,Pu0.v0)
    Puu = rho*((pi1+pi0)-abs(pi1-pi0))/2 + (1-rho)*pi1*pi0
    Puc = pi1 - Puu
    Pcu = pi0 - Puu
    Pcc = 1 - Puu - Puc - Pcu
    Suu1 = round((Puc[2]*pi1[1]*Su1.v1 - Puc[1]*pi1[2]*Su1.v0) / (Puc[2]*pi1[1] - Puc[1]*pi1[2]), 8)
    Suu0 = round((Pcu[2]*pi0[1]*Su0.v1 - Pcu[1]*pi0[2]*Su0.v0) / (Pcu[2]*pi0[1] - Pcu[1]*pi0[2]), 8)
    sumPuu = Puu[1]*Pv+Puu[2]*(1-Pv)
    sumPuc = Puc[1]*Pv+Puc[2]*(1-Pv)
    sumPcu = Pcu[1]*Pv+Pcu[2]*(1-Pv)
    if (sum(is.na(Suu1))>0) Suu1 = (Su1.v1*Puu[1]*Pv/sumPuu)+(Su1.v0*Puu[2]*(1-Pv)/sumPuu)
    if (sum(is.na(Suu0))>0) Suu0 = (Su0.v1*Puu[1]*Pv/sumPuu)+(Su0.v0*Puu[2]*(1-Pv)/sumPuu)
    
    Suc1 = round((Puu[2]*pi1[1]*Su1.v1 - Puu[1]*pi1[2]*Su1.v0) / (Puu[2]*pi1[1] - Puu[1]*pi1[2]), 8)
    Scu0 = round((Puu[2]*pi0[1]*Su0.v1 - Puu[1]*pi0[2]*Su0.v0) / (Puu[2]*pi0[1] - Puu[1]*pi0[2]), 8)
    eps = abs(Puu[2]*pi1[1]-Puu[1]*pi1[2])
    eps = max(eps, 0.001)
    uunull = abs(Puu[2]*pi1[1]-Puu[1]*pi1[2])<=eps
    Suc1[uunull] = (Su1.v1*Puc[1]*Pv/sumPuc)[uunull]+(Su1.v0*Puc[2]*(1-Pv)/sumPuc)[uunull]
    Scu0[uunull] = (Su0.v1*Pcu[1]*Pv/sumPcu)[uunull]+(Su0.v0*Pcu[2]*(1-Pv)/sumPcu)[uunull]
    
    Suu1 = c(1, Suu1)
    Suu0 = c(1, Suu0)
    Suu1 = Suu1[od1]
    Suu0 = Suu0[od0]
    Puu = Puu[2-newdata[V]]
    Puc = Puc[2-newdata[V]]
    Pcu = Pcu[2-newdata[V]]
    Pcc = Pcc[2-newdata[V]]
    
    Suc1 = c(1, Suc1)
    Scu0 = c(1, Scu0)
    Suc1 = Suc1[od1]
    Scu0 = Scu0[od0]
  }
  #tmin = min(max(time1),max(time0))
  #n1 = sum(time1<tmin)
  #n0 = sum(time0<tmin)
  #RMST1 = sum(diff(c(0,time1[1:n1]))*c(1,Suu1[1:(n1-1)])) + (tmin-time1[n1])*Suu1[n1]
  #RMST0 = sum(diff(c(0,time0[1:n0]))*c(1,Suu0[1:(n0-1)])) + (tmin-time0[n0])*Suu0[n0]
  RMST1 = sum(diff(c(0,time1))*c(1,Suu1[-length(Suu1)]))
  RMST0 = sum(diff(c(0,time0))*c(1,Suu0[-length(Suu0)]))
  TDU = RMST1 - RMST0
  return(list(time1=time1,Suu1=Suu1,time0=time0,Suu0=Suu0,rho=rho,
              Puu=Puu,Puc=Puc,Pcu=Pcu,Pcc=Pcc,TDU=TDU,
              Suc1=Suc1, Scu0=Scu0))
}

# identification formula based
surv_uu1 <- function(Time,status,X,V,W=NULL,TRT,data,rho=1,newdata=NULL){
  timeform = formula(paste0('Surv(',Time,',',status,')~',paste0(c(X,V,W),collapse='+')))
  cureform = formula(paste0('~',paste0(c(X,V),collapse='+')))
  treated = data[,TRT] == 1
  n = nrow(data)
  Xm = as.matrix(data[,X])
  XWm = cbind(1,as.matrix(data[,c(X,W)]))
  Vm = as.matrix(data[,V])
  #betavx = solve(t(XWm)%*%XWm)%*%t(XWm)%*%Vm
  #Vm = as.matrix(Vm - XWm%*%betavx)
  #colnames(Vm) = V
  #data[,V] = Vm
  cure1 = smcure(timeform, cureform, data=data[treated,], model='ph', Var=FALSE)
  names(cure1$b) = cure1$bnm
  betax1 = cure1$b[X]
  betav1 = cure1$b[V]
  cure0 = smcure(timeform, cureform, data=data[!treated,], model='ph', Var=FALSE)
  names(cure0$b) = cure0$bnm
  betax0 = cure0$b[X]
  betav0 = cure0$b[V]
  
  Su1 = predictsmcure(cure1, newX=data[,c(X,V,W)], newZ=data[,c(X,V)], model='ph')
  Su0 = predictsmcure(cure0, newX=data[,c(X,V,W)], newZ=data[,c(X,V)], model='ph')
  time1 = Su1$prediction[,1+n]
  time0 = Su0$prediction[,1+n]
  
  if (is.null(newdata)){
    pi1 = cure1$b[1] + (Xm%*%betax1)%*%t(rep(1,n)) + rep(1,n)%*%t(Vm%*%betav1)
    pi1 = expit(pi1)
    pi0 = cure0$b[1] + (Xm%*%betax0)%*%t(rep(1,n)) + rep(1,n)%*%t(Vm%*%betav0)
    pi0 = expit(pi0)
    pi.uu = rho*((pi1+pi0)-abs(pi1-pi0))/2 + (1-rho)*pi1*pi0
    pi.uc = pi1 - pi.uu
    pi.cu = pi0 - pi.uu
    pi.cc = 1 - pi.uu - pi.uc - pi.cu
    pi.uu = round(pi.uu,6)
    pi.uc = round(pi.uc,6)
    pi.cu = round(pi.cu,6)
    pi.cc = round(pi.cc,6)
    betavx = solve(t(XWm)%*%XWm)%*%t(XWm)%*%Vm
    h1 = rowMeans(as.matrix(Vm))
    h2 = rowMeans(as.matrix(XWm%*%betavx))
    h = rep(1,n)%*%t(h1) - h2%*%t(rep(1,n))
    g1 = (pi1+0.0001)/(pi.uc+0.0001)*h
    g0 = (pi0+0.0001)/(pi.cu+0.0001)*h
    Su1 = (t(Su1$prediction[,1:n]) - 1) / diag(pi1) + 1
    Su0 = (t(Su0$prediction[,1:n]) - 1) / diag(pi0) + 1
    Suu1 = as.numeric(colMeans(diag(g1)/rowMeans(g1)*Su1*rowMeans(pi.uu))/mean(pi.uu))
    Suu0 = as.numeric(colMeans(diag(g0)/rowMeans(g0)*Su0*rowMeans(pi.uu))/mean(pi.uu))
    #if (max(diag(g1))>10000) Suu1 = as.numeric(colMeans(Su1[diag(g1)>10000,]))
    #if (max(diag(g0))>10000) Suu0 = as.numeric(colMeans(Su0[diag(g0)>10000,]))
    Suu1 = c(1,sort(Suu1,decreasing=TRUE))
    Suu0 = c(1,sort(Suu0,decreasing=TRUE))
    Suu1[Suu1>1] = 1
    Suu1[Suu1<0] = 0
    Suu0[Suu0>1] = 1
    Suu0[Suu0<0] = 0
    time1 = sort(time1)
    time0 = sort(time0)
    Puu = mean(diag(pi.uu))
    Puc = mean(diag(pi.uc))
    Pcu = mean(diag(pi.cu))
    Pcc = mean(diag(pi.cc))
    RMST1 = sum(diff(c(0,time1,max(time1)))*Suu1)
    RMST0 = sum(diff(c(0,time0,max(time0)))*Suu0)
    TDU = RMST1 - RMST0
  } else {
    m = nrow(newdata)
    Xm = newdata[,X]
    XWm = cbind(1,newdata[,c(X,W)])
    if(!is.null(newdata[,V])){
      if (nrow(newdata)==1) {
        newV = t(newdata[,V])
        XWm = cbind(1,t(newdata[,c(X,W)]))
      }
      #newV = newV - XWm%*%betavx
      pi1 = as.numeric(expit(cure1$b[1] + Xm%*%betax1 + newV%*%betav1))
      pi0 = as.numeric(expit(cure0$b[1] + Xm%*%betax0 + newV%*%betav0))
      Puu = rho*((pi1+pi0)-abs(pi1-pi0))/2 + (1-rho)*pi1*pi0
      Puc = pi1 - Puu
      Pcu = pi0 - Puu
      Pcc = 1 - Puu - Puc - Pcu
    } else {
      Puu = Puc = Pcu = Pcc = NULL
    }
    pi1 = cure1$b[1] + (Xm%*%betax1)%*%t(rep(1,n)) + rep(1,m)%*%t(Vm%*%betav1)
    pi1 = expit(pi1)
    pi0 = cure0$b[1] + (Xm%*%betax0)%*%t(rep(1,n)) + rep(1,m)%*%t(Vm%*%betav0)
    pi0 = expit(pi0)
    pi.uu = rho*((pi1+pi0)-abs(pi1-pi0))/2 + (1-rho)*pi1*pi0
    pi.uc = pi1 - pi.uu
    pi.cu = pi0 - pi.uu
    pi.cc = 1 - pi.uu - pi.uc - pi.cu
    
    h1 = rowMeans(as.matrix(Vm))
    h = rep(1,m)%*%t(h1) - mean(h1)
    g1 = (pi1+0.0001)/(pi.uc+0.0001)*h
    g0 = (pi0+0.0001)/(pi.cu+0.0001)*h
    CSuu1 = CSuu0 = NULL
    newX = data[,c(X,V,W)]
    newZ = data[,c(X,V)]
    for (i in 1:m){
      newX[,c(X,W)] = newdata[i,c(X,W)]
      newZ[,X] = newdata[i,X]
      g1x = g1[i,]
      g0x = g0[i,]
      Su1 = predictsmcure(cure1, newX=newX, newZ=newZ, model='ph')
      Su0 = predictsmcure(cure0, newX=newX, newZ=newZ, model='ph')
      Su1 = (t(Su1$prediction[,1:n]) - 1) / Su1$newuncureprob[1,] + 1
      Su0 = (t(Su0$prediction[,1:n]) - 1) / Su0$newuncureprob[1,] + 1
      CSuu1 = rbind(CSuu1, as.numeric(colMeans(g1x*Su1/mean(g1x))))
      CSuu0 = rbind(CSuu0, as.numeric(colMeans(g0x*Su0/mean(g0x))))
    }
    Suu1 = as.matrix(CSuu1[,order(time1)])
    Suu0 = as.matrix(CSuu0[,order(time0)])
    if (ncol(Suu1)==1) {
      Suu1 = t(Suu1)
      Suu0 = t(Suu0)
    }
    Suu1 = cbind(1,Suu1)
    Suu0 = cbind(1,Suu0)
    time1 = sort(time1)
    time0 = sort(time0)
    RMST1 = colSums(diff(c(0,time1,max(time1)))*t(Suu1))
    RMST0 = colSums(diff(c(0,time0,max(time0)))*t(Suu0))
    TDU = RMST1 - RMST0
  }
  
  Suu1[Suu1>1] = 1
  Suu1[Suu1<0] = 0
  Suu0[Suu0>1] = 1
  Suu0[Suu0<0] = 0
  return(list(time1=c(0,time1),Suu1=Suu1,
              time0=c(0,time0),Suu0=Suu0,
              Puu=Puu,Puc=Puc,Pcu=Pcu,Pcc=Pcc,TDU=TDU))
}

# principal ignorability
surv_pi <- function(Time,status,X,V=NULL,W=NULL,TRT,data,rho=1,newdata=NULL){
  timeform = formula(paste0('Surv(',Time,',',status,')~',paste0(c(X,W),collapse='+')))
  cureform = formula(paste0('~',paste0(c(X,V),collapse='+')))
  treated = data[,TRT] == 1
  n = nrow(data)
  Xm = as.matrix(data[,X])
  XWm = cbind(1,as.matrix(data[,c(X,W)]))
  Vm = as.matrix(data[,V])
  # Fit mixture cure model
  cure1 = smcure(timeform, cureform, data=data[treated,], model='ph', Var=FALSE)
  cure0 = smcure(timeform, cureform, data=data[!treated,], model='ph', Var=FALSE)
  names(cure1$b) = cure1$bnm
  betax1 = cure1$b[X]
  betav1 = cure1$b[V]
  names(cure0$b) = cure0$bnm
  betax0 = cure0$b[X]
  betav0 = cure0$b[V]
  # Estimate survival probability in UU
  if (is.null(newdata)){
    Pu1 = as.numeric(expit(cure1$b[1] + cbind(Xm,Vm)%*%cure1$b[-1]))
    Pu0 = as.numeric(expit(cure0$b[1] + cbind(Xm,Vm)%*%cure0$b[-1]))
    pi.uu = rho*((Pu1+Pu0)-abs(Pu1-Pu0))/2 + (1-rho)*Pu1*Pu0
    pi.uc = Pu1 - pi.uu
    pi.cu = Pu0 - pi.uu
    pi.cc = 1 - pi.uu - pi.uc - pi.cu
    Su1 = predictsmcure(cure1, newX=data[,c(X,W)], newZ=data[,c(X,V)], model='ph')
    Su0 = predictsmcure(cure0, newX=data[,c(X,W)], newZ=data[,c(X,V)], model='ph')
    t1 =  c(0,Su1$prediction[,1+n])
    t0 =  c(0,Su0$prediction[,1+n])
    time1 = sort(t1)
    time0 = sort(t0)
    Suu1 = (t(Su1$prediction[,1:n]) - 1) / Su1$newuncureprob[1,] + 1
    Suu0 = (t(Su0$prediction[,1:n]) - 1) / Su0$newuncureprob[1,] + 1
    Suu1 = rbind(cbind(1,Suu1),t1)
    Suu0 = rbind(cbind(1,Suu0),t0)
    Suu1 = Suu1[1:n,order(t1)]
    Suu0 = Suu0[1:n,order(t0)]
    Puu = mean(pi.uu)
    Puc = mean(pi.uc)
    Pcu = mean(pi.cu)
    Pcc = mean(pi.cc)
    Suu1 = colMeans(Suu1*pi.uu)/Puu
    Suu0 = colMeans(Suu0*pi.uu)/Puu
  } else {
    Su1 = predictsmcure(cure1, newX=data.frame(t(newdata[c(X,W)])), newZ=data.frame(t(newdata[c(X,V)])), model='ph')
    Su0 = predictsmcure(cure0, newX=data.frame(t(newdata[c(X,W)])), newZ=data.frame(t(newdata[c(X,V)])), model='ph')
    time1 = c(0,sort(Su1$prediction[,2]))
    time0 = c(0,sort(Su0$prediction[,2]))
    Pu1 = Su1$newuncureprob[1,]
    Pu0 = Su0$newuncureprob[1,]
    Suu1 = (Su1$prediction[,1] - 1) / Pu1 + 1
    Suu0 = (Su0$prediction[,1] - 1) / Pu0 + 1
    Puu = rho*((Pu1+Pu0)-abs(Pu1-Pu0))/2 + (1-rho)*Pu1*Pu0
    Puc = Pu1 - Puu
    Pcu = Pu0 - Puu
    Pcc = 1 - Puu - Puc - Pcu
    Suu1 = c(1, sort(Suu1, decreasing=TRUE))
    Suu0 = c(1, sort(Suu0, decreasing=TRUE))
  }
  RMST1 = sum(diff(c(0,time1))*Suu1)
  RMST0 = sum(diff(c(0,time0))*Suu0)
  TDU = RMST1 - RMST0
  return(list(time1=time1,Suu1=Suu1,
              time0=time0,Suu0=Suu0,
              Puu=Puu,Puc=Puc,Pcu=Pcu,Pcc=Pcc,TDU=TDU))
}

# bootstrap to obtain confidence interval
surv_uu_ci <- function(Time,status,X,V,W=NULL,TRT,data,rho=1,newdata=NULL,
                       method='uu',conf.int=.95,nboot=100,seed=2024){
  set.seed(seed)
  if (method=='uu') {
    fit = surv_uu(Time,status,X,V,W,TRT,data,rho,newdata)
  } else {
    fit = surv_pi(Time,status,X,V,W,TRT,data,rho,newdata)
  }
  t1 = fit$time1
  t0 = fit$time0
  Suu1 = as.numeric(fit$Suu1)
  Suu0 = as.numeric(fit$Suu0)
  S1 = S0 = NULL
  Puu = Puc = Pcu = Pcc = TDU = NULL
  for (b in 1:nboot){
    subs = sample(nrow(data),replace=TRUE)
    if (method=='uu') {
      fit.b = try(surv_uu(Time,status,X,V,W,TRT,data[subs,],rho,newdata), silent=TRUE)
    } else {
      fit.b = try(surv_pi(Time,status,X,V,W,TRT,data[subs,],rho,newdata), silent=TRUE)
    }
    if ('try-error' %in% class(fit.b)) next
    S1 = rbind(S1,matchy(as.numeric(fit.b$Suu1),fit.b$time1,t1))
    S0 = rbind(S0,matchy(as.numeric(fit.b$Suu0),fit.b$time0,t0))
    Puu = append(Puu,fit.b$Puu)
    Puc = append(Puc,fit.b$Puc)
    Pcu = append(Pcu,fit.b$Pcu)
    Pcc = append(Pcc,fit.b$Pcc)
    TDU = append(TDU,fit.b$TDU)
  }
  a = (1-conf.int)/2
  se.uu = quantile(Puu,c(a,1-a))
  se.uc = quantile(Puc,c(a,1-a))
  se.cu = quantile(Pcu,c(a,1-a))
  se.cc = quantile(Pcc,c(a,1-a))
  se.tdu = quantile(TDU,c(a,1-a))
  if (se.uc[1]<1e-3) se.uc = quantile(Puc,c(0,conf.int))
  if (se.cu[1]<1e-3) se.cu = quantile(Pcu,c(0,conf.int))
  S1[is.na(S1)] = 1
  S0[is.na(S0)] = 1
  S1[S1>1] = 1
  S0[S0>1] = 1
  S1[S1<0] = 0
  S0[S0<0] = 0
  ci1_u = as.numeric(apply(S1,2,function(j) quantile(j,1-a)))
  ci1_l = as.numeric(apply(S1,2,function(j) quantile(j,a)))
  ci0_u = as.numeric(apply(S0,2,function(j) quantile(j,1-a)))
  ci0_l = as.numeric(apply(S0,2,function(j) quantile(j,a)))
  return(list(time1=t1,Suu1=Suu1,ci1_u=ci1_u,ci1_l=ci1_l,
              time0=t0,Suu0=Suu0,ci0_u=ci0_u,ci0_l=ci0_l,
              Puu=fit$Puu,Puc=fit$Puc,Pcu=fit$Pcu,Pcc=fit$Pcc,
              se.uu=se.uu,se.uc=se.uc,se.cu=se.cu,se.cc=se.cc,
              TDU=fit$TDU,se.tdu=se.tdu))
}
