rm(list=ls())
set.seed(100)

# Model:
# beta.psi=c(beta20,beta21,beta22,beta23,beta24,beta10,beta11,beta12,beta13,psi0,psi1,psi2,psi3)
# Q2=beta20+beta21*Q2+beta22*S2+beta23*P2+beta24*A1 + A2*(psi0+psi1*Q2+psi2*S2)
# Q1=beta10+beta11*Q1+beta12*S1+beta13*P1           + A1*(psi0+psi1*Q1+psi2*S1+psi3*P1)

g_shared=function(dta,tmt,ini){
  # function of shared g-estimation.
  
  # dta: data of interest
  # tmt: either 1 or 2 for different treatment modeling
  # ini: initial value used
  
  Q1=dta$Q1
  Q2=dta$Q2
  A1=dta$A1
  A2=dta$A2
  S1=dta$S1
  S2=dta$S2
  P1=dta$P1
  P2=dta$P2
  Y =dta$Y
  
  # Keep only participants having P2 and S2 for stage 2, S2 and P2 are key variables indicate whether someone got treated at stage 2.
  index2=!is.na(dta$P2) & !is.na(dta$S2) & is.finite(dta$S2) 
  temp2=dta[index2,] 

  # Specify the blip model
  H2_psi=cbind(1,temp2$Q2,temp2$S2,0)
  H1_psi=cbind(1,Q1,S1,P1)
  
  # Specify the treatment free outcome model
  H2_beta=cbind(1,temp2$Q2,temp2$S2,temp2$P2,temp2$A1)
  H1_beta=cbind(1,Q1,S1,P1)
  
  # Specify the treatment model
  if(tmt==1){
    stage1=glm(A1~Age+Sex+Q1+S1+P1-1,family = "binomial",data=dta)
    stage2=glm(A2~Age+Sex+Q2+S2+P2-1,family = "binomial",data=temp2)
    alpha1=stage1$coefficients;alpha1
    alpha2=stage2$coefficients;alpha2
    EA2=stage2$fitted.values
    EA1=stage1$fitted.values
  }
  else if(tmt==2){
    stage1=glm(A1~Age+Sex+Race+School+Emplcat+Privins+Q1+S1+P1-1,family = "binomial",data=dta)
    stage2=glm(A2~Age+Sex+Race+School+Emplcat+Privins+Q2+S2+P2-1,family = "binomial",data=temp2)
    alpha1=stage1$coefficients;alpha1
    alpha2=stage2$coefficients;alpha2
    EA2=stage2$fitted.values
    EA1=stage1$fitted.values
  }
  
  # Estimation of psi using the following matrix
  
  # X_omega=[H_beta,H_psi*(A-E[A])]
  # X_delta=[H_beta,A*H_psi]
  # delta=[beta,psi]
  # delta_hat=solve(t(X_omega)*X_delta)*t(X_omega)*Y_hat
  
  error=1e-6
  diff=1e3
  psi=as.matrix(ini,nrow=4,ncol=1) #psi0, psi1, psi2, psi3
  temp_betapsi=rep(0,ncol(H1_beta)+ncol(H2_beta)+ncol(H1_psi)) #5+3+5 beta1, beta2, psi
  j=0
  dta$Y.S1=dta$Y
  temp2$Y.S2=temp2$Y
  
  # Iterate until converge
  while(diff > error){
    
    j=j+1
    
    # Construct the outcome matrix Ypseu
    Y2=temp2$Y.S2 
    dta[index2,]$Y.S1=dta[index2,]$Y-A2[index2]*(H2_psi%*%psi)+0.5*(H2_psi%*%psi+abs(H2_psi%*%psi))
    Y1=dta$Y.S1 #Y1 will change along with psi
    Ypseu=c(Y1,Y2)
    
    # Construct the covariate matrix X_omega, X_delta.
    Xbeta1=cbind(H1_beta,matrix(0,nrow=nrow(H1_beta),ncol=ncol(H2_beta)))
    Xbeta2=cbind(matrix(0,nrow=nrow(H2_beta),ncol=ncol(H1_beta)),H2_beta)
    Xbeta=rbind(Xbeta1,Xbeta2);dim(Xbeta)
    
    AXpsi=rbind(A1*H1_psi,A2[index2]*H2_psi);dim(AXpsi)
    Xdelta=cbind(Xbeta,AXpsi);dim(Xdelta)

    XpsiAE=rbind(H1_psi*(A1-EA1), H2_psi*(A2[index2]-EA2));dim(XpsiAE)
    
    Xomega=cbind(Xbeta,XpsiAE);dim(Xomega)

    # Store betapsi from previous iteration.
    old_betapsi=temp_betapsi
    
    # Solve psi for the current iteration.
    temp_betapsi=solve(t(Xomega)%*%Xdelta)%*%t(Xomega)%*%Ypseu
    psi=temp_betapsi[(ncol(H1_beta)+ncol(H2_beta)+1):length(temp_betapsi)]
    
    diff=sum(abs(temp_betapsi-old_betapsi))
  }
  return(list(psi=psi,iter=j))
}


boot_g=function(dta,tmt,ini,bootsize){
  # function of bootstrapping shared g-estimation.
  
  # dta: data of interest
  # tmt: 1 or 2 for different treatment models
  # ini: initial value used for g-est
  # bootsize: total number of bootstrapping
  
  index2=!is.na(dta$P2) & !is.na(dta$S2) & is.finite(dta$S2)
  s1index=c(1:nrow(dta))
  bootpsi=matrix(NA,bootsize,4)
  
  for(i in 1:bootsize){
    set.seed(i)
    # re-sampling
    new_s1index=sample(s1index,replace = TRUE,size=length(s1index))
    news1=dta[new_s1index,]
    bootpsi[i,]=g_shared(news1,tmt,ini)$psi
  }
  return(list(bootpsi=bootpsi))
}

# Calculate bootstrapping statistics with 1000 samples
ini.zero=c(0,0,0,0)
boot_result=boot_g(dta,tmt=1,ini.zero,bootsize = 1000)
boot_result=boot_g(dta,tmt=2,ini.zero,bootsize = 1000)

# Collect relevent results
bootmean=c(mean(boot_result$bootpsi[,1]),mean(boot_result$bootpsi[,2]),mean(boot_result$bootpsi[,3]),mean(boot_result$bootpsi[,4]))
bootsd=c(sd(boot_result$bootpsi[,1]),sd(boot_result$bootpsi[,2]),sd(boot_result$bootpsi[,3]),sd(boot_result$bootpsi[,4]))

