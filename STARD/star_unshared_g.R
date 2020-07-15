rm(list=ls())
set.seed(100)

# Model:
# beta.psi=c(beta20,beta21,beta22,beta23,beta24,beta10,beta11,beta12,beta13,psi0,psi1,psi2,psi3)
# Q2=beta20+beta21*Q2+beta22*S2+beta23*P2+beta24*A1 + A2*(psi0+psi1*Q2+psi2*S2)
# Q1=beta10+beta11*Q1+beta12*S1+beta13*P1           + A1*(psi0+psi1*Q1+psi2*S1+psi3*P1)


g_unshared=function(dta,tmt){
  # function of unshared g-estimation.
  
  # dta: data of interest
  # tmt: either 1 or 2 for different treatment modeling

  # Keep only participants having P2 and S2 for stage 2, S2 and P2 are key variables indicating whether someone got treated at stage 2.
  index2=!is.na(dta$P2) & !is.na(dta$S2) & is.finite(dta$S2) 
  temp2=dta[index2,]

  # Specify the blip model
  H2_psi=cbind(1,temp2$Q2,temp2$S2)
  H1_psi=cbind(1,dta$Q1,dta$S1,dta$P1)

  # Specify the treatment free outcome model
  H2_beta=cbind(1,temp2$Q2,temp2$S2,temp2$P2,temp2$A1)
  H1_beta=cbind(1,dta$Q1,dta$S1,dta$P1)

  # Specify the treatment model
  if(tmt==1){
    stage1=glm(A1~Age+Sex+Q1+S1+P1-1,family = "binomial",data=dta)
    stage2=glm(A2~Age+Sex+Q2+S2+P2-1,family = "binomial",data=temp2)
    alpha1=stage1$coefficients;alpha1
    alpha2=stage2$coefficients;alpha2
    EA2=stage2$fitted.values
    EA1=stage1$fitted.values
  }else if(tmt==2){
    stage1=glm(A1~Age+Sex+Race+School+Emplcat+Privins+Q1+S1+P1-1,family = "binomial",data=dta)
    stage2=glm(A2~Age+Sex+Race+School+Emplcat+Privins+Q2+S2+P2-1,family = "binomial",data=temp2)
    alpha1=stage1$coefficients;alpha1
    alpha2=stage2$coefficients;alpha2
    EA2=stage2$fitted.values
    EA1=stage1$fitted.values
  }
    
  # Construct covariate matrix for stage 2
  X2_omega=cbind(H2_beta,H2_psi*(temp2$A2-EA2))
  X2_delta=cbind(H2_beta,temp2$A2*H2_psi)
  
  # Solve beta psi for stage 2
  beta.psi2=solve(t(X2_omega)%*%X2_delta)%*%t(X2_omega)%*%temp2$Y
  psi.s2=beta.psi2[6:8]
  
  # Construct covariate matrix for stage 2
  X1_omega=cbind(H1_beta,H1_psi*(dta$A1-EA1))
  X1_delta=cbind(H1_beta,dta$A1*H1_psi)

  # Calculate pseudo outcome and construct outcome matrix Y
  Y1.pseu=dta$Y
  Y1.pseu[index2]=dta[index2,]$Y-temp2$A2*(H2_psi%*%psi.s2)+0.5*(H2_psi%*%psi.s2+abs(H2_psi%*%psi.s2))

  # Solve for beta psi for stage 1
  beta.psi1=solve(t(X1_omega)%*%X1_delta)%*%t(X1_omega)%*%Y1.pseu
  psi.s1=beta.psi1[5:8]

  # Simple average estimate
  psi1_sa=(psi.s1[1]+psi.s2[1])/2
  psi2_sa=(psi.s1[2]+psi.s2[2])/2
  psi3_sa=(psi.s1[3]+psi.s2[3])/2
  psi4_sa= psi.s1[4]

  psi_sa=c(psi1_sa,psi2_sa,psi3_sa,psi4_sa);psi_sa

  return(list(psi.s1=psi.s1,psi.s2=psi.s2,psi_SA=psi_sa))
}


boot_g=function(dta,tmt,bootsize){
  # function of bootstrapping unshared g-estimation.
  
  # dta: data of interest
  # tmt: set to be 1 by default
  # bootsize: total number of bootstrapping
  
  index2=!is.na(dta$P2) & !is.na(dta$S2) & is.finite(dta$S2)
  s1index=c(1:nrow(dta))
  bootpsi=matrix(NA,bootsize,4)
  
  for(i in 1:bootsize){
    set.seed(i)
    new_s1index=sample(s1index,replace = TRUE,size=length(s1index))
    news1=dta[new_s1index,]
    bootpsi[i,]=g_unshared(news1,tmt)$psi_SA
  }
  return(list(bootpsi=bootpsi))
}

# Calculate bootstrapping statistics with 1000 samples
boot_result=boot_g(dta,tmt=1,bootsize = 1000)
boot_result=boot_g(dta,tmt=2,bootsize = 1000)

# Collect relevent results
bootmean=c(mean(boot_result$bootpsi[,1]),mean(boot_result$bootpsi[,2]),mean(boot_result$bootpsi[,3]),mean(boot_result$bootpsi[,4]))
bootsd=c(sd(boot_result$bootpsi[,1]),sd(boot_result$bootpsi[,2]),sd(boot_result$bootpsi[,3]),sd(boot_result$bootpsi[,4]))
