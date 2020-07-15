rm(list=ls())
set.seed(100)

# Model:
# beta.psi=c(beta20,beta21,beta22,beta23,beta24,beta10,beta11,beta12,beta13,psi0,psi1,psi2,psi3)
# Q2=beta20+beta21*Q2+beta22*S2+beta23*P2+beta24*A1 + A2*(psi0+psi1*Q2+psi2*S2)
# Q1=beta10+beta11*Q1+beta12*S1+beta13*P1           + A1*(psi0+psi1*Q1+psi2*S1+psi3*P1)

q_unshared=function(dta){
  # function of unshared Q-learning.
  
  # dta: data of interest
  
  Q1=dta$Q1
  Q2=dta$Q2
  A1=dta$A1
  A2=dta$A2
  S1=dta$S1
  S2=dta$S2
  P1=dta$P1
  P2=dta$P2
  Y =dta$Y
  
  dta$A2Q2=dta$A2*dta$Q2
  dta$A2S1=dta$A2*dta$S1
  
  dta$A1Q1=dta$A1*dta$Q1
  dta$A1S1=dta$A1*dta$S1
  dta$A1P1=dta$A1*dta$P1
  
  # Keep only participants having P2 and S2 for stage 2, S2 and P2 are key variables indicating whether someone got treated at stage 2.
  index2=!is.na(dta$P2) & !is.na(dta$S2) & is.finite(dta$S2) #discard NA/Inf S2/S1, choose who entered Stage2
  temp2=dta[index2,]

  # Derive beta psi of stage 2 using regression  
  Q2mod=lm(Y~1+Q2+S2+P2+A1+A2+A2Q2+A2S1,data=temp2)
  betapsi.2=Q2mod$coefficients
  
  # Calculate Stage 1 pseudo-outcome
  Y1.pseu=betapsi.2[1]+betapsi.2[2]*Q2+betapsi.2[3]*S2+betapsi.2[4]*P2+betapsi.2[5]*A1+
    0.5*(betapsi.2[6]+betapsi.2[7]*Q2+betapsi.2[8]*S2+abs(betapsi.2[6]+betapsi.2[7]*Q2+betapsi.2[8]*S2))
  Y1.pseu[!index2]=dta[!index2,]$Y
  
  # Derive beta psi of stage 1 using regression
  Q1mod=lm(Y1.pseu~1+Q1+S1+P1+A1+A1Q1+A1S1+A1P1,data=dta)
  betapsi.1=Q1mod$coefficients
  
  # Collect results
  results=list()
  results$psi.s2=betapsi.2[6:8]
  results$psi.s1=betapsi.1[5:8]
  
  # Simple average estimates
  psi0_SA=(betapsi.1[5]+betapsi.2[6])/2
  psi1_SA=(betapsi.1[6]+betapsi.2[7])/2
  psi2_SA=(betapsi.1[7]+betapsi.2[8])/2
  psi3_SA= betapsi.1[8]
  
  results$psi_SA=as.numeric(c(psi0_SA,psi1_SA,psi2_SA,psi3_SA))
  
  return(results)
}


boot_g=function(dta,bootsize){
  # function of bootstrapping unshared Q-learning.
  
  index2=!is.na(dta$P2) & !is.na(dta$S2) & is.finite(dta$S2)
  s1index=c(1:nrow(dta))
  bootpsi=matrix(NA,bootsize,4)
  
  for(i in 1:bootsize){
    set.seed(i)
    new_s1index=sample(s1index,replace = TRUE,size=length(s1index))
    news1=dta[new_s1index,]
    bootpsi[i,]=q_unshared(news1)$psi_SA
  }
  return(list(bootpsi=bootpsi))
}

# Calculate bootstrapping statistics with 1000 samples
boot_result=boot_g(dta,bootsize = 1000)

# Collect relevent results
bootmean=c(mean(boot_result$bootpsi[,1]),mean(boot_result$bootpsi[,2]),mean(boot_result$bootpsi[,3]),mean(boot_result$bootpsi[,4]))
bootsd=c(sd(boot_result$bootpsi[,1]),sd(boot_result$bootpsi[,2]),sd(boot_result$bootpsi[,3]),sd(boot_result$bootpsi[,4]))

