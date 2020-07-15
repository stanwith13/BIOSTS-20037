rm(list=ls())
set.seed(100)
library(corpcor) # pseudoinverse

# Model:
#beta.psi=c(beta20,beta21,beta22,beta23,beta24,beta10,beta11,beta12,beta13,psi0,psi1,psi2,psi3)
#Q2=beta20+beta21*Q2+beta22*S2+beta23*P2+beta24*A1 + A2*(psi0+psi1*Q2+psi2*S2)
#Q1=beta10+beta11*Q1+beta12*S1+beta13*P1           + A1*(psi0+psi1*Q1+psi2*S1+psi3*P1)

q_pseudo=function(dta,beta.psi=NULL){
  # function calculating the pseudo outcome for Q-learning
  
  # dta: data of interest
  # beta.psi: the estimated beta and psi used for the outcome, set to be NULL by default
  
  Q1=dta$Q1
  Q2=dta$Q2
  A1=dta$A1
  A2=dta$A2
  S1=dta$S1
  S2=dta$S2
  P1=dta$P1
  P2=dta$P2
  Y =dta$Y
  
  # Keep only participants having P2 and S2 for stage 2
  index2=!is.na(dta$P2) & !is.na(dta$S2) & is.finite(dta$S2) #discard NA/Inf S2/S1, choose who entered Stage2
  temp2=dta[index2,] 
  Y.S2=dta[index2,]$Y 
  
  betapsi.2 <- NULL
  
  # Use beta.psi if provided, otherwise estimate from as regression
  if(!is.null(beta.psi)){
    betapsi.2=beta.psi[c(1:5,10:12)] #5+3
  }else{
    betapsi.2=as.numeric(lm(Y.S2~1+Q2+S2+P2+A1+ A2+A2*Q2+A2*S2,data=temp2)$coefficients)
  } 
  
  # Pseudo outcome
  Y1.pseudo=betapsi.2[1]+betapsi.2[2]*Q2+betapsi.2[3]*S2+betapsi.2[4]*P2+betapsi.2[5]*A1+
    0.5*(betapsi.2[6]+betapsi.2[7]*Q2+betapsi.2[8]*S2+abs(betapsi.2[6]+betapsi.2[7]*Q2+betapsi.2[8]*S2))
  
  Y1.pseudo[!index2]=Y[!index2]
  result=list()
  result$Y1.pseudo=Y1.pseudo
  return(result)
}


q_shared<-function(dta,ini="ZERO"){ 
  # function of shared Q-learning.
  
  # dta: data of interest
  # ini: intial value used, set to be zero by default
  
  Q1=dta$Q1
  Q2=dta$Q2
  A1=dta$A1
  A2=dta$A2
  S1=dta$S1
  S2=dta$S2
  P1=dta$P1
  P2=dta$P2
  Y =dta$Y
  
  # Keep only participants having P2 and S2 for stage 2, S2 and P2 are key variables indicating whether someone got treated at stage 2.
  index2=!is.na(dta$P2) & !is.na(dta$S2) & is.finite(dta$S2)
  temp2=dta[index2,] #273, data for stage 2

  # Initialize the outcome matrix
  Y1.pseudo=q_pseudo(dta,ini)$Y1.pseudo;Y1.pseudo
  Y.S2=dta[index2,]$Y
  Y.star=matrix(c(Y.S2,Y1.pseudo),,1);Y.star
  
  # Specify blip model
  H2_psi=cbind(1,temp2$Q2,temp2$S2,0)
  H1_psi=cbind(1,Q1,S1,P1)
  
  # Specify expected treatment free outcome model
  H2_beta=cbind(1,temp2$Q2,temp2$S2,temp2$P2,temp2$A1)
  H1_beta=cbind(1,Q1,S1,P1)
  
  # Combine blip and expected tmt free to get matrix of Q-model
  H2=cbind(H2_beta,matrix(0,nrow=nrow(H2_beta),ncol=ncol(H1_beta)),H2_psi)
  H1=cbind(matrix(0,nrow=nrow(H1_beta),ncol=ncol(H2_beta)),H1_beta,H1_psi)
  Z=rbind(H2,H1)
  
  # Use Moore-Penrose for pseudoinverse to solve the betapsi
  beta.psi=pseudoinverse(t(Z)%*%Z)%*%(t(Z)%*%Y.star)#use
  psi=beta.psi[10:13]
  
  # Initialize variables for iteration
  temp.beta.psi=rep(0,13)
  diff=1e3
  error=1e-5
  iter=1 #since we already ran the first iter. hence 1 not 0
  
  # Iterate until converge
  while (diff>error & (iter<=100) ){
    iter=iter+1
    
    # Store previous beta psi
    temp.beta.psi=beta.psi
    
    # Update pseudo outcome 
    Y1.pseudo <- q_pseudo(dta,beta.psi)$Y1.pseudo
    Y.star=matrix(c(Y.S2,Y1.pseudo),,1)
    
    # Solve beta psi for current iteration
    beta.psi<-pseudoinverse(t(Z)%*%Z)%*%(t(Z)%*%Y.star) #update
    
    # Update psi
    psi=beta.psi[10:13]
    
    diff=sum(abs(temp.beta.psi-beta.psi))
  }
  result=list()
  result$beta.psi <- beta.psi
  result$iter<-iter
  result$psi<-psi
  return(result)
}


boot_g=function(dta,ini,bootsize){
  # function of bootstrapping shared Q-learning.
  
  # dta: data of interest
  # ini: initial value used for g-est
  # bootsize: total number of bootstrapping
  
  index2=!is.na(dta$P2) & !is.na(dta$S2) & is.finite(dta$S2)
  s1index=c(1:nrow(dta))
  bootpsi=matrix(NA,bootsize,4)
  
  for(i in 1:bootsize){
    set.seed(i)
    new_s1index=sample(s1index,replace = TRUE,size=length(s1index))
    news1=dta[new_s1index,]
    bootpsi[i,]=q_shared(news1,ini)$psi
  }
  return(list(bootpsi=bootpsi))
}

# Calculate bootstrapping statistics with 1000 samples
ini.zero=rep(0,13)
boot_result=boot_g(dta,ini.zero,bootsize = 1000)

# Collect relevent results
bootmean=c(mean(boot_result$bootpsi[,1]),mean(boot_result$bootpsi[,2]),mean(boot_result$bootpsi[,3]),mean(boot_result$bootpsi[,4]))
bootsd=c(sd(boot_result$bootpsi[,1]),sd(boot_result$bootpsi[,2]),sd(boot_result$bootpsi[,3]),sd(boot_result$bootpsi[,4]))
