rm(list=ls())
  
# Install "DTRreg" package for unshared Q-learning and G-estimation
if(!"DTRreg" %in% installed.packages()){
  install.packages("DTRreg")
}
library(DTRreg)

# Define expit function
expit=function(x){1/(1+exp(-x))}


# Data Generation ---------------------------------------------------------

dta_gen=function(n){
  # n: Simulation sample size.
  
  # Covariates
  X1=rnorm(n,10,5)
  X2=rnorm(n,1.25*X1,5)
  
  # Treatment
  A1=rbinom(n,1,expit(0.05*X1));mean(A1)
  A2=rbinom(n,1,expit(-0.05*X2));mean(A2)
  
  # Psi, parameter of interest
  psi0=8
  psi1=-1.2
  psi2=8
  
  # Model Specification and Outcome generation
  
  # Blip
  gamma2=A2*(psi0+psi1*X2+psi2*A1)
  gamma1=A1*(psi0+psi1*X1)
  
  # Regret 
  maxgamma2=0.5*(abs(psi0+psi1*X2+psi2*A1)+psi0+psi1*X2+psi2*A1)
  mu2=maxgamma2-gamma2
  maxgamma1=0.5*(abs(psi0+psi1*X1)+psi0+psi1*X1)
  mu1=maxgamma1-gamma1
  
  # Beta, parameter of treatment free outcome
  beta10=30
  beta11=3
  
  # Outcome
  Y=rnorm(n,beta10+beta11*X1,60)-mu2-mu1 #G1=Y+mu1+mu2, hence Y=G1-mu1-mu2.
  
  # Data of Simulation
  sim.data=data.frame(X1,X2,A1,A2,Y)
  
  return(sim.data)
}



# Unshared Q-learning -----------------------------------------------------

output_unsharedq=function(scene,sample,n){
  # A function estimate mean, sd, bias, MSE and matching rate using unshared q-learning.
  
  # scene: when scene = 1 use a regular non-flexible model specification, when scene = 2 use a flexible model specification.
  # sample: size of total samples used.
  # n: Simulated sample size.
      
  psi10=rep(0,sample)
  psi11=rep(0,sample)
  psi20=rep(0,sample)
  psi21=rep(0,sample)
  psi22=rep(0,sample)
  
  # Loop through all the samples, and calculate estimate for each sample.
  for(i in 1:sample){

    # Using a pre-fixed sequence for the random seeds to ensure the same estimates each time running this script.
    set.seed(i)

    # Generate the simulated data with size n.
    sim.data=dta_gen(n)
    
    # Extract variables from simluated data
    X1=sim.data$X1
    X2=sim.data$X2
    A1=sim.data$A1
    A2=sim.data$A2
    Y =sim.data$Y
    X1sq=X1^2
    A1X1sq=A1*X1^2
    X2sq=X2^2
    
    # Specify blip model, treatment model and treatment free outcome model for the 1st and 2nd simulation respectively.
    if(scene==1){ #non-flex model
      blip.mod=list(~X1,~X2+A1)
      treat.mod=list(A1~X1,A2~X2)
      tf.mod=list(~X1,~X2)
      H2=cbind(1,X2,A2,A2*X2,A2*A1)
    }else if(scene==2){ #flex model
      blip.mod=list(~X1,~X2+A1)
      treat.mod=list(A1~X1,A2~X2)
      tf.mod=list(~X1,~X1sq+A1X1sq+X2+X2sq)
      H2=cbind(1,X1sq,A1X1sq,X2,X2sq,A2,A2*X2,A2*A1)
    }
    
    # Fit the model using DTRreg::DTRreg with q-learning.
    mod1 <- DTRreg(Y, blip.mod, treat.mod, tf.mod, method = "qlearn")
    
    # Extract the estimates of Psi.
    psi10[i]=sapply(mod1$psi[1], function(x){as.numeric(x[1])})
    psi11[i]=sapply(mod1$psi[1], function(x){as.numeric(x[2])})
    
    psi20[i]=sapply(mod1$psi[2], function(x){as.numeric(x[1])})
    psi21[i]=sapply(mod1$psi[2], function(x){as.numeric(x[2])})
    psi22[i]=sapply(mod1$psi[2], function(x){as.numeric(x[3])})
    
    beta.psi=c(mod1$beta[[2]],mod1$psi[[2]])
  }
  
  # Summarizing estimates over all the samples
  psi0=8
  psi1=-1.2
  psi2=8
  
  hatpsi10=mean(psi10)
  hatpsi11=mean(psi11)
  hatpsi20=mean(psi20)
  hatpsi21=mean(psi21)
  hatpsi22=mean(psi22)
  
  # Simple average estimates.
  hatpsi0_sa=(hatpsi10+hatpsi20)/2;hatpsi0_sa
  hatpsi1_sa=(hatpsi11+hatpsi21)/2;hatpsi1_sa
  hatpsi2_sa=hatpsi22;hatpsi2_sa
  varpsi0_sa=0.25*(var(psi10)+var(psi20));varpsi0_sa
  varpsi1_sa=0.25*(var(psi11)+var(psi21));varpsi1_sa
  varpsi2_sa=var(psi22);varpsi2_sa
  hatpsi_sa=c(hatpsi0_sa,hatpsi1_sa,hatpsi2_sa)
  varpsi_sa=c(varpsi0_sa,varpsi1_sa,varpsi2_sa)
  bias_sa=c(hatpsi0_sa-psi0,hatpsi1_sa-psi1,hatpsi2_sa-psi2);bias_sa
  
  # Inverse variance weighted average estimates.
  hatpsi0_ivwa=((1/var(psi10))*hatpsi10+(1/var(psi20))*hatpsi20)/(1/var(psi10)+1/var(psi20));hatpsi0_ivwa
  hatpsi1_ivwa=((1/var(psi11))*hatpsi11+(1/var(psi21))*hatpsi21)/(1/var(psi11)+1/var(psi21));hatpsi1_ivwa
  hatpsi2_ivwa=hatpsi22;hatpsi2_ivwa
  varpsi0_ivwa=1/(1/var(psi10)+1/var(psi20));varpsi0_ivwa
  varpsi1_ivwa=1/(1/var(psi11)+1/var(psi21));varpsi1_ivwa
  varpsi2_ivwa=var(psi22);varpsi2_ivwa
  hatpsi_ivwa=c(hatpsi0_ivwa,hatpsi1_ivwa,hatpsi2_ivwa)
  varpsi_ivwa=c(varpsi0_ivwa,varpsi1_ivwa,varpsi2_ivwa)
  bias_ivwa=c(hatpsi0_ivwa-psi0,hatpsi1_ivwa-psi1,hatpsi2_ivwa-psi2);bias_ivwa
  
  # Collect all relevent results.
  results=list()
  results$hatpsi_s1=c(hatpsi10,hatpsi11)
  results$sdhatpsi_s1=c(sqrt(var(psi10)),sqrt(var(psi11)))
  results$hatpsi_s2=c(hatpsi20,hatpsi21,hatpsi22)
  results$sdhatpsi_s2=c(sqrt(var(psi20)),sqrt(var(psi21)),sqrt(var(psi22)))
  
  results$hatpsi_sa=hatpsi_sa
  results$sdpsi_sa=sqrt(varpsi_sa)
  results$bias_sa=bias_sa
  results$MSE_sa[1:3]<-varpsi_sa+bias_sa^2
  
  results$hatpsi_ivwa=hatpsi_ivwa
  results$sdpsi_ivwa=sqrt(varpsi_ivwa)
  results$bias_ivwa=bias_ivwa
  
  
  # Stage specific and overall matching rate.
  M_sa=rep(NA,sample)
  M_sa_overall=rep(NA,sample)
  
  M_ivwa=rep(NA,sample)
  M_ivwa_overall=rep(NA,sample)
  
  # Iterate through all samples, and calculate the correct matching rate.
  for(i in 1:sample){
    set.seed(i)
    
    X1=dta_gen(n)$X1
    X2=dta_gen(n)$X2
    A1=dta_gen(n)$A1
    A2=dta_gen(n)$A2
    
    s1_sa=(psi0+psi1*X1>0)==(hatpsi0_sa+hatpsi1_sa*X1>0)
    s2_sa=(psi0+psi1*X2+psi2*A1>0)==(hatpsi0_sa+hatpsi1_sa*X2+hatpsi2_sa*A1>0)
    
    s1_ivwa=(psi0+psi1*X1>0)==(hatpsi0_ivwa+hatpsi1_ivwa*X1>0)
    s2_ivwa=(psi0+psi1*X2+psi2*A1>0)==(hatpsi0_ivwa+hatpsi1_ivwa*X2+hatpsi2_ivwa*A1>0)
    
    s_sa=s1_sa*s2_sa
    s_ivwa=s1_ivwa*s2_ivwa
    
    M_sa[i]=(sum(s1_sa)+sum(s2_sa))/(2*n)
    M_sa_overall[i]=sum(s_sa)/n
    
    M_ivwa[i]=(sum(s1_ivwa)+sum(s2_ivwa))/(2*n)
    M_ivwa_overall[i]=sum(s_ivwa)/n
  }
  
  results$M_sa_stage<-mean(M_sa)
  results$M_sa_overall<-mean(M_sa_overall)
  
  results$M_ivwa_stage<-mean(M_ivwa)
  results$M_ivwa_overall<-mean(M_ivwa_overall)
  
  return(results)
}

# Unshared G-estimation ---------------------------------------------------

output_unsharedg=function(scene,sample,n){
  # A function estimate mean, sd, bias, MSE and matching rate using unshared g-estimation.
  
  # scene: when scene = 1 use a regular non-flexible model specification, when scene = 2 use a flexible model specification.
  # sample: Size of sample used in total.
  # n: Simulated sample size.
  
  psi10=rep(0,sample)
  psi11=rep(0,sample)
  psi20=rep(0,sample)
  psi21=rep(0,sample)
  psi22=rep(0,sample)
  
  # Loop through all samples, and calculate estimate for each sample.
  for(i in 1:sample){
    
    # Using a pre-fixed sequence for the random seeds to ensure the same estimates each time running this script.
    set.seed(i)
    
    # Generate the simulated data with size n.
    sim.data=dta_gen(n)
    
    # Extract variables from simluated data.
    X1=sim.data$X1
    X2=sim.data$X2
    A1=sim.data$A1
    A2=sim.data$A2
    Y =sim.data$Y
    X1sq=X1^2
    A1X1sq=A1*X1^2
    X2sq=X2^2
    
    # Specify blip model, treatment model and treatment free outcome model for the first and second simulation respectively.
    if(scene==1){ #non-flex model
      blip.mod=list(~X1,~X2+A1)
      treat.mod=list(A1~X1,A2~X2)
      tf.mod=list(~X1,~X2)
      H2=cbind(1,X2,A2,A2*X2,A2*A1)
    }else if(scene==2){ #flex model
      blip.mod=list(~X1,~X2+A1)
      treat.mod=list(A1~X1,A2~X2)
      tf.mod=list(~X1,~X1sq+A1X1sq+X2+X2sq)
      H2=cbind(1,X1sq,A1X1sq,X2,X2sq,A2,A2*X2,A2*A1)
    }
    
    # Fit the model using DTRreg::DTRreg with q-learning.
    mod1 <- DTRreg(Y, blip.mod, treat.mod, tf.mod, method = "gest")
    
    # Extract the estimates of Psi.
    psi10[i]=sapply(mod1$psi[1], function(x){as.numeric(x[1])})
    psi11[i]=sapply(mod1$psi[1], function(x){as.numeric(x[2])})
    psi20[i]=sapply(mod1$psi[2], function(x){as.numeric(x[1])})
    psi21[i]=sapply(mod1$psi[2], function(x){as.numeric(x[2])})
    psi22[i]=sapply(mod1$psi[2], function(x){as.numeric(x[3])})
    
    beta.psi=c(mod1$beta[[2]],mod1$psi[[2]])
  }

  # Summarizing estimates over all the samples
  psi0=8
  psi1=-1.2
  psi2=8
  
  hatpsi10=mean(psi10)
  hatpsi11=mean(psi11)
  hatpsi20=mean(psi20)
  hatpsi21=mean(psi21)
  hatpsi22=mean(psi22)
  
  # Simple average estimates.
  hatpsi0_sa=(hatpsi10+hatpsi20)/2;hatpsi0_sa
  hatpsi1_sa=(hatpsi11+hatpsi21)/2;hatpsi1_sa
  hatpsi2_sa=hatpsi22;hatpsi2_sa
  varpsi0_sa=0.25*(var(psi10)+var(psi20));varpsi0_sa
  varpsi1_sa=0.25*(var(psi11)+var(psi21));varpsi1_sa
  varpsi2_sa=var(psi22);varpsi2_sa
  hatpsi_sa=c(hatpsi0_sa,hatpsi1_sa,hatpsi2_sa)
  varpsi_sa=c(varpsi0_sa,varpsi1_sa,varpsi2_sa)
  bias_sa=c(hatpsi0_sa-psi0,hatpsi1_sa-psi1,hatpsi2_sa-psi2);bias_sa
  
  # Inverse variance weighted average estimates.
  hatpsi0_ivwa=((1/var(psi10))*hatpsi10+(1/var(psi20))*hatpsi20)/(1/var(psi10)+1/var(psi20));hatpsi0_ivwa
  hatpsi1_ivwa=((1/var(psi11))*hatpsi11+(1/var(psi21))*hatpsi21)/(1/var(psi11)+1/var(psi21));hatpsi1_ivwa
  hatpsi2_ivwa=hatpsi22;hatpsi2_ivwa
  varpsi0_ivwa=1/(1/var(psi10)+1/var(psi20));varpsi0_ivwa
  varpsi1_ivwa=1/(1/var(psi11)+1/var(psi21));varpsi1_ivwa
  varpsi2_ivwa=var(psi22);varpsi2_ivwa
  hatpsi_ivwa=c(hatpsi0_ivwa,hatpsi1_ivwa,hatpsi2_ivwa)
  varpsi_ivwa=c(varpsi0_ivwa,varpsi1_ivwa,varpsi2_ivwa)
  bias_ivwa=c(hatpsi0_ivwa-psi0,hatpsi1_ivwa-psi1,hatpsi2_ivwa-psi2);bias_ivwa
  
  # Collect all relevent results.
  results=list()
  
  results$hatpsi_s1=c(hatpsi10,hatpsi11)
  results$sdhatpsi_s1=c(sqrt(var(psi10)),sqrt(var(psi11)))
  results$hatpsi_s2=c(hatpsi20,hatpsi21,hatpsi22)
  results$sdhatpsi_s2=c(sqrt(var(psi20)),sqrt(var(psi21)),sqrt(var(psi22)))
  
  results$hatpsi_sa=hatpsi_sa
  results$sdpsi_sa=sqrt(varpsi_sa)
  results$bias_sa=bias_sa
  results$MSE_sa <- varpsi_sa+bias_sa^2
  
  results$hatpsi_ivwa=hatpsi_ivwa
  results$sd_ivwa=sqrt(varpsi_ivwa)
  results$bias_ivwa=bias_ivwa
  
  # Stage specific and overall matching rate.
  M_sa=rep(NA,sample)
  M_sa_overall=rep(NA,sample)
  
  M_ivwa=rep(NA,sample)
  M_ivwa_overall=rep(NA,sample)
  
  # Iterate through all samples, and calculate the correct matching rate.
  for(i in 1:sample){
    set.seed(i)
    
    X1=dta_gen(n)$X1
    X2=dta_gen(n)$X2
    A1=dta_gen(n)$A1
    A2=dta_gen(n)$A2
    
    s1_sa=(psi0+psi1*X1>0)==(hatpsi0_sa+hatpsi1_sa*X1>0)
    s2_sa=(psi0+psi1*X2+psi2*A1>0)==(hatpsi0_sa+hatpsi1_sa*X2+hatpsi2_sa*A1>0)
    
    s1_ivwa=(psi0+psi1*X1>0)==(hatpsi0_ivwa+hatpsi1_ivwa*X1>0)
    s2_ivwa=(psi0+psi1*X2+psi2*A1>0)==(hatpsi0_ivwa+hatpsi1_ivwa*X2+hatpsi2_ivwa*A1>0)
    
    s_sa=s1_sa*s2_sa
    s_ivwa=s1_ivwa*s2_ivwa
    
    M_sa[i]=(sum(s1_sa)+sum(s2_sa))/(2*n)
    M_sa_overall[i]=sum(s_sa)/n
    
    M_ivwa[i]=(sum(s1_ivwa)+sum(s2_ivwa))/(2*n)
    M_ivwa_overall[i]=sum(s_ivwa)/n
  }
  
  results$M_sa_stage<-mean(M_sa)
  results$M_sa_overall<-mean(M_sa_overall)
  
  results$M_ivwa_stage<-mean(M_ivwa)
  results$M_ivwa_overall<-mean(M_ivwa_overall)
  
  return(results)
}

# Shared Q-learning ------------------------------------------------------

# Pseudo Outcome of Q-learning function
q_pseudo=function(dta,scene,beta.psi=NULL){
  # A function calculate the pseudo outcome of shared Q-learning, treatment free outcome + optimal blip.
  
  # dta: a data.frame of interest
  # scene: 1 or 2 for different model specifications, corresponds to simulation table 1 and 2.
  # beta.psi: Initial value, set to Null by default.
  
  Y=dta$Y 
  A1=dta$A1
  A2=dta$A2
  X1=dta$X1
  X2=dta$X2
  
  X1sq=X1^2
  A1X1sq=A1*X1^2
  X2sq=X2^2

  A2X2=A2*X2
  A2A1=A2*A1
  
  betapsi.2 <- NULL

  # If initial value beta.psi is provided, use that, otherwise calculate from linear regression.
  if (scene == 1) {
    if (!is.null(beta.psi)) { 
      betapsi.2 = beta.psi[c(1:2, 5:7)] #1:2 are beta, 5:7 are shared psi0, psi1, psi2.
    } else{
      betapsi.2 = lm(Y ~ 1 + X2 + A2 + A2X2 + A2A1)$coefficients
    }
    
    Y1.pseudo = betapsi.2[1] + betapsi.2[2] * X2 +
      0.5 * ((betapsi.2[3] + betapsi.2[4] * X2 + betapsi.2[5] * A1) + 
               abs(betapsi.2[3] +betapsi.2[4] * X2 + betapsi.2[5] * A1))
  } else if (scene == 2) {
    if (!is.null(beta.psi)) {
      betapsi.2 = beta.psi[c(1:5, 8:10)] #1:5 are beta, 8:10 are shared psi0, psi1, psi2.
    } else{
      betapsi.2 = lm(Y ~ 1 + X1sq + A1X1sq + X2 + X2sq + A2 + A2X2 + A2A1)$coefficients
    }
    
    Y1.pseudo = betapsi.2[1] + betapsi.2[2] * X1sq + betapsi.2[3] * A1X1sq +
      betapsi.2[4] * X2 + betapsi.2[5] * X2sq +
      0.5 * ((betapsi.2[6] + betapsi.2[7] * X2 + betapsi.2[8] * A1) + 
               abs(betapsi.2[6] +betapsi.2[7] * X2 + betapsi.2[8] * A1))
  }
  
  result=list()
  result$Y1.pseudo=Y1.pseudo
  return(result)
}

# Shared Q-learning function
sharedq<-function(dta,scene,ini="ZERO"){ 
  # Function of shared Q-learning
  
  # dta: data of interest.
  # scene: 1 or 2, for different model specifications, corresponds to simulation 1 and 2, respectively.
  # ini: Initial value, set to be a zero valued vector by default.
  
  Y=dta$Y 
  A1=dta$A1
  A2=dta$A2
  X1=dta$X1
  X2=dta$X2
  
  # Initialize with different scenario and initial value.
  initial=list()
  if (scene == 1) {
    if (ini == "ZERO") {
      initial$beta.psi[1:7] = rep(0, 7) 
    } else{
      initial$beta.psi[1:7] = ini
    }
  } else if (scene == 2) {
    if (ini == "ZERO") {
      initial$beta.psi[1:10] = rep(0, 10)
    } else{
      initial$beta.psi[1:10] = ini
    }
  }
  
  # Call q_pseudo function to calculate pseudo outcome to start up.
  Y1.pseudo=q_pseudo(dta,scene,beta.psi=initial$beta.psi)$Y1.pseudo;Y1.pseudo
  
  # Construct overall outcome matrix Y.star to start up, by stacking the real outcome and pseudo outcome.
  Y.star=matrix(c(Y,Y1.pseudo),,1);Y.star

  # Construct covariates matrix Z.
  if(scene==1){
    H1=cbind(0,0,1,X1,A1,A1*X1,0)
    H2=cbind(1,X2,0,0,A2,A2*X2,A2*A1)
    psi.start=5
  } else if(scene==2){
    H1=cbind(0,0,0,0,0,1,X1,A1,A1*X1,0)
    H2=cbind(1,X1^2,A1*X1^2,X2,X2^2,0,0,A2,A2*X2,A2*A1) #flexible terms added in H2, increase estimate
    psi.start=8
  }
  Z=rbind(H2,H1)
  
  # Solve for a initial beta.psi.
  beta.psi=qr.solve(Z,Y.star)
  
  # Extract initial Psi.
  psi=beta.psi[psi.start:length(beta.psi)]
  
  # Initialize variables.
  psi.prev=c(0,0,0)
  temp.beta.psi=rep(0,length(beta.psi))
  diff=1e3
  error=1e-5
  iter=1 #since we already ran the first iter. hence 1 not 0
  
  # Repeatedly solve for beta.psi until converge.
  while (diff>error & (iter<=1000) ){
    # Log the number of iteration takes to converge.
    iter=iter+1
    # Log psi and beta.psi from previous iteration.
    psi.prev=psi
    temp.beta.psi<-beta.psi
    # Calculate pseudo outcome
    Y1.pseudo <- q_pseudo(dta,scene,beta.psi)$Y1.pseudo
    # Construct the new outcome matrix for this iteration
    Y.star=matrix(c(Y,Y1.pseudo),,1)
    # Solve beta.psi for this iteration
    beta.psi <- qr.solve(Z,Y.star) #update
    # Update Psi
    psi=beta.psi[psi.start:length(beta.psi)] 
    # Calculate the difference of beta.psi between two iterations.
    diff=sum(abs(beta.psi-temp.beta.psi))
  }
  
  result=list()
  result$beta.psi <- beta.psi
  result$iter<-iter
  result$psi<-psi
  return(result)
}

# Shared Q-learning estimation
output_sharedq=function(scene,sample,n,ini){
  # A function estimate mean, sd, bias, MSE and matching rate using shared Q-learning.
  
  # scene: 1 or 2, for different model specifications, corresponds to simulation 1 and 2, respectively.
  # sample: Size of total sample used.
  # n: Simulated sample size.
  # ini: Initial value used.
  
  psi0=8
  psi1=-1.2
  psi2=8
  
  if(scene==1){
    betapsi=matrix(NA,sample,7)
  }else{ 
    betapsi=matrix(NA,sample,10)
    }
  
  iter=rep(NA,sample)
  psi=matrix(NA,sample,3)
  M=rep(NA,sample)
  M_overall=rep(NA,sample)
  
  # Iterate through all samples and estimate psi.
  for(i in 1:sample){
    set.seed(i)
    dta=dta_gen(n)
    
    temp=sharedq(dta,scene,ini) #temp store results
    iter[i]=temp$iter
    betapsi[i,]=temp$beta.psi
    psi[i,]=temp$psi
  }
  
  meaniter=mean(iter)
  
  hatpsi0=mean(psi[,1]);hatpsi0
  hatpsi1=mean(psi[,2]);hatpsi1
  hatpsi2=mean(psi[,3]);hatpsi2
  
  varpsi0=var(psi[,1]);varpsi0
  varpsi1=var(psi[,2]);varpsi1
  varpsi2=var(psi[,3]);varpsi2
  
  bias0=hatpsi0-psi0
  bias1=hatpsi1-psi1
  bias2=hatpsi2-psi2
  
  hatpsi=c(hatpsi0,hatpsi1,hatpsi2)
  varpsi=c(varpsi0,varpsi1,varpsi2)
  bias=c(bias0,bias1,bias2)
  
  # Calculate matching rate.
  for(i in 1:sample){
    set.seed(i)
    dta=dta_gen(n)
    X1=dta$X1
    X2=dta$X2
    A1=dta$A1
    A2=dta$A2
    
    s1=(psi0+psi1*X1>0)==(hatpsi0+hatpsi1*X1>0)
    s2=(psi0+psi1*X2+psi2*A1>0)==(hatpsi0+hatpsi1*X2+hatpsi2*A1>0)
    
    s=s1*s2
    M[i]=(sum(s1)+sum(s2))/(2*n)
    M_overall[i]=sum(s)/n
  }
  
  # Collect all relevent results
  results=list()
  results$hatpsi=hatpsi
  results$sdpsi=sqrt(varpsi)
  results$bias=bias
  results$M_stage=mean(M)
  results$M_overall=mean(M_overall)
  results$MSE_p=varpsi+bias^2

  return(results)
}


# Shared G-estimation -----------------------------------------------------

sharedg=function(data,scene,ini){
  # Function of shared G-estimation
  
  # dta: data of interest.
  # scene: 1 or 2, for different model specifications, corresponds to simulation 1 and 2, respectively.
  # ini: Initial value, set to be a zero valued vector by default.
  
  X1=data$X1
  X2=data$X2
  A1=data$A1
  A2=data$A2
  Y =data$Y
  n=length(Y)
  
  # Correct blip
  H2_psi=cbind(1,X2,A1)
  H1_psi=cbind(1,X1,0) 
  
  # Specify model for scene 1 and 2.
  if (scene==1)
  {
    alpha1=glm(A1~X1-1,family = "binomial")$coefficients;alpha1
    alpha2=glm(A2~X2-1,family = "binomial")$coefficients;alpha2
    
    EA1=expit(alpha1*X1) 
    EA2=expit(alpha2*X2)
    
    H1_beta=cbind(1,X1)
    H2_beta=cbind(1,X2) 
  } else if (scene==2)
  {
    alpha1=glm(A1~X1-1,family = "binomial")$coefficients;alpha1
    alpha2=glm(A2~X2-1,family = "binomial")$coefficients;alpha2
    
    EA1=expit(alpha1*X1) 
    EA2=expit(alpha2*X2)
    
    H2_beta=cbind(1,X1^2,A1*X1^2,X2,X2^2)
    H1_beta=cbind(1,X1)
  }
  
  # Converge track of psi.
  psi0_track=c()
  psi1_track=c()
  psi2_track=c()
  
  # Initialize variables.
  ini_psi=ini
  error=1e-6
  diff=1e3
  psi=ini_psi
  j=0
  
  if(scene==1){
    temp_betapsi=rep(0,ncol(H1_beta)+ncol(H2_psi))
  }else if(scene==2){
    temp_betapsi=rep(0,ncol(H1_beta)+ncol(H2_beta)+ncol(H2_psi))
  }
  
  # Repeat until converge
  while(diff > error){
    
    j=j+1
    
    # Construct outcome variables
    Y2=Y
    Y1=Y-A2*(H2_psi%*%psi)+0.5*(H2_psi%*%psi+abs(H2_psi%*%psi))
    Ypseu=c(Y1,Y2)
    
    # Construct covariates
    if(scene==1){
      Xbeta=rbind(H1_beta,H2_beta)
    }else if(scene==2){
      Xbeta=rbind(cbind(H1_beta,matrix(0,nrow=nrow(H2_beta),ncol=ncol(H2_beta))),
                  cbind(matrix(0,nrow=nrow(H1_beta),ncol=ncol(H1_beta)),H2_beta))
    }

    AXpsi=rbind(A1*H1_psi,A2*H2_psi);dim(AXpsi)
    Xdelta=cbind(Xbeta,AXpsi);dim(Xdelta)
    
    XpsiAE=rbind(H1_psi*(A1-EA1), H2_psi*(A2-EA2));dim(XpsiAE)
    
    Xomega=cbind(Xbeta,XpsiAE);dim(Xomega)
    
    # Store the temporary beta psi
    old_betapsi=temp_betapsi
    # Solve beta psi for the current iteration
    temp_betapsi=solve(t(Xomega)%*%Xdelta)%*%t(Xomega)%*%Ypseu
  
    # Extract psi and beta 
    if(scene==1){
      psi=temp_betapsi[(ncol(H1_beta)+1):length(temp_betapsi)]
      beta2=temp_betapsi[1:ncol(H2_beta)]
    }else if(scene==2){
      psi=temp_betapsi[(ncol(H1_beta)+ncol(H2_beta)+1):length(temp_betapsi)]
      beta2=temp_betapsi[(ncol(H1_beta)+1):(ncol(H1_beta)+ncol(H2_beta))]
    }

    psi0_track=c(psi0_track,psi[1])
    psi1_track=c(psi1_track,psi[2])
    psi2_track=c(psi2_track,psi[3])
    
    # Calculate the difference between betapsi from current and previous iteration.
    diff=sum(abs(temp_betapsi-old_betapsi))
  }
  
  return(list(iter=j,betapsi=temp_betapsi,psi0_track=psi0_track,psi1_track=psi1_track,psi2_track=psi2_track))
}

output_shareg=function(scene,sample,n,ini){
  # A function estimate mean, sd, bias, MSE and matching rate using shared G-estimation.
  
  # scene: 1 or 2, for different model specifications, corresponds to simulation 1 and 2, respectively.
  # sample: Size of total sample used.
  # n: Simulated sample size.
  # ini: Initial value used.
  
  psi0=8
  psi1=-1.2
  psi2=8
  psi=c(psi0,psi1,psi2)
  
  if(scene==1){
    betapsi=matrix(NA,sample,5)
  } else if(scene==2){
    betapsi=matrix(NA,sample,10)
  }

  # Iterate through all sample, and call sharedg to estimate psi for each sample. 
  for(i in 1:sample){
    set.seed(i)
    sim.data=dta_gen(n)
    temp=sharedg(sim.data,scene,ini)
    betapsi[i,]=temp$betapsi
  }

  # Extract psi
  if(scene==1){
    hatpsi0=mean(betapsi[,3]);hatpsi0
    hatpsi1=mean(betapsi[,4]);hatpsi1
    hatpsi2=mean(betapsi[,5]);hatpsi2
    
    varpsi0=var(betapsi[,3]);varpsi0
    varpsi1=var(betapsi[,4]);varpsi1
    varpsi2=var(betapsi[,5]);varpsi2
  } else if(scene==2){
    hatpsi0=mean(betapsi[,8]);hatpsi0
    hatpsi1=mean(betapsi[,9]);hatpsi1
    hatpsi2=mean(betapsi[,10]);hatpsi2
    
    varpsi0=var(betapsi[,8]);varpsi0
    varpsi1=var(betapsi[,9]);varpsi1
    varpsi2=var(betapsi[,10]);varpsi2
  }
  
  hatpsi=c(hatpsi0,hatpsi1,hatpsi2)
  varpsi=c(varpsi0,varpsi1,varpsi2)
  biaspsi=hatpsi-psi
  
  # Iterate through all samples, derive matching rate and calculate the average over all samples.
  M=rep(NA,100)
  M_overall=rep(NA,100)
  for(i in 1:100){
    set.seed(i)
    
    X1=dta_gen(n)$X1
    X2=dta_gen(n)$X2
    A1=dta_gen(n)$A1
    A2=dta_gen(n)$A2
    
    s1=(psi0+psi1*X1>0)==(hatpsi0+hatpsi1*X1>0)
    s2=(psi0+psi1*X2+psi2*A1>0)==(hatpsi0+hatpsi1*X2+hatpsi2*A1>0)
    
    s=s1*s2
    M[i]=(sum(s1)+sum(s2))/(2*n)
    M_overall[i]=sum(s)/n
  }
  
  # Collect all relevent results.
  results=list()
  results$hatpsi=hatpsi
  results$sdpsi=sqrt(varpsi)
  results$bias=biaspsi
  results$M_stage=mean(M)
  results$M_overall=mean(M_overall)
  results$MSE_p=biaspsi^2+varpsi
  
  return(results)
}
