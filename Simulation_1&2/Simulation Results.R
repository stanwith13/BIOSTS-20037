rm(list=ls())
source('Simulation Functions.R')
set.seed(100)
sample.size=1000

#' Call functions from "Simulation Functions" and simulate for different scenario and sample size.

##################################################################
##                          Section 1:                          ##
##                          Unshared Q                          ##
##################################################################

output_unsharedq_s1_200=output_unsharedq(scene = 1,sample=sample.size,n=200)
output_unsharedq_s1_500=output_unsharedq(scene = 1,sample=sample.size,n=500)
output_unsharedq_s1_2000=output_unsharedq(scene = 1,sample=sample.size,n=2000)
output_unsharedq_s1_5000=output_unsharedq(scene = 1,sample=sample.size,n=5000)

output_unsharedq_s2_200=output_unsharedq(scene = 2,sample=sample.size,n=200)
output_unsharedq_s2_500=output_unsharedq(scene = 2,sample=sample.size,n=500)
output_unsharedq_s2_2000=output_unsharedq(scene = 2,sample=sample.size,n=2000)
output_unsharedq_s2_5000=output_unsharedq(scene = 2,sample=sample.size,n=5000)

##################################################################
##                          Section 2:                          ##
##                          Unshared G                          ##
##################################################################

output_unsharedg_s1_200=output_unsharedg(scene = 1,sample=sample.size,n=200)
output_unsharedg_s1_500=output_unsharedg(scene = 1,sample=sample.size,n=500)
output_unsharedg_s1_2000=output_unsharedg(scene = 1,sample=sample.size,n=2000)
output_unsharedg_s1_5000=output_unsharedg(scene = 1,sample=sample.size,n=5000)

output_unsharedg_s2_200=output_unsharedg(scene = 2,sample=sample.size,n=200)
output_unsharedg_s2_500=output_unsharedg(scene = 2,sample=sample.size,n=500)
output_unsharedg_s2_2000=output_unsharedg(scene = 2,sample=sample.size,n=2000)
output_unsharedg_s2_5000=output_unsharedg(scene = 2,sample=sample.size,n=5000)

##################################################################
##                          Section 3:                          ##
##                           Shared Q                           ##
##################################################################
output_sharedq_s1_200=output_sharedq(scene = 1,sample=sample.size,n=200,ini="ZERO")
output_sharedq_s1_500=output_sharedq(scene = 1,sample=sample.size,n=500,ini="ZERO")
output_sharedq_s1_2000=output_sharedq(scene = 1,sample=sample.size,n=2000,ini="ZERO")
output_sharedq_s1_5000=output_sharedq(scene = 1,sample=sample.size,n=5000,ini="ZERO")

output_sharedq_s2_200=output_sharedq(scene = 2,sample=sample.size,n=200,ini="ZERO")
output_sharedq_s2_500=output_sharedq(scene = 2,sample=sample.size,n=500,ini="ZERO")
output_sharedq_s2_2000=output_sharedq(scene = 2,sample=sample.size,n=2000,ini="ZERO")
output_sharedq_s2_5000=output_sharedq(scene = 2,sample=sample.size,n=5000,ini="ZERO")


##################################################################
##                          Section 4:                          ##
##                           Shared G                           ##
##################################################################
output_sharedg_s1_200=output_shareg(scene = 1,sample=sample.size,n=200,ini=c(0,0,0))
output_sharedg_s1_500=output_shareg(scene = 1,sample=sample.size,n=500,ini=c(0,0,0))
output_sharedg_s1_2000=output_shareg(scene = 1,sample=sample.size,n=2000,ini=c(0,0,0))
output_sharedg_s1_5000=output_shareg(scene = 1,sample=sample.size,n=5000,ini=c(0,0,0))

output_sharedg_s2_200=output_shareg(scene = 2,sample=sample.size,n=200,ini=c(0,0,0))
output_sharedg_s2_500=output_shareg(scene = 2,sample=sample.size,n=500,ini=c(0,0,0))
output_sharedg_s2_2000=output_shareg(scene = 2,sample=sample.size,n=2000,ini=c(0,0,0))
output_sharedg_s2_5000=output_shareg(scene = 2,sample=sample.size,n=5000,ini=c(0,0,0))


# Save the results in an R object.
save.image("sim_results.RData")
load(here::here("sim_results.RData"))


# Summary of Results ------------------------------------------------------

fn.summary <- function(unshare_q, share_q, unshare_g, share_g){
  df_est <- data.frame(
    psi0 <- c(unshare_q$hatpsi_sa[1], 
              share_q$hatpsi[1], 
              unshare_g$hatpsi_sa[1], 
              share_g$hatpsi[1]),
    
    sd_psi0 <- c(unshare_q$sdpsi_sa[1], 
                 share_q$sdpsi[1], 
                 unshare_g$sdpsi_sa[1], 
                 share_g$sdpsi[1]),
    
    bias_psi0 <- c(unshare_q$bias_sa[1],
                   share_q$bias[1],
                   unshare_g$bias_sa[1],
                   share_g$bias[1]),
    
    MSE_psi0 <- c(unshare_q$MSE_sa[1], 
                  share_q$MSE_p[1],
                  unshare_g$MSE_sa[1], 
                  share_g$MSE_p[1]),
    
    psi1 <- c(unshare_q$hatpsi_sa[2],
               share_q$hatpsi[2], 
               unshare_g$hatpsi_sa[2], 
               share_g$hatpsi[2]),
    
    sd_psi1 <- c(unshare_q$sdpsi_sa[2], 
                  share_q$sdpsi[2], 
                  unshare_g$sdpsi_sa[2], 
                  share_g$sdpsi[2]),
    
    bias_psi1 <- c(unshare_q$bias_sa[2], 
                    share_q$bias[2],
                    unshare_g$bias_sa[2], 
                    share_g$bias[2]),
    
    MSE_psi1 <- c(unshare_q$MSE_sa[2],
                   share_q$MSE_p[2],
                   unshare_g$MSE_sa[2], 
                   share_g$MSE_p[2]),
    
    psi2 <- c(unshare_q$hatpsi_sa[3], 
               share_q$hatpsi[3], 
               unshare_g$hatpsi_sa[3], 
               share_g$hatpsi[3]),
    
    sd_psi2 <- c(unshare_q$sdpsi_sa[3], 
                  share_q$sdpsi[3], 
                  unshare_g$sdpsi_sa[3], 
                  share_g$sdpsi[3]),
    
    bias_psi2 <- c(unshare_q$bias_sa[3], 
                    share_q$bias[3],
                    unshare_g$bias_sa[3],
                    share_g$bias[3]),
    
    MSE_psi2 <- c(unshare_q$MSE_sa[3],
                   share_q$MSE_p[3],
                   unshare_g$MSE_sa[3], 
                   share_g$MSE_p[3]),
    
    M_bar <- c(unshare_q$M_sa_stage, 
               share_q$M_stage,
               unshare_g$M_sa_stage, 
               share_g$M_stage),
    
    M_tilde <- c(unshare_q$M_sa_overall,
                 share_q$M_overall,
                 unshare_g$M_sa_overall,
                 share_g$M_overall)
  )
  
  
  names(df_est) <- c("psi0", "sd_psi0", "bias_psi0", "MSE_psi0",
                     "psi1", "sd_psi1", "bias_psi1", "MSE_psi1",
                     "psi2", "sd_psi2", "bias_psi2", "MSE_psi2",
                     "M_bar", "M_tilde")
  row.names(df_est) <- c("unshared-Q SA", "shared-Q", "unshared-G SA", "shared-G")
  
  df_est <- round(df_est, 2)
  
  return(df_est)
}




# n = 200, scene = 1 --------------------------------------------------------
fn.summary(output_unsharedq_s1_200, output_sharedq_s1_200, output_unsharedg_s1_200, output_sharedg_s1_200)

# n = 500, scene = 1 --------------------------------------------------------
fn.summary(output_unsharedq_s1_500, output_sharedq_s1_500, output_unsharedg_s1_500, output_sharedg_s1_500)

# n = 2000, scene = 1 --------------------------------------------------------
fn.summary(output_unsharedq_s1_2000, output_sharedq_s1_2000, output_unsharedg_s1_2000, output_sharedg_s1_2000)

# n = 5000, scene = 1 --------------------------------------------------------
fn.summary(output_unsharedq_s1_5000, output_sharedq_s1_5000, output_unsharedg_s1_5000, output_sharedg_s1_5000)



# n = 200, scene = 2 --------------------------------------------------------
fn.summary(output_unsharedq_s2_200, output_sharedq_s2_200, output_unsharedg_s2_200, output_sharedg_s2_200)

# n = 500, scene = 2 --------------------------------------------------------
fn.summary(output_unsharedq_s2_500, output_sharedq_s2_500, output_unsharedg_s2_500, output_sharedg_s2_500)

# n = 2000, scene = 2 --------------------------------------------------------
fn.summary(output_unsharedq_s2_2000, output_sharedq_s2_2000, output_unsharedg_s2_2000, output_sharedg_s2_2000)

# n = 5000, scene = 2 --------------------------------------------------------
fn.summary(output_unsharedq_s2_5000, output_sharedq_s2_5000, output_unsharedg_s2_5000, output_sharedg_s2_5000)