get_gcomp <- function(input_data, version = 0){
  # fit linear outcome regression
  Cs = paste0('C', 1:5, collapse  = ' + ')
  if(version == 0){
    or_fit <- glm(paste0("Y ~ A*M1 + A*M2 +", Cs), 
                  family = gaussian(), data = input_data)
    # fit linear mediator regression
    M1_fit <- glm(paste0("M1 ~ A + ", Cs), family = gaussian(), data = input_data)
    M2_fit <- glm(paste0("M2 ~ A + ", Cs), family = gaussian(), data = input_data)
    
    # obtain coefficients
    thetas <- or_fit$coefficients
    beta1s <- M1_fit$coefficients
    beta2s <- M2_fit$coefficients
    # interventional direct effect
    direct_gcomp_est <- thetas["A"] + thetas["A:M1"]*(beta1s[1] + beta1s["A"]*0 + sum(beta1s[3:7]*apply(input_data[,1:5], 2, mean))) + 
      thetas["A:M2"]*(beta2s[1] + beta2s["A"]*0 + sum(beta2s[3:7]*apply(input_data[,1:5], 2, mean)))
    # interventional indirect effect via M1
    indirectM1_gcomp_est <- (thetas["M1"] +thetas["A:M1"])*beta1s["A"]
    # interventional indirect effect via M2
    indirectM2_gcomp_est <- (thetas["M2"] + thetas["A:M2"])*beta2s["A"]
    # covarM1M2 
    covarM1M2_gcomp_est <- 0
    # total effect
    total_gcomp_est <- direct_gcomp_est + indirectM1_gcomp_est + indirectM2_gcomp_est 
  }else if(version == 1){
    or_fit <- glm(paste0("Y ~ A*M1 + A*M2 + M1*M2 +", Cs), 
          family = gaussian(), data = input_data)
    # fit linear mediator regression
    M1_fit <- glm(paste0("M1 ~ A + ", Cs), family = gaussian(), data = input_data)
    M2_fit <- glm(paste0("M2 ~ A + ", Cs), family = gaussian(), data = input_data)
    
    # obtain coefficients
    thetas <- or_fit$coefficients
    beta1s <- M1_fit$coefficients
    beta2s <- M2_fit$coefficients
    # interventional direct effect
    direct_gcomp_est <- thetas["A"] + thetas["A:M1"]*(beta1s[1] + beta1s["A"]*0 + sum(beta1s[3:7]*apply(input_data[,1:5], 2, mean))) + 
      thetas["A:M2"]*(beta2s[1] + beta2s["A"]*0 + sum(beta2s[3:7]*apply(input_data[,1:5], 2, mean)))
    # interventional indirect effect via M1
    indirectM1_gcomp_est <- (thetas["M1"] + thetas["M1:M2"]*(beta2s[1] + beta2s["A"]*0 + sum(beta2s[3:7]*apply(input_data[,1:5], 2, mean))) + 
                               thetas["A:M1"])*beta1s["A"]
    # interventional indirect effect via M2
    indirectM2_gcomp_est <- (thetas["M2"] + thetas["M1:M2"]*(beta1s[1] + beta1s["A"]*1 + sum(beta1s[3:7]*apply(input_data[,1:5], 2, mean))) + 
                               thetas["A:M2"])*beta2s["A"]
    # covarM1M2 
    covarM1M2_gcomp_est <- 0
    # total effect
    total_gcomp_est <- direct_gcomp_est + indirectM1_gcomp_est + indirectM2_gcomp_est 
  }else if(version == 2){
    or_fit <- lm(paste0("Y ~ A*M1 + A*M2 + M1*M2 +", Cs), data = input_data)
    # fit linear mediator regression
    M1_fit <- lm(paste0("M1 ~ A + ", Cs), data = input_data)
    M2_fit <- lm(paste0("M2 ~ A + ", Cs, "+ M1 + A*M1"), data = input_data)
    # obtain coefficients
    thetas <- or_fit$coefficients
    beta1s <- M1_fit$coefficients
    beta2s <- M2_fit$coefficients
    # interventional direct effect
    direct_gcomp_est <- thetas["A"] + thetas["A:M1"]*(beta1s[1] + beta1s["A"]*0 + sum(beta1s[3:7]*apply(input_data[,1:5], 2, mean))) + 
      thetas["A:M2"]*(beta2s[1] + beta2s["A"]*0 + sum(beta2s[3:7]*apply(input_data[,1:5], 2, mean)))
    # interventional indirect effect via M1
    indirectM1_gcomp_est <- (thetas["M1"] + thetas["M1:M2"]*(beta2s[1] + beta2s["A"]*0 + sum(beta2s[3:7]*apply(input_data[,1:5], 2, mean))) + 
                               thetas["A:M1"])*beta1s["A"]
    # interventional indirect effect via M2
    indirectM2_gcomp_est <- (thetas["M2"] + thetas["M1:M2"]*(beta1s[1] + beta1s["A"]*1 + sum(beta1s[3:7]*apply(input_data[,1:5], 2, mean))) + 
                               thetas["A:M2"])*beta2s["A"]
    # covarM1M2 
    covarM1M2_gcomp_est <- summary(M1_fit)$sigma^2 * thetas["M1:M2"] * beta2s["A:M1"]
    # total effect
    total_gcomp_est <- direct_gcomp_est + indirectM1_gcomp_est + indirectM2_gcomp_est + covarM1M2_gcomp_est
    
  }
  return(list(total_gcomp_est = total_gcomp_est,
              direct_gcomp_est = direct_gcomp_est,
              indirectM1_gcomp_est = indirectM1_gcomp_est, 
              indirectM2_gcomp_est = indirectM2_gcomp_est,
              covarM1M2_gcomp_est = covarM1M2_gcomp_est))
}

get_gcomp_boot <- function(full_data, version = 0){
  # g-comp confidence interval(bootstrap)
  M <- 1000
  # empty vectors to hold results
  boot_gcomp_direct <- boot_gcomp_indirectM1 <-boot_gcomp_indirectM2 <- boot_gcomp_total <- boot_gcomp_covarM1M2 <- rep(NA, M)
  # set a seed to ensure reproducibility
  set.seed(parm$seed[i])
  for(m in 1:M){
    # sample id's with replacement
    boot_ids <- sample(1:parm$n[i], replace = TRUE)
    boot_data <- full_data[boot_ids,]
    # get estimates based on bootstrap data
    boot_gcomp_est <- get_gcomp(boot_data, version = version)
    boot_gcomp_direct[m] <- boot_gcomp_est$direct_gcomp_est
    boot_gcomp_indirectM1[m] <- boot_gcomp_est$indirectM1_gcomp_est
    boot_gcomp_indirectM2[m] <- boot_gcomp_est$indirectM2_gcomp_est
    boot_gcomp_covarM1M2[m] <- boot_gcomp_est$covarM1M2_gcomp_est
    boot_gcomp_total[m] <- boot_gcomp_est$total_gcomp_est
  }
  
  gcomp_cil <- c(apply(cbind(boot_gcomp_total, boot_gcomp_direct, boot_gcomp_indirectM1, boot_gcomp_indirectM2, boot_gcomp_covarM1M2), 2, quantile, p = 0.025))
  gcomp_ciu <- c(apply(cbind(boot_gcomp_total, boot_gcomp_direct, boot_gcomp_indirectM1, boot_gcomp_indirectM2, boot_gcomp_covarM1M2), 2, quantile, p = 0.975))
  return(list(cil = gcomp_cil, ciu = gcomp_ciu))
}

