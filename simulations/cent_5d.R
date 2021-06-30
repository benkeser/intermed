#! /usr/bin/env Rscript

# get environment variables
MYSCRATCH <- Sys.getenv('MYSCRATCH')
RESULTDIR <- Sys.getenv('RESULTDIR')
STEPSIZE <- as.numeric(Sys.getenv('STEPSIZE'))
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# set defaults if nothing comes from environment variables
MYSCRATCH[is.na(MYSCRATCH)] <- '.'
RESULTDIR[is.na(RESULTDIR)] <- '.'
STEPSIZE[is.na(STEPSIZE)] <- 1
TASKID[is.na(TASKID)] <- 0

# get command lines arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Not enough arguments. Please use args 'listsize', 'prepare', 'run <itemsize>' or 'merge'")
}

# simulation parameters
ns <- c(250, 500, 1000, 2000)
seed <- 1:1000
parm <- expand.grid(n = ns, seed = seed)

# save directory
save_dir <- "~/intermed/out/"
code_dir <- "~/intermed/"

# load packages
library(SuperLearner)
library(intermed)
# a function to simulate data 
# here's a 5-dimensional example, where only 2 variables matter
make_data <- function(n, 
                      M1_int = -1, M2_int = -1,
                      M1_C1 = -0.125, M2_C1 = 0.075,
                      M1_A = 0.25, M2_A = -0.25){
  C <- data.frame(C1 = runif(n), C2 = runif(n), C5 = rbinom(n, 1, 0.5),
                  C4 = rbinom(n, 1, 0.25), C3 = runif(n))
  g0 <- plogis(-1 + 0.25*C$C1 + 0.5*C$C2)
  A <- rbinom(n, 1, g0)
  # A increases values of M1 and M2
  success_prob_M1 <- plogis(M1_int + M1_C1 * C$C1 + M1_A * A)
  success_prob_M2 <- plogis(M2_int + M2_C1 * C$C1 + M2_A * A)

  # for checking joint distribution for extreme values
  # C1_val <- c(0, 1)
  # A_val <- c(0, 1)
  # cov_grid <- expand.grid(C1 = C1_val, A = A_val)
  # success_prob_M1_t <- plogis(-1 + 0.25 * cov_grid$C1 + 0.25 * cov_grid$A)
  # success_prob_M2_t <- plogis(-1 + 0.25 * cov_grid$C1 + 0.35 * cov_grid$A)

  # marg_M1 <- c(dgeom(0:4, success_prob_M1_t), 1 - sum(dgeom(0:4, success_prob_M1_t)))
  # marg_M2 <- c(dgeom(0:4, success_prob_M2_t), 1 - sum(dgeom(0:4, success_prob_M2_t)))
  # combn_grid <- expand.grid(marg_M1 = marg_M1, marg_M2 = marg_M2)
  # combn_grid$joint <- combn_grid$marg_M1 * combn_grid$marg_M2

  # probably should switch this to geometric distribution, so that we 
  # can be sure that we are nailing the mediator distribution estimation
  M1 <- rgeom(n, success_prob_M1)
  M2 <- rgeom(n, success_prob_M2)
  M1[M1 > 4] <- 5
  M2[M2 > 4] <- 5
  # higher values of M1 and M2 lead to higher probability that Y = 1
  Qbar0 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1 + 0.5 * M2 + 0.25 * A)

  Y <- rbinom(n, 1, Qbar0)
  return(list(C = C, A = A, M1 = M1, M2 = M2, Y = Y))
}

get_truth <- function(big_n = 1e6,
                      M1_int = -1, M2_int = -1,
                      M1_C1 = -0.125, M2_C1 = 0.075,
                      M1_A = 0.25, M2_A = -0.25){

  C <- data.frame(C1 = runif(big_n), C2 = runif(big_n))

  success_prob_M1_A0 <- plogis(M1_int + M1_C1 * C$C1 + M1_A * 0)
  success_prob_M2_A0 <- plogis(M2_int + M2_C1 * C$C1 + M2_A * 0)

  success_prob_M1_A1 <- plogis(M1_int + M1_C1 * C$C1 + M1_A * 1)
  success_prob_M2_A1 <- plogis(M2_int + M2_C1 * C$C1 + M2_A * 1)

  M1_A0 <- rgeom(big_n, success_prob_M1_A0)
  M2_A0 <- rgeom(big_n, success_prob_M2_A0)
  M1_A1 <- rgeom(big_n, success_prob_M1_A1)
  M2_A1 <- rgeom(big_n, success_prob_M2_A1)

  M1_A0[M1_A0 > 4] <- 5
  M2_A0[M2_A0 > 4] <- 5
  M1_A1[M1_A1 > 4] <- 5
  M2_A1[M2_A1 > 4] <- 5

  # total effect
  Qbar0_A1 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A1 + 0.5 * M2_A1 + 0.25 * 1)
  Qbar0_A0 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A0 + 0.5 * M2_A0 + 0.25 * 0)
  total_effect <- mean(Qbar0_A1 - Qbar0_A0)

  # direct effect 
  Qbar0_A1 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A0 + 0.5 * M2_A0 + 0.25 * 1)
  Qbar0_A0 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A0 + 0.5 * M2_A0 + 0.25 * 0)
  direct_effect <- mean(Qbar0_A1 - Qbar0_A0)

  # indirect effect through M1
  Qbar0_A1 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A1 + 0.5 * M2_A0 + 0.25 * 1)
  Qbar0_A0 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A0 + 0.5 * M2_A0 + 0.25 * 1)
  indirect_effect_M1 <- mean(Qbar0_A1 - Qbar0_A0)

  # indirect effect through M2
  Qbar0_A1 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A1 + 0.5 * M2_A1 + 0.25 * 1)
  Qbar0_A0 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A1 + 0.5 * M2_A0 + 0.25 * 1)
  indirect_effect_M2 <- mean(Qbar0_A1 - Qbar0_A0)

  return(list(total = total_effect,
        direct = direct_effect,
        indirect_M1 = indirect_effect_M1, 
        indirect_M2 = indirect_effect_M2))
}

# functions to compute to use truth as estimators
SL_Qbar_true <- function(...){
  fit <- list(NULL)
  class(fit) <- "SL_Qbar_true"
  return(list(fit = fit, pred = NULL))
}

predict.SL_Qbar_true <- function(object, newdata, ...){
  plogis(-1 + newdata$C1 - newdata$C2 + 0.5 * newdata$M1 + 0.5 * newdata$M2 + newdata$A)
}

SL_Q_M1_true <- function(...){
  fit <- list(NULL)
  class(fit) <- "SL_Q_M1_true"
  return(list(fit = fit, pred = NULL)) 
}
SL_Q_M2_true <- function(...){
  fit <- list(NULL)
  class(fit) <- "SL_Q_M2_true"
  return(list(fit = fit, pred = NULL)) 
}

predict.SL_Q_M1_true <- function(object, newdata, 
                                 M1_int = -1, 
                                 M1_C1 = 0.25, 
                                 M1_A = 0.25, ...){
  plogis(M1_int + M1_C1 * newdata$C1 + M1_A * newdata$A)
}

predict.SL_Q_M2_true <- function(object, newdata, 
                                 M2_int = -1, 
                                 M2_C1 = 0.25, 
                                 M2_A = 0.35, ...){
  plogis(M2_int + M2_C1 * newdata$C1 + M2_A * newdata$A)
}

SL_g_true <- function(newX, ...){
  fit <- list(NULL)
  class(fit) <- "SL_g_true"
  return(list(fit = fit, pred = 1 - plogis(-1 + newX$C1 + newX$C2))) 
}

predict.SL_g_true <- function(object, newdata, ...){
  1 - plogis(-1 + newdata$C1 + newdata$C2)
}

# get the list size #########
if (args[1] == 'listsize') {
  cat(nrow(parm))
}

# execute prepare job ##################
if (args[1] == 'prepare') {
  print(paste0('initial datasets saved to: ~/drinf/scratch/dataList ... .RData'))
}

# execute parallel job #################################################
if (args[1] == 'run') {
  if (length(args) < 2) {
    stop("Not enough arguments. 'run' needs a second argument 'id'")
  }
  id <- as.numeric(args[2])
  print(paste(Sys.time(), "arrid:" , id, "TASKID:",
              TASKID, "STEPSIZE:", STEPSIZE))
  for (i in (id+TASKID):(id+TASKID+STEPSIZE-1)) {
    print(paste(Sys.time(), "i:" , i))
    print(parm[i,])
    
    # set seed
    set.seed(parm$seed[i])

    # make data
    data <- make_data(n = parm$n[i])

    rslt <- intermed(Y = data$Y, C = data$C, M1 = data$M1, M2 = data$M2, A = data$A, 
                     a = 1, 
                     a_star = 0,
                     n_SL = 1,
                     SL_Qbar = c("SL.glm", "SL.earth", "SL.ranger"),
                     SL_g = c("SL.glm", "SL.earth", "SL.ranger"),
                     SL_Q_M = list(M1 = c("SL.glm", "SL.earth", "SL.ranger"), 
                                   M2 = c("SL.glm", "SL.earth", "SL.ranger")),
                     tolg = 1e-3, 
                     targeted_se = FALSE, 
                     return_models = FALSE,
                     verbose = FALSE,
                     stratify = FALSE,
                     max_iter = 10)

    # get truth
    set.seed(1234)
    truth <- c(unlist(get_truth()), 0)

    # get confidence intervals
    all_ci <- ci(rslt, est = c("tmle", "aiptw"))
    
    # check for truth in CIs
    in_ci <- function(truth, ci){
      truth > min(ci) & truth < max(ci)
    }
    tmle_cover <- aiptw_cover <- rep(NA, 5)
    for(j in seq_len(5)){
      tmle_cover[j] <- in_ci(truth[j], all_ci$tmle[j, c(2,3)])
      aiptw_cover[j] <- in_ci(truth[j], all_ci$aiptw[j, c(2,3)])
    }
    out <- c(parm$n[i], parm$seed[i], t(all_ci$tmle), t(all_ci$aiptw), tmle_cover, aiptw_cover,
             truth)
    names(out) <- c("n", "seed", 
                    "total_tmle_est", "total_tmle_cil", "total_tmle_ciu",
                    "direct_tmle_est","direct_tmle_cil","direct_tmle_ciu",
                    "indirectM1_tmle_est","indirectM1_tmle_cil","indirectM1_tmle_ciu",
                    "indirectM2_tmle_est","indirectM2_tmle_cil","indirectM2_tmle_ciu",
                    "covarM1M2_tmle_est","covarM1M2_tmle_cil","covarM1M2_tmle_ciu",
                    "total_aiptw_est", "total_aiptw_cil", "total_aiptw_ciu",
                    "direct_aiptw_est","direct_aiptw_cil","direct_aiptw_ciu",
                    "indirectM1_aiptw_est","indirectM1_aiptw_cil","indirectM1_aiptw_ciu",
                    "indirectM2_aiptw_est","indirectM2_aiptw_cil","indirectM2_aiptw_ciu",
                    "covarM1M2_aiptw_est","covarM1M2_aiptw_cil","covarM1M2_aiptw_ciu",
                    "total_tmle_cover", "direct_tmle_cover", "indirectM1_tmle_cover", "indirectM2_tmle_cover", "covarM1M2_tmle_cover",
                    "total_aiptw_cover", "direct_aiptw_cover", "indirectM1_aiptw_cover", "indirectM2_aiptw_cover", "covarM1M2_aiptw_cover",
                    "true_total", "true_direct", "true_indirectM1", "true_indirectM2", "true_covarM1M2"
                    )
    save(out, file = paste0("~/intermed/out/fitsl_", i, ".RData"))
  }
}

# merge job ###########################
if (args[1] == 'merge') {   
  # parameters
  ns <- c(250, 500, 1000, 2000)
  seed <- 1:1000
  parm <- expand.grid(n = ns, seed = seed)

  save_dir <- "~/intermed/out/"
  n_out <- 47
  
  rslt <- matrix(NA, nrow = nrow(parm), ncol = n_out)
  for(i in 1:nrow(parm)){
    tmp <- tryCatch({
      load(paste0(save_dir, "fituflmsl3_", i, ".RData"))
      out
    }, error = function(e){
      rep(NA, n_out)
    })
    rslt[i, ] <- tmp
  }

  out <- data.frame(rslt)
  colnames(out) <- c("n", "seed", 
                    "total_tmle_est", "total_tmle_cil", "total_tmle_ciu",
                    "direct_tmle_est","direct_tmle_cil","direct_tmle_ciu",
                    "indirectM1_tmle_est","indirectM1_tmle_cil","indirectM1_tmle_ciu",
                    "indirectM2_tmle_est","indirectM2_tmle_cil","indirectM2_tmle_ciu",
                    "covarM1M2_tmle_est","covarM1M2_tmle_cil","covarM1M2_tmle_ciu",
                    "total_aiptw_est", "total_aiptw_cil", "total_aiptw_ciu",
                    "direct_aiptw_est","direct_aiptw_cil","direct_aiptw_ciu",
                    "indirectM1_aiptw_est","indirectM1_aiptw_cil","indirectM1_aiptw_ciu",
                    "indirectM2_aiptw_est","indirectM2_aiptw_cil","indirectM2_aiptw_ciu",
                    "covarM1M2_aiptw_est","covarM1M2_aiptw_cil","covarM1M2_aiptw_ciu",
                    "total_tmle_cover", "direct_tmle_cover", "indirectM1_tmle_cover", "indirectM2_tmle_cover", "covarM1M2_tmle_cover",
                    "total_aiptw_cover", "direct_aiptw_cover", "indirectM1_aiptw_cover", "indirectM2_aiptw_cover", "covarM1M2_aiptw_cover",
                    "truth_total", "truth_direct", "truth_indirectM1", "truth_indirectM2", "truth_covarM1M2"
                    )

  save(out, file = "~/intermed/all_out_intermed_uflmsl3.RData")

  get_bias <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- mean(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      aiptw <- mean(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      return(c(tmle, aiptw))
    })
  }

  get_rootnbias <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
      aiptw <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
      return(c(tmle, aiptw))
    })
  }

  get_mse <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- mean((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      aiptw <- mean((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      return(c(tmle, aiptw))
    })
  }
  get_nmse <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      aiptw <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      return(c(tmle, aiptw))
    })
  }
  get_sd <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- sd(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))], na.rm = TRUE)
      aiptw <- sd(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] , na.rm = TRUE)
      return(c(tmle, aiptw))
    })
  }
  get_rootnsd <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))], na.rm = TRUE)
      aiptw <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] , na.rm = TRUE)
      return(c(tmle, aiptw))
    })
  }
  get_se_est <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- mean(x[ , grep(paste0(which_eff, "_tmle_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_tmle_cil"), colnames(x))]) / (2 * 1.96)
      aiptw <- mean(x[ , grep(paste0(which_eff, "_aiptw_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_aiptw_cil"), colnames(x))]) / (2 * 1.96)
      return(c(tmle, aiptw))
    })
  }
  get_cover <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- mean(x[ , grep(paste0(which_eff, "_tmle_cover"), colnames(x))], na.rm = TRUE)
      aiptw <- mean(x[ , grep(paste0(which_eff, "_aiptw_cover"), colnames(x))], na.rm = TRUE)
      return(c(tmle, aiptw))
    })
  }
  
  get_density_orc <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- density((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / sd(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))]), na.rm = TRUE)
      aiptw <- density((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / sd(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))]), na.rm = TRUE)
      return(list(tmle, aiptw))
    })
  }
  get_density_est <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle_se <- (x[ , grep(paste0(which_eff, "_tmle_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_tmle_cil"), colnames(x))]) / (2 * 1.96)
      aiptw_se <- (x[ , grep(paste0(which_eff, "_aiptw_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_aiptw_cil"), colnames(x))]) / (2 * 1.96)
      tmle <- density((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / tmle_se, na.rm = TRUE)
      aiptw <- density((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / aiptw_se, na.rm = TRUE)
      return(list(tmle, aiptw))
    })
  }

  # check oracle CI coverage
  check_oracle_ci <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
    oracle_se_tmle <- sd(x[, paste0(which_eff, "_tmle_est")], na.rm = TRUE)
    oracle_se_aiptw <- sd(x[, paste0(which_eff, "_aiptw_est")], na.rm = TRUE)
    tmle_cil <- x[, paste0(which_eff, "_tmle_est")] - 1.96 * oracle_se_tmle
    tmle_ciu <- x[, paste0(which_eff, "_tmle_est")] + 1.96 * oracle_se_tmle
    aiptw_cil <- x[, paste0(which_eff, "_aiptw_est")] - 1.96 * oracle_se_aiptw
    aiptw_ciu <- x[, paste0(which_eff, "_aiptw_est")] + 1.96 * oracle_se_aiptw
    tmle_cover <- tmle_cil < x[, paste0("truth_", which_eff)] & tmle_ciu > x[, paste0("truth_", which_eff)]
    aiptw_cover <- aiptw_cil < x[, paste0("truth_", which_eff)] & aiptw_ciu > x[, paste0("truth_", which_eff)]
    return(c(mean(tmle_cover, na.rm = TRUE), mean(aiptw_cover, na.rm = TRUE)))
    })
  }
  format_out <- function(out, summary_fn,
                         all_eff = c("total", "direct", "indirectM1", 
                                     "indirectM2", "covarM1M2")){
    if(!grepl("get_density", summary_fn)){    
      all_out <- sapply(all_eff, function(x){
        suppressWarnings(a <- data.frame(Reduce(rbind, 
                          do.call(summary_fn, 
                                  args = list(out = out, 
                                              which_eff = x))),
                   n = c(250, 500, 1000, 2000),
                   stringsAsFactors = FALSE))
        colnames(a) <- c("tmle", "os", "n"); row.names(a) <- NULL; a
      }, simplify = FALSE)
    }else{
      all_out <- sapply(all_eff, function(x){
        do.call(summary_fn, args = list(out = out, which_eff = x))
      })
    }
  }
  # # fix coverage for indirect M2
  # out$indirectM2_tmle_cover <- out$indirectM2_tmle_cil < out$truth_indirectM2 & out$indirectM2_tmle_ciu > out$truth_indirectM2
  # out$indirectM2_aiptw_cover <- out$indirectM2_aiptw_cil < out$truth_indirectM2 & out$indirectM2_aiptw_ciu > out$truth_indirectM2

  load("~/Dropbox/Emory/Flu/inter_med/sim_rslt/all_out_intermed_sl.RData")
  
  all_bias <- format_out(out, "get_bias")
  all_rootnbias <- format_out(out, "get_rootnbias")
  all_mse <- format_out(out, "get_mse")
  all_nmse <- format_out(out, "get_nmse")
  # !!! could get true asymptotic variance by simulating large data set
  # !!! using true versions of nuisance estimators and evaluating variance of
  # !!! EIF. 
  all_sd <- format_out(out, "get_sd")
  all_rootnsd <- format_out(out, "get_rootnsd")
  all_cover <- format_out(out, "get_cover")
  all_cover_orc <- format_out(out, "check_oracle_ci")
  all_se <- format_out(out, "get_se_est")

  all_density_orc <- format_out(out, "get_density_orc")
  all_density_est <- format_out(out, "get_density_est")


  blank_plot <- function(...){
    plot(1e-10, 1e-10, pch = "", bty = "n", xaxt = "n", ...)
  }

  plot_one_est_row <- function(which_eff = "total",
                all_bias, all_mse, all_sd, 
                all_cover, all_cover_orc, all_se,
                bias_ylim = NULL){
    # six panel plot, two rows
    # top row = bias, sd, MSE
    # bottom row = samp dist. TMLE, samp. dist AIPTW, coverage ()
    # bias plot
    bias_range <- range(c(all_bias[[which_eff]]$tmle, all_bias[[which_eff]]$os))
    max_bias_val <- max(abs(bias_range))

    if(!is.null(bias_ylim)){
      yl <- c(-1.05 * max_bias_val, 1.05 * max_bias_val)
    }
    blank_plot(xlim = c(0.5, 4.5), ylim = bias_ylim,
               xlab = "n", ylab = "Bias")
    axis(side = 1, at = 1:4, labels = all_bias[[which_eff]]$n)
    abline(h = 0, lty = 3)
    points(y = all_bias[[which_eff]]$tmle, x = 1:4, type = "b")
    points(y = all_bias[[which_eff]]$os, x = 1:4, type = "b", pch = 2)

    # sd plot
    sd_range <- range(c(all_sd[[which_eff]]$tmle, all_sd[[which_eff]]$os))
    max_sd_val <- max(abs(sd_range))
    blank_plot(xlim = c(0.5, 4.5), ylim = c(0, max_sd_val * 1.1),
               xlab = "n", ylab = "Standard deviation")
    axis(side = 1, at = 1:4, labels = all_sd[[which_eff]]$n)
    # abline(h = 0, lty = 3)
    points(y = all_sd[[which_eff]]$tmle, x = 1:4, type = "b")
    points(y = all_sd[[which_eff]]$os, x = 1:4, type = "b", pch = 2)

    # mean squared error plot
    mse_range <- range(c(all_mse[[which_eff]]$tmle, all_mse[[which_eff]]$os))
    max_mse_val <- max(abs(mse_range))
    min_mse_val <- min(abs(mse_range))
    blank_plot(xlim = c(0.5, 4.5), ylim = c(min_mse_val * 0.9, max_mse_val * 1.1), log = "y", 
               xlab = "n", ylab = "Mean squared error")
    axis(side = 1, at = 1:4, labels = all_mse[[which_eff]]$n)
    abline(h = 0, lty = 3)
    points(y = all_mse[[which_eff]]$tmle, x = 1:4, type = "b")
    points(y = all_mse[[which_eff]]$os, x = 1:4, type = "b", pch = 2)
  }

  # plot results for point estimates
  pdf("~/Dropbox/Emory/Flu/inter_med/estimation.pdf",
    height = 7, width = 7)
  layout(matrix(c(1, 1, 1, 2:13), nrow =5, byrow = TRUE),
         heights = c(0.25, 1, 1, 1, 1))
  par(oma = c(0, 2.1, 0, 0), mar = c(2.9, 2.9, 0.1, 0.1), mgp = c(1.8, 0.5, 0))
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend(x = 0.37, y = 0.5, ncol = 2, pch = c(2, 1), title = "Estimator",
       legend = c("One-step", "TMLE"), xpd = NA)
  at_loc <- seq(0, 1, length = 6)[-c(1,6)] + c(0.052, 0.0125, -0.022, -0.055)[4:1]
  plot_one_est_row(which_eff = "direct", all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                   bias_ylim = c(-0.02, 0.02))
  mtext(side = 2, outer = TRUE, line = 0, expression(psi[A]), at = at_loc[4])
  plot_one_est_row(which_eff = "indirectM1", all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                   bias_ylim = c(-0.02, 0.02))
  mtext(side = 2, outer = TRUE, line = 0, expression(psi[M[1]]), at = at_loc[3])
  plot_one_est_row(which_eff = "indirectM2", all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                   bias_ylim = c(-0.02, 0.02))
  mtext(side = 2, outer = TRUE, line = 0, expression(psi[M[2]]), at = at_loc[2])
  plot_one_est_row(which_eff = "covarM1M2", all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                   bias_ylim = c(-0.02, 0.02))
  mtext(side = 2, outer = TRUE, line = 0, expression(psi[M[1]*","*M[2]]), at = at_loc[1])
  dev.off()
  # plot results for sampling distributions
  plot_one_samp_dist <- function(which_eff = "total",
                                 which_est = c("tmle", "aiptw"),
                all_density_orc, all_density_est,
                all_cover, all_cover_orc){

    if(which_eff == "direct"){
      xlab1 = expression("("*psi["n,A"]*" - "*psi[A]*")/se("*psi["n,A"]*")")
      xlab2 = expression(n^{1/2}*"("*psi["n,A"]*" - "*psi[A]*")/"*sigma["n,A"])
    }
    
    if(which_eff == "indirectM1"){
      xlab1 = expression("("*psi["n,"*M[1]]*" - "*psi[M[1]]*")/se("*psi["n,"*M[1]]*")")
      xlab2 = expression(n^{1/2}*"("*psi["n,"*M[1]]*" - "*psi[M[1]]*")/"*sigma["n,"*M[1]])
    }

    if(which_eff == "indirectM2"){
      xlab1 = expression("("*psi["n,"*M[2]]*" - "*psi[M[2]]*")/se("*psi["n,"*M[2]]*")")
      xlab2 = expression(n^{1/2}*"("*psi["n,"*M[2]]*" - "*psi[M[2]]*")/"*sigma["n,"*M[2]])
    }

    if(which_eff == "covarM1M2"){
      xlab1 = expression("("*psi["n,"*M[1]*","*M[2]]*" - "*psi[M[1]*","*M[2]]*")/se("*psi["n,"*M[1]*","*M[2]]*")")
      xlab2 = expression(n^{1/2}*"("*psi["n,"*M[1]*","*M[2]]*" - "*psi[M[1]*","*M[2]]*")/"*sigma["n,"*M[1]*","*M[2]])
    }
    
    # layout(matrix(1:(3 * length(which_est)), byrow = TRUE, nrow = 2))

    if("tmle" %in% which_est){
      # plot sampling distribution of TMLE by sample size
      blank_plot(xlim = c(-3.5, 3.5), ylim = c(0, 0.55),
                 xlab = xlab1, ylab = "Density")
      axis(side = 1, at = seq(-3, 3, by = 1))
      x_seq <- seq(-4, 4, length = 2000)
      lines(x = x_seq, y = dnorm(x_seq), lty = 1, lwd = 2)
      grays <- paste0("gray", c(80, 60, 40, 20))
      for(i in 1:4){
        lines(all_density_orc[,which_eff][[i]][[1]], col = grays[i], lty = 2)
      }

      # plot sampling distribution of TMLE by sample size
      blank_plot(xlim = c(-3.5, 3.5), ylim = c(0, 0.55),
                 xlab = xlab2, ylab = "Density")
      axis(side = 1, at = seq(-3, 3, by = 1))
      x_seq <- seq(-4, 4, length = 2000)
      lines(x = x_seq, y = dnorm(x_seq), lty = 1, lwd = 2)
      grays <- paste0("gray", c(80, 60, 40, 20))
      for(i in 1:4){
        lines(all_density_est[,which_eff][[i]][[1]], col = grays[i], lty = 2)
      }

      # plot oracle confidence interval coverage
      blank_plot(xlim = c(0.5, 4.5), ylim = c(0.5, 1),
                 xlab = "n", ylab = "Coverage probability")
      axis(side = 1, at = 1:4, labels = all_mse[[which_eff]]$n)
      abline(h = 0.95, lty = 3)
      points(y = all_cover_orc[[which_eff]]$tmle, x = 1:4, type = "b", pch = 1)
      points(y = all_cover[[which_eff]]$tmle, x = 1:4, type = "b", pch = 19)
    }

    if("aiptw" %in% which_est){
      # plot sampling distribution of aiptw by sample size
      blank_plot(xlim = c(-3.5, 3.5), ylim = c(0, 0.55),
                 xlab = xlab1, ylab = "Density")
      axis(side = 1, at = seq(-3, 3, by = 1))
      x_seq <- seq(-4, 4, length = 2000)
      lines(x = x_seq, y = dnorm(x_seq), lty = 1, lwd = 2)
      grays <- paste0("gray", c(80, 60, 40, 20))
      for(i in 1:4){
        lines(all_density_orc[,which_eff][[i]][[2]], col = grays[i], lty = 2)
      }

      # plot sampling distribution of TMLE by sample size
      blank_plot(xlim = c(-3.5, 3.5), ylim = c(0, 0.55),
                 xlab = xlab2, ylab = "Density")
      axis(side = 1, at = seq(-3, 3, by = 1))
      x_seq <- seq(-4, 4, length = 2000)
      lines(x = x_seq, y = dnorm(x_seq), lty = 1, lwd = 2)
      grays <- paste0("gray", c(80, 60, 40, 20))
      for(i in 1:4){
        lines(all_density_est[,which_eff][[i]][[2]], col = grays[i], lty = 2)
      }

      # plot oracle confidence interval coverage
      blank_plot(xlim = c(0.5, 4.5), ylim = c(0.5, 1),
                 xlab = "n", ylab = "Coverage probability")
      axis(side = 1, at = 1:4, labels = all_mse[[which_eff]]$n)
      abline(h = 0.95, lty = 3)
      points(y = all_cover_orc[[which_eff]]$os, x = 1:4, type = "b", pch = 2)
      points(y = all_cover[[which_eff]]$os, x = 1:4, type = "b", pch = 17)
    }

  }



pdf("~/Dropbox/Emory/Flu/inter_med/inference_tmle.pdf",
    height = 7, width = 7)
layout(matrix(c(1, 1, 2, 3:14), nrow = 5, byrow = TRUE),
       heights = c(0.25, rep(1, 4)))

par(mar = c(1.5, 0, 1.5, 0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", ncol = 5, lty = c(rep(2,4), 1), col = c(paste0("gray", c(80, 60, 40, 20)), 1),
       legend = c(250, 500, 1000, 2000, expression(infinity)), title = "n", xpd = TRUE)
par(mar = c(1.5, 0, 1.5, 0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", ncol = 2, pch = c(19, 1), title = "Standard error",
       legend = c("Estimated", "Oracle"), xpd = TRUE)

for(eff in c("direct", "indirectM1", "indirectM2", "covarM1M2")){
  par(mar = c(2.9, 2.9, 0.05, 0.05), mgp = c(1.8, 0.5, 0))
  plot_one_samp_dist(which_eff = eff, 
                     which_est = "tmle",
                     all_density_orc = all_density_orc, 
                     all_density_est = all_density_est, 
                     all_cover = all_cover, 
                     all_cover_orc = all_cover_orc)
}
dev.off()

pdf("~/Dropbox/Emory/Flu/inter_med/inference_aiptw.pdf",
    height = 7, width = 7)
layout(matrix(c(1, 1, 2, 3:14), nrow = 5, byrow = TRUE),
       heights = c(0.25, rep(1, 4)))

par(mar = c(1.5, 0, 1.5, 0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", ncol = 5, lty = c(rep(2,4), 1), col = c(paste0("gray", c(80, 60, 40, 20)), 1),
       legend = c(250, 500, 1000, 2000, expression(infinity)), title = "n", xpd = TRUE)
par(mar = c(1.5, 0, 1.5, 0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", ncol = 2, pch = c(17, 2), title = "Standard error",
       legend = c("Estimated", "Oracle"), xpd = TRUE)

for(eff in c("direct", "indirectM1", "indirectM2", "covarM1M2")){
  par(mar = c(2.9, 2.9, 0.05, 0.05), mgp = c(1.8, 0.5, 0))
  plot_one_samp_dist(which_eff = eff, 
                     which_est = "aiptw",
                     all_density_orc = all_density_orc, 
                     all_density_est = all_density_est, 
                     all_cover = all_cover, 
                     all_cover_orc = all_cover_orc)
}
dev.off()




par(mar = c(3.2, 3.2, 0.25, 0.25), mgp = c(2.2, 0.5, 0))
layout(matrix(1:12, nrow = 4, byrow = TRUE))
for(eff in c("direct", "indirectM1", "indirectM2", "covarM1M2")){
  plot_one_samp_dist(which_eff = eff, 
                     which_est = "aiptw",
                     all_density_orc = all_density_orc, 
                     all_density_est = all_density_est, 
                     all_cover = all_cover, 
                     all_cover_orc = all_cover_orc)
}





## need to figure out what's going on with covariant TMLE
