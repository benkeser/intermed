#! /usr/bin/env Rscript

# install.packages("~/scratch/intermed", type = "source", repos = NULL)

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

make_data <- function(n = 1e2, 
                      A_success = "plogis(-1 + 0.125 * C$C1 + 0.25*C$C2)",
                      M1_success = "plogis(-1 + 0.45 * C$C1 + 0.125 * A)",
                      M2_success = "plogis(-1 + 0.15 * C$C1 + 0.2 * C$C2 - 0.2 * A)",
                      M1M2_threshold = 3,
                      Y_success = "-1 + C$C1 - C$C2 + 0.25 * M1 + 0.25 * M2 + 0.25 * A",
                      get_truth = FALSE){
  C <- data.frame(C1 = runif(n), C2 = runif(n), C5 = rbinom(n, 1, 0.5),
                  C4 = rbinom(n, 1, 0.25), C3 = runif(n))
  if(!get_truth){
    g0 <- eval(parse(text = A_success))
    A <- rbinom(n, 1, g0)
    success_prob_M1 <- eval(parse(text = M1_success))
    success_prob_M2 <- eval(parse(text = M2_success))
    M1 <- rgeom(n, success_prob_M1)
    M2 <- rgeom(n, success_prob_M2)
    M1[M1 > (M1M2_threshold - 1)] <- M1M2_threshold
    M2[M2 > (M1M2_threshold - 1)] <- M1M2_threshold
    Qbar0 <- eval(parse(text = Y_success))
    Y <- Qbar0 + rnorm(n)
    return(list(C = C, A = A, M1 = M1, M2 = M2, Y = Y))
  }else{
    success_prob_M1_A1 <- eval(parse(text = gsub("A", "1", M1_success)))
    success_prob_M1_A0 <- eval(parse(text = gsub("A", "0", M1_success)))
    success_prob_M2_A1 <- eval(parse(text = gsub("A", "1", M2_success)))
    success_prob_M2_A0 <- eval(parse(text = gsub("A", "0", M2_success)))
    M1_A0 <- rgeom(n, success_prob_M1_A0)
    M2_A0 <- rgeom(n, success_prob_M2_A0)
    M1_A1 <- rgeom(n, success_prob_M1_A1)
    M2_A1 <- rgeom(n, success_prob_M2_A1)
    M1_A0[M1_A0 > (M1M2_threshold - 1)] <- M1M2_threshold
    M2_A0[M2_A0 > (M1M2_threshold - 1)] <- M1M2_threshold
    M1_A1[M1_A1 > (M1M2_threshold - 1)] <- M1M2_threshold
    M2_A1[M2_A1 > (M1M2_threshold - 1)] <- M1M2_threshold

    # total effect
    Qbar0_A1 <- eval(parse(text = gsub("M2", "M2_A1", gsub("M1", "M1_A1", gsub("A", "1", Y_success)))))
    Qbar0_A0 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A0", gsub("A", "0", Y_success)))))
    total_effect <- mean(Qbar0_A1 - Qbar0_A0)
    # direct effect
    Qbar0_A1 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A0", gsub("A", "1", Y_success)))))
    Qbar0_A0 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A0", gsub("A", "0", Y_success)))))
    direct_effect <- mean(Qbar0_A1 - Qbar0_A0)
    
    # indirect effect through M1
    Qbar0_A1 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A1", gsub("A", "1", Y_success)))))
    Qbar0_A0 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A0", gsub("A", "1", Y_success)))))
    indirect_effect_M1 <- mean(Qbar0_A1 - Qbar0_A0)

    # indirect effect through M2
    Qbar0_A1 <- eval(parse(text = gsub("M2", "M2_A1", gsub("M1", "M1_A1", gsub("A", "1", Y_success)))))
    Qbar0_A0 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A1", gsub("A", "1", Y_success)))))
    indirect_effect_M2 <- mean(Qbar0_A1 - Qbar0_A0)

    return(list(total = total_effect,
          direct = direct_effect,
          indirect_M1 = indirect_effect_M1, 
          indirect_M2 = indirect_effect_M2))
  }
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
                     targeted_se = TRUE, 
                     return_models = FALSE,
                     verbose = FALSE,
                     stratify = FALSE,
                     max_iter = 10)

    rslt_parametric <- intermed(Y = data$Y, C = data$C, 
                                M1 = data$M1, M2 = data$M2, A = data$A, 
                               a = 1, 
                               a_star = 0,
                               n_SL = 1,
                               SL_Qbar = "SL.glm",
                               SL_g = "SL.glm",
                               SL_Q_M = list(M1 = "SL.glm", 
                                             M2 = "SL.glm"),
                               tolg = 1e-3, 
                               targeted_se = TRUE, 
                               return_models = FALSE,
                               verbose = FALSE,
                               stratify = FALSE,
                               max_iter = 10)

    # Jialu add code for parametric G-comp estimator
    # ...

    # get truth
    set.seed(1234)
    truth <- c(unlist(make_data(n = 1e6, get_truth = TRUE)), 0)

    # get confidence intervals
    all_ci <- ci(rslt, est = c("tmle", "aiptw"))
    all_ci_parametric <- ci(rslt_parametric, est = c("tmle", "aiptw"))

    # check for truth in CIs
    in_ci <- function(truth, ci){
      truth > min(ci) & truth < max(ci)
    }
    tmle_cover <- aiptw_cover <- rep(NA, 5)
    for(j in seq_len(5)){
      tmle_cover[j] <- in_ci(truth[j], all_ci$tmle[j, c(2,3)])
      aiptw_cover[j] <- in_ci(truth[j], all_ci$aiptw[j, c(2,3)])
    }
    tmle_cover_parametric <- aiptw_cover_parametric <- rep(NA, 5)
    for(j in seq_len(5)){
      tmle_cover_parametric[j] <- in_ci(truth[j], all_ci_parametric$tmle[j, c(2,3)])
      aiptw_cover_parametric[j] <- in_ci(truth[j], all_ci_parametric$aiptw[j, c(2,3)])
    }
    

    out <- c(parm$n[i], parm$seed[i], 
             rslt$plugin, 
             rslt_parametric$plugin,
             t(all_ci$tmle), 
             t(all_ci$aiptw), 
             tmle_cover, 
             aiptw_cover,
             t(all_ci_parametric$tmle), 
             t(all_ci_parametric$aiptw), 
             tmle_cover_parametric, 
             aiptw_cover_parametric,             
             truth) # add parametric G-comp output


    names(out) <- c("n", "seed", 
                    "total_plugin_est", "direct_plugin_est", "indirectM1_plugin_est", "indirectM2_plugin_est", "covarM1M2_plugin_est",
                    "total_parametric_plugin_est", "direct_parametric_plugin_est", "indirectM1_parametric_plugin_est", "indirectM2_parametric_plugin_est", "covarM1M2_parametric_plugin_est",
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
                    "total_parametric_tmle_est", "total_parametric_tmle_cil", "total_parametric_tmle_ciu",
                    "direct_parametric_tmle_est","direct_parametric_tmle_cil","direct_parametric_tmle_ciu",
                    "indirectM1_parametric_tmle_est","indirectM1_parametric_tmle_cil","indirectM1_parametric_tmle_ciu",
                    "indirectM2_parametric_tmle_est","indirectM2_parametric_tmle_cil","indirectM2_parametric_tmle_ciu",
                    "covarM1M2_parametric_tmle_est","covarM1M2_parametric_tmle_cil","covarM1M2_parametric_tmle_ciu",
                    "total_parametric_aiptw_est", "total_parametric_aiptw_cil", "total_parametric_aiptw_ciu",
                    "direct_parametric_aiptw_est","direct_parametric_aiptw_cil","direct_parametric_aiptw_ciu",
                    "indirectM1_parametric_aiptw_est","indirectM1_parametric_aiptw_cil","indirectM1_parametric_aiptw_ciu",
                    "indirectM2_parametric_aiptw_est","indirectM2_parametric_aiptw_cil","indirectM2_parametric_aiptw_ciu",
                    "covarM1M2_parametric_aiptw_est","covarM1M2_parametric_aiptw_cil","covarM1M2_parametric_aiptw_ciu",
                    "total_parametric_tmle_cover", "direct_parametric_tmle_cover", "indirectM1_parametric_tmle_cover", "indirectM2_parametric_tmle_cover", "covarM1M2_parametric_tmle_cover",
                    "total_parametric_aiptw_cover", "direct_parametric_aiptw_cover", "indirectM1_parametric_aiptw_cover", "indirectM2_parametric_aiptw_cover", "covarM1M2_parametric_aiptw_cover",
                    "true_total", "true_direct", "true_indirectM1", "true_indirectM2", "true_covarM1M2"
                    ) # add names of parametric G-comp output, e.g., total_gcomp_est, total_gcomp_cil, total_gcomp_ciu, total_gcomp_cover... 
    save(out, file = paste0("~/intermed/out/newfit6_", i, ".RData"))
  }
}

# merge job ###########################
if (args[1] == 'merge') {   
  # parameters
  ns <- c(250, 500, 1000, 2000)
  seed <- 1:1000
  parm <- expand.grid(n = ns, seed = seed)

  save_dir <- "~/intermed/out/"
  n_out <- 97
  
  rslt <- matrix(NA, nrow = nrow(parm), ncol = n_out)
  for(i in 1:nrow(parm)){
    tmp <- tryCatch({
      load(paste0(save_dir, "newfit1_", i, ".RData"))
      out
    }, error = function(e){
      rep(NA, n_out)
    })
    rslt[i, ] <- tmp
  }

  out <- data.frame(rslt)
  colnames(out) <- c("n", "seed", 
                    "total_plugin_est", "direct_plugin_est", "indirectM1_plugin_est", "indirectM2_plugin_est", "covarM1M2_plugin_est",
                    "total_parametric_plugin_est", "direct_parametric_plugin_est", "indirectM1_parametric_plugin_est", "indirectM2_parametric_plugin_est", "covarM1M2_parametric_plugin_est",
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
                    "total_parametric_tmle_est", "total_parametric_tmle_cil", "total_parametric_tmle_ciu",
                    "direct_parametric_tmle_est","direct_parametric_tmle_cil","direct_parametric_tmle_ciu",
                    "indirectM1_parametric_tmle_est","indirectM1_parametric_tmle_cil","indirectM1_parametric_tmle_ciu",
                    "indirectM2_parametric_tmle_est","indirectM2_parametric_tmle_cil","indirectM2_parametric_tmle_ciu",
                    "covarM1M2_parametric_tmle_est","covarM1M2_parametric_tmle_cil","covarM1M2_parametric_tmle_ciu",
                    "total_parametric_aiptw_est", "total_parametric_aiptw_cil", "total_parametric_aiptw_ciu",
                    "direct_parametric_aiptw_est","direct_parametric_aiptw_cil","direct_parametric_aiptw_ciu",
                    "indirectM1_parametric_aiptw_est","indirectM1_parametric_aiptw_cil","indirectM1_parametric_aiptw_ciu",
                    "indirectM2_parametric_aiptw_est","indirectM2_parametric_aiptw_cil","indirectM2_parametric_aiptw_ciu",
                    "covarM1M2_parametric_aiptw_est","covarM1M2_parametric_aiptw_cil","covarM1M2_parametric_aiptw_ciu",
                    "total_parametric_tmle_cover", "direct_parametric_tmle_cover", "indirectM1_parametric_tmle_cover", "indirectM2_parametric_tmle_cover", "covarM1M2_parametric_tmle_cover",
                    "total_parametric_aiptw_cover", "direct_parametric_aiptw_cover", "indirectM1_parametric_aiptw_cover", "indirectM2_parametric_aiptw_cover", "covarM1M2_parametric_aiptw_cover",
                    "truth_total", "truth_direct", "truth_indirectM1", "truth_indirectM2", "truth_covarM1M2"
                    )

  save(out, file = "~/intermed/all_out_intermed_new1.RData")

  get_bias <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- mean(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      aiptw <- mean(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      tmle2 <- mean(x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      aiptw2 <- mean(x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      plugin <- mean(x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      plugin2 <- mean(x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2))
    })
  }

  get_rootnbias <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
      aiptw <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
      tmle2 <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
      aiptw2 <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
      plugin <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
      plugin2 <- mean(x$n[1]^(1/2)*(x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)]), na.rm = TRUE)
      return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2))
    })
  }

  get_mse <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- mean((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      aiptw <- mean((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      tmle2 <- mean((x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      aiptw2 <- mean((x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      plugin <- mean((x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      plugin2 <- mean((x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2))
    })
  }

  get_nmse <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      aiptw <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      tmle2 <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      aiptw2 <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      plugin <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      plugin2 <- x$n[1]*mean((x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
      return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2))
    })
  }
  get_sd <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- sd(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))], na.rm = TRUE)
      aiptw <- sd(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] , na.rm = TRUE)
      tmle2 <- sd(x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))], na.rm = TRUE)
      aiptw2 <- sd(x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] , na.rm = TRUE)
      plugin <- sd(x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))], na.rm = TRUE)
      plugin2 <- sd(x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] , na.rm = TRUE)
      return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2))
    })
  }
  get_rootnsd <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))], na.rm = TRUE)
      aiptw <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] , na.rm = TRUE)
      tmle2 <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))], na.rm = TRUE)
      aiptw2 <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] , na.rm = TRUE)
      plugin <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_plugin_est"), colnames(x))], na.rm = TRUE)
      plugin2 <- x$n[1]^(1/2)*sd(x[ , grep(paste0(which_eff, "_parametric_plugin_est"), colnames(x))] , na.rm = TRUE)
      return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2))
    })
  }
  get_se_est <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- mean(x[ , grep(paste0(which_eff, "_tmle_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_tmle_cil"), colnames(x))]) / (2 * 1.96)
      aiptw <- mean(x[ , grep(paste0(which_eff, "_aiptw_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_aiptw_cil"), colnames(x))]) / (2 * 1.96)
      tmle2 <- mean(x[ , grep(paste0(which_eff, "_parametric_tmle_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_parametric_tmle_cil"), colnames(x))]) / (2 * 1.96)
      aiptw2 <- mean(x[ , grep(paste0(which_eff, "_parametric_aiptw_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_parametric_aiptw_cil"), colnames(x))]) / (2 * 1.96)
      plugin <- mean(x[ , grep(paste0(which_eff, "_plugin_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_plugin_cil"), colnames(x))]) / (2 * 1.96)
      plugin2 <- mean(x[ , grep(paste0(which_eff, "_parametric_plugin_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_parametric_plugin_cil"), colnames(x))]) / (2 * 1.96)
      return(c(tmle, aiptw, tmle2, aiptw2, plugin, plugin2))
    })
  }
  get_cover <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- mean(x[ , grep(paste0(which_eff, "_tmle_cover"), colnames(x))], na.rm = TRUE)
      aiptw <- mean(x[ , grep(paste0(which_eff, "_aiptw_cover"), colnames(x))], na.rm = TRUE)
      tmle2 <- mean(x[ , grep(paste0(which_eff, "_parametric_tmle_cover"), colnames(x))], na.rm = TRUE)
      aiptw2 <- mean(x[ , grep(paste0(which_eff, "_parametric_aiptw_cover"), colnames(x))], na.rm = TRUE)
      return(c(tmle, aiptw, tmle2, aiptw2))
    })
  }
  
  get_density_orc <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle <- density((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / sd(x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))]), na.rm = TRUE)
      aiptw <- density((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / sd(x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))]), na.rm = TRUE)
      tmle2 <- density((x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / sd(x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))]), na.rm = TRUE)
      aiptw2 <- density((x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / sd(x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))]), na.rm = TRUE)
      return(list(tmle, aiptw, tmle2, aiptw2))
    })
  }
  get_density_est <- function(out, which_eff = "total"){
    by(out, out$n, function(x){
      tmle_se <- (x[ , grep(paste0(which_eff, "_tmle_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_tmle_cil"), colnames(x))]) / (2 * 1.96)
      aiptw_se <- (x[ , grep(paste0(which_eff, "_aiptw_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_aiptw_cil"), colnames(x))]) / (2 * 1.96)
      tmle <- density((x[ , grep(paste0(which_eff, "_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / tmle_se, na.rm = TRUE)
      aiptw <- density((x[ , grep(paste0(which_eff, "_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / aiptw_se, na.rm = TRUE)
      tmle_se2 <- (x[ , grep(paste0(which_eff, "_parametric_tmle_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_parametric_tmle_cil"), colnames(x))]) / (2 * 1.96)
      aiptw_se2 <- (x[ , grep(paste0(which_eff, "_parametric_aiptw_ciu"), colnames(x))] - x[ , grep(paste0(which_eff, "_parametric_aiptw_cil"), colnames(x))]) / (2 * 1.96)
      tmle2 <- density((x[ , grep(paste0(which_eff, "_parametric_tmle_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / tmle_se2, na.rm = TRUE)
      aiptw2 <- density((x[ , grep(paste0(which_eff, "_parametric_aiptw_est"), colnames(x))] - x[,paste0("truth_", which_eff)]) / aiptw_se2, na.rm = TRUE)
      return(list(tmle, aiptw, tmle2, aiptw2))
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
      oracle_se_tmle2 <- sd(x[, paste0(which_eff, "_parametric_tmle_est")], na.rm = TRUE)
      oracle_se_aiptw2 <- sd(x[, paste0(which_eff, "_parametric_aiptw_est")], na.rm = TRUE)
      tmle_cil2 <- x[, paste0(which_eff, "_parametric_tmle_est")] - 1.96 * oracle_se_tmle2
      tmle_ciu2 <- x[, paste0(which_eff, "_parametric_tmle_est")] + 1.96 * oracle_se_tmle2
      aiptw_cil2 <- x[, paste0(which_eff, "_parametric_aiptw_est")] - 1.96 * oracle_se_aiptw2
      aiptw_ciu2 <- x[, paste0(which_eff, "_parametric_aiptw_est")] + 1.96 * oracle_se_aiptw2
      tmle_cover2 <- tmle_cil2 < x[, paste0("truth_", which_eff)] & tmle_ciu2 > x[, paste0("truth_", which_eff)]
      aiptw_cover2 <- aiptw_cil2 < x[, paste0("truth_", which_eff)] & aiptw_ciu2 > x[, paste0("truth_", which_eff)]
      return(c(mean(tmle_cover, na.rm = TRUE), mean(aiptw_cover, na.rm = TRUE),
               mean(tmle_cover2, na.rm = TRUE), mean(aiptw_cover2, na.rm = TRUE)))
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
        # 2 means parametric
        num_cols <- ncol(a)
        if(ncol(a) == 5){
          colnames(a) <- c("tmle", "os", "tmle2", "os2", "n")
        }else{
          colnames(a) <- c("tmle", "os", "tmle2", "os2", "plugin", "plugin2", "n")
        }
        row.names(a) <- NULL
        return(a)
      }, simplify = FALSE)
    }else{
      all_out <- sapply(all_eff, function(x){
        do.call(summary_fn, args = list(out = out, which_eff = x))
      })
      return(all_out)
    }
  }
  # # fix coverage for indirect M2
  # out$indirectM2_tmle_cover <- out$indirectM2_tmle_cil < out$truth_indirectM2 & out$indirectM2_tmle_ciu > out$truth_indirectM2
  # out$indirectM2_aiptw_cover <- out$indirectM2_aiptw_cil < out$truth_indirectM2 & out$indirectM2_aiptw_ciu > out$truth_indirectM2

  load("~/Dropbox/Emory/Flu/inter_med/sim_rslt/all_out_intermed_new1.RData")
  
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
  pdf("~/Dropbox/Emory/Flu/inter_med/estimation_new1.pdf",
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
                                 parametric = FALSE,
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
        lines(all_density_orc[,which_eff][[i]][[1 + ifelse(parametric, 2, 0)]],
              col = grays[i], lty = 2)
      }

      # plot sampling distribution of TMLE by sample size
      blank_plot(xlim = c(-3.5, 3.5), ylim = c(0, 0.55),
                 xlab = xlab2, ylab = "Density")
      axis(side = 1, at = seq(-3, 3, by = 1))
      x_seq <- seq(-4, 4, length = 2000)
      lines(x = x_seq, y = dnorm(x_seq), lty = 1, lwd = 2)
      grays <- paste0("gray", c(80, 60, 40, 20))
      for(i in 1:4){
        lines(all_density_est[,which_eff][[i]][[1 + ifelse(parametric, 2, 0)]],
              col = grays[i], lty = 2)
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
        lines(all_density_orc[,which_eff][[i]][[2 + ifelse(parametric, 2, 0)]], col = grays[i], lty = 2)
      }

      # plot sampling distribution of TMLE by sample size
      blank_plot(xlim = c(-3.5, 3.5), ylim = c(0, 0.55),
                 xlab = xlab2, ylab = "Density")
      axis(side = 1, at = seq(-3, 3, by = 1))
      x_seq <- seq(-4, 4, length = 2000)
      lines(x = x_seq, y = dnorm(x_seq), lty = 1, lwd = 2)
      grays <- paste0("gray", c(80, 60, 40, 20))
      for(i in 1:4){
        lines(all_density_est[,which_eff][[i]][[2 + ifelse(parametric, 2, 0)]], col = grays[i], lty = 2)
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



pdf("~/Dropbox/Emory/Flu/inter_med/inference_tmle_new1.pdf",
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

pdf("~/Dropbox/Emory/Flu/inter_med/inference_aiptw_new1.pdf",
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


pdf("~/Dropbox/Emory/Flu/inter_med/inference_parametric_tmle_new1.pdf",
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
                     parametric = TRUE,
                     all_density_orc = all_density_orc, 
                     all_density_est = all_density_est, 
                     all_cover = all_cover, 
                     all_cover_orc = all_cover_orc)
}
dev.off()

pdf("~/Dropbox/Emory/Flu/inter_med/inference_parametric_aiptw_new1.pdf",
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
                     parametric = TRUE,
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
