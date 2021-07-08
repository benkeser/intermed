make_data_discrete <- function(n = 1e2, 
                      A_success = "plogis(-1 + 0.125 * C$C1 + 0.25*C$C2)",
                      M1_success = "plogis(-1.1 + 0.45 * C$C1 + 0.125 * A)",
                      M2_success = "plogis(-1.1 + 0.15 * C$C1 + 0.2 * C$C2 - 0.2 * A)",
                      M1M2_threshold = 6,
                      Y_success = "plogis(-1 + C$C1 - C$C2 + 0.25 * M1 + 0.25 * M2 + 0.25 * A)",
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
    Y <- rbinom(n, 1, Qbar0)
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

make_data_discrete2 <- function(..., M1M2_threshold = 3){
  make_data_discrete(..., M1M2_threshold = M1M2_threshold)
}

make_data_discrete3 <- function(..., M1M2_threshold = 12){
  make_data_discrete(..., M1M2_threshold = M1M2_threshold)
}

do_one_discrete <- function(i, parm, j = NULL){
  # set seed
  set.seed(parm$seed[i])

  # make data
  data <- do.call(paste0('make_data_discrete', j), args = list(n = parm$n[i]))

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
                             targeted_se = FALSE, 
                             return_models = FALSE,
                             verbose = FALSE,
                             stratify = FALSE,
                             max_iter = 10)
  set.seed(1234)
  truth <- c(unlist(make_data(n = 1e6, get_truth = TRUE)), 0)

  all_ci <- ci(rslt, est = c("tmle", "aiptw"))
  all_ci_parametric <- ci(rslt_parametric, est = c("tmle", "aiptw"))
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
           truth)


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
                  )
  return(out)
}

make_data_continuous <- function(n = 1e2, 
                      A_success = "plogis(-1 + 0.25 * C$C1 - 0.5 * C$C2 * C$C1)",
                      M1_success = "C$C1 - 0.5 * A * C$C1",
                      M2_success = "C$C2 - C$C2 * A",
                      M1M2_threshold = 6,
                      Y_success = "plogis(-1 + C$C1 - C$C2 + 0.5 * M2 + 0.25 * A - 0.5 * M1)",
                      get_truth = FALSE){
  C <- data.frame(C1 = runif(n), C2 = runif(n), C5 = rbinom(n, 1, 0.5),
                  C4 = rbinom(n, 1, 0.25), C3 = runif(n))
  if(!get_truth){
    g0 <- eval(parse(text = A_success))
    A <- rbinom(n, 1, g0)
    success_prob_M1 <- eval(parse(text = M1_success))
    success_prob_M2 <- eval(parse(text = M2_success))
    M1 <- rnorm(n, success_prob_M1, 1)
    M2 <- rnorm(n, success_prob_M2, 1)
    # M1[M1 > (M1M2_threshold - 1)] <- M1M2_threshold
    # M2[M2 > (M1M2_threshold - 1)] <- M1M2_threshold
    Qbar0 <- eval(parse(text = Y_success))
    Y <- rbinom(n, 1, Qbar0)
    return(list(C = C, A = A, M1 = M1, M2 = M2, Y = Y))
  }else{
    success_prob_M1_A1 <- eval(parse(text = gsub("A", "1", M1_success)))
    success_prob_M1_A0 <- eval(parse(text = gsub("A", "0", M1_success)))
    success_prob_M2_A1 <- eval(parse(text = gsub("A", "1", M2_success)))
    success_prob_M2_A0 <- eval(parse(text = gsub("A", "0", M2_success)))
    M1_A0 <- rnorm(n, success_prob_M1_A0, 1)
    M2_A0 <- rnorm(n, success_prob_M2_A0, 1)
    M1_A1 <- rnorm(n, success_prob_M1_A1, 1)
    M2_A1 <- rnorm(n, success_prob_M2_A1, 1)

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

do_one_continuous <- function(i, parm, ...){
  # set seed
  set.seed(parm$seed[i])

  # make data
  data <- make_data_continuous(n = parm$n[i])

  # let's do 5, 10, 20, 40 bins
  cut_M <- function(n_breaks = 6, M){
    as.numeric(cut(M, breaks = quantile(M, p = seq(0, 1, length = n_breaks)), include.lowest = TRUE))
  }

  M1_cut <- cut_M(M = data$M1, n_breaks = 6)
  M2_cut <- cut_M(M = data$M2, n_breaks = 6)
  rslt_1 <- intermed(Y = data$Y, C = data$C, 
                              M1 = M1_cut, M2 = M2_cut, A = data$A, 
                             a = 1, 
                             a_star = 0,
                             n_SL = 1,
                             SL_Qbar = c("SL.earth","SL.step.interaction","SL.ranger"),
                             SL_g = c("SL.earth","SL.step.interaction","SL.ranger"),
                             SL_Q_M = list(M1 = c("SL.earth","SL.step.interaction","SL.ranger"), 
                                           M2 = c("SL.earth","SL.step.interaction","SL.ranger")),
                             tolg = 1e-3, 
                             targeted_se = FALSE, 
                             return_models = FALSE,
                             verbose = FALSE,
                             stratify = FALSE,
                             max_iter = 10)

  M1_cut <- cut_M(M = data$M1, n_breaks = 11)
  M2_cut <- cut_M(M = data$M2, n_breaks = 11)
  rslt_2 <- intermed(Y = data$Y, C = data$C, 
                     M1 = M1_cut, M2 = M2_cut, A = data$A, 
                     a = 1, 
                     a_star = 0,
                     n_SL = 1,
                     SL_Qbar = c("SL.earth","SL.step.interaction","SL.ranger"),
                     SL_g = c("SL.earth","SL.step.interaction","SL.ranger"),
                     SL_Q_M = list(M1 = c("SL.earth","SL.step.interaction","SL.ranger"), 
                                   M2 = c("SL.earth","SL.step.interaction","SL.ranger")),
                     tolg = 1e-3, 
                     targeted_se = FALSE, 
                     return_models = FALSE,
                     verbose = FALSE,
                     stratify = FALSE,
                     max_iter = 10)

  # get truth
  set.seed(1234)
  truth <- c(unlist(make_data(n = 1e6, get_truth = TRUE)), 0)
  # get confidence intervals
  all_ci_1 <- ci(rslt_1, est = c("tmle", "aiptw"))
  all_ci_2 <- ci(rslt_2, est = c("tmle", "aiptw"))

  # check for truth in CIs
  in_ci <- function(truth, ci){
    truth > min(ci) & truth < max(ci)
  }

  tmle_cover <- aiptw_cover <- rep(NA, 5)
  for(j in seq_len(5)){
    tmle_cover[j] <- in_ci(truth[j], all_ci_1$tmle[j, c(2,3)])
    aiptw_cover[j] <- in_ci(truth[j], all_ci_1$aiptw[j, c(2,3)])
  }

  tmle_cover2 <- aiptw_cover2 <- rep(NA, 5)
  for(j in seq_len(5)){
    tmle_cover2[j] <- in_ci(truth[j], all_ci_2$tmle[j, c(2,3)])
    aiptw_cover2[j] <- in_ci(truth[j], all_ci_2$aiptw[j, c(2,3)])
  }
  

  out <- c(parm$n[i], parm$seed[i], 
           rslt_1$plugin, 
           rslt_2$plugin,
           t(all_ci_1$tmle), 
           t(all_ci_1$aiptw), 
           tmle_cover, 
           aiptw_cover,
           t(all_ci_2$tmle), 
           t(all_ci_2$aiptw), 
           tmle_cover2, 
           aiptw_cover2,             
           truth)

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
                  )
  return(out)
}


make_data_robustness <- function(n = 1e2, 
                      A_success = "plogis(-1 + 0.25 * C$C1 - 0.5 * C$C2 * C$C1)",
                      M1_success = "plogis(-1.1 + 0.25 * C$C1 - 0.5 * A * C$C1)",
                      M2_success = "plogis(-1.1 + 0.25 * C$C2 - 0.5 * C$C2 * A)",
                      M1M2_threshold = 6,
                      Y_success = "plogis(-1 + C$C1 - C$C2 + 0.25 * M2 + 0.25 * M1 * A - 0.5 * M1 * A)",
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
    Y <- rbinom(n, 1, Qbar0)
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

do_one_robustness <- function(i, parm, ...){
  # set seed
  set.seed(parm$seed[i])

  # make data
  data <- make_data_robustness(n = parm$n[i])

  rslt_parametric_1 <- intermed(Y = data$Y, C = data$C, 
                              M1 = data$M1, M2 = data$M2, A = data$A, 
                             a = 1, 
                             a_star = 0,
                             n_SL = 1,
                             SL_Qbar = "SL.step.interaction",
                             SL_g = "SL.glm",
                             SL_Q_M = list(M1 = "SL.step.interaction", 
                                           M2 = "SL.step.interaction"),
                             tolg = 1e-3, 
                             targeted_se = FALSE, 
                             return_models = FALSE,
                             verbose = FALSE,
                             stratify = FALSE,
                             max_iter = 10)
  rslt_parametric_2 <- intermed(Y = data$Y, C = data$C, 
                              M1 = data$M1, M2 = data$M2, A = data$A, 
                             a = 1, 
                             a_star = 0,
                             n_SL = 1,
                             SL_Qbar = "SL.glm",
                             SL_g = "SL.step.interaction",
                             SL_Q_M = list(M1 = "SL.step.interaction", 
                                           M2 = "SL.step.interaction"),
                             tolg = 1e-3, 
                             targeted_se = FALSE, 
                             return_models = FALSE,
                             verbose = FALSE,
                             stratify = FALSE,
                             max_iter = 10)

  rslt_parametric_3 <- intermed(Y = data$Y, C = data$C, 
                              M1 = data$M1, M2 = data$M2, A = data$A, 
                             a = 1, 
                             a_star = 0,
                             n_SL = 1,
                             SL_Qbar = "SL.step.interaction",
                             SL_g = "SL.step.interaction",
                             SL_Q_M = list(M1 = "SL.glm", 
                                           M2 = "SL.step.interaction"),
                             tolg = 1e-3, 
                             targeted_se = FALSE, 
                             return_models = FALSE,
                             verbose = FALSE,
                             stratify = FALSE,
                             max_iter = 10)

    rslt_parametric_4 <- intermed(Y = data$Y, C = data$C, 
                            M1 = data$M1, M2 = data$M2, A = data$A, 
                           a = 1, 
                           a_star = 0,
                           n_SL = 1,
                           SL_Qbar = "SL.step.interaction",
                           SL_g = "SL.step.interaction",
                           SL_Q_M = list(M1 = "SL.step.interaction", 
                                         M2 = "SL.glm"),
                           tolg = 1e-3, 
                           targeted_se = FALSE, 
                           return_models = FALSE,
                           verbose = FALSE,
                           stratify = FALSE,
                           max_iter = 10)
  
  # get truth
  set.seed(1234)
  truth <- c(unlist(make_data(n = 1e6, get_truth = TRUE)), 0)

  # get confidence intervals
  all_ci_parametric_1 <- ci(rslt_parametric_1, est = c("tmle", "aiptw"))
  all_ci_parametric_2 <- ci(rslt_parametric_2, est = c("tmle", "aiptw"))
  all_ci_parametric_3 <- ci(rslt_parametric_3, est = c("tmle", "aiptw"))
  all_ci_parametric_4 <- ci(rslt_parametric_4, est = c("tmle", "aiptw"))

  # check for truth in CIs
  in_ci <- function(truth, ci){
    truth > min(ci) & truth < max(ci)
  }
  tmle_cover <- c(
    in_ci(truth[3], all_ci_parametric_1$tmle[3, c(2,3)]),
    in_ci(truth[3], all_ci_parametric_2$tmle[3, c(2,3)]),
    in_ci(truth[3], all_ci_parametric_3$tmle[3, c(2,3)]),
    in_ci(truth[3], all_ci_parametric_4$tmle[3, c(2,3)])
  )

  aiptw_cover <- c(
    in_ci(truth[3], all_ci_parametric_1$aiptw[3, c(2,3)]),
    in_ci(truth[3], all_ci_parametric_2$aiptw[3, c(2,3)]),
    in_ci(truth[3], all_ci_parametric_3$aiptw[3, c(2,3)]),
    in_ci(truth[3], all_ci_parametric_4$aiptw[3, c(2,3)])
  )


  out <- c(parm$n[i], parm$seed[i], 
           rslt_parametric_1$plugin[3], 
           rslt_parametric_2$plugin[3], 
           rslt_parametric_3$plugin[3], 
           rslt_parametric_4$plugin[3], 
           t(all_ci_parametric_1$tmle)[,3], 
           t(all_ci_parametric_2$tmle)[,3], 
           t(all_ci_parametric_3$tmle)[,3], 
           t(all_ci_parametric_4$tmle)[,3], 
           t(all_ci_parametric_1$aiptw)[,3], 
           t(all_ci_parametric_2$aiptw)[,3], 
           t(all_ci_parametric_3$aiptw)[,3], 
           t(all_ci_parametric_4$aiptw)[,3], 
           tmle_cover, 
           aiptw_cover,
           truth[3])


  names(out) <- c("n", "seed", 
                  "indirectM1_plugin1_est",
                  "indirectM1_plugin2_est",
                  "indirectM1_plugin3_est",
                  "indirectM1_plugin4_est",
                  "indirectM1_tmle1_est","indirectM1_tmle1_cil","indirectM1_tmle1_ciu",
                  "indirectM1_tmle2_est","indirectM1_tmle2_cil","indirectM1_tmle2_ciu",
                  "indirectM1_tmle3_est","indirectM1_tmle3_cil","indirectM1_tmle3_ciu",
                  "indirectM1_tmle4_est","indirectM1_tmle4_cil","indirectM1_tmle4_ciu",
                  "indirectM1_aiptw1_est","indirectM1_aiptw1_cil","indirectM1_aiptw1_ciu",
                  "indirectM1_aiptw2_est","indirectM1_aiptw2_cil","indirectM1_aiptw2_ciu",
                  "indirectM1_aiptw3_est","indirectM1_aiptw3_cil","indirectM1_aiptw3_ciu",
                  "indirectM1_aiptw4_est","indirectM1_aiptw4_cil","indirectM1_tmle4_ciu",
                  "indirectM1_tmle1_cover", 
                  "indirectM1_tmle2_cover", 
                  "indirectM1_tmle3_cover", 
                  "indirectM1_tmle4_cover", 
                  "indirectM1_aiptw1_cover", 
                  "indirectM1_aiptw2_cover", 
                  "indirectM1_aiptw3_cover", 
                  "indirectM1_aiptw4_cover",
                  "true_indirectM1"
                  )
  return(out)
}

make_data_gcomp <- function(n = 1e2, 
                        A_success = "plogis(-1 + 0.125 * C$C1 + 0.25*C$C2)",
                        M1_success = "plogis(-0.25 + 0.45 * C$C1 + 0.125 * A + 0.5 * C$C1 * A)",
                        M2_success = "plogis(-0.25 + 0.2 * C$C2 - 0.2 * A + 0.25 * C$C2 * A)",
                        M1M2_threshold = 5,
                        Y_success = "-1 + C$C1 - C$C2 + 0.05 * M1 + 0.8 * (M1 - 3)^2 - 0.25 * M2 - 0.8*(M2 - 3)^2 + 0.125 * A + 0.15 * M1*A - 0.15*M2*A",
                        get_truth = FALSE){
  C <- data.frame(C1 = runif(n), C2 = runif(n), C5 = rbinom(n, 1, 0.5),
                  C4 = rbinom(n, 1, 0.25), C3 = runif(n))
  if(!get_truth){
    g0 <- eval(parse(text = A_success))
    A <- rbinom(n, 1, g0)
    success_prob_M1 <- eval(parse(text = M1_success))
    success_prob_M2 <- eval(parse(text = M2_success))
    M1 <- rbinom(n, 5, success_prob_M1)
    M2 <- rbinom(n, 5, success_prob_M2)
    M1[M1 > (M1M2_threshold - 1)] <- M1M2_threshold
    M2[M2 > (M1M2_threshold - 1)] <- M1M2_threshold
    Qbar0 <- eval(parse(text = Y_success))
    Y <- Qbar0 + rnorm(n, 1)
    return(list(C = C, A = A, M1 = M1, M2 = M2, Y = Y))
  }else{
    success_prob_M1_A1 <- eval(parse(text = gsub("A", "1", M1_success)))
    success_prob_M1_A0 <- eval(parse(text = gsub("A", "0", M1_success)))
    success_prob_M2_A1 <- eval(parse(text = gsub("A", "1", M2_success)))
    success_prob_M2_A0 <- eval(parse(text = gsub("A", "0", M2_success)))
    M1_A0 <- rbinom(n, 5, success_prob_M1_A0)
    M2_A0 <- rbinom(n, 5, success_prob_M2_A0)
    M1_A1 <- rbinom(n, 5, success_prob_M1_A1)
    M2_A1 <- rbinom(n, 5, success_prob_M2_A1)
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
                indirect_M2 = indirect_effect_M2,
                covar_M1M2 = total_effect - (direct_effect + indirect_effect_M1 + indirect_effect_M2
          )))
  }
}

do_one_gcomp <- function(i, parm, ...){
  # set seed
  set.seed(parm$seed[i])

  # make data
  data <- make_data_gcomp(n = parm$n[i], M1M2_threshold = 5)

  rslt <- intermed(Y = data$Y, C = data$C, M1 = data$M1, M2 = data$M2, A = data$A, 
                   a = 1, 
                   a_star = 0,
                   n_SL = 1,
                   SL_Qbar = c("SL.glm", "SL.earth", "SL.ranger", "SL.gam"),
                   SL_g = c("SL.glm", "SL.earth", "SL.ranger"),
                   SL_Q_M = list(M1 = c("SL.glm", "SL.earth", "SL.ranger", "SL.gam"), 
                                 M2 = c("SL.glm", "SL.earth", "SL.ranger", "SL.gam")),
                   tolg = 1e-2, 
                   targeted_se = TRUE, 
                   return_models = FALSE,
                   verbose = FALSE,
                   stratify = FALSE,
                   max_iter = 2)

  rslt_parametric <- intermed(Y = data$Y, C = data$C, 
                              M1 = data$M1, M2 = data$M2, A = data$A, 
                             a = 1, 
                             a_star = 0,
                             n_SL = 1,
                             SL_Qbar = NULL,
                             glm_Qbar = "M1*A + I(M1^2) + M2*A + I(M2^2) + C1 + C2",
                             glm_g = "C1 + C2",
                             SL_g = NULL,
                             glm_Q_M = list(M1 = "C1 + A", 
                                           M2 = "C2 + A"),
                             SL_Q_M = NULL,
                             tolg = 1e-2, 
                             targeted_se = TRUE, 
                             return_models = FALSE,
                             verbose = FALSE,
                             stratify = FALSE,
                             max_iter = 0)

  # create a data frame
  full_data <- Reduce(cbind, data)
  colnames(full_data) <- c(colnames(data$C), names(data)[-1])
  full_data = full_data[,c(1:2, 5:3, 6:9)]
  # get g-comp on original data
  gcomp_est0 <- unlist(get_gcomp(full_data, version = 0))
  gcomp_est1 <- unlist(get_gcomp(full_data, version = 1))
  gcomp_est2 <- unlist(get_gcomp(full_data, version = 2))
  # parametric g-comp estimator
  gcomp_ci0 <- get_gcomp_boot(full_data, 0)
  gcomp_cil0 <- gcomp_ci0$cil
  gcomp_ciu0 <- gcomp_ci0$ciu

  gcomp_ci1 <- get_gcomp_boot(full_data, 1)
  gcomp_cil1 <- gcomp_ci1$cil
  gcomp_ciu1 <- gcomp_ci1$ciu

  gcomp_ci2 <- get_gcomp_boot(full_data, 2)
  gcomp_cil2 <- gcomp_ci2$cil
  gcomp_ciu2 <- gcomp_ci2$ciu

  # get truth
  set.seed(1234)
  truth <- c(unlist(make_data(n = 1e6, get_truth = TRUE, M1M2_threshold = 5)))

  # get confidence intervals
  all_ci <- ci(rslt, est = c("tmle", "aiptw"))
  all_ci_parametric <- ci(rslt_parametric, est = c("tmle", "aiptw"))

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
  gcomp_cover0 <- rep(NA, 5)
  gcomp_cover1 <- rep(NA, 5)
  gcomp_cover2 <- rep(NA, 5)
  for(j in seq_len(5)){
    gcomp_cover0[j] <- in_ci(truth[j], c(gcomp_cil0[j], gcomp_ciu0[j]))
    gcomp_cover1[j] <- in_ci(truth[j], c(gcomp_cil1[j], gcomp_ciu1[j]))
    gcomp_cover2[j] <- in_ci(truth[j], c(gcomp_cil2[j], gcomp_ciu2[j]))
  }

  out <- c(parm$n[i], parm$seed[i], parm$version[i],
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
           gcomp_est0, gcomp_cil0, gcomp_ciu0,
           gcomp_cover0,
           gcomp_est1,
           gcomp_cil1, 
           gcomp_ciu1,
           gcomp_cover1,
           gcomp_est2,
           gcomp_cil2, 
           gcomp_ciu2,
           gcomp_cover2,
           truth) 

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
                  "total_gcomp_est0", "direct_gcomp_est0", "indirectM1_gcomp_est0", "indirectM2_gcomp_est0", "covarM1M2_gcomp_est0",
                  "total_gcomp_cil0", "direct_gcomp_cil0", "indirectM1_gcomp_cil0", "indirectM2_gcomp_cil0", "covarM1M2_gcomp_cil0",
                  "total_gcomp_ciu0", "direct_gcomp_ciu0", "indirectM1_gcomp_ciu0", "indirectM2_gcomp_ciu0", "covarM1M2_gcomp_ciu0",
                  "total_gcomp_cover0", "direct_gcomp_cover0", "indirectM1_gcomp_cover0", "indirectM2_gcomp_cover0", "covarM1M2_gcomp_cover0",
                  "total_gcomp_est1", "direct_gcomp_est1", "indirectM1_gcomp_est1", "indirectM2_gcomp_est1", "covarM1M2_gcomp_est1",
                  "total_gcomp_cil1", "direct_gcomp_cil1", "indirectM1_gcomp_cil1", "indirectM2_gcomp_cil1", "covarM1M2_gcomp_cil1",
                  "total_gcomp_ciu1", "direct_gcomp_ciu1", "indirectM1_gcomp_ciu1", "indirectM2_gcomp_ciu1", "covarM1M2_gcomp_ciu1",
                  "total_gcomp_cover1", "direct_gcomp_cover1", "indirectM1_gcomp_cover1", "indirectM2_gcomp_cover1", "covarM1M2_gcomp_cover1",
                  "total_gcomp_est2", "direct_gcomp_est2", "indirectM1_gcomp_est2", "indirectM2_gcomp_est2", "covarM1M2_gcomp_est2",
                  "total_gcomp_cil2", "direct_gcomp_cil2", "indirectM1_gcomp_cil2", "indirectM2_gcomp_cil2", "covarM1M2_gcomp_cil2",
                  "total_gcomp_ciu2", "direct_gcomp_ciu2", "indirectM1_gcomp_ciu2", "indirectM2_gcomp_ciu2", "covarM1M2_gcomp_ciu2",
                  "total_gcomp_cover2", "direct_gcomp_cover2", "indirectM1_gcomp_cover2", "indirectM2_gcomp_cover2", "covarM1M2_gcomp_cover2",
                  "true_total", "true_direct", "true_indirectM1", "true_indirectM2", "true_covarM1M2"
  ) 
  return(out)
}
