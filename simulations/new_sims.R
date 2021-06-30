make_data <- function(n = 1e2, 
                      A_success = "plogis(-1 + 0.25*C$C1 + 0.5*C$C2)",
                      M1_success = "plogis(-1 + 0.25 * C$C1 + 0.25 * A)",
                      M2_success = "plogis(-1 + 0.25 * C$C1 + 0.35 * A)",
                      M1M2_threshold = 5,
                      Y_success = "plogis(-1 + C$C1 - C$C2 + 0.5 * M1 + 0.5 * M2 + 1 * A)",
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
