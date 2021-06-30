# simple simulation
make_data <- function(n){
	C <- data.frame(C1 = runif(n), C2 = runif(n))
	g0 <- plogis(-1 + C$C1 + C$C2)
	A <- rbinom(n, 1, g0)
	# A increases values of M1 and M2
	success_prob_M1 <- plogis(-1 + C$C1 + 1.5 * A)
	success_prob_M2 <- plogis(-1 + C$C1 + 0.5 * A)

	# probably should switch this to geometric distribution, so that we 
	# can be sure that we are nailing the mediator distribution estimation
	M1 <- rgeom(n, success_prob_M1)
	M2 <- rgeom(n, success_prob_M2)
	M1[M1 > 4] <- 5
	M2[M2 > 4] <- 5
	# higher values of M1 and M2 lead to higher probability that Y = 1
	Qbar0 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1 + 0.5 * M2 + A)

	Y <- rbinom(n, 1, Qbar0)
	return(list(C = C, A = A, M1 = M1, M2 = M2, Y = Y))
}

# figure out truth
get_truth <- function(big_n = 1e6){

	C <- data.frame(C1 = runif(big_n), C2 = runif(big_n))

	success_prob_M1_A0 <- plogis(-1 + C$C1 + 1.5 * 0)
	success_prob_M2_A0 <- plogis(-1 + C$C1 + 0.5 * 0)

	success_prob_M1_A1 <- plogis(-1 + C$C1 + 1.5 * 1)
	success_prob_M2_A1 <- plogis(-1 + C$C1 + 0.5 * 1)

	M1_A0 <- rgeom(big_n, success_prob_M1_A0)
	M2_A0 <- rgeom(big_n, success_prob_M2_A0)
	M1_A1 <- rgeom(big_n, success_prob_M1_A1)
	M2_A1 <- rgeom(big_n, success_prob_M2_A1)

	M1_A0[M1_A0 > 4] <- 5
	M2_A0[M2_A0 > 4] <- 5
	M1_A1[M1_A1 > 4] <- 5
	M2_A1[M2_A1 > 4] <- 5

	# total effect
	Qbar0_A1 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A1 + 0.5 * M2_A1 + 1)
	Qbar0_A0 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A0 + 0.5 * M2_A0 + 0)
	total_effect <- mean(Qbar0_A1 - Qbar0_A0)

	# direct effect	
	Qbar0_A1 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A0 + 0.5 * M2_A0 + 1)
	Qbar0_A0 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A0 + 0.5 * M2_A0 + 0)
	direct_effect <- mean(Qbar0_A1 - Qbar0_A0)

	# indirect effect through M1
	Qbar0_A1 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A1 + 0.5 * M2_A0 + 1)
	Qbar0_A0 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A0 + 0.5 * M2_A0 + 1)
	indirect_effect_M1 <- mean(Qbar0_A1 - Qbar0_A0)

	# indirect effect through M2
	Qbar0_A1 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A1 + 0.5 * M2_A1 + 1)
	Qbar0_A0 <- plogis(-1 + C$C1 - C$C2 + 0.5 * M1_A1 + 0.5 * M2_A0 + 1)
	indirect_effect_M2 <- mean(Qbar0_A1 - Qbar0_A0)

	return(list(total = total_effect,
				direct = direct_effect,
				indirect_M1 = indirect_effect_M1, 
				indirect_M2 = indirect_effect_M2))
}

# try out our estimation routine
devtools::load_all("~/Dropbox/R/intermed")
data <- make_data(n = 1000)
debug(intermed)
rslt <- intermed(Y = data$Y, C = data$C, M1 = data$M1, M2 = data$M2, A = data$A, 
                 a = 1, 
                 a_star = 0,
                 n_SL = 1,
                 glm_Qbar = "C1 + C2 + M1 + M2 + A", 
                 glm_g = "C1 + C2",
                 glm_Q_M = list(M1 = "C1*I(bin_id == 6) + A", M2 = "C1*I(bin_id == 6) + A"),
                 tolg = 1e-3, 
                 targeted_se = FALSE, 
                 return_models = TRUE,
                 verbose = FALSE,
                 stratify = FALSE,
                 max_iter = 3)
truth <- c(unlist(get_truth()), 0)
cbind(truth, rslt$tmle$est, rslt$aiptw$est)
