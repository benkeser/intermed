#' Target the outcome regression
#' 
#' @inheritParams intermed
#' 
#' @param all_mediator_values Grid of all combinations of M1 and M2
#' @param which_effects Which effects to include in the targeted of the outcome regression
#' @param bound_pred Should predictions be bounded?
#' @param epsilon_threshold To avoid extreme values of fluctuation parameters (indicating 
#' likely numerical instability), we truncate the value this parameter can take.
#' @return A list containing all outcome regression evaluations needed
#' for downstream calculations.

target_Qbar <- function(Y, A, M1, M2, a, a_star,
                        all_mediator_values,
                        Qbar_n,
                        gn,
                        Q_M_n, 
                        which_effects = c("direct", "indirectM1", "indirectM2"), 
                        bound_pred = FALSE, 
                        epsilon_threshold = 5, # truncate large epsilon values
                        ...){
	# create clever covariates
	# I(A_i = a)/g(a | C_i) * (Q_M1(M1_i | a, C_i) - Q_M1(M1_i | a_star, C_i)) * Q_M2(M2_i | a_star, C_i)/Q_{M1,M2}(M1_i, M2_i | a, C_i)
	# I(A_i = a)/g(a | C_i) Q_{M1,M2}(M1_i, M2_i | a_star, C_i) / Q_{M1,M2}(M1_i, M2_i | a, C_i)
	# I(A_i = a_star)/g(a_star | C_i)
	unique_M1_values <- unique(M1)
	unique_M2_values <- unique(M2)

	# get an n-length vector of M1 marginal at M1_i and C_i under a_star and a
	Q_M1_a_star <- unlist(mapply(Q_M_n_i = Q_M_n, M_i = M1, FUN = extract_marginal,
	                    MoreArgs = list(a_val = "a_star", 
	                                    unique_M_values = unique_M1_values,
	                                    mediator = "M1"),
	                    SIMPLIFY = FALSE), use.names = FALSE)
	Q_M1_a <- unlist(mapply(Q_M_n_i = Q_M_n, M_i = M1, FUN = extract_marginal,
	                    MoreArgs = list(a_val = "a", 
	                                    unique_M_values = unique_M1_values,
	                                    mediator = "M1"),
	                    SIMPLIFY = FALSE), use.names = FALSE)
	# get an n-length vector of M2 marginal at M2_i and C_i under a_star 
	Q_M2_a_star <- unlist(mapply(Q_M_n_i = Q_M_n, M_i = M2, FUN = extract_marginal,
	                    MoreArgs = list(a_val = "a_star", 
	                                    unique_M_values = unique_M2_values,
	                                    mediator = "M2"),
	                    SIMPLIFY = FALSE), use.names = FALSE)	
	Q_M2_a <- unlist(mapply(Q_M_n_i = Q_M_n, M_i = M2, FUN = extract_marginal,
	                    MoreArgs = list(a_val = "a", 
	                                    unique_M_values = unique_M2_values,
	                                    mediator = "M2"),
	                    SIMPLIFY = FALSE), use.names = FALSE)
	# get an n-length vector of joint dist. at M1_i, M2_i and C_i under a
	Q_M1M2_a <- unlist(mapply(Q_M_n_i = Q_M_n, M1_i = M1, M2_i = M2, FUN = extract_joint,
	                   MoreArgs = list(a_val = "a",
	                                   all_mediator_values = all_mediator_values),
	                   SIMPLIFY = FALSE), use.names = FALSE)	
    # get an n-length vector of joint dist. at M1_i, M2_i and C_i under a_star
	Q_M1M2_a_star <- unlist(mapply(Q_M_n_i = Q_M_n, M1_i = M1, M2_i = M2, FUN = extract_joint,
	                   MoreArgs = list(a_val = "a_star",
	                                   all_mediator_values = all_mediator_values),
	                   SIMPLIFY = FALSE), use.names = FALSE)

	# make "clever covariates"; note the first entry of gn corresponds to the probability 
	# that A = a_star | C = C_i; the second entry of gn corresponds to the probability that 
	# A = a | C = C_i
	# for indirect effect through M1
	n <- length(Y)
	H_bound <- sqrt(n) # try it out
	H_indirectM1_1_obs <- as.numeric(A == a)/gn[[2]] * Q_M1_a * Q_M2_a_star / Q_M1M2_a
	H_indirectM1_1_obs[H_indirectM1_1_obs > H_bound] <- H_bound
	H_indirectM1_2_obs <- - as.numeric(A == a)/gn[[2]] * Q_M1_a_star * Q_M2_a_star / Q_M1M2_a
	H_indirectM1_2_obs[H_indirectM1_2_obs < -H_bound] <- -H_bound
	# for indirect through M2
	H_indirectM2_1_obs <- as.numeric(A == a)/gn[[2]] * (Q_M2_a * Q_M1_a) / Q_M1M2_a
	H_indirectM2_1_obs[H_indirectM2_1_obs > H_bound] <- H_bound
	H_indirectM2_2_obs <- - as.numeric(A == a)/gn[[2]] * (Q_M2_a_star * Q_M1_a) / Q_M1M2_a # -H_indirectM1_1_obs
	H_indirectM2_2_obs[H_indirectM2_2_obs < -H_bound] <- -H_bound
	# for direct
	H_direct_1_obs <- as.numeric(A == a)/gn[[2]] * Q_M1M2_a_star / Q_M1M2_a
	H_direct_1_obs[H_direct_1_obs > H_bound] <- H_bound
	H_direct_2_obs <- as.numeric(A == a_star)/gn[[1]]
	H_direct_2_obs[H_direct_2_obs > H_bound] <- H_bound

	covariate_values <- NULL
	covariate_names <- NULL
	if("indirectM1" %in% which_effects){
		covariate_names <- c(covariate_names, "H_indirectM1_1_obs", "H_indirectM1_2_obs")
		covariate_values <- cbind(covariate_values, H_indirectM1_1_obs, H_indirectM1_2_obs)
	}
	if("indirectM2" %in% which_effects){
		covariate_names <- c(covariate_names, "H_indirectM2_1_obs")
		covariate_values <- cbind(covariate_values, H_indirectM2_1_obs)
		if(!("indirectM1" %in% which_effects)){
			covariate_names <- c(covariate_names, "H_indirectM2_2_obs")
			covariate_values <- cbind(covariate_values, H_indirectM2_2_obs)
		}
	}
	if("direct" %in% which_effects){
		covariate_names <- c(covariate_names, "H_direct_1_obs", "H_direct_2_obs")
		covariate_values <- cbind(covariate_values, H_direct_1_obs, H_direct_2_obs)
	}
	num_covariates <- dim(covariate_values)[2]

	# H1_obs <- as.numeric(A == a)/gn[[2]] * (Q_M1_a - Q_M1_a_star) * Q_M2_a_star / Q_M1M2_a
	# this is for indirect effect through M2
	# note that there is a shared component of H1 and H2, so one could also produce
	# three covariates to knock out the relevant terms. for now, we stick to two
	# H2_obs <- as.numeric(A == a)/gn[[2]] * (Q_M2_a - Q_M2_a_star) * Q_M1_a / Q_M1M2_a
	# this is for direct effect
	# H3_obs <- as.numeric(A == a)/gn[[2]] * Q_M1M2_a_star / Q_M1M2_a
	# this is for direct effect
	# H4_obs <- as.numeric(A == a_star)/gn[[1]]

	# make offset term
	Qbar_M1M2C_obs <- unlist(mapply(Qbar_n_i = Qbar_n, A_i = A, M1_i = M1, M2_i = M2,
	                         FUN = function(Qbar_n_i, A_i, M1_i, M2_i){
	                         	extract_Qbar_obs(
	                         	  Qbar_n_i = Qbar_n_i, 
	                         	  M1_i = M1_i, M2_i = M2_i, 
	                         	  all_mediator_values = all_mediator_values, 
                         		  a_val = ifelse(A_i == a_star, "a_star", "a")
             	                )
	                         }, SIMPLIFY = FALSE), use.names = FALSE)

	# try scaling by min and max of initial Qbar
	if(bound_pred){
		allQbar <- Reduce(c, lapply(Qbar_n, function(x){unlist(x[[1]], use.names = FALSE)}))
		ell_scale <- 0
		u_scale <- max(allQbar)
	}else{
		ell_scale <- min(Y)
		u_scale <- max(Y)
	}
	scaled_Y <- (Y - ell_scale)/(u_scale - ell_scale)
	scaled_offset <- SuperLearner::trimLogit((Qbar_M1M2C_obs - ell_scale) / (u_scale - ell_scale))

	target_data <- data.frame(covariate_values,
	                          scaled_Y = scaled_Y, 
	                          scaled_offset = scaled_offset)
	colnames(target_data)[1:num_covariates] <- paste0("H", 1:num_covariates)

	if(TRUE){
		llik <- function(eps, data){
			p <- plogis(data$scaled_offset + as.matrix(data[,1]) %*% eps)
			return(-sum(data$scaled_Y * log(p) + (1 - data$scaled_Y) * log(1 - p)))
		}
		# fluc_mod <- optim(par = rep(0, num_covariates), fn = llik, data = target_data)
		# epsilon <- fluc_mod$par
		# try doing it one by one
		epsilon <- NULL
		for(i in 1:num_covariates){			
			fluc_modi <- optim(par = rep(0, 1), fn = llik, data = target_data[,c(i,6,7)],
			                   method = "Brent", lower = -10, upper = 10)
			epsilon <- c(epsilon, fluc_modi$par)
			target_data$scaled_offset <- target_data$scaled_offset + epsilon[i] * target_data[,i]
		}


	}else{
		fluc_mod <- suppressWarnings(glm(as.formula(paste0("scaled_Y ~ offset(scaled_offset) -1 + ",
		                                                   paste0("H", 1:num_covariates, collapse = "+"))), 
		                data = target_data, family = binomial(), start = rep(0, num_covariates)))
		epsilon <- coef(fluc_mod)
	}
	# truncate epsilon
	epsilon[epsilon > epsilon_threshold] <- epsilon_threshold
	epsilon[epsilon < -epsilon_threshold] <- -epsilon_threshold

	# we'll actually need to evaluate the targeted fit under each value of (M1, M2)
	# so go through and replace Qbar_n_i$Qbar_a_0[[1]] and [[2]] with respective values
	# reminder that the structure of Q_M_n is 
	# outer [[1]] = distributions under a_star
	# outer [[2]] = distributions under a
	# inner [[1]] = joint (M1, M2)
	# inner [[2]] = marginal M1
	# inner [[3]] = marginal M2
	# max_Y <- max(Y); min_Y <- min(Y)
	max_Y <- u_scale; min_Y <- ell_scale
	Qbar_n_tmle <- mapply(Q_M_n_i = Q_M_n, Qbar_n_i = Qbar_n, A = A,
	       gn_a_star_i = gn[[1]], gn_a_i = gn[[2]], 
	    FUN = function(Q_M_n_i, Qbar_n_i, A, gn_a_i, gn_a_star_i){
	    	Q_M1_a_frame <- data.frame(M1 = unique(M1), 
	    	                            Q_M1_a = Q_M_n_i[[2]][[2]])
	    	Q_M1_a_star_frame <- data.frame(M1 = unique(M1), 
	    	                            Q_M1_a_star = Q_M_n_i[[1]][[2]])	    	
	    	Q_M2_a_star_frame <- data.frame(M2 = unique(M2), 
	    	                            Q_M2_a_star = Q_M_n_i[[1]][[3]])	    	
	    	Q_M2_a_frame <- data.frame(M2 = unique(M2), 
	    	                            Q_M2_a = Q_M_n_i[[2]][[3]])
	    	Q_M1M2_a_frame <- data.frame(all_mediator_values, 
	    	                            Q_M1M2_a = Q_M_n_i[[2]][[1]])	    	
	    	Q_M1M2_a_star_frame <- data.frame(all_mediator_values, 
	    	                            Q_M1M2_a_star = Q_M_n_i[[1]][[1]])
	    	Q_M1M2_a_frame$id  <- seq_len(nrow(Q_M1M2_a_frame))

	    	nuisance_frame <- Reduce("reduce_merge", list(Q_M1M2_a_frame, Q_M1M2_a_star_frame,
	    	                                              Q_M2_a_star_frame, Q_M1_a_star_frame, 
	    	                                              Q_M1_a_frame, Q_M2_a_frame))
	    	nuisance_frame <- nuisance_frame[order(nuisance_frame$id), ]
	    	nuisance_frame$logit_Qbar_a_star <- SuperLearner::trimLogit((Qbar_n_i$Qbar_a_0[[1]] - min_Y)/(max_Y - min_Y))
	    	nuisance_frame$logit_Qbar_a <- SuperLearner::trimLogit((Qbar_n_i$Qbar_a_0[[2]] - min_Y)/(max_Y - min_Y))
			nuisance_frame$Qbar_a_star <- Qbar_n_i$Qbar_a_0[[1]] 
	    	nuisance_frame$Qbar_a <- Qbar_n_i$Qbar_a_0[[2]] 
	    	nuisance_frame$gn_a_star <- gn_a_star_i
	    	nuisance_frame$gn_a <- gn_a_i
	    	# browser()
			n_obs <- length(nuisance_frame$gn_a)
			H_a_matrix <- matrix(NA, ncol = num_covariates, nrow = n_obs)
			H_a_star_matrix <- matrix(0, ncol = num_covariates, nrow = n_obs)
			ct <- 0
			if("indirectM1" %in% which_effects){
				ct <- ct + 1
				H_a_matrix[,ct] <- with(nuisance_frame, 1 / gn_a * Q_M1_a * Q_M2_a_star / Q_M1M2_a)
				ct <- ct + 1
				H_a_matrix[,ct] <- with(nuisance_frame, - 1 / gn_a * Q_M1_a_star * Q_M2_a_star / Q_M1M2_a)
			}
			if("indirectM2" %in% which_effects){
				ct <- ct + 1
				H_a_matrix[,ct] <- with(nuisance_frame, 1 / gn_a * Q_M2_a * Q_M1_a / Q_M1M2_a)
				if(!("indirectM1" %in% which_effects)){
					ct <- ct + 1
					H_a_matrix[,ct] <- with(nuisance_frame, - 1 / gn_a * Q_M2_a_star * Q_M1_a / Q_M1M2_a)
				}
			}
			if("direct" %in% which_effects){
				ct <- ct + 1
				H_a_matrix[,ct] <- with(nuisance_frame, 1 / gn_a * Q_M1M2_a_star / Q_M1M2_a)
				ct <- ct + 1
				H_a_matrix[,ct] <- 0
				H_a_star_matrix[,dim(H_a_star_matrix[2])] <- with(nuisance_frame, 1 / gn_a_star)
			}
   	
	    	Qbar_tmle_a <- with(nuisance_frame, 
                plogis(logit_Qbar_a + H_a_matrix %*% epsilon)*(max_Y - min_Y) + min_Y
           	)

           	Qbar_tmle_a_star <- with(nuisance_frame, 
                plogis(logit_Qbar_a_star + H_a_star_matrix %*% epsilon)*(max_Y - min_Y) + min_Y
           	)

           	out <- Qbar_n_i
           	out$Qbar_a_0[[1]] <- as.numeric(Qbar_tmle_a_star)
           	out$Qbar_a_0[[2]] <- as.numeric(Qbar_tmle_a)
           	return(out)
		}, SIMPLIFY = FALSE)
	return(Qbar_n_tmle)
}


#' Target the iterated regressions
#' 
#' @inheritParams intermed
#' @inheritParams target_Qbar
#' 
#' @param Qbar Outcome regression estimates
#' @param Qbarbar Iterated mean estimates
#' @return A list containing all outcome regression evaluations needed
#' for downstream calculations.

target_Qbarbar <- function(Qbar, Qbarbar, Y, A, a, a_star, gn, M1, M2,
                           all_mediator_values, max_iter, 
                           tol = 1/(sqrt(length(Y)) * log(length(Y)))){
	
	# have confirmed this is correct
	M1_times_M2_star_a <- target_Qbarbar_M1_times_M2_star_a(
		Qbarbar = Qbarbar, Y = Y, A = A, a = a, a_star = a_star, tol = tol, gn = gn,
		max_iter = max_iter
    )

	# have confirmed this is correct
    M1_times_M2_a <- target_Qbarbar_M1_times_M2_a(
		Qbarbar = Qbarbar, Y = Y, A = A, a = a, a_star = a_star, tol = tol, gn = gn,
		max_iter = max_iter
    )
    # have confirmed this is correct
    M1_star_times_M2_star_a <- target_Qbarbar_M1_star_times_M2_star_a(
		Qbarbar = Qbarbar, Y = Y, A = A, a = a, a_star = a_star, tol = tol, gn = gn
    )

	conditional_direct_effect <- target_conditional_direct_effect(
        Qbarbar = Qbarbar, Qbar = Qbar, Y = Y, A = A, a = a, a_star = a_star, gn = gn, M1 = M1, M2 = M2,
        all_mediator_values = all_mediator_values
    )
    total_effect <- target_conditional_total_effect(
        Qbarbar = Qbarbar, Qbar = Qbar, Y = Y, A = A, a = a, a_star = a_star, gn = gn
    )

    Qbarbar$M1_times_M2_star_a <- M1_times_M2_star_a
    Qbarbar$M1_times_M2_a <- M1_times_M2_a
    Qbarbar$M1_star_times_M2_star_a <- M1_star_times_M2_star_a
    Qbarbar$conditional_direct_effect <- conditional_direct_effect
    Qbarbar$M1_star_M2_star_a_star <- total_effect[[1]]
    Qbarbar$M1_M2_a <- total_effect[[2]]

    return(Qbarbar)
}



#' Function for targeting Qbarbar_M1_times_M2_star_a
#' 
#' There are two interesting features of this 
#' targeting problem. First, we see that the nuisance parameter Qbarbar_M1_times_M2_star_a
#' can be viewed in two ways: (1) the conditional mean of Qbarbar_M1_a given C with respect
#' to the marginal of M_2 given A = a_star, C; (2) the conditional mean of Qbarbar_M2_star_a given C 
#' with respect to the marginal of M_1 given A = a, C. The natural inclination then 
#' is to use a sum loss function. It seems that we're able to do that here.
#' However, to generate the proper score, we would need to 
#' consider a submodel for the conditional mean of Qbarbar_M1_a/Qbarbar_M2_star_a given A and C; 
#' since, the inverse probability of treatment weight is a function of A. We also cannot
#' include the IPTW as part of the loss function since we need one of the weights to be 
#' negative. So, we resort to an iterative approach, where we define a submodel and loss
#' for the conditional mean of Qbarbar_M1_a and then a loss for the conditional mean of
#' Qbarbar_M2_star_a. We iterate until the empirical mean of this portion of the
#' canonical gradient is smaller than \code{tol}. 
#' 
#' @inheritParams intermed
#' @inheritParams target_Qbar
#' @inheritParams target_Qbarbar
#' @inheritParams target_Qbarbar_M1_times_M2_star_a
#' 
#' @importFrom SuperLearner trimLogit
#' 
target_Qbarbar_M1_times_M2_star_a <- function(Qbarbar, Y, A, a, a_star, gn, 
                                        tol = 1/(sqrt(length(Y)) * log(length(Y))), 
                                        max_iter = 25, iterative = FALSE, 
                                        ...){
	# outcome of regression = { Qbarbar_M2_star_a if A = a
	# 						  { Qbarbar_M1_star_a if A = a_star
	# covariate for regression = { I(A = a)/(gn[[2]]) if A = a
	# covariate for regression = { I(A = a_star)/(gn[[1]]) if A = a_star
	# offset for regression = logit(Qbarbar_M1_m2_star_a)
	n <- length(Qbarbar[[1]])
	min_Y <- min(Y); max_Y <- max(Y)
	scaled_outcome <- rep(0, n)
	scaled_outcome[A == a] <- (Qbarbar$M2_star_a[A == a] - min_Y)/(max_Y - min_Y)
	scaled_outcome[A == a_star] <- (Qbarbar$M1_a[A == a_star] - min_Y)/(max_Y - min_Y)
	
	H_A <- rep(0, n)
	H_A[A == a] <- 1 / gn[[2]][A == a]
	H_A[A == a_star] <- 1 / gn[[1]][A == a_star]

	if(iterative){
		Qbarbar_M1_times_M2_star_a_tmle <- Qbarbar$M1_times_M2_star_a
		Pn_D <- Inf
		iter <- 0
		while(abs(Pn_D) > tol & iter <= max_iter){
			iter <- iter + 1
			# target via submodel through Qbarbar, 
			scaled_offset <- rep(0, n)
			scaled_offset <- (Qbarbar_M1_times_M2_star_a_tmle - min_Y)/(max_Y - min_Y)

			fluc_data <- data.frame(
			  scaled_outcome = scaled_outcome,
			  H_A = H_A,
			  scaled_offset = SuperLearner::trimLogit(scaled_offset)
		    )[A == a, ]

		    fluc_mod <- suppressWarnings(glm(scaled_outcome ~ H_A -1 + offset(scaled_offset),
		                    family = stats::binomial(), data = fluc_data, start = 0))

			pred_data_a <- data.frame(
			  scaled_outcome = scaled_outcome,
			  H_A = 1 / gn[[2]],
			  scaled_offset = SuperLearner::trimLogit(scaled_offset)	 
	        )
		    Qbarbar_M1_times_M2_star_a_tmle <- predict(fluc_mod, 
		                                         newdata = pred_data_a,
		                                         type = "response")*(max_Y - min_Y) + min_Y

		    # now for a_star folks
		    scaled_offset <- rep(0, n)
			scaled_offset <- (Qbarbar_M1_times_M2_star_a_tmle - min_Y)/(max_Y - min_Y)

			fluc_data <- data.frame(
			  scaled_outcome = scaled_outcome,
			  H_A = H_A,
			  scaled_offset = SuperLearner::trimLogit(scaled_offset)	                  	 
		    )[A == a_star, ]

		    fluc_mod <- suppressWarnings(glm(scaled_outcome ~ H_A -1 + offset(scaled_offset),
		                    family = binomial(), data = fluc_data, start = 0))

			pred_data_a_star <- data.frame(
			  scaled_outcome = scaled_outcome,
			  H_A = 1 / gn[[1]],
			  scaled_offset = SuperLearner::trimLogit(scaled_offset)	 
	        )
		    Qbarbar_M1_times_M2_star_a_tmle <- predict(fluc_mod, 
		                                         newdata = pred_data_a_star,
		                                         type = "response")*(max_Y - min_Y) + min_Y

		    D_i <- (2*as.numeric(A == a) - as.numeric(A %in% c(a,a_star))) / ifelse(A == a, gn[[2]], gn[[1]]) * 
		    			(ifelse(A == a, Qbarbar$M2_star_a, Qbarbar$M1_a) - Qbarbar_M1_times_M2_star_a_tmle)
			Pn_D <- mean(D_i)
			cat(mean(Pn_D), "\n")
		}
	}else{
		scaled_offset <- (Qbarbar$M1_times_M2_star_a - min_Y)/(max_Y - min_Y)

		fluc_data <- data.frame(
			  scaled_outcome = scaled_outcome,
			  H_A = H_A,
			  scaled_offset = SuperLearner::trimLogit(scaled_offset)
	    )

	    fluc_mod <- suppressWarnings(glm(scaled_outcome ~ offset(scaled_offset),
                    family = stats::binomial(), weights = fluc_data$H_A, 
                    data = fluc_data, start = 0))
	    Qbarbar_M1_times_M2_star_a_tmle <- fitted(fluc_mod)*(max_Y - min_Y) + min_Y
	}
	return(Qbarbar_M1_times_M2_star_a_tmle)
}

#' Function for targeting Qbarbar_M1_times_M2_a
#' 
#' There are two interesting features of this 
#' targeting problem. First, we see that the nuisance parameter Qbarbar_M1_times_M2_a
#' can be viewed in two ways: (1) the conditional mean of Qbarbar_M1_a given C with respect
#' to the marginal of M_2 given A = a, C; (2) the conditional mean of Qbarbar_M2_a given C 
#' with respect to the marginal of M_1 given A = a, C. The natural inclination then 
#' is to use a sum loss function, which it seems we can do here. 
#' 
#' @inheritParams intermed
#' @inheritParams target_Qbar
#' @inheritParams target_Qbarbar
#' @inheritParams target_Qbarbar_M1_times_M2_star_a
#' @inheritParams intermed
#' @inheritParams target_Qbar
#' @inheritParams target_Qbarbar
#' @importFrom SuperLearner trimLogit
#' 
target_Qbarbar_M1_times_M2_a <- function(Qbarbar, Y, A, a, a_star, gn, 
                                        tol = 1/(sqrt(length(Y)) * log(length(Y))),
                                        max_iter = 25,  iterative = FALSE,
                                        ...){
	# outcome of regression = { Qbarbar_M2_a if A = a_star
	# 						  { Qbarbar_M1_a if A = a
	# covariate for regression = { I(A = a)/(gn[[2]]) if A = a
	# covariate for regression = { I(A = a_star)/(gn[[1]]) if A = a_star
	# offset for regression = logit(Qbarbar_M1_times_M2_star_a)
	n <- length(Qbarbar[[1]])
	min_Y <- min(Y); max_Y <- max(Y)
	scaled_outcome <- c((Qbarbar$M1_a[A == a] - min_Y)/(max_Y - min_Y),
	                    (Qbarbar$M2_a[A == a] - min_Y)/(max_Y - min_Y))
	# scaled_outcome[A == a] <- (Qbarbar$M1_a[A == a] - min_Y)/(max_Y - min_Y)
	# scaled_outcome[A == a_star] <- (Qbarbar$M2_a[A == a_star] - min_Y)/(max_Y - min_Y)
	
	# H_A <- rep(0, n)
	H_A <- rep(1 / gn[[2]][A == a], 2)
	# H_A[A == a_star] <- 1 / gn[[1]][A == a_star]

	if(iterative){
		Qbarbar_M1_times_M2_a_tmle <- Qbarbar$M1_star_times_M2_a
		Pn_D <- Inf
		iter <- 0
		while(abs(Pn_D) > tol & iter <= max_iter){
			iter <- iter + 1
			# target via submodel through Qbarbar, 
			scaled_offset <- rep(0, n)
			scaled_offset <- (Qbarbar_M1_times_M2_a_tmle - min_Y)/(max_Y - min_Y)

			fluc_data <- data.frame(
			  scaled_outcome = scaled_outcome,
			  H_A = H_A,
			  scaled_offset = SuperLearner::trimLogit(scaled_offset)
		    )[A == a, ]

		    fluc_mod <- suppressWarnings(glm(scaled_outcome ~ H_A -1 + offset(scaled_offset),
		                    family = stats::binomial(), data = fluc_data, start = 0))

			pred_data_a <- data.frame(
			  scaled_outcome = scaled_outcome,
			  H_A = 1 / gn[[2]],
			  scaled_offset = SuperLearner::trimLogit(scaled_offset)	 
	        )
		    Qbarbar_M1_times_M2_a_tmle <- predict(fluc_mod, 
		                                         newdata = pred_data_a,
		                                         type = "response")*(max_Y - min_Y) + min_Y

		    # now for a_star folks
			scaled_offset <- (Qbarbar_M1_times_M2_a_tmle - min_Y)/(max_Y - min_Y)

			fluc_data <- data.frame(
			  scaled_outcome = scaled_outcome,
			  H_A = H_A,
			  scaled_offset = SuperLearner::trimLogit(scaled_offset)	                  	 
		    )[A == a_star, ]

		    fluc_mod <- suppressWarnings(glm(scaled_outcome ~ H_A -1 + offset(scaled_offset),
		                    family = binomial(), data = fluc_data, start = 0))

			pred_data_a_star <- data.frame(
			  scaled_outcome = scaled_outcome,
			  H_A = 1 / gn[[1]],
			  scaled_offset = SuperLearner::trimLogit(scaled_offset)	 
	        )
		    Qbarbar_M1_times_M2_a_tmle <- predict(fluc_mod, 
		                                         newdata = pred_data_a_star,
		                                         type = "response")*(max_Y - min_Y) + min_Y

		    D_i <- (2*as.numeric(A == a) - as.numeric(A %in% c(a,a_star))) / ifelse(A == a, gn[[2]], gn[[1]]) * 
		    			(ifelse(A == a, Qbarbar$M1_star_a, Qbarbar$M2_a) - Qbarbar_M1_times_M2_a_tmle)
			Pn_D <- mean(D_i)
			cat(mean(Pn_D), "\n")
		}
	}else{
		scaled_offset <- rep((Qbarbar$M1_times_M2_a[A == a] - min_Y)/(max_Y - min_Y), 2)
		fluc_data <- data.frame(
		  scaled_outcome = scaled_outcome,
		  H_A = H_A,
		  scaled_offset = SuperLearner::trimLogit(scaled_offset)
	    )
	    fluc_mod <- suppressWarnings(glm(scaled_outcome ~ offset(scaled_offset),
                    family = stats::binomial(), weights = fluc_data$H_A, 
                    data = fluc_data, start = 0))

	    Qbarbar_M1_times_M2_a_tmle <- predict(fluc_mod, type = "response", newdata = data.frame(
            scaled_outcome = Qbarbar$M1_a, # doesn't matter
            H_A = 1 / gn[[2]],
            scaled_offset = SuperLearner::trimLogit((Qbarbar$M1_times_M2_a - min_Y)/(max_Y - min_Y))                                                                        
     	))*(max_Y - min_Y) + min_Y
	}
	return(Qbarbar_M1_times_M2_a_tmle)
}

#' Function for targeting Qbarbar_M1_star_times_M2_star_a
#' 
#' There are two interesting features of this 
#' targeting problem. First, we see that the nuisance parameter Qbarbar_M1_star_times_M2_star_a
#' can be viewed in two ways: (1) the conditional mean of Qbarbar_M1_star_a given C with respect
#' to the marginal of M_2 given A = a_star, C; (2) the conditional mean of Qbarbar_M2_star_a given C 
#' with respect to the marginal of M_1 given A = a_star, C. The natural inclination then 
#' is to use a sum loss function. Here it looks like we actually can use a sum 
#' loss approach, so long as the IPTW are incorporated into the loss function. We make
#' two copies of each observation with A = a_star; assign Qbarbar_M1_star_a as outcome in
#' half and Qbarbar_M2_star_a in the other half; then do one-shot targeting
#' @inheritParams intermed
#' @inheritParams target_Qbar
#' @inheritParams target_Qbarbar
#' @inheritParams target_Qbarbar_M1_times_M2_star_a
#' @importFrom SuperLearner trimLogit
#' 
target_Qbarbar_M1_star_times_M2_star_a <- function(Qbarbar, Y, A, a, a_star, gn,
                                        tol = 1/(sqrt(length(Y)) * log(length(Y))), 
                                        ...){
	# outcome of regression = { Qbarbar_M2_star_a if A = a
	# 						  { Qbarbar_M1_star_a if A = a_star
	# covariate for regression = { I(A = a)/(gn[[2]]) if A = a
	# covariate for regression = { I(A = a_star)/(gn[[1]]) if A = a_star
	# offset for regression = logit(Qbarbar_M1_m2_star_a)
	n <- length(Qbarbar[[1]])
	min_Y <- min(Y); max_Y <- max(Y)
	scaled_outcome <- c((Qbarbar$M1_star_a[A == a_star] - min_Y)/(max_Y - min_Y),
	                    (Qbarbar$M2_star_a[A == a_star] - min_Y)/(max_Y - min_Y))
	H_A <- rep(1 / gn[[1]][A == a_star], 2)
	
	scaled_offset <- rep((Qbarbar$M1_star_times_M2_star_a[A == a_star] - min_Y)/(max_Y - min_Y),2)
	# scaled_offset <- rep((Qbarbar_M1_star_times_M2_star_a[A == a_star] - min_Y)/(max_Y - min_Y),2)

	fluc_data <- data.frame(
	  scaled_outcome = scaled_outcome,
	  H_A = H_A,
	  scaled_offset = SuperLearner::trimLogit(scaled_offset)
    )

	fluc_mod <- suppressWarnings(glm(scaled_outcome ~ offset(scaled_offset),
                    family = stats::binomial(), weights = fluc_data$H_A, 
                    data = fluc_data, start = 0))

	scaled_offset <- (Qbarbar$M1_star_times_M2_star_a - min_Y)/(max_Y - min_Y)

	pred_data_a <- data.frame(
	  scaled_outcome = Qbarbar$M1_star_a, # doesn't matter
	  H_A = 1 / gn[[1]],
	  scaled_offset = SuperLearner::trimLogit(scaled_offset)	 
    )
    Qbarbar_M1_star_times_M2_star_a_tmle <- predict(fluc_mod, 
                                         newdata = pred_data_a,
                                         type = "response")*(max_Y - min_Y) + min_Y
	    
    D_i <- (as.numeric(A == a_star)/gn[[1]] * (Qbarbar$M1_star_a - Qbarbar_M1_star_times_M2_star_a_tmle)
    + as.numeric(A == a_star)/gn[[1]] * (Qbarbar$M2_star_a - Qbarbar_M1_star_times_M2_star_a_tmle))
	Pn_D <- mean(D_i)
	
	return(Qbarbar_M1_star_times_M2_star_a_tmle)
}



#' The naming convention is a bit different here, because we're actually
#' able to go after the conditional direct effect, well, directly. That is, 
#' we can define a loss function whose minimizer defines the conditional mean of the
#' difference between Qbar_a
#' and Qbar_a_star with respect to the joint distribution of M1 and M2 given C and A = a_star. 
#' We can then define a submodel through this conditional mean difference (which is 
#' exactly the conditional direct effect) and target this quantity directly. 
#' @inheritParams intermed
#' @inheritParams target_Qbar
#' @inheritParams target_Qbarbar
#' @inheritParams target_Qbarbar_M1_times_M2_star_a

target_conditional_direct_effect <- function(Qbarbar, all_mediator_values, gn,
                                             Qbar, # should come in as TMLE
                                             Y, A, a, a_star, M1, M2, 
                                             target_conditional = TRUE,
                                             epsilon_threshold = 5, 
                                             bound_pred = FALSE, 
                                             universal = TRUE, 
                                             deps = 1e-5, 
                                             max_iter = 10000, ...){
	n <- length(Qbarbar[[1]])
	
	Qbar_M1M2_a <- unlist(mapply(Qbar_n_i = Qbar, M1_i = M1, M2_i = M2,
	                             FUN = extract_Qbar_obs,
	                             MoreArgs = list( 
	                         	  all_mediator_values = all_mediator_values, 
                         		  a_val = "a"
             	                ), SIMPLIFY = FALSE), use.names = FALSE)	
	Qbar_M1M2_a_star <- unlist(mapply(Qbar_n_i = Qbar, M1_i = M1, M2_i = M2,
	                             FUN = extract_Qbar_obs,
	                             MoreArgs = list( 
	                         	  all_mediator_values = all_mediator_values, 
                         		  a_val = "a_star"
             	                ), SIMPLIFY = FALSE), use.names = FALSE)

	# enforce that targeted estimates don't stray out of range of initial estimates

	if(target_conditional){
		min_Y <- min(Y); max_Y <- max(Y)
		ell <- min_Y - max_Y; u <- max_Y - min_Y

		scaled_outcome <- (Qbar_M1M2_a - Qbar_M1M2_a_star - ell)/(u - ell)

		H_A <- as.numeric(A == a_star) / gn[[1]]
		
		unscaled_offset <- Qbarbar$M1_star_M2_star_a - Qbarbar$M1_star_M2_star_a_star
		scaled_offset_out <- (unscaled_offset - ell)/(u - ell)

		fluc_data <- data.frame(
		  scaled_outcome = scaled_outcome,
		  H_A = H_A,
		  scaled_offset = SuperLearner::trimLogit(scaled_offset_out)
	    )

		if(!bound_pred){
			if(!universal){
				pred_data_a <- data.frame(
					  scaled_outcome = 0, # doesn't matter
					  H_A = 1 / gn[[1]],
					  scaled_offset = SuperLearner::trimLogit(scaled_offset_out)	 
			    )
				fluc_mod <- suppressWarnings(glm(scaled_outcome ~ -1 + offset(scaled_offset) + H_A,
			                    family = stats::binomial(), data = fluc_data, start = 0))

			    cond_direct_effect <- predict(fluc_mod, 
			                                  newdata = pred_data_a,
			                                  type = "response")*(u - ell) + ell
			}else{
				init_deriv <- mean(as.numeric(A == a_star)/gn[[1]] * (Qbar_M1M2_a - Qbar_M1M2_a_star - (Qbarbar$M1_star_M2_star_a - Qbarbar$M1_star_M2_star_a_star)))
			    if(init_deriv < 0){
			    	deps <- -deps
			    }
			    nllik <- function(scaled_conditional_direct_effect, scaled_outcome, A, a_star){
			        -sum(as.numeric(A == a_star) * 
			             	(scaled_outcome * log(scaled_conditional_direct_effect) + 
			             	 (1 - scaled_outcome) * log(1 - scaled_conditional_direct_effect)))
    			}
    			current_loss <- nllik(scaled_offset_out, scaled_outcome, A, a_star)
    			previous_loss <- Inf
    			old_scaled_conditional_direct_effect <- scaled_offset_out
    			iter <- 0
    			all_loss <- all_deriv <- rep(NA, max_iter)
    			cur_deriv <- init_deriv
    			cond_direct_effect <- scaled_offset_out * (u - ell) + ell
    			while(previous_loss >= current_loss & iter < max_iter & 
    			      abs(cur_deriv) > 1/(sqrt(n) * log(n))){
    				iter <- iter + 1
    				previous_loss <- current_loss
    				new_scaled_conditional_direct_effect <- plogis(SuperLearner::trimLogit(old_scaled_conditional_direct_effect) + deps * 1 / gn[[1]])
    				current_loss <- nllik(new_scaled_conditional_direct_effect, scaled_outcome, A, a_star)
    				all_loss[iter] <- current_loss
    				old_scaled_conditional_direct_effect <- new_scaled_conditional_direct_effect
	    			cond_direct_effect <- new_scaled_conditional_direct_effect * (u - ell) + ell
					cur_deriv <- mean(as.numeric(A == a_star)/gn[[1]] * (Qbar_M1M2_a - Qbar_M1M2_a_star - cond_direct_effect))
					all_deriv[iter] <- cur_deriv
    			}
				final_deriv <- cur_deriv
			}
		}else{
			stop("nothing here")
		}
	}else{
		# target QbarbarM1_star_m2_star_a
		if(!all(Y %in% c(0,1))){
			ell <- min(Qbar_M1M2_a); u <- max(Qbar_M1M2_a)
		}else{
			ell <- 0; u <- 1
		}
		scaled_outcome <- (Qbar_M1M2_a - ell)/(u - ell)

		H_A <- as.numeric(A == a_star) / gn[[1]]
		
		scaled_offset_out <- (Qbarbar$M1_star_M2_star_a  - ell)/(u - ell)

		fluc_data <- data.frame(
		  scaled_outcome = scaled_outcome,
		  H_A = H_A,
		  scaled_offset = SuperLearner::trimLogit(scaled_offset_out)
	    )
	    fluc_mod <- suppressWarnings(glm(scaled_outcome ~ -1 + offset(scaled_offset) + H_A,
		                    family = stats::binomial(), data = fluc_data, start = 0))

		pred_data_a <- data.frame(
		  scaled_outcome = 0, # doesn't matter
		  H_A = 1 / gn[[1]],
		  scaled_offset = SuperLearner::trimLogit(scaled_offset_out)	 
	    )
	    Qbar_M1M2_a_tmle <- predict(fluc_mod, newdata = pred_data_a, type = "response")*(u - ell) + ell

	 	# target QbarbarM1_star_m2_star_a_star
		if(!all(Y %in% c(0,1))){
			ell <- min(Qbar_M1M2_a_star); u <- max(Qbar_M1M2_a_star)
		}
		scaled_outcome <- (Qbar_M1M2_a_star - ell)/(u - ell)

		H_A <- as.numeric(A == a_star) / gn[[1]]
		
		scaled_offset_out <- (Qbarbar$M1_star_M2_star_a_star  - ell)/(u - ell)

		fluc_data <- data.frame(
		  scaled_outcome = scaled_outcome,
		  H_A = H_A,
		  scaled_offset = SuperLearner::trimLogit(scaled_offset_out)
	    )
	    fluc_mod <- suppressWarnings(glm(scaled_outcome ~ -1 + offset(scaled_offset) + H_A,
		                    family = stats::binomial(), data = fluc_data, start = 0))

		pred_data_a <- data.frame(
		  scaled_outcome = 0, # doesn't matter
		  H_A = 1 / gn[[1]],
		  scaled_offset = SuperLearner::trimLogit(scaled_offset_out)	 
	    )
	    Qbar_M1M2_a_star_tmle <- predict(fluc_mod, newdata = pred_data_a, type = "response")*(u - ell) + ell
		cond_direct_effect <- Qbar_M1M2_a_tmle - Qbar_M1M2_a_star_tmle 
	}
	return(cond_direct_effect)
}


#' The naming convention is a bit different here, because we're actually
#' able to go after the conditional total effect directly. That is, 
#' we can define a loss function whose minimizer defines the conditional mean of the
#' difference between Qbar_a
#' and Qbar_a_star with respect to the joint distribution of M1 and M2 given C and A.
#' We can then define a submodel through this conditional mean difference (which is 
#' exactly the conditional total effect) and target this quantity directly. 
#' @inheritParams intermed
#' @inheritParams target_Qbar
#' @inheritParams target_Qbarbar
#' @inheritParams target_Qbarbar_M1_times_M2_star_a

target_conditional_total_effect <- function(Qbarbar, gn, 
                                            Qbar, # should come in as TMLE
                                            Y, A, a, a_star, ...){
	n <- length(Qbarbar[[1]])
	min_Y <- min(Y); max_Y <- max(Y)
	scaled_outcome <- (Y - min_Y)/(max_Y - min_Y)

	H_A <- (2*as.numeric(A == a) - as.numeric(A %in% c(a,a_star)))/ifelse(A == a, gn[[2]], gn[[1]])
	
	unscaled_offset <- ifelse(A == a, Qbarbar$M1_M2_a, Qbarbar$M1_star_M2_star_a_star)
	scaled_offset <- (unscaled_offset - min_Y)/(max_Y - min_Y)

	fluc_data <- data.frame(
	  scaled_outcome = scaled_outcome,
	  H_A = H_A,
	  scaled_offset = SuperLearner::trimLogit(scaled_offset)
    )

	fluc_mod <- suppressWarnings(glm(scaled_outcome ~ -1 + H_A + offset(scaled_offset),
                    family = stats::binomial(), 
                    data = fluc_data, start = 0))

	unscaled_offset <- Qbarbar$M1_M2_a
	scaled_offset <- (unscaled_offset - min_Y)/(max_Y - min_Y)
	pred_data_a <- data.frame(
	  scaled_outcome = Qbarbar$M1_star_a, # doesn't matter
	  H_A = 1 / gn[[2]],
	  scaled_offset = SuperLearner::trimLogit(scaled_offset)	 
    )	
    unscaled_offset <- Qbarbar$M1_star_M2_star_a_star
	scaled_offset <- (unscaled_offset - min_Y)/(max_Y - min_Y)
    pred_data_a_star <- data.frame(
	  scaled_outcome = Qbarbar$M1_star_a, # doesn't matter
	  H_A = -1 / gn[[1]],
	  scaled_offset = SuperLearner::trimLogit(scaled_offset)	 
    )
    Qbarbar_a <- predict(fluc_mod, newdata = pred_data_a,
                         type = "response")*(max_Y - min_Y) + min_Y
	Qbarbar_a_star <- predict(fluc_mod, newdata = pred_data_a_star,
                         type = "response")*(max_Y - min_Y) + min_Y
	    
    D_i <- (2*as.numeric(A == a) - as.numeric(A %in% c(a, a_star)))/ifelse(A == a, gn[[2]], gn[[1]]) * 
    			(Y - ifelse(A == a, Qbarbar_a, Qbarbar_a_star))
	Pn_D <- mean(D_i)
	
	return(list(Qbarbar_a_star, Qbarbar_a))
}

