#' Function to evaluate the canonical gradient of 
#' the total effect at a given set of estimated nuisance parameters
#' and for each observation.
#' 
#' @param Y The outcome
#' @param A The treatment
#' @param Qbarbar A list; needs to have entries named M1_M2_a and M1_star_M2_star_a_star corresponding
#' to, respectively, the outcome regression under A = a, marginalized with respect to joint 
#' mediator distribution given C and A = a, and the outcome regression under A = a_star, marginalized
#' with respect to the joint mediator distribution given C and A = a_star. 
#' @param gn List of mediator values
#' @param a Treatment value
#' @param a_star Referent treatment value


evaluate_eif_total <- function(Y, A, Qbarbar, gn, a, a_star, ...){
	n <- length(Qbarbar[[1]])
	min_Y <- min(Y); max_Y <- max(Y)

	eif <- as.numeric(A == a)/gn[[2]] * (Y - Qbarbar$M1_M2_a) -
			 as.numeric(A == a_star)/gn[[1]] * (Y - Qbarbar$M1_star_M2_star_a_star) + 
			 	 Qbarbar$M1_M2_a - Qbarbar$M1_star_M2_star_a_star - 
			 	 	mean(Qbarbar$M1_M2_a - Qbarbar$M1_star_M2_star_a_star)
	return(eif)
}

#' Helper function to evaluate the total effect
evaluate_total_effect <- function(Qbarbar){
	return(mean(Qbarbar$M1_M2_a - Qbarbar$M1_star_M2_star_a_star))
}

#' Function to evaluate the canonical gradient of 
#' the direct effect at a given set of estimated nuisance parameters
#' and for each observation.
#' 
#' @param Y The outcome
#' @param A The treatment
#' @param Qbarbar A list; needs to have entries named M1_M2_a and M1_star_M2_star_a_star corresponding
#' to, respectively, the outcome regression under A = a, marginalized with respect to joint 
#' mediator distribution given C and A = a, and the outcome regression under A = a_star, marginalized
#' with respect to the joint mediator distribution given C and A = a_star. 
#' @param Qbar Outcome regression list
#' @param Q_M Mediator distribution list
#' @param gn List of propensity scores
#' @param a Treatment value
#' @param a_star Referent treatment value
#' @param use_conditional Boolean. Should the conditional_total_effect be used (generally used for when evaluating
#' EIF of targeted nuisance parameters, since there we target directly the difference)

evaluate_eif_direct <- function(Y, A, M1, M2, Qbar, Q_M, Qbarbar, gn, a, a_star, use_conditional, all_mediator_values, 
                                ...){
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
	# get an n-length vector of joint dist. at M1_i, M2_i and C_i under a
	Q_M1M2_a <- unlist(mapply(Q_M_n_i = Q_M, M1_i = M1, M2_i = M2, FUN = extract_joint,
	                   MoreArgs = list(a_val = "a",
	                                   all_mediator_values = all_mediator_values),
	                   SIMPLIFY = FALSE), use.names = FALSE)	
    # get an n-length vector of joint dist. at M1_i, M2_i and C_i under a_star
	Q_M1M2_a_star <- unlist(mapply(Q_M_n_i = Q_M, M1_i = M1, M2_i = M2, FUN = extract_joint,
	                   MoreArgs = list(a_val = "a_star",
	                                   all_mediator_values = all_mediator_values),
	                   SIMPLIFY = FALSE), use.names = FALSE)

	if(!use_conditional){
		conditional_effect <- (Qbarbar$M1_star_M2_star_a - Qbarbar$M1_star_M2_star_a_star)
	}else{
		conditional_effect <- Qbarbar$conditional_direct_effect
	}
	# unscaled_offset <- Qbarbar$M1_star_M2_star_a - Qbarbar$M1_star_M2_star_a_star
	eif <- as.numeric(A == a)/gn[[2]] * Q_M1M2_a_star / Q_M1M2_a * (Y - Qbar_M1M2_a) -
			 as.numeric(A == a_star)/gn[[1]] * (Y - Qbar_M1M2_a_star) +
			 	as.numeric(A == a_star)/gn[[1]] * (Qbar_M1M2_a - Qbar_M1M2_a_star - 
			 	                                   	conditional_effect) + 
			 		conditional_effect - mean(conditional_effect)
	return(eif)
}

#' Helper function to evaluate the direct effect
evaluate_direct_effect <- function(Qbarbar, use_conditional = TRUE){
	if(use_conditional){
		return(mean(Qbarbar$conditional_direct_effect))
	}else{
		return(mean(Qbarbar$M1_star_M2_star_a - Qbarbar$M1_star_M2_star_a_star))
	}
}

#' Function to evaluate the canonical gradient of 
#' the indirect effect through M1 at a given set of estimated nuisance parameters
#' and for each observation.
#' 
#' @param Y The outcome
#' @param A The treatment
#' @param M1 First mediator
#' @param M2 Second mediator
#' @param Qbarbar A list; needs to have entries named M1_M2_a and M1_star_M2_star_a_star corresponding
#' to, respectively, the outcome regression under A = a, marginalized with respect to joint 
#' mediator distribution given C and A = a, and the outcome regression under A = a_star, marginalized
#' with respect to the joint mediator distribution given C and A = a_star. 
#' @param gn List of mediator values
#' @param Q_M Mediator distribution list
#' @param Qbar Outcome regression list 
#' @param a Treatment value
#' @param a_star Referent treatment value
#' 

evaluate_eif_indirect_M1 <- function(Y, A, M1, M2, Qbar, Q_M, Qbarbar, gn, a, a_star, all_mediator_values, ...){
	n <- length(Qbarbar[[1]])
	unique_M1_values <- unique(M1)
	unique_M2_values <- unique(M2)
	Qbar_M1M2_a <- unlist(mapply(Qbar_n_i = Qbar, M1_i = M1, M2_i = M2,
                             FUN = extract_Qbar_obs,
                             MoreArgs = list( 
                         	  all_mediator_values = all_mediator_values, 
                     		  a_val = "a"
         	                ), SIMPLIFY = FALSE), use.names = FALSE)	
	# get an n-length vector of M1 marginal at M1_i and C_i under a_star and a
	Q_M1_a_star <- unlist(mapply(Q_M_n_i = Q_M, M_i = M1, FUN = extract_marginal,
	                    MoreArgs = list(a_val = "a_star", 
	                                    unique_M_values = unique_M1_values,
	                                    mediator = "M1"),
	                    SIMPLIFY = FALSE), use.names = FALSE)
	Q_M1_a <- unlist(mapply(Q_M_n_i = Q_M, M_i = M1, FUN = extract_marginal,
	                    MoreArgs = list(a_val = "a", 
	                                    unique_M_values = unique_M1_values,
	                                    mediator = "M1"),
	                    SIMPLIFY = FALSE), use.names = FALSE)
	Q_M2_a_star <- unlist(mapply(Q_M_n_i = Q_M, M_i = M2, FUN = extract_marginal,
	                    MoreArgs = list(a_val = "a_star", 
	                                    unique_M_values = unique_M2_values,
	                                    mediator = "M2"),
	                    SIMPLIFY = FALSE), use.names = FALSE)	
	
	# get an n-length vector of joint dist. at M1_i, M2_i and C_i under a
	Q_M1M2_a <- unlist(mapply(Q_M_n_i = Q_M, M1_i = M1, M2_i = M2, FUN = extract_joint,
	                   MoreArgs = list(a_val = "a",
	                                   all_mediator_values = all_mediator_values),
	                   SIMPLIFY = FALSE), use.names = FALSE)	
    
	eif <- as.numeric(A == a)/gn[[2]] * (Q_M1_a - Q_M1_a_star) * Q_M2_a_star / Q_M1M2_a * (Y - Qbar_M1M2_a) + 
			as.numeric(A == a)/gn[[2]] * (Qbarbar$M2_star_a - Qbarbar$M1_times_M2_star_a) - 
			as.numeric(A == a_star)/gn[[1]] * (Qbarbar$M2_star_a - Qbarbar$M1_star_times_M2_star_a) + 
			as.numeric(A == a_star)/gn[[1]] * ((Qbarbar$M1_a - Qbarbar$M1_star_a) - (Qbarbar$M1_times_M2_star_a - Qbarbar$M1_star_times_M2_star_a)) + 
				Qbarbar$M1_times_M2_star_a - Qbarbar$M1_star_times_M2_star_a - 
					mean(Qbarbar$M1_times_M2_star_a - Qbarbar$M1_star_times_M2_star_a)
	return(eif)
}

evaluate_indirect_effect_M1 <- function(Qbarbar){
	return(mean(Qbarbar$M1_times_M2_star_a - Qbarbar$M1_star_times_M2_star_a))
}

#' Function to evaluate the canonical gradient of 
#' the indirect effect through M2 at a given set of estimated nuisance parameters
#' and for each observation.
#' 
#' @param Y The outcome
#' @param A The treatment
#' @param M1 First mediator
#' @param M2 Second mediator
#' @param Qbarbar A list; needs to have entries named M1_M2_a and M1_star_M2_star_a_star corresponding
#' to, respectively, the outcome regression under A = a, marginalized with respect to joint 
#' mediator distribution given C and A = a, and the outcome regression under A = a_star, marginalized
#' with respect to the joint mediator distribution given C and A = a_star. 
#' @param gn List of mediator values
#' @param Q_M Mediator distribution list
#' @param Qbar Outcome regression list 
#' @param a Treatment value
#' @param a_star Referent treatment value
#' 

evaluate_eif_indirect_M2 <- function(Y, A, M1, M2, Qbar, Q_M, Qbarbar, gn, a, a_star, all_mediator_values, ...){
	n <- length(Qbarbar[[1]])
	unique_M2_values <- unique(M2)
	unique_M1_values <- unique(M1)
	Qbar_M1M2_a <- unlist(mapply(Qbar_n_i = Qbar, M1_i = M1, M2_i = M2,
                             FUN = extract_Qbar_obs,
                             MoreArgs = list( 
                         	  all_mediator_values = all_mediator_values, 
                     		  a_val = "a"
         	                ), SIMPLIFY = FALSE), use.names = FALSE)	
	Q_M2_a <- unlist(mapply(Q_M_n_i = Q_M, M_i = M2, FUN = extract_marginal,
	                    MoreArgs = list(a_val = "a", 
	                                    unique_M_values = unique_M2_values,
	                                    mediator = "M2"),
	                    SIMPLIFY = FALSE), use.names = FALSE)
	Q_M2_a_star <- unlist(mapply(Q_M_n_i = Q_M, M_i = M2, FUN = extract_marginal,
	                    MoreArgs = list(a_val = "a_star", 
	                                    unique_M_values = unique_M2_values,
	                                    mediator = "M2"),
	                    SIMPLIFY = FALSE), use.names = FALSE)		
	Q_M1_a <- unlist(mapply(Q_M_n_i = Q_M, M_i = M1, FUN = extract_marginal,
	                    MoreArgs = list(a_val = "a", 
	                                    unique_M_values = unique_M1_values,
	                                    mediator = "M1"),
	                    SIMPLIFY = FALSE), use.names = FALSE)

	# get an n-length vector of joint dist. at M1_i, M2_i and C_i under a
	Q_M1M2_a <- unlist(mapply(Q_M_n_i = Q_M, M1_i = M1, M2_i = M2, FUN = extract_joint,
	                   MoreArgs = list(a_val = "a",
	                                   all_mediator_values = all_mediator_values),
	                   SIMPLIFY = FALSE), use.names = FALSE)	
    
	eif <- as.numeric(A == a)/gn[[2]] * (Q_M2_a - Q_M2_a_star) * Q_M1_a / Q_M1M2_a * (Y - Qbar_M1M2_a) + 
			as.numeric(A == a)/gn[[2]] * (Qbarbar$M1_a - Qbarbar$M1_times_M2_a) - 
			as.numeric(A == a_star)/gn[[1]] * (Qbarbar$M1_a - Qbarbar$M1_times_M2_star_a) + 
			as.numeric(A == a)/gn[[2]] * ((Qbarbar$M2_a - Qbarbar$M2_star_a) - (Qbarbar$M1_times_M2_a - Qbarbar$M1_times_M2_star_a)) + 
				Qbarbar$M1_times_M2_a - Qbarbar$M1_times_M2_star_a - 
					mean(Qbarbar$M1_times_M2_a - Qbarbar$M1_times_M2_star_a)
	return(eif)
}

evaluate_indirect_effect_M2 <- function(Qbarbar){
	return(mean(Qbarbar$M1_times_M2_a - Qbarbar$M1_times_M2_star_a))
}

evaluate_covariance_effect_M1M2 <- function(Qbarbar, use_conditional){
	total_eff <- evaluate_total_effect(Qbarbar = Qbarbar)
	direct_eff <- evaluate_direct_effect(Qbarbar = Qbarbar, use_conditional = use_conditional)
	indirect_eff_M1 <- evaluate_indirect_effect_M1(Qbarbar = Qbarbar)
	indirect_eff_M2 <- evaluate_indirect_effect_M2(Qbarbar = Qbarbar)
	return(total_eff - direct_eff - indirect_eff_M1 - indirect_eff_M2)
}

