#' AIPTW and TMLE estimates of interventional mediation effects
#' 
#' This function computes multiply robust estimates of interventional mediation
#' effects with two mediators. Two such estimates are produced, an augmented
#' inverse probability of treatment weighted (AIPTW) estimator and a targeted minimum loss
#' estimator (TMLE). The AIPTW was found to perform better in most simulations.
#' 
#' Under the hood, the function fits several nuisance regressions and combines them
#' into the final estimate. Users have several options for fitting these regressions
#' as described below. In particular, \code{SuperLearner} can be used to flexibly 
#' estimate the regressions.
#' 
#' @param Y A vector of continuous or binary outcomes.
#' @param A A vector of binary treatment assignment (assumed to be equal to 0 or
#'  1).
#' @param C A \code{data.frame} of named covariates.
#' @param M1 A \code{vector} of mediators.
#' @param M2 A \code{vector} of mediators.
#' @param all_mediator_values All combinations of M1 and M2
#' @param DeltaM Indicator of missing outcome (assumed to be equal to 0 if
#'  missing 1 if observed).
#' @param DeltaA Indicator of missing treatment (assumed to be equal to 0 if
#'  missing 1 if observed).
#' @param SL_Qbar A vector of characters or a list describing the Super Learner
#'  library to be used for the outcome regression.
#' @param SL_g A vector of characters or a list describing the Super Learner
#'  library to be used for the propensity score.
#' @param SL_Q_M A list with two named entries \code{M1} and \code{M2}, specifying
#'  super learner libraries for estimation of the pooled hazard regression that is
#'  mapped into an estimate of the mediator distributions.
#' @param a The label for the treatment. The effects estimates returned pertain
#' to estimation of interventional effects of \code{a} versus \code{a_star}.
#' @param a_star The label for the treatment. The effects estimates returned pertain
#' to estimation of interventional effects of \code{a} versus \code{a_star}.
#' @param n_SL The number of repeated super learner runs to execute for each
#' regression.
#' @param glm_Qbar A character describing a formula to be used in the call to
#'  \code{glm} for the outcome regression. The formula may include \code{colnames(C)},
#' \code{"A"} (if \code{stratify = FALSE}), \code{"M1"}, and \code{"M2"} as terms in 
#' the formula.
#' @param glm_g A character describing a formula to be used in the call to
#'  \code{glm} for the outcome regression. The formula may include \code{colnames(C)},
#'  in the formula.
#' @param glm_Q_M A list with two named entries \code{M1} and \code{M2}, specifying
#'  the regression formula for estimation of the pooled hazard regression that is
#'  mapped into an estimate of the mediator distributions. The formula for \code{M1}
#'  can include \code{"A"} (if \code{stratify = FALSE}) and \code{colnames(C)}, 
#' while the formula for \code{M2} can additionally include \code{M1}.
#' @param tolg The truncation level for the propensity score
#' @param tol The tolerance for stopping the iterative targeting procedure.
#' @param targeted_se A boolean indicating whether to return the standard error estimates
#' based on targeted nuisance parameters or the initial estimates of nuisance parameters.
#' @param return_models A boolean indicating whether to return the fitted models for 
#' each of the nuisance regressions. If \code{TRUE} then the output will include a list 
#' with named entries \code{g} (propensity score fit), \code{Qbar} (outcome regression fit),
#' and \code{Q_M} (mediator distribution fits).
#' @param Qbar_n Power users may wish to pass in their own properly formatted list of the
#' outcome regression so that
#' nuisance parameters can be fitted outside of \code{intermed}.
#' @param Q_M_n Power users may wish to pass in their own properly formatted list of the
#' mediator distributions so that nuisance parameters can be fitted outside 
#' of \code{intermed}.
#' @param gn Power users may wish to pass in their own properly formatted list of the
#' propensity score so that
#' nuisance parameters can be fitted outside of \code{intermed}.
#' @param use_future A boolean indicating whether to use the \code{future} package 
#' to parallelize computations
#' @param cvFolds Number of cross-validation folds to use if CVTMLE and CV-one step are
#' desired
#' @param max_iter The maximum number of iterations for the TMLE
#' @param verbose A boolean indicating whether to print status updates.
#' @export
#' 
#' @return An object of class \code{"intermed"}.
#' \describe{
#'  \item{\code{aiptw}}{A \code{list} of point estimates and
#'        estimated covariance matrix from the one-step estimator}
#'  \item{\code{tmle}}{A \code{list} of point estimates and
#'        estimated covariance matrix from the TMLE}
#'  \item{\code{plugin}}{Plugin estimates of the mediation effects}
#'  \item{\code{fm}}{Fitted models if \code{return_models = TRUE}, \code{NULL} otherwise}
#' }

intermed <- function(Y, C, M1, M2, A, 
                     DeltaA = as.numeric(!is.na(A)),
                     DeltaM = as.numeric(!is.na(M1)),
                     a = 1, 
                     a_star = 0,
                     SL_Qbar = NULL, 
                     SL_g = NULL, 
                     SL_Q_M = NULL, 
                     n_SL = 1,
                     glm_Qbar = NULL, 
                     glm_g = NULL,
                     glm_Q_M = NULL, 
                     tolg = 1e-2, 
                     tol = 1 / (sqrt(length(Y))*log(length(Y))),
                     targeted_se = FALSE, 
                     return_models = FALSE, 
                     cvFolds = 1, 
                     use_future = FALSE, 
                     Qbar_n = NULL, 
                     Q_M_n = NULL,
                     gn = NULL, 
                     max_iter = 50,
                     verbose = FALSE,
                     stratify = FALSE, 
                     ...){
	n <- length(Y)
	a_0 <- c(a_star, a)

	#! TO DO: don't know if this condition will work in general or not 
	stopifnot(all(which(is.na(M1)) == which(is.na(M2))))

    # make a data.frame of all mediator values
    all_mediator_values <- expand.grid(M1 = unique(M1),
                                       M2 = unique(M2))
    n_mediator_values <- nrow(all_mediator_values)

	# check for user input Q and g 
	Qbar_n_user <- !is.null(Qbar_n)
  gn_user <- !is.null(gn)
  Q_M_n_user <- !is.null(Q_M_n)

    # make validation row list for super learner
	valid_rows <- drtmle:::make_validRows(cvFolds, n = n, n_SL = n_SL)
    if (n_SL > 1) {
        valid_rows <- rep(valid_rows, n_SL)
    }

    if (is.null(gn)) {
        if (use_future) {
              gnOut <- future.apply::future_lapply(X = valid_rows, 
                FUN = estimateG, A = A, W = C, DeltaA = DeltaA, 
                DeltaY = DeltaM, tolg = tolg, verbose = verbose, 
                stratify = stratify, returnModels = return_models, 
                SL_g = SL_g, glm_g = glm_g, a_0 = a_0)
        }
        else {
            gnOut <- lapply(X = valid_rows, FUN = drtmle:::estimateG, A = A, 
                W = C, DeltaA = DeltaA, DeltaY = DeltaM, tolg = tolg, 
                verbose = verbose, stratify = stratify, returnModels = return_models, 
                SL_g = SL_g, glm_g = glm_g, a_0 = a_0)
        }
        gn <- gnOut[[1]]$est
        #reorder_list(gnOut, a_0 = a_0, valid_rows = valid_rows, n_SL = n_SL, n = n,
        #                  which_nuisance = "PS")
        gnMod <- drtmle:::extract_models(gnOut)
    }else {
        gn <- lapply(gn, function(g) {
            g[g < tolg] <- tolg
            g
        })
    }
    # right now it's a bit weird that gn is structured as a 2-length list with n elements
    # whereas the other nuisances are n-length lists. Probably not too hard to keep doing
    # this bookkeeping, but something to keep in mind. 

	# -------------------------------
	# estimate outcome regression
	# -------------------------------
	if (is.null(Qbar_n)) {
		Qfam <- if(all(Y %in% c(0,1))){ binomial() }else{ gaussian() }
		if(use_future){
		  Qbar_n_out <- future.apply::future_lapply(
		    X = valid_rows, FUN = estimate_Qbar,
		    Y = Y, A = A, C = C, M1 = M1, M2 = M2, 
            all_mediator_values = all_mediator_values,
		    DeltaA = DeltaA, DeltaM = DeltaM,
		    verbose = verbose,
		    return_models = return_models,
		    SL_Qbar = SL_Qbar, a_0 = a_0,
		    stratify = stratify,
		    glm_Qbar = glm_Qbar,
		    family = Qfam
		  )
		}else{
		  Qbar_n_out <- lapply(
		    valid_rows, FUN = estimate_Qbar,
		    Y = Y, A = A, C = C, M1 = M1, M2 = M2,
		    DeltaA = DeltaA, DeltaM = DeltaM,
        all_mediator_values = all_mediator_values,
		    verbose = verbose,
		    return_models = return_models,
		    SL_Qbar = SL_Qbar, a_0 = a_0,
		    stratify = stratify,
		    glm_Qbar = glm_Qbar,
		    family = Qfam
		  )
		}
		# re-order predictions
		Qbar_n <- reorder_list(Qbar_n_out, a_0 = a_0, valid_rows = valid_rows, 
		                   n_SL = n_SL, n = n, which_nuisance = "OR",
                           n_mediator_values = n_mediator_values)

		# obtain list of outcome regression fits
		Qbar_n_models <- drtmle:::extract_models(Qbar_n_out)
	}

    # -------------------------------
    # estimate mediator distributions
    # -------------------------------
    if (is.null(Q_M_n)) {
        if(use_future){
          Q_M_out <- future.apply::future_lapply(
            X = valid_rows, FUN = estimate_Q_M,
            A = A, C = C, M1 = M1, M2 = M2, 
            all_mediator_values = all_mediator_values,
            DeltaA = DeltaA, DeltaM = DeltaM,
            verbose = verbose,
            return_models = return_models,
            SL_Q = SL_Q_M, a_0 = a_0,
            stratify = stratify,
            glm_Q_M = glm_Q_M
          )
        }else{
          Q_M_out <- lapply(
            valid_rows, FUN = estimate_Q_M,
            Y = Y, A = A, C = C, M1 = M1, M2 = M2,
            DeltaA = DeltaA, DeltaM = DeltaM,
            all_mediator_values = all_mediator_values,
            verbose = verbose,
            return_models = return_models,
            SL_Q_M = SL_Q_M, a_0 = a_0,
            stratify = stratify,
            glm_Q_M = glm_Q_M,
          )
        }
        # re-order predictions
        Q_M_n <- reorder_list(Q_M_out, a_0 = a_0, valid_rows = valid_rows, 
                              n_SL = n_SL, n = n, which_nuisance = "MED",
                              n_mediator_values = n_mediator_values,
                              n_M1_values = length(unique(M1)),
                              n_M2_values = length(unique(M2)))
        #! TO DO: Get rid of the stupid names that come out after doing 
        #! predict. 

        # obtain list of outcome regression fits
        Q_M_n_models <- drtmle:::extract_models(Q_M_out)
    }
	
    # get marginalized parameters

    # first we target Qbar
    # then we marginalize the resultant Qbar over relevant distributions
    #    - Qbar a_star over joint of M1 M2 under a_star
    #    - Qbar a over joint of M1 M2 under a
    #    - Qbar a over marginal of M2 under a_star
    #    - Qbar a over marginal of M1 under a_star
    #    - Qbar a over marginal of M1 under a

    Qbar_n_tmle <- target_Qbar(Y = Y, A = A, M1 = M1, M2 = M2, 
                                a = a, a_star = a_star,
                                all_mediator_values = all_mediator_values,
                                Qbar_n = Qbar_n, gn = gn, Q_M_n = Q_M_n)

    # marginalize targeted outcome regressions
    Qbarbar_n <- get_Qbarbar(Qbar_n = Qbar_n, Q_M_n = Q_M_n,
                             all_mediator_values = all_mediator_values,
                             unique_M1_values = unique(M1),
                             unique_M2_values = unique(M2))    

    # marginalize targeted outcome regressions
    Qbarbar_n_star <- get_Qbarbar(Qbar_n = Qbar_n_tmle, Q_M_n = Q_M_n,
                             all_mediator_values = all_mediator_values,
                             unique_M1_values = unique(M1),
                             unique_M2_values = unique(M2))

    # targeted them
    Qbarbar_n_tmle <- target_Qbarbar(Qbar = Qbar_n_tmle, Qbarbar = Qbarbar_n_star, 
                                     Y = Y, A = A, a = a, a_star = a_star,
                                     M1 = M1, M2 = M2, all_mediator_values = all_mediator_values,
                                     tol = tol, gn = gn, max_iter = max_iter) 

    # get EIFs for TMLEs (sanity check)
    eif_total_tmle <- evaluate_eif_total(Y = Y, A = A, Qbar = Qbar_n_tmle, 
                                    Qbarbar = Qbarbar_n_tmle, 
                                    gn = gn, a = a, a_star = a_star)
    eif_direct_tmle <- evaluate_eif_direct(Y = Y, A = A, M1 = M1, M2 = M2, 
                                      Qbar = Qbar_n_tmle, Q_M = Q_M_n, 
                                      Qbarbar = Qbarbar_n_tmle, gn = gn, all_mediator_values = all_mediator_values,
                                      a = a, a_star = a_star, use_conditional = TRUE)
    eif_indirect_M1_tmle <- evaluate_eif_indirect_M1(Y = Y, A = A, M1 = M1, M2 = M2, 
                                      Qbar = Qbar_n_tmle, Q_M = Q_M_n, 
                                      Qbarbar = Qbarbar_n_tmle, gn = gn, 
                                      a = a, a_star = a_star, all_mediator_values = all_mediator_values)
    eif_indirect_M2_tmle <- evaluate_eif_indirect_M2(Y = Y, A = A, M1 = M1, M2 = M2, 
                                      Qbar = Qbar_n_tmle, Q_M = Q_M_n, 
                                      Qbarbar = Qbarbar_n_tmle, gn = gn, 
                                      a = a, a_star = a_star, all_mediator_values = all_mediator_values)

    eif_total <- evaluate_eif_total(Y = Y, A = A, Qbar = Qbar_n, 
                                    Qbarbar = Qbarbar_n, 
                                    gn = gn, a = a, a_star = a_star)
    eif_direct <- evaluate_eif_direct(Y = Y, A = A, M1 = M1, M2 = M2, 
                                      Qbar = Qbar_n, Q_M = Q_M_n, 
                                      Qbarbar = Qbarbar_n, gn = gn, all_mediator_values = all_mediator_values,
                                      a = a, a_star = a_star, use_conditional = FALSE)
    eif_indirect_M1 <- evaluate_eif_indirect_M1(Y = Y, A = A, M1 = M1, M2 = M2, 
                                      Qbar = Qbar_n, Q_M = Q_M_n, 
                                      Qbarbar = Qbarbar_n, gn = gn, 
                                      a = a, a_star = a_star, all_mediator_values = all_mediator_values)
    eif_indirect_M2 <- evaluate_eif_indirect_M2(Y = Y, A = A, M1 = M1, M2 = M2, 
                                      Qbar = Qbar_n, Q_M = Q_M_n, 
                                      Qbarbar = Qbarbar_n, gn = gn, 
                                      a = a, a_star = a_star, all_mediator_values = all_mediator_values)

    # point estimates for TMLE
    tmle_est <- c(
      evaluate_total_effect(Qbarbar = Qbarbar_n_tmle),
      evaluate_direct_effect(Qbarbar = Qbarbar_n_tmle, use_conditional = TRUE),
      evaluate_indirect_effect_M1(Qbarbar = Qbarbar_n_tmle),
      evaluate_indirect_effect_M2(Qbarbar = Qbarbar_n_tmle),
      evaluate_covariance_effect_M1M2(Qbarbar = Qbarbar_n_tmle, use_conditional = TRUE)
    )

    # point estimates for one-step
    plugins <- c(
      evaluate_total_effect(Qbarbar = Qbarbar_n),
      evaluate_direct_effect(Qbarbar = Qbarbar_n, use_conditional = FALSE),
      evaluate_indirect_effect_M1(Qbarbar = Qbarbar_n),
      evaluate_indirect_effect_M2(Qbarbar = Qbarbar_n),
      evaluate_covariance_effect_M1M2(Qbarbar = Qbarbar_n, use_conditional = FALSE)
    )
    # browser()
    corrections <- c(
      mean(eif_total), mean(eif_direct), mean(eif_indirect_M1), mean(eif_indirect_M2), 
      mean(eif_total - eif_direct - eif_indirect_M1 - eif_indirect_M2)     
    )

    onestep_est <- plugins + corrections

    # covariance
    eif_matrix <- cbind(eif_total, eif_direct, eif_indirect_M1, eif_indirect_M2, 
                        eif_total - eif_direct - eif_indirect_M1 - eif_indirect_M2)
    plugin_cov <- cov(eif_matrix) / n

    if(targeted_se){
      eif_matrix_tmle <- cbind(eif_total_tmle, eif_direct_tmle, eif_indirect_M1_tmle, eif_indirect_M2_tmle, 
                          eif_total_tmle - eif_direct_tmle - eif_indirect_M1_tmle - eif_indirect_M2_tmle)
      tmle_cov <- cov(eif_matrix_tmle) / n
    }else{
      tmle_cov <- plugin_cov
    }

    out <- list()
    out$tmle <- list(est = tmle_est, cov = tmle_cov)
    out$aiptw <- list(est = onestep_est, cov = plugin_cov)
    out$plugin <- plugins
    if(return_models){
      out$fm <- list(g = gnMod, Qbar = Qbar_n_models, Q_M = Q_M_n_models)
    }else{
      out$fm <- NULL
    }
    class(out) <- "intermed"
    return(out)
}

#' Compute confidence intervals for drtmle and adaptive_iptw@
#' @param ... Arguments to be passed to method
#' @export
ci <- function(...) {
  UseMethod("ci")
}

#' Confidence intervals for interventional mediation effects
#' 
#' @param object An object of class \code{"intermed"}
#' @param est The estimate to obtain CI around. Options are \code{"tmle"} and \code{"aiptw"}
#' @param level Level of the confidence interval
#' @param ... Other options (ignored)
#' @export
ci.intermed <- function(object, est = "aiptw", level = 0.95, ...){
  out <- vector(mode = "list", length = length(est))
  names(out) <- est
  for (i in seq_along(est)) {
    out[[i]] <- matrix(NA, nrow = 5, ncol = 3)
    for (j in seq_len(5)) {
      out[[i]][j, ] <-
        rep(object[[est[i]]]$est[j], 3) +
        stats::qnorm(c(0.5, (1 - level) / 2, (1 + level) / 2)) *
          rep(sqrt(object[[est[i]]]$cov[j, j]), 3)
    }
    row.names(out[[i]]) <- c("Total", "Direct", "Indirect M1", "Indirect M2", "Covar. M1/M2")
    colnames(out[[i]]) <- c("est", "cil", "ciu")
  }
  return(out)
}

#' Print the AIPTW results 
#' @export
#' @method print intermed

print.intermed <- function(x, ...) {
  tmp <- list(
    est = cbind(x$aiptw$est),
    cov = x$aiptw$cov
  )
  row.names(tmp$est) <- c("Total", "Direct", "Indirect M1", "Indirect M2", "Covar. M1/M2")
  colnames(tmp$est) <- ""
  row.names(tmp$cov) <- colnames(tmp$cov) <- c("Total", "Direct", "Indirect M1", "Indirect M2", "Covar. M1/M2")
  print(tmp)
  invisible(tmp)
}