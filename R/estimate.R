#' estimate_Qbar
#'
#' Function to estimate outcome regression as a function of \code{A, C, M1,} 
#' and \code{M2}. Because later we will need to marginalize these estimates over
#' estimated distributions of \code{M1} and \code{M2}, the output includes the predicted 
#' value for each C_i, i = 1, ..., n and for every value of \code{A} in 
#' \code{a_0}. The output is formatted as an n-length list where there is one entry
#' for each observation. This entry includes a list of predicted values under each treatment
#' \code{Qbar_a_0}, which is itself a list with a vector of predictions for each value of 
#' \code{a_0}. Also included is an entry called \code{which_M1_obs}, which indicates rows of
#' the \code{all_mediator_values} that correspond to this observation's observed value of 
#' \code{M1} and \code{M2}. Similarly, there is a vector \code{which_M2_obs}, and also a vector
#' \code{which_M1_M2_obs}, which indicates the row of \code{all_mediator_values} that 
#' corresponds to this observation's observed value of BOTH \code{M1} and \code{M2}. 
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
#' @param verbose A boolean indicating whether to print status updates.
#' @param return_models A boolean indicating whether to return model fits for the
#'  outcome regression, propensity score, and reduced-dimension regressions.
#' @param glm_Qbar A character describing a formula to be used in the call to
#'  \code{glm} for the outcome regression.
#' @param a_0 A list of fixed treatment values
#' @param family A character passed to \code{SuperLearner}
#' @param stratify A \code{boolean} indicating whether to estimate the outcome
#'  regression separately for observations with \code{A} equal to 0/1 (if
#'  \code{TRUE}) or to pool across \code{A} (if \code{FALSE}).
#' @param valid_rows A \code{list} of length \code{cvFolds} containing the row
#'  indexes of observations to include in validation fold.
#' @param ... Additional arguments (not currently used)
#'
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#
estimate_Qbar <- function(Y, A, M1, M2, C, DeltaA, DeltaM, SL_Qbar, glm_Qbar = NULL, a_0, stratify,
                      family, verbose = FALSE, return_models = FALSE,
                      valid_rows, all_mediator_values, ...) {
  if (is.null(SL_Qbar) & is.null(glm_Qbar)) {
    stop("Specify Super Learner library or GLM formula for Q")
  }
  if (!is.null(SL_Qbar) & !is.null(glm_Qbar)) {
    warning(paste0(
      "Super Learner library and GLM formula specified.",
      " Proceeding with Super Learner only."
    ))
    glm_Qbar <- NULL
  }
  # subset data into training and validation sets
  if (length(valid_rows) != length(Y)) {
    trainY <- Y[-valid_rows]
    trainA <- A[-valid_rows]
    trainC <- C[-valid_rows, , drop = FALSE]
    trainM1 <- M1[-valid_rows]
    trainM2 <- M2[-valid_rows]
    trainDeltaA <- DeltaA[-valid_rows]
    trainDeltaM <- DeltaM[-valid_rows]
    validC <- C[valid_rows, , drop = FALSE]
    validA <- A[valid_rows]
    validM1 <- M1[valid_rows]
    validM2 <- M2[valid_rows]
    validY <- Y[valid_rows]
    validDeltaM <- DeltaM[valid_rows]
    validDeltaA <- DeltaA[valid_rows]
    valid_n <- length(valid_rows)
  } else {
    trainA <- validA <- A
    trainC <- validC <- C
    trainY <- validY <- Y
    trainM1 <- validM1 <- M1
    trainM2 <- validM2 <- M2
    trainDeltaA <- validDeltaA <- DeltaA
    trainDeltaM <- validDeltaM <- DeltaM
    valid_n <- length(A)
  }

  # include only DeltaA = 1 and DeltaM = 1 folks
  include <- (trainDeltaA == 1) & (trainDeltaM == 1)

  # Super Learner
  if (!is.null(SL_Qbar)) {
    if (!stratify) {
      if (length(SL_Qbar) > 1 | is.list(SL_Qbar)) {
        fm <- SuperLearner::SuperLearner(
          Y = trainY[include],
          X = data.frame(A = trainA, trainC, M1 = trainM1, M2 = trainM2)[include, , drop = FALSE],
          verbose = verbose, family = family, SL.library = SL_Qbar,
          method = if(family$family == "binomial"){
            drtmle:::tmp_method.CC_nloglik()
          }else{
              drtmle:::tmp_method.CC_LS()
          }
          )

        # prediction is done here
        Qbar_list <- vector(mode = "list", length = valid_n)
        for(i in seq_len(valid_n)){
          Qn <- sapply(a_0, function(x) {
            as.numeric(stats::predict(
              fm,
              newdata = suppressWarnings(data.frame(A = x, M1 = all_mediator_values$M1, 
                                   M2 = all_mediator_values$M2, validC[i, , drop = FALSE])),
              onlySL = TRUE
            )[[1]])
          }, simplify = FALSE)
          idx_M1 <- which(all_mediator_values$M1 == validM1[i])
          idx_M2 <- which(all_mediator_values$M2 == validM2[i])
          idx_M1_M2 <- intersect(idx_M1, idx_M2)
          Qbar_list[[i]] <- list(
            Qbar_a_0 = Qn, which_M1_obs = idx_M1, 
            which_M2_obs = idx_M2, which_M1_M2_obs = idx_M1_M2
          )          
        }    
      } else if (length(SL_Qbar) == 1) {
        fm <- do.call(SL_Qbar, args = list(
          Y = trainY[include],
          X = data.frame(A = trainA, trainC, M1 = trainM1, M2 = trainM2)[include, , drop = FALSE],
          verbose = verbose, newX = data.frame(A = validA, validC, M1 = validM1, M2 = validM2),
          obsWeights = rep(1, length(trainA[include])),
          family = family
        ))            
        Qbar_list <- vector(mode = "list", length = valid_n)
        for(i in seq_len(valid_n)){
          Qn <- sapply(a_0, function(x) {
            as.numeric(stats::predict(object = fm$fit, newdata = data.frame(A = x, M1 = all_mediator_values$M1, 
                                   M2 = all_mediator_values$M2, validC[i, , drop = FALSE])
            ))  
          }, simplify = FALSE)
          idx_M1 <- which(all_mediator_values$M1 == validM1[i])
          idx_M2 <- which(all_mediator_values$M2 == validM2[i])
          idx_M1_M2 <- intersect(idx_M1, idx_M2)
          Qbar_list[[i]] <- list(
            Qbar_a_0 = Qn, which_M1_obs = idx_M1, 
            which_M2_obs = idx_M2, which_M1_M2_obs = idx_M1_M2
          )
        }
      }
    } else { # stratified regressions
      if (length(SL_Qbar) > 1 | is.list(SL_Qbar)) {
        Qbar_list <- vector(mode = "list", length = valid_n)        
        fm <- sapply(a_0, function(x) {
          include2 <- trainA == x
          # handle NAs properly
          include2[is.na(include2)] <- FALSE
          fm <- SuperLearner::SuperLearner(
            Y = trainY[include2 & include],
            X = data.frame(trainC, M1 = trainM1, M2 = trainM2)[include2 & include, , drop = FALSE],
            newX = data.frame(validC, M1 = validM1, M2 = validM2), 
            verbose = verbose, family = family,
            SL.library = SL_Qbar, method = if(family$family ==
              "binomial"){
              tmp_method.CC_nloglik()
              }else{
                tmp_method.CC_LS()
              })
          return(fm)
        }, simplify = FALSE)
        for(i in seq_len(valid_n)){
          Qbar_list[[i]] <- vector(mode = "list", length = 4)
          names(Qbar_list[[i]]) <- c("Qbar_a_0", "which_M1_obs",
                                     "which_M2_obs", "which_M1_M2_obs")
          Qbar_list[[i]]$Qbar_a_0[[1]] <- predict(fm[[1]], newdata = data.frame(validC[i, , drop = FALSE],
                                                                                    M1 = all_mediator_values$M1,
                                                                                    M2 = all_mediator_values$M2))         
          Qbar_list[[i]]$Qbar_a_0[[2]] <- predict(fm[[2]], newdata = data.frame(validC[i, , drop = FALSE],
                                                                                    M1 = all_mediator_values$M1,
                                                                                    M2 = all_mediator_values$M2))
          Qbar_list[[i]]$which_M1_obs <- which(all_mediator_values$M1 == validM1[i])
          Qbar_list[[i]]$which_M2_obs <- which(all_mediator_values$M2 == validM2[i])
          Qbar_list[[i]]$which_M1_M2_obs <- intersect(Qbar_list[[i]]$which_M1_obs, Qbar_list[[i]]$which_M2_obs)        
        }
      } else if (length(SL_Qbar) == 1) {
        Qbar_list <- vector(mode = "list", length = valid_n)        
        fm <- sapply(a_0, function(x) {
          include2 <- trainA == x
          # handle NAs properly
          include2[is.na(include2)] <- FALSE
          # call function
          fm <- do.call(SL_Qbar, args = list(
            Y = trainY[include2 & include],
            X = data.frame(trainC, M1 = trainM1, M2 = trainM2)[include2 & include, , drop = FALSE],
            newX = data.frame(validC, M1 = validM1, M2 = validM2), verbose = verbose,
            obsWeights = rep(1, sum(include2 & include)),
            family = family
          ))
          return(fm)
        }, simplify = FALSE)
        for(i in seq_len(valid_n)){
          Qbar_list[[i]] <- vector(mode = "list", length = 4)
          names(Qbar_list[[i]]) <- c("Qbar_a_0", "which_M1_obs",
                                     "which_M2_obs", "which_M1_M2_obs")
          Qbar_list[[i]]$Qbar_a_0[[1]] <- predict(fm[[1]], newdata = data.frame(validC[i, , drop = FALSE],
                                                                                M1 = all_mediator_values$M1,
                                                                                M2 = all_mediator_values$M2))         
          Qbar_list[[i]]$Qbar_a_0[[2]] <- predict(fm[[2]], newdata = data.frame(validC[i, , drop = FALSE],
                                                                                M1 = all_mediator_values$M1,
                                                                                M2 = all_mediator_values$M2))
          Qbar_list[[i]]$which_M1_obs <- which(all_mediator_values$M1 == validM1[i])
          Qbar_list[[i]]$which_M2_obs <- which(all_mediator_values$M2 == validM2[i])
          Qbar_list[[i]]$which_M1_M2_obs <- intersect(Qbar_list[[i]]$which_M1_obs, Qbar_list[[i]]$which_M2_obs)        
        }
      } # end else if length(SL_Q) == 1
    } # end stratify
  } # end super learner

  # GLM
  if (!is.null(glm_Qbar)) {
    if (!stratify) {
      fm <- stats::glm(
        stats::as.formula(paste0("Y~", glm_Qbar)),
        data = data.frame(Y = trainY, A = trainA, trainC, M1 = trainM1, M2 = trainM2)[
          include, ,
          drop = FALSE
        ], family = family
      )
      Qbar_list <- vector(mode = "list", length = valid_n)        
      for(i in seq_len(valid_n)){
        Qbar_list[[i]] <- vector(mode = "list", length = 4)
        names(Qbar_list[[i]]) <- c("Qbar_a_0", "which_M1_obs",
                                   "which_M2_obs", "which_M1_M2_obs")        
        Qbar_list[[i]]$Qbar_a_0 <- sapply(a_0, function(a, fm) {
          stats::predict(
            fm,
            newdata = suppressWarnings(data.frame(A = a, validC[i, , drop = FALSE], M1 = all_mediator_values$M1,
                                 M2 = all_mediator_values$M2)),
            type = "response"
          )
        }, fm = fm, simplify = FALSE)
        Qbar_list[[i]]$which_M1_obs <- which(all_mediator_values$M1 == validM1[i])
        Qbar_list[[i]]$which_M2_obs <- which(all_mediator_values$M2 == validM2[i])
        Qbar_list[[i]]$which_M1_M2_obs <- intersect(Qbar_list[[i]]$which_M1_obs, Qbar_list[[i]]$which_M2_obs)        
      }
    } else {
      fm <- sapply(a_0, function(a) {
        include2 <- trainA == a
        # handle NAs properly
        include2[is.na(include2)] <- FALSE
        fm <- stats::glm(
          stats::as.formula(paste0(
            "trainY[include2 & include] ~ ", glm_Qbar
          )),
          data = trainC[include2 & include, , drop = FALSE],
          family = family
        )
        return(fm)
      }, simplify = FALSE)
      Qbar_list <- vector(mode = "list", length = valid_n)      
      for(i in seq_len(valid_n)){
        Qbar_list[[i]] <- vector(mode = "list", length = 4)
        names(Qbar_list[[i]]) <- c("Qbar_a_0", "which_M1_obs",
                                   "which_M2_obs", "which_M1_M2_obs")        
        Qbar_list[[i]]$Qbar_a_0 <- lapply(fm, function(fm) {
          stats::predict(
            fm,
            newdata = suppressWarnings(data.frame(validC[i, , drop = FALSE], M1 = all_mediator_values$M1,
                                 M2 = all_mediator_values$M2)),
            type = "response"
          )
        }, fm = fm, simplify = FALSE)
        Qbar_list[[i]]$which_M1_obs <- which(all_mediator_values$M1 == validM1[i])
        Qbar_list[[i]]$which_M2_obs <- which(all_mediator_values$M2 == validM2[i])
        Qbar_list[[i]]$which_M1_M2_obs <- intersect(Qbar_list[[i]]$which_M1_obs, Qbar_list[[i]]$which_M2_obs)        
      }
    }
  }
  out <- list(Qbar_list = Qbar_list, fm = NULL)
  if (return_models) {
    out$fm <- fm
  }
  return(out)
}

#' Helper function to reorder lists according to cvFolds
#' 
#' @param a_list Structured list of nuisance parameters
#' @param a_0 Treatment levels
#' @param valid_rows List of rows of data in validation data for
#' each split.
#' @param which_nuisance Structure of reordering is slightly different for 
#' each relevant nuisance parameter ("OR", "PS", "MED")
#' @param n_SL Number of super learners. If >1, then predictions
#' are averaged over repeated super learners
#' @param n Sample size
#' @param n_mediator_values Number of unique combinations of mediators
reorder_list <- function(a_list, 
                         a_0, 
                         valid_rows,
                         n_SL = 1, 
                         which_nuisance,
                         n,
                         n_mediator_values,
                         n_M1_values = 0,
                         n_M2_values = 0){

  n_cvFolds <- length(valid_rows) / n_SL

  if(which_nuisance == "PS"){
    reduced_out <- vector(mode = "list", length = length(a_0))
  } else {
    reduced_out <- vector(mode = "list", length = n)
  }

  for(j in seq_along(reduced_out)){
    if(which_nuisance == "PS"){ 
      reduced_out[[j]] <- rep(0, n)
    }else if(which_nuisance == "OR"){
      reduced_out[[j]] <- list(rep(0, n_mediator_values), rep(0, n_mediator_values))  
    }else if(which_nuisance == "MED"){
      reduced_out[[j]] <- list(list(rep(0, n_mediator_values), rep(0, n_M1_values), rep(0, n_M2_values)),  
                               list(rep(0, n_mediator_values), rep(0, n_M1_values), rep(0, n_M2_values)))
    }
  }

  # re-order predictions
  for(v in seq_len(n_SL)){
    if(which_nuisance == "PS"){
      outListValid <- unlist(a_list[(n_cvFolds * (v-1) + 1):(v*n_cvFolds)], 
                             recursive = FALSE, use.names = FALSE)
      # this is in 0/1 format 
      outListUnOrd <- do.call(Map, c(c, outListValid[seq(1, length(outListValid), 2)]))
      outList <- vector(mode = "list", length = length(a_0))
      for (i in seq_along(a_0)) {
        outList[[i]] <- rep(NA, n)
        # works because valid_rows are the same across repeated SLs
        outList[[i]][unlist(valid_rows)[1:n]] <- outListUnOrd[[i]]
      }
      reduced_out <- mapply(x = reduced_out, y = outList, function(x,y){
        x + y
      }, SIMPLIFY = FALSE)
    }else if(which_nuisance == "OR"){
      # first subset to output from this particular super learner
      # then subset to predictions, as opposed to fm's
      # then get only predictions
      first_subset <- a_list[(n_cvFolds * (v-1) + 1):(v*n_cvFolds)]
      second_subset <- unlist(lapply(first_subset, "[[" , 1), recursive = FALSE)
      out_list <- vector(mode = "list", length = n)
      # 
      ct <- 0
      for (i in unlist(valid_rows)[1:n]){
        ct <- ct + 1
        out_list[[i]] <- second_subset[[ct]]
      }
      reduced_out <- mapply(x = reduced_out, y = out_list, function(x,y){
        list(x[[1]] + y$Qbar_a_0[[1]],
             x[[2]] + y$Qbar_a_0[[2]])
      }, SIMPLIFY = FALSE)
    }else if(which_nuisance == "MED"){
      first_subset <- a_list[(n_cvFolds * (v-1) + 1):(v*n_cvFolds)]
      second_subset <- unlist(lapply(first_subset, "[[" , 1), recursive = FALSE)
      out_list <- vector(mode = "list", length = n)
      ct <- 0
      for (i in unlist(valid_rows)[1:n]){
        ct <- ct + 1
        out_list[[i]] <- second_subset[[ct]]
      }
      reduced_out <- mapply(x = reduced_out, y = out_list, function(x,y){
        list(list(x[[1]][[1]] + y[[1]][[1]],
                  x[[1]][[2]] + y[[1]][[2]],
                  x[[1]][[3]] + y[[1]][[3]]),
             list(x[[2]][[1]] + y[[2]][[1]],
                  x[[2]][[2]] + y[[2]][[2]],
                  x[[2]][[3]] + y[[2]][[3]]))
      }, SIMPLIFY = FALSE)
    }
  }
  # list(rep(0, n_mediator_values), rep(0, n_M1_values), rep(0, n_M2_values))
  if(which_nuisance == "PS"){
    out <- lapply(reduced_out, function(x){ x / n_SL })
  }else if(which_nuisance == "OR"){
    out <- out_list
    for(i in seq_len(n)){
      out[[i]]$Qbar_a_0[[1]] <- reduced_out[[i]][[1]] / n_SL
      out[[i]]$Qbar_a_0[[2]] <- reduced_out[[i]][[2]] / n_SL
    }
  }else if(which_nuisance == "MED"){
    out <- out_list
    for(i in seq_len(n)){
      out[[i]][[1]][[1]] <- reduced_out[[i]][[1]][[1]] / n_SL
      out[[i]][[1]][[2]] <- reduced_out[[i]][[1]][[2]] / n_SL
      out[[i]][[1]][[3]] <- reduced_out[[i]][[1]][[3]] / n_SL
      out[[i]][[2]][[1]] <- reduced_out[[i]][[2]][[1]] / n_SL
      out[[i]][[2]][[2]] <- reduced_out[[i]][[2]][[2]] / n_SL
      out[[i]][[2]][[3]] <- reduced_out[[i]][[2]][[3]] / n_SL
    }
  }
  return(out)
}



#' estimate_Q_M
#'
#' Helper function to estimate the mediator distribution. 
#' Returns an n-length list, where each entry is a 2-length list
#' corresponding to mediator distributions under each treatment assignment.
#' Within each of these is another list where there are three entries corresponding
#' to the bivariate distribution, and each marginal distribution. 
#' 
#' The bivariate distribution is estimated by estimating the conditional
#' distribution of M1 given A, C, and M2 and the marginal distribution of M1 given
#' A and C. In each case, we use a hazard-based estimation approach for estimating
#' these distributions. The 
#' 
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
#' @param SL_Q A vector of characters or a list describing the Super Learner
#'  library to be used for the outcome regression.
#' @param verbose A boolean indicating whether to print status updates.
#' @param return_models A boolean indicating whether to return model fits for the
#'  outcome regression, propensity score, and reduced-dimension regressions.
#' @param glm_Q_M A character describing a formula to be used in the call to
#'  \code{glm} for the outcome regression.
#' @param a_0 A list of fixed treatment values
#' @param family A character passed to \code{SuperLearner}
#' @param stratify A \code{boolean} indicating whether to estimate the outcome
#'  regression separately for observations with \code{A} equal to 0/1 (if
#'  \code{TRUE}) or to pool across \code{A} (if \code{FALSE}).
#' @param valid_rows A \code{list} of length \code{cvFolds} containing the row
#'  indexes of observations to include in validation fold.
#' @param ... Additional arguments (not currently used)
#'
#' @importFrom SuperLearner SuperLearner trimLogit
#' @importFrom stats predict glm as.formula
#
estimate_Q_M <- function(A, M1, M2, C, DeltaA, DeltaM, SL_Q_M, glm_Q_M = NULL, a_0, stratify,
                      verbose = FALSE, return_models = FALSE,
                      valid_rows, all_mediator_values, ...) {
  if (is.null(SL_Q_M) & is.null(glm_Q_M)) {
    stop("Specify Super Learner library or GLM formula for Q_M")
  }
  if (!is.null(SL_Q_M) & !is.null(glm_Q_M)) {
    warning(paste0(
      "Super Learner library and GLM formula specified.",
      " Proceeding with Super Learner only."
    ))
    glm_Q_M <- NULL
  }
  # subset data into training and validation sets
  if (length(valid_rows) != length(M1)) {
    trainA <- A[-valid_rows]
    trainC <- C[-valid_rows, , drop = FALSE]
    trainM1 <- M1[-valid_rows]
    trainM2 <- M2[-valid_rows]
    trainDeltaA <- DeltaA[-valid_rows]
    trainDeltaM <- DeltaM[-valid_rows]
    validC <- C[valid_rows, , drop = FALSE]
    validA <- A[valid_rows]
    validM1 <- M1[valid_rows]
    validM2 <- M2[valid_rows]
    validDeltaM <- DeltaM[valid_rows]
    validDeltaA <- DeltaA[valid_rows]
    valid_n <- length(valid_rows)
  } else {
    trainA <- validA <- A
    trainC <- validC <- C
    trainM1 <- validM1 <- M1
    trainM2 <- validM2 <- M2
    trainDeltaA <- validDeltaA <- DeltaA
    trainDeltaM <- validDeltaM <- DeltaM
    valid_n <- length(A)
  }

  # include only DeltaA = 1 and DeltaM = 1 folks
  include <- (trainDeltaA == 1) & (trainDeltaM == 1)

  # make a long formatted data set for the conditional hazard estimation
  # of M1 given A, C, M2
  long_train_data_M1 <- format_long_hazards(A = trainM1[include], W = data.frame(trainC, M2 = trainM2, A = trainA)[include,, drop = FALSE], 
                      wts = rep(1, length(M1)), type = "equal_range",
                      n_bins = length(unique(M1)), breaks = NULL)
  # remove weights column
  long_train_data_M1[[1]] <- long_train_data_M1[[1]][,-ncol(long_train_data_M1[[1]])]
  
  #~~~~~
  # remove max value of M1 rows
  max_bin_id <- max(long_train_data_M1$data$bin_id)
  long_train_data_M1$data <- long_train_data_M1$data[-which(long_train_data_M1$data$bin_id == max_bin_id), ]
  #~~~~~

  # make a long formatted data set for the conditional hazard estimation
  # of M2 given A, C
  long_train_data_M2 <- format_long_hazards(A = trainM2[include], 
                      W = data.frame(trainC, A = trainA)[include,, drop = FALSE], 
                      wts = rep(1, length(M2)), type = "equal_range",
                      n_bins = length(unique(M2)), breaks = NULL)
  # remove weights column
  long_train_data_M2[[1]] <- long_train_data_M2[[1]][,-ncol(long_train_data_M2[[1]])]

  #~~~~~
  # remove max value of M1 rows
  max_bin_id <- max(long_train_data_M2$data$bin_id)
  long_train_data_M2$data <- long_train_data_M2$data[-which(long_train_data_M2$data$bin_id == max_bin_id), ]
  #~~~~~

  # Super Learner
  if (!is.null(SL_Q_M)) {
    if (!stratify) {
      if (length(SL_Q_M) > 1 | is.list(SL_Q_M)) {
        fm_M1_given_M2 <- SuperLearner::SuperLearner(
          Y = long_train_data_M1$data$in_bin,
          X = long_train_data_M1$data[ , 3:ncol(long_train_data_M1$data)],
          id = long_train_data_M1$data$obs_id,
          verbose = verbose, family = binomial(), SL.library = SL_Q_M,
          method = drtmle:::tmp_method.CC_nloglik()
        )
        
        fm_M2 <- SuperLearner::SuperLearner(
          Y = long_train_data_M2$data$in_bin,
          X = long_train_data_M2$data[ , 3:ncol(long_train_data_M2$data)],
          id = long_train_data_M2$data$obs_id,
          verbose = verbose, family = binomial(), SL.library = SL_Q_M,
          method = drtmle:::tmp_method.CC_nloglik()
        )
        hp <- "SuperLearner"
      } else if (length(SL_Q_M) == 1) {
        fm_M1_given_M2 <- do.call(SL_Q_M, args = list(
          Y = long_train_data_M1$data$in_bin,
          X = long_train_data_M1$data[ , 3:ncol(long_train_data_M1$data)],
          verbose = verbose, newX = long_train_data_M1$data[ , 3:ncol(long_train_data_M1$data)],
          obsWeights = rep(1, length(long_train_data_M1$data[ , 1])),
          family = binomial()
        ))

        fm_M2 <- do.call(SL_Q_M, args = list(
          Y = long_train_data_M2$data$in_bin,
          X = long_train_data_M2$data[ , 3:ncol(long_train_data_M2$data)],
          verbose = verbose, newX = long_train_data_M2$data[ , 3:ncol(long_train_data_M2$data)],
          obsWeights = rep(1, length(long_train_data_M2$data[ , 1])),
          family = binomial()
        ))
        hp <- "single_algo"
      }
      # get predicted value back... 
      list_by_a_0 <- sapply(a_0, FUN = predict_density, 
                            sl_fit_conditional = fm_M1_given_M2, 
                            sl_fit_marginal = fm_M2,
                            all_mediator_values = all_mediator_values, 
                            how_predict = hp,
                            validC = validC, 
                            M1 = M1, M2 = M2, valid_n = valid_n,
                            stratify = FALSE,
                            simplify = FALSE)
      # re-format 
      Q_M_list <- mapply(x = list_by_a_0[[1]], y = list_by_a_0[[2]], FUN = function(x,y){ list(x, y) },
                    SIMPLIFY = FALSE)

    } else { # stratified regressions
      if (length(SL_Q) > 1 | is.list(SL_Q)) {
        fm <- sapply(a_0, function(x) {
          include2 <- long_train_data_M1$data$A == x
          
          fm_M1_given_M2 <- SuperLearner::SuperLearner(
            Y = long_train_data_M1$data$in_bin[include2],
            # remove A column
            X = long_train_data_M1$data[include2 , 3:(ncol(long_train_data_M1$data)-1)],
            id = long_train_data_M1$data$obs_id[include2],
            verbose = verbose, family = binomial(),
            SL.library = SL_Q_M, method = drtmle:::tmp_method.CC_nloglik()
          )

          fm_M2 <- SuperLearner::SuperLearner(
            Y = long_train_data_M2$data$in_bin[include2],
            # remove A column
            X = long_train_data_M2$data[include2 , 3:(ncol(long_train_data_M2$data)-1)],
            id = long_train_data_M2$data$obs_id[include2],
            verbose = verbose, family = binomial(),
            SL.library = SL_Q_M, method = drtmle:::tmp_method.CC_nloglik()
          )
          
          return(list(fm_M1_given_M2, fm_M2))
        }, simplify = FALSE)
        hp <- "SuperLearner"
      } else if (length(SL_Q) == 1) {
        fm <- sapply(a_0, function(x) {
          include2 <- long_train_data_M1$data$A == x

          # call function
          fm_M1_given_M2 <- do.call(SL_Q_M, args = list(
            Y = long_train_data_M1$data$in_bin[include2],
            # remove A column
            X = long_train_data_M1$data[include2 , 3:(ncol(long_train_data_M1$data)-1)],
            newX = long_train_data_M1$data[include2 , 3:(ncol(long_train_data_M1$data)-1)], 
            verbose = verbose,
            obsWeights = rep(1, length(long_train_data_M1$data[include2 , 3])),
            family = binomial()
          ))
          fm_M2 <- do.call(SL_Q_M, args = list(
            Y = long_train_data_M2$data$in_bin[include2],
            # remove A column
            X = long_train_data_M2$data[include2 , 3:(ncol(long_train_data_M2$data)-1)],
            newX = long_train_data_M2$data[include2 , 3:(ncol(long_train_data_M2$data)-1)], 
            verbose = verbose,
            obsWeights = rep(1, length(long_train_data_M2$data[include2 , 3])),
            family = binomial()
          ))
          return(list(fm_M1_given_M2, fm_M2))
        }, simplify = FALSE)   
        hp <- "single_algo"       
      } # end else if length(SL_Q) == 1

      list_by_a_0 <- mapply(a_val = a_0, fm = fm, FUN = function(a_val, fm){
        predict_density(sl_fit_conditional = fm[[1]], 
                        sl_fit_marginal = fm[[2]],
                        all_mediator_values = all_mediator_values, 
                        validC = validC, 
                        how_predict = hp,
                        M1 = M1, M2 = M2, valid_n = valid_n,
                        stratify = TRUE)
      }, SIMPLIFY = FALSE)
      # re-format 
      Q_M_list <- mapply(x = list_by_a_0[[1]], y = list_by_a_0[[2]], FUN = function(x,y){ list(x, y) },
                         SIMPLIFY = FALSE)
    } # end stratify
  } # end super learner

  # GLM
  if (!is.null(glm_Q_M)) {
    if (!stratify) {
      fm_M1_given_M2 <- stats::glm(
        stats::as.formula(paste0("in_bin ~ ", glm_Q_M$M1)),
        data = long_train_data_M1$data[ , 2:ncol(long_train_data_M1$data)], family = binomial()
      )
      fm_M2 <- stats::glm(
        stats::as.formula(paste0("in_bin ~ ", glm_Q_M$M2)),
        data = long_train_data_M2$data[ , 2:ncol(long_train_data_M2$data)], family = binomial()
      )
      
      # get predicted value back... 
      list_by_a_0 <- sapply(a_0, FUN = predict_density, 
                            sl_fit_conditional = fm_M1_given_M2, 
                            sl_fit_marginal = fm_M2,
                            all_mediator_values = all_mediator_values, 
                            validC = validC, 
                            how_predict = "glm",
                            M1 = M1, M2 = M2, valid_n = valid_n,
                            stratify = FALSE,
                            simplify = FALSE)
      # re-format 
      Q_M_list <- mapply(x = list_by_a_0[[1]], y = list_by_a_0[[2]], FUN = function(x,y){ list(x, y) },
                    SIMPLIFY = FALSE)
    } else {
      fm <- sapply(a_0, function(a) {
        include2 <- long_train_data_M1$data$A == x
        fm_M1_given_M2 <- stats::glm(
          stats::as.formula(paste0(
            "in_bin ~ ", glm_Q_M$M1
          )),
          data = long_train_data_M1[include2, 2:(ncol(long_train_data_M1$data)-1)],
          family = binomial()
        )
        fm_M2 <- stats::glm(
          stats::as.formula(paste0(
            "in_bin ~ ", glm_Q_M$M2
          )),
          data = long_train_data_M2[include2, 2:(ncol(long_train_data_M1$data)-1)],
          family = binomial()
        )        
        return(list(fm_M1_given_M2, fm_M2))
      }, simplify = FALSE)
      list_by_a_0 <- mapply(a_val = a_0, fm = fm, FUN = function(a_val, fm){
      predict_density(a_val = a_val,
                      sl_fit_conditional = fm[[1]], 
                      sl_fit_marginal = fm[[2]],
                      all_mediator_values = all_mediator_values, 
                      validC = validC, 
                      how_predict = "glm",
                      M1 = M1, M2 = M2, valid_n = valid_n,
                      stratify = TRUE)
      }, SIMPLIFY = FALSE)
      # re-format 
      Q_M_list <- mapply(x = list_by_a_0[[1]], y = list_by_a_0[[2]], FUN = function(x,y){ list(x, y) },
                         SIMPLIFY = FALSE)    
    }
  }
  out <- list(Q_M_list = Q_M_list, fm = NULL)

  if (return_models) {
    if(!stratify){
      out$fm <- list(fm_M1_given_M2, fm_M2)
    }else{
      out$fm <- fm
    }
  }
  return(out)
}


#' Helper function to marginalize outcome regression over mediator 
#' distributions. Used to get each of the relevant nuisance parameters needed
#' to evaluate the total, direct, and indirect effects. 

get_Qbarbar <- function(Qbar_n, Q_M_n, unique_M1_values, unique_M2_values, all_mediator_values, ...){
 
  # get Qbarbar_M1_times_M2_star_a (indirect)
  M1_times_M2_star_a <- mapply(FUN = get_M1_times_M2_star_a,
                               Qbar_n_i = Qbar_n,
                               Q_M_n_i = Q_M_n,
                               MoreArgs = list(all_mediator_values = all_mediator_values,
                                               unique_M1_values = unique_M1_values, 
                                               unique_M2_values = unique_M2_values),
                               SIMPLIFY = TRUE)
  # get Qbarbar_M1_star_times_M2_star_a (indirect)
  M1_star_times_M2_star_a <- mapply(FUN = get_M1_star_times_M2_star_a,
                             Qbar_n_i = Qbar_n,
                             Q_M_n_i = Q_M_n,
                             MoreArgs = list(all_mediator_values = all_mediator_values,
                                             unique_M1_values = unique_M1_values, 
                                             unique_M2_values = unique_M2_values),
                             SIMPLIFY = TRUE)
  # get Qbarbar_M1_star_times_M2_a (indirect)
  M1_star_times_M2_a <- mapply(FUN = get_M1_star_times_M2_a,
                             Qbar_n_i = Qbar_n,
                             Q_M_n_i = Q_M_n,
                             MoreArgs = list(all_mediator_values = all_mediator_values,
                                             unique_M1_values = unique_M1_values, 
                                             unique_M2_values = unique_M2_values),
                             SIMPLIFY = TRUE)
  # get Qbarbar_M1_star_M2_star_a_star (total + direct, iterative?)
  M1_star_M2_star_a_star <- mapply(FUN = get_M1_star_M2_star_a_star,
                             Qbar_n_i = Qbar_n,
                             Q_M_n_i = Q_M_n,
                             MoreArgs = list(all_mediator_values = all_mediator_values,
                                             unique_M1_values = unique_M1_values, 
                                             unique_M2_values = unique_M2_values),
                             SIMPLIFY = TRUE)
  
  # get Qbarbar_M1_star_M2_star_a (indirect)
  M1_star_M2_star_a <- mapply(FUN = get_M1_star_M2_star_a,
                             Qbar_n_i = Qbar_n,
                             Q_M_n_i = Q_M_n,
                             MoreArgs = list(all_mediator_values = all_mediator_values,
                                             unique_M1_values = unique_M1_values, 
                                             unique_M2_values = unique_M2_values),
                             SIMPLIFY = TRUE)

  # get Qbarbar_M1_M2_a (total)
  M1_M2_a <- mapply(FUN = get_M1_M2_a,
                     Qbar_n_i = Qbar_n,
                     Q_M_n_i = Q_M_n,
                     MoreArgs = list(all_mediator_values = all_mediator_values,
                                     unique_M1_values = unique_M1_values, 
                                     unique_M2_values = unique_M2_values),
                     SIMPLIFY = TRUE)

  # get Qbarbar_M1_star_a
  M1_star_a <- mapply(FUN = get_M1_star_a,
                   Qbar_n_i = Qbar_n,
                   Q_M_n_i = Q_M_n,
                   MoreArgs = list(all_mediator_values = all_mediator_values,
                                   unique_M1_values = unique_M1_values, 
                                   unique_M2_values = unique_M2_values),
                   SIMPLIFY = TRUE)
  # get Qbarbar_M1_a
  M1_a <- mapply(FUN = get_M1_a,
                   Qbar_n_i = Qbar_n,
                   Q_M_n_i = Q_M_n,
                   MoreArgs = list(all_mediator_values = all_mediator_values,
                                   unique_M1_values = unique_M1_values, 
                                   unique_M2_values = unique_M2_values),
                   SIMPLIFY = TRUE)

  # get Qbarbar_M2_star_a
  M2_star_a <- mapply(FUN = get_M2_star_a,
                   Qbar_n_i = Qbar_n,
                   Q_M_n_i = Q_M_n,
                   MoreArgs = list(all_mediator_values = all_mediator_values,
                                   unique_M1_values = unique_M1_values, 
                                   unique_M2_values = unique_M2_values),
                   SIMPLIFY = TRUE)
  # get Qbarbar_M2_a
  M2_a <- mapply(FUN = get_M2_a,
                   Qbar_n_i = Qbar_n,
                   Q_M_n_i = Q_M_n,
                   MoreArgs = list(all_mediator_values = all_mediator_values,
                                   unique_M1_values = unique_M1_values, 
                                   unique_M2_values = unique_M2_values),
                   SIMPLIFY = TRUE)

  return(list(M1_times_M2_star_a = M1_times_M2_star_a,
              M1_star_times_M2_star_a = M1_star_times_M2_star_a,
              M1_star_times_M2_a = M1_star_times_M2_a, 
              M1_star_M2_star_a_star = M1_star_M2_star_a_star,
              M1_star_M2_star_a = M1_star_M2_star_a,
              M1_M2_a = M1_M2_a,
              M1_a = M1_a, M1_star_a = M1_star_a,
              M2_a = M2_a, M2_star_a = M2_star_a))
}

get_M2_star_a <- function(Qbar_n_i, Q_M_n_i, all_mediator_values, unique_M1_values, unique_M2_values){
  Qbarn_a_frame <- data.frame(all_mediator_values, Qbar_n_a_i = Qbar_n_i$Qbar_a_0[[2]])
  Qbarn_a_frame <- Qbarn_a_frame[Qbar_n_i$which_M1_obs, ]
  Q_M2_a_star_frame <- data.frame(M2 = unique_M2_values, 
                              Q_M2_a_star = Q_M_n_i[[1]][[3]])
  nuisance_frame <- reduce_merge(Qbarn_a_frame, Q_M2_a_star_frame)
  out <- sum(with(nuisance_frame,
                  Qbar_n_a_i * Q_M2_a_star))
  return(out)
}

get_M2_a <- function(Qbar_n_i, Q_M_n_i, all_mediator_values, unique_M1_values, unique_M2_values){
  Qbarn_a_frame <- data.frame(all_mediator_values, Qbar_n_a_i = Qbar_n_i$Qbar_a_0[[2]])
  Qbarn_a_frame <- Qbarn_a_frame[Qbar_n_i$which_M1_obs, ]
  Q_M2_a_frame <- data.frame(M2 = unique_M2_values, 
                              Q_M2_a = Q_M_n_i[[2]][[3]])
  nuisance_frame <- reduce_merge(Qbarn_a_frame, Q_M2_a_frame)
  out <- sum(with(nuisance_frame,
                  Qbar_n_a_i * Q_M2_a))
  return(out)
}

get_M1_star_a <- function(Qbar_n_i, Q_M_n_i, all_mediator_values, unique_M1_values, unique_M2_values){
  Qbarn_a_frame <- data.frame(all_mediator_values, Qbar_n_a_i = Qbar_n_i$Qbar_a_0[[2]])
  Qbarn_a_frame <- Qbarn_a_frame[Qbar_n_i$which_M2_obs, ]
  Q_M1_a_star_frame <- data.frame(M1 = unique_M1_values, 
                              Q_M1_a_star = Q_M_n_i[[1]][[2]])
  nuisance_frame <- reduce_merge(Qbarn_a_frame, Q_M1_a_star_frame)
  out <- sum(with(nuisance_frame,
                  Qbar_n_a_i * Q_M1_a_star))
  return(out)
}

get_M1_a <- function(Qbar_n_i, Q_M_n_i, all_mediator_values, unique_M1_values, unique_M2_values){
  Qbarn_a_frame <- data.frame(all_mediator_values, Qbar_n_a_i = Qbar_n_i$Qbar_a_0[[2]])
  Qbarn_a_frame <- Qbarn_a_frame[Qbar_n_i$which_M2_obs, ]
  Q_M1_a_frame <- data.frame(M1 = unique_M1_values, 
                              Q_M1_a = Q_M_n_i[[2]][[2]])
  nuisance_frame <- reduce_merge(Qbarn_a_frame, Q_M1_a_frame)
  out <- sum(with(nuisance_frame,
                  Qbar_n_a_i * Q_M1_a))
  return(out)
}

get_M1_star_times_M2_star_a <- function(Qbar_n_i, Q_M_n_i, all_mediator_values, unique_M1_values, unique_M2_values){
  Qbarn_a_frame <- data.frame(all_mediator_values, Qbar_n_a_i = Qbar_n_i$Qbar_a_0[[2]],
                              id = seq_along(Qbar_n_i$Qbar_a_0[[2]]))
  Q_M1_a_star_frame <- data.frame(M1 = unique_M1_values, 
                              Q_M1_a_star = Q_M_n_i[[1]][[2]])        
  Q_M2_a_star_frame <- data.frame(M2 = unique_M2_values, 
                              Q_M2_a_star = Q_M_n_i[[1]][[3]])


  nuisance_frame <- Reduce("reduce_merge", list(Qbarn_a_frame, Q_M2_a_star_frame, Q_M1_a_star_frame))
  nuisance_frame <- nuisance_frame[order(nuisance_frame$id), ]
  out <- sum(with(nuisance_frame,
                  Qbar_n_a_i * Q_M1_a_star * Q_M2_a_star))
  return(out)
}

get_M1_times_M2_star_a <- function(Qbar_n_i, Q_M_n_i, all_mediator_values, unique_M1_values, unique_M2_values){
  Qbarn_a_frame <- data.frame(all_mediator_values, Qbar_n_a_i = Qbar_n_i$Qbar_a_0[[2]],
                              id = seq_along(Qbar_n_i$Qbar_a_0[[2]]))
  Q_M1_a_frame <- data.frame(M1 = unique_M1_values, 
                            Q_M1_a = Q_M_n_i[[2]][[2]])        
  Q_M2_a_star_frame <- data.frame(M2 = unique_M2_values, 
                              Q_M2_a_star = Q_M_n_i[[1]][[3]])

  nuisance_frame <- Reduce("reduce_merge", list(Qbarn_a_frame, Q_M2_a_star_frame, Q_M1_a_frame))
  nuisance_frame <- nuisance_frame[order(nuisance_frame$id), ]
  out <- sum(with(nuisance_frame,
                  Qbar_n_a_i * Q_M1_a * Q_M2_a_star))
  return(out)
}

get_M1_star_times_M2_a <- function(Qbar_n_i, Q_M_n_i, all_mediator_values, unique_M1_values, unique_M2_values){
  Qbarn_a_frame <- data.frame(all_mediator_values, Qbar_n_a_i = Qbar_n_i$Qbar_a_0[[2]],
                              id = seq_along(Qbar_n_i$Qbar_a_0[[2]]))
  Q_M1_a_star_frame <- data.frame(M1 = unique_M1_values, 
                            Q_M1_a_star = Q_M_n_i[[1]][[2]])        
  Q_M2_a_frame <- data.frame(M2 = unique_M2_values, 
                            Q_M2_a = Q_M_n_i[[2]][[3]])

  nuisance_frame <- Reduce("reduce_merge", list(Qbarn_a_frame, Q_M2_a_frame, Q_M1_a_star_frame))
  nuisance_frame <- nuisance_frame[order(nuisance_frame$id), ]
  out <- sum(with(nuisance_frame,
                  Qbar_n_a_i * Q_M1_a_star * Q_M2_a))
  return(out)
}


get_M1_star_M2_star_a_star <- function(Qbar_n_i, Q_M_n_i, all_mediator_values, ...){
  out <- sum(Qbar_n_i$Qbar_a_0[[1]] * Q_M_n_i[[1]][[1]])
  return(out)
}

get_M1_star_M2_star_a <- function(Qbar_n_i, Q_M_n_i, all_mediator_values, ...){
  out <- sum(Qbar_n_i$Qbar_a_0[[2]] * Q_M_n_i[[1]][[1]])
  return(out)
}

get_M1_M2_a <- function(Qbar_n_i, Q_M_n_i, all_mediator_values, ...){
  out <- sum(Qbar_n_i$Qbar_a_0[[2]] * Q_M_n_i[[2]][[1]])
  return(out)
}
