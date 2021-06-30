#' Generate long format hazards data for conditional density estimation
#'
#' @param A The \code{numeric} vector or similar of the observed values of an
#'  intervention for a group of observational units of interest.
#' @param W A \code{data.frame}, \code{matrix}, or similar giving the values of
#'  baseline covariates (potential confounders) for the observed units whose
#'  observed intervention values are provided in the previous argument.
#' @param wts A \code{numeric} vector of observation-level weights. The default
#'  is to weight all observations equally.
#' @param type A \code{character} indicating the strategy to be used in creating
#'  bins along the observed support of the intervention \code{A}. For bins of
#'  equal range, use "equal_range" and consider consulting the documentation of
#'  \code{ggplot2::cut_interval} for more information. To ensure each bins has
#'  the same number of points, use "equal_mass" and consult the documentation of
#'  \code{ggplot2::cut_number} for details.
#' @param n_bins Only used if \code{type} is set to \code{"equal_range"} or
#'  \code{"equal_mass"}. This \code{numeric} value indicates the number of bins
#'  that the support of the intervention \code{A} is to be divided into.
#' @param breaks A \code{numeric} vector of break points to be used in dividing
#'  up the support of \code{A}. This is passed as a \code{...} argument to
#'  \code{base::cut.default} by either \code{cut_interval} or \code{cut_number}.
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom ggplot2 cut_interval cut_number
#' @importFrom future.apply future_lapply
#
format_long_hazards <- function(A, W, wts = rep(1, length(A)),
                                type = c(
                                  "equal_range", "equal_mass"
                                ),
                                n_bins = NULL, breaks = NULL) {
  # clean up arguments
  type <- match.arg(type)

  # set grid along A and find interval membership of observations along grid
  if (is.null(breaks) & !is.null(n_bins)) {
    if (type == "equal_range") {
      bins <- ggplot2::cut_interval(A, n_bins,
        right = FALSE,
        ordered_result = TRUE, dig.lab = 12
      )
    } else if (type == "equal_mass") {
      bins <- ggplot2::cut_number(A, n_bins,
        right = FALSE,
        ordered_result = TRUE, dig.lab = 12
      )
    }
    # https://stackoverflow.com/questions/36581075/extract-the-breakpoints-from-cut
    breaks_left <- as.numeric(sub(".(.+),.+", "\\1", levels(bins)))
    breaks_right <- as.numeric(sub(".+,(.+).", "\\1", levels(bins)))
    bin_length <- round(breaks_right - breaks_left, 3)
    bin_id <- as.numeric(bins)
    all_bins <- matrix(seq_len(n_bins), ncol = 1)
    # for predict method, only need to assign observations to existing intervals
  } else if (!is.null(breaks)) {
    # NOTE: findInterval() and cut() might return slightly different results...
    bin_id <- findInterval(A, breaks, all.inside = TRUE)
    all_bins <- matrix(seq_along(breaks), ncol = 1)
  } else {
    stop("Combination of arguments `breaks`, `n_bins` incorrectly specified.")
  }

  # loop over observations to create expanded set of records for each
  reformat_each_obs <- lapply(seq_along(A), function(i) {
    # browser()
    # create indicator and "turn on" indicator for interval membership
    bin_indicator <- rep(0, nrow(all_bins))
    bin_indicator[bin_id[i]] <- 1
    id <- rep(i, nrow(all_bins))

    # get correct value of baseline variables and repeat along intervals
    if (is.null(dim(W))) {
      # assume vector
      obs_w <- rep(W[i], nrow(all_bins))
      names_w <- "W"
    } else {
      # assume two-dimensional array
      obs_w <- rep(as.numeric(W[i, ]), nrow(all_bins))
      obs_w <- matrix(obs_w, ncol = ncol(W), byrow = TRUE)

      # use names from array if present
      if (is.null(names(W))) {
        names_w <- paste("W", seq_len(ncol(W)), sep = "_")
      } else {
        names_w <- names(W)
      }
    }

    # get correct value of weights and repeat along intervals
    # NOTE: the weights are always a vector
    obs_wts <- rep(wts[i], nrow(all_bins))

    # create data table with membership indicator and interval limits
    suppressWarnings(
      hazards_df <- data.frame(cbind(
        id, bin_indicator,
        all_bins, obs_w,
        obs_wts
      ))
    )

    # trim records to simply end at the failure time for a given observation
    hazards_df_reduced <- hazards_df[seq_len(bin_id[i]), ]

    # give explicit names and add to appropriate position in list
    colnames(hazards_df_reduced) <- c("obs_id", "in_bin", "bin_id", names_w, "wts")

    return(hazards_df_reduced)
  })

  # combine observation-level hazards data into larger structure
  reformatted_data <- do.call(rbind, reformat_each_obs)
  out <- list(
    data = reformatted_data,
    breaks =
      if (exists("breaks_left")) {
        breaks_left
      } else {
        NULL
      },
    bin_length =
      if (exists("bin_length")) {
        bin_length
      } else {
        NULL
      }
  )
  return(out)
}

################################################################################

#' Map a predicted hazard to a predicted density for a single observation
#'
#' For a single observation, map a predicted hazard of failure (occurrence in a
#' particular bin, under a given partitioning of the support) to a density.
#'
#' @param hazard_pred_single_obs A \code{numeric} vector of the predicted hazard
#'  of failure in a given bin (under a given partitioning of the support) for a
#'  single observational unit based on a long format data structure (as produced
#'  by \code{\link{format_long_hazards}}). This is simply the probability that
#'  the observed value falls in a corresponding bin, given that it has not yet
#'  failed (fallen in a previous bin).
#'
#
map_hazard_to_density <- function(hazard_pred_single_obs) {
  # number of records for the given observation
  n_records <- length(hazard_pred_single_obs)

  # NOTE: pred_hazard = (1 - pred) if 0 in this bin * pred if 1 in this bin
  if (n_records > 1) {
    hazard_prefailure <- 1 - hazard_pred_single_obs[-n_records]
    hazard_at_failure <- hazard_pred_single_obs[n_records]
    hazard_predicted <- c(hazard_prefailure, hazard_at_failure)
  } else {
    hazard_predicted <- hazard_pred_single_obs
  }

  # multiply hazards across rows to construct the individual-level density
  density_pred_from_hazards <- prod(hazard_predicted)
  return(density_pred_from_hazards)
}

#' Helper function to predict from a super learner fit on the hazard scale
#' to obtain a prediction on the density scale
#' @param sl_fit_conditional A fitted SuperLearner object for the M1 | M2, C, A regression
#' @param sl_fit_marginal A fitted SuperLearner object for the M2 | C, A regression
#' @param all_mediator_values data.frame of all observed mediator values
#' @param a_val Relevant treatment level
#' @param validC Confounders in the validation set
#' @param valid_n Number of obs in validation set
#' @param M1 All M1 values
#' @param M2 All M2 values
#' @param stratify Were sl_fit's stratified?
#' @param how_predict Needed to figure out how to predict from model. Possible choices
#' are "SuperLearner", "single_algo", and "glm"
#' 
predict_density <- function(sl_fit_conditional, 
                            sl_fit_marginal,
                            how_predict, 
                            all_mediator_values, a_val, validC,
                            M1, M2, valid_n, stratify){
  n_mediator_values <- nrow(all_mediator_values)
  valid_n <- nrow(validC)
  
  # list_by_a_0 <- vector(mode = "list", length = 2)
  
  
  
  #### GET PREDICTION OF CONDITIONAL DENSITY OF M1 | M2, C, A = a_val
  # here we copy each observed C a number of times so that
  # each observation has one row for each value of M2
  if(!stratify){
    conf_data <- data.frame(validC[sort(rep(seq_len(valid_n), n_mediator_values)),], 
                            M2 = all_mediator_values$M2, A = a_val)
  }else{
    conf_data <- data.frame(validC[sort(rep(seq_len(valid_n), n_mediator_values)),], 
                            M2 = all_mediator_values$M2)
  }
  newdata <- format_long_hazards(A = rep(all_mediator_values$M1, valid_n),                                      
                                 W = conf_data,
                                 wts = rep(1, length(M1)), type = "equal_range",
                                 n_bins = length(unique(M1)))
  # remove weights column
  newdata$data <- newdata$data[,-ncol(newdata$data)]

  # get predictions, cast to numeric
  if(how_predict == "SuperLearner"){
    newdata$data$haz_pred <- as.numeric(predict(sl_fit_conditional, newdata = newdata$data[ , 3:ncol(newdata$data)])$pred)
  }else if(how_predict == "single_algo"){
    newdata$data$haz_pred <- as.numeric(predict(sl_fit_conditional$fit, newdata = newdata$data[ , 3:ncol(newdata$data)]))
  }else if(how_predict == "glm"){
    newdata$data$haz_pred <- as.numeric(predict(sl_fit_conditional, newdata = newdata$data[ , 2:ncol(newdata$data)],
                                                type = "response"))
  }
  # get density estimates 
  tmp <- c(by(newdata$data, newdata$data$obs_id, function(x){
    map_hazard_to_density(x$haz_pred)
  }))
  # split up by observation
  long_length <- length(conf_data[,1])
  split_densities <- split(tmp, sort(rep(seq_len(valid_n), n_mediator_values)))
  # # now normalize by dividing each conditional density by sum
  # normalized_densities_conditional <- lapply(split_densities, function(x){
  #   # normalize for each unique value of M2
  #   dens_this_id <- unlist(by(data.frame(all_mediator_values, dens = x), 
  #                             all_mediator_values$M2, function(y){
  #                               y$dens/ sum(y$dens)
  #                             }), 
  #                         use.names = FALSE)
  #   return(dens_this_id)
  # })  

  # normalize by replacing each value for the max observed value of M1 
  # by 1 - sum(other density values)
  normalized_densities_conditional <- lapply(split_densities, function(x){
    # normalize for each unique value of M2
    dens_this_id <- unlist(by(data.frame(all_mediator_values, dens = x), 
                              all_mediator_values$M2, function(y){
                                # y$dens/ sum(y$dens)
                                y$dens[y$M1 == max(y$M1)] <- 1 - sum(y$dens[y$M1 != max(y$M1)])
                                return(y$dens)
                              }), 
                          use.names = FALSE)
    return(dens_this_id)
  })

  #### GET PREDICTION OF CONDITIONAL DENSITY OF M2 | C, A = a_val
  # here we copy each observed C a number of times so that
  # each observation has multiple rows for each value of M2
  # While this is a bit redundant, it makes for convenient formatting
  # when multiplying by the conditional density estimates
  if(!stratify){
    conf_data <- data.frame(validC[sort(rep(seq_len(valid_n), n_mediator_values)),], 
                            A = a_val)
  }else{
    conf_data <- data.frame(validC[sort(rep(seq_len(valid_n), n_mediator_values)),])
  }
  newdata <- format_long_hazards(A = rep(all_mediator_values$M2, valid_n),                                      
                                 W = conf_data,
                                 wts = rep(1, length(M2)), type = "equal_range",
                                 n_bins = length(unique(M2)))
  # remove weights column
  newdata$data <- newdata$data[,-ncol(newdata$data)]
  # get predictions, cast to numeric
  if(how_predict == "SuperLearner"){
    newdata$data$haz_pred <- as.numeric(predict(sl_fit_marginal, newdata = newdata$data[ , 3:ncol(newdata$data)])$pred)
  }else if(how_predict == "single_algo"){
    newdata$data$haz_pred <- as.numeric(predict(sl_fit_marginal$fit, newdata = newdata$data[ , 3:ncol(newdata$data)]))
  }else if(how_predict == "glm"){
    newdata$data$haz_pred <- as.numeric(predict(sl_fit_marginal, newdata = newdata$data[ , 2:ncol(newdata$data)],
                                                type = "response"))
  }  
  # get density estimates 
  tmp <- c(by(newdata$data, newdata$data$obs_id, function(x){
    map_hazard_to_density(x$haz_pred)
  }))
  # split up by observation
  long_length <- length(conf_data[,1])
  split_densities <- split(tmp, sort(rep(seq_len(valid_n), n_mediator_values)))

  # # now normalize by dividing each conditional density by sum
  # normalized_densities_marginal <- lapply(split_densities, function(x){
  #   # normalize 
  #   # note that each marginal density value for M2 = m2 is copied
  #   # length(unique(all_mediator_values$M1)) times. So we divide the sum
  #   # by this number to get a properly normalized density
  #   dens_this_id <- x / (sum(x) / length(unique(all_mediator_values$M1)))        
  #   return(dens_this_id)
  # })
  
  #~~~~~ 
  max_M2 <- max(M2)
  # now normalize by dividing each conditional density by sum
  normalized_densities_marginal <- lapply(split_densities, function(x){
    # normalize 
    # note that each marginal density value for M2 = m2 is copied
    # length(unique(all_mediator_values$M1)) times. So we divide the sum
    # by this number to get a properly normalized density
    x[all_mediator_values$M2 == max_M2] <- 1 - sum(x[all_mediator_values$M2 != max_M2]) / length(unique(all_mediator_values$M1))
    return(x)
  })

  # joint density estimates
  joint_densities <- mapply(normalized_densities_marginal, normalized_densities_conditional, FUN = "*",
                            SIMPLIFY = FALSE)
  
  # reduce marginal estimates of M2 to only unique values of M2 
  # note these are ordered according to unique(M2)
  marginal_M2_densities <- lapply(normalized_densities_marginal, function(x){
    x[!duplicated(all_mediator_values$M2)]
  })

  # get marginal estimates of M1 from joint
  # note these are ordered according to unique(M1)
  marginal_M1_densities <- lapply(joint_densities, function(x){
    sapply(unique(M1), function(y){ sum(x[all_mediator_values$M1 == y]) })      
  })

  out <- mapply(x = joint_densities, y = marginal_M1_densities, z = marginal_M2_densities, 
                                FUN = function(x,y,z){list(x,y,z)}, SIMPLIFY = FALSE)
  return(out)
}


#' Helper function for extracting marginal distributions at particular values of 
#' the mediator. 
#' @param Q_M_n_i An entry in the \code{Q_M_n} list, where the first entry in the list
#' corresponds to mediator distributions under \code{a_star}; the second entry in the list
#' corresponds to mediator distributions under \code{a}. In each of these lists, we have
#' the joint of M1,M2; the marginal of M1; and the marginal of M2. 
#' @param mediator Which mediator are you interested in extracting the marginal of? Possible
#' values are \code{"M1"} and \code{"M2"}. 
#' @param M_i The particular value at which you are interested in evaluating the marginal density.
#' Note that the general idea is to \code{mapply} over this function to extract the marginal
#' under each observed value of \code{C} and extract the marginal distribution at the corresponding
#' value of e.g., \code{M1}. 
#' @param unique_M_values The unique values that the selected \code{mediator} can assume. Recall that
#' the ordering of the marginal entries in \code{Q_M_n_i} is in this order. 
#' @param a_val Which treatment are you interested in extracting the marginal under? Possible values
#' are \code{"a_star"} and \code{"a"}
#' @param ... Not used

extract_marginal <- function(Q_M_n_i, mediator = "M1", M_i, unique_M_values, 
                             a_val = "a_star", ...){
  # which entry of the marginal to grab
  M_idx <- which(unique_M_values %in% M_i)
  # under which treatment
  a_idx <- ifelse(a_val == "a_star", 1, 2)
  # and which marginal
  marg_idx <- ifelse(mediator == "M1", 2, 3)
  return(Q_M_n_i[[a_idx]][[marg_idx]][M_idx])
}

#' Helper function for extracting joint distributions at particular values of 
#' the mediator. 
#' @param Q_M_n_i An entry in the \code{Q_M_n} list, where the first entry in the list
#' corresponds to mediator distributions under \code{a_star}; the second entry in the list
#' corresponds to mediator distributions under \code{a}. In each of these lists, we have
#' the joint of M1,M2; the marginal of M1; and the marginal of M2. 
#' @param all_mediator_values All possible combinations of \code{M1} and \code{M2}. Used to identify
#' which entry in \code{Q_M_n_i} that corresponds to values specified in \code{M1_i} and \code{M2_i}.
#' @param M1_i The particular value of \code{M1} at which you are interested in evaluating the joint density.
#' Note that the general idea is to \code{mapply} over this function to extract the joint
#' under each observed value of \code{C} and extract the marginal distribution at the corresponding
#' value of \code{M1} and \code{M2}.
#' @param M1_i The particular value of \code{M1} at which you are interested in evaluating the joint density.
#' Note that the general idea is to \code{mapply} over this function to extract the joint
#' under each observed value of \code{C} and extract the marginal distribution at the corresponding
#' value of \code{M1} and \code{M2}.
#' @param a_val Which treatment are you interested in extracting the marginal under? Possible values
#' are \code{"a_star"} and \code{"a"}
#' @param ... Not used

extract_joint <- function(Q_M_n_i, M1_i, M2_i, all_mediator_values,
                          a_val = "a_star", ...){
  # which entry of the marginal to grab
  M1_idx <- which(all_mediator_values$M1 %in% M1_i)
  M2_idx <- which(all_mediator_values$M2 %in% M2_i)
  M1_M2_idx <- intersect(M1_idx, M2_idx)

  # under which treatment
  a_idx <- ifelse(a_val == "a_star", 1, 2)
  
  return(Q_M_n_i[[a_idx]][[1]][M1_M2_idx])
}

#' Helper function for extracting Qbar at the observed values of confounders and
#' mediators under \code{a} or \code{a_star}
extract_Qbar_obs <- function(Qbar_n_i, M1_i, M2_i, all_mediator_values, 
                         a_val = "a_star"){
  a_idx <- ifelse(a_val == "a_star", 1, 2)
  return(Qbar_n_i$Qbar_a_0[[a_idx]][Qbar_n_i$which_M1_M2_obs])
}

#' Helper to merge covariates in targeting
#' 
#' @param ... Arguments passed in 
reduce_merge <- function(...){ merge(..., all.x = TRUE, sort = FALSE) } 
