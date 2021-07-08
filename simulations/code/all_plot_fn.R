# check for truth in CIs
in_ci <- function(truth, ci){
  truth >= min(ci) & truth <= max(ci)
}

get_bias_mr <- function(out, which_eff = "indirectM1"){
    by(out, out$n, function(x){
      tmle1 <- mean(x[ , grep(paste0(which_eff, "_tmle1_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      aiptw1 <- mean(x[ , grep(paste0(which_eff, "_aiptw1_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      tmle2 <- mean(x[ , grep(paste0(which_eff, "_tmle2_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      aiptw2 <- mean(x[ , grep(paste0(which_eff, "_aiptw2_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      tmle3 <- mean(x[ , grep(paste0(which_eff, "_tmle3_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      aiptw3 <- mean(x[ , grep(paste0(which_eff, "_aiptw3_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      tmle4 <- mean(x[ , grep(paste0(which_eff, "_tmle4_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      aiptw4 <- mean(x[ , grep(paste0(which_eff, "_aiptw4_est"), colnames(x))] - x[,paste0("truth_", which_eff)], na.rm = TRUE)
      return(c(tmle1, aiptw1,tmle2, aiptw2,tmle3, aiptw3,tmle4, aiptw4))
    })
  }

get_sd_mr <- function(out, which_eff = "indirectM1"){
  by(out, out$n, function(x){
    tmle1 <- sd(x[ , grep(paste0(which_eff, "_tmle1_est"), colnames(x))], na.rm = TRUE)
    aiptw1 <- sd(x[ , grep(paste0(which_eff, "_aiptw1_est"), colnames(x))] , na.rm = TRUE)
    tmle2 <- sd(x[ , grep(paste0(which_eff, "_tmle2_est"), colnames(x))], na.rm = TRUE)
    aiptw2 <- sd(x[ , grep(paste0(which_eff, "_aiptw2_est"), colnames(x))] , na.rm = TRUE)
    tmle3 <- sd(x[ , grep(paste0(which_eff, "_tmle3_est"), colnames(x))], na.rm = TRUE)
    aiptw3 <- sd(x[ , grep(paste0(which_eff, "_aiptw3_est"), colnames(x))] , na.rm = TRUE)
    tmle4 <- sd(x[ , grep(paste0(which_eff, "_tmle4_est"), colnames(x))], na.rm = TRUE)
    aiptw4 <- sd(x[ , grep(paste0(which_eff, "_aiptw4_est"), colnames(x))] , na.rm = TRUE)
    return(c(tmle1, aiptw1,tmle2, aiptw2,tmle3, aiptw3,tmle4, aiptw4))
  })
}
get_mse_mr <- function(out, which_eff = "indirectM1"){
  by(out, out$n, function(x){
    tmle1 <- mean((x[ , grep(paste0(which_eff, "_tmle1_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    aiptw1 <- mean((x[ , grep(paste0(which_eff, "_aiptw1_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    tmle2 <- mean((x[ , grep(paste0(which_eff, "_tmle2_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    aiptw2 <- mean((x[ , grep(paste0(which_eff, "_aiptw2_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    tmle3 <- mean((x[ , grep(paste0(which_eff, "_tmle3_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    aiptw3 <- mean((x[ , grep(paste0(which_eff, "_aiptw3_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    tmle4 <- mean((x[ , grep(paste0(which_eff, "_tmle4_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    aiptw4 <- mean((x[ , grep(paste0(which_eff, "_aiptw4_est"), colnames(x))] - x[,paste0("truth_", which_eff)])^2, na.rm = TRUE)
    return(c(tmle1, aiptw1,tmle2, aiptw2,tmle3, aiptw3,tmle4, aiptw4))
  })
}

format_out_mr <- function(out, summary_fn,
                       all_eff = c("indirectM1")){
  if(!grepl("get_density", summary_fn)){    
    all_out <- sapply(all_eff, function(x){
      suppressWarnings(a <- data.frame(Reduce(rbind, 
                        do.call(summary_fn, 
                                args = list(out = out, 
                                            which_eff = x))),
                 n = c(250, 500, 1000, 2000),
                 stringsAsFactors = FALSE))
      colnames(a) <- c("tmle1","os1","tmle2","os2","tmle3","os3","tmle4","os4")
      row.names(a) <- NULL
      a
    }, simplify = FALSE)
  }else{
    all_out <- sapply(all_eff, function(x){
      do.call(summary_fn, args = list(out = out, which_eff = x))
    })
  }
  return(all_out)
}


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


blank_plot <- function(...){
  plot(1e-10, 1e-10, pch = "", bty = "n", xaxt = "n", ...)
}




plot_one_est_row <- function(which_eff = "total",
              all_bias, all_mse, all_sd, 
              all_cover, all_cover_orc, all_se,
              tmle_name = "tmle",
              os_name = "os",
              bias_ylim = NULL){
  # six panel plot, two rows
  # top row = bias, sd, MSE
  # bottom row = samp dist. TMLE, samp. dist AIPTW, coverage ()
  # bias plot
  bias_range <- range(c(all_bias[[which_eff]][[tmle_name]], all_bias[[which_eff]][[os_name]]))
  max_bias_val <- max(abs(bias_range))

  if(!is.null(bias_ylim)){
    yl <- c(-1.05 * max_bias_val, 1.05 * max_bias_val)
  }
  blank_plot(xlim = c(0.5, 4.5), ylim = bias_ylim,
             xlab = "n", ylab = "Bias")
  axis(side = 1, at = 1:4, labels = all_bias[[which_eff]]$n)
  abline(h = 0, lty = 3)
  points(y = all_bias[[which_eff]][[tmle_name]], x = 1:4, type = "b")
  points(y = all_bias[[which_eff]][[os_name]], x = 1:4, type = "b", pch = 2)

  # sd plot
  sd_range <- range(c(all_sd[[which_eff]][[tmle_name]], all_sd[[which_eff]][[os_name]]))
  max_sd_val <- max(abs(sd_range))
  blank_plot(xlim = c(0.5, 4.5), ylim = c(0, max_sd_val * 1.1),
             xlab = "n", ylab = "Standard deviation")
  axis(side = 1, at = 1:4, labels = all_sd[[which_eff]]$n)
  # abline(h = 0, lty = 3)
  points(y = all_sd[[which_eff]][[tmle_name]], x = 1:4, type = "b")
  points(y = all_sd[[which_eff]][[os_name]], x = 1:4, type = "b", pch = 2)

  # mean squared error plot
  mse_range <- range(c(all_mse[[which_eff]][[tmle_name]], all_mse[[which_eff]][[os_name]]))
  max_mse_val <- max(abs(mse_range))
  min_mse_val <- min(abs(mse_range))
  blank_plot(xlim = c(0.5, 4.5), ylim = c(min_mse_val * 0.9, max_mse_val * 1.1), log = "y", 
             xlab = "n", ylab = "Mean squared error")
  axis(side = 1, at = 1:4, labels = all_mse[[which_eff]]$n)
  abline(h = 0, lty = 3)
  points(y = all_mse[[which_eff]][[tmle_name]], x = 1:4, type = "b")
  points(y = all_mse[[which_eff]][[os_name]], x = 1:4, type = "b", pch = 2)
}


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



# heat map for mediator distributions
# 4x4 tile of 6x6 heat map
# cols = A = 0/1
# rows = C 0.1, 0.9  
make_heat_plot <- function(
                      M1_success = "plogis(-1.1 + 0.45 * C$C1 + 0.125 * A)",
                      M2_success = "plogis(-1.1 + 0.15 * C$C1 + 0.2 * C$C2 - 0.2 * A)",
                      M1M2_threshold = 6, m_breaks = seq(0, 0.1, length = 11)
              ){
  library(ggplot2)
  C <- data.frame(C1 = c(0.1, 0.9, 0.1, 0.9), C2 = c(0.1, 0.9, 0.1, 0.9))
  A <- c(0, 0, 1, 1)
  n <- length(A)
  success_prob_M1 <- eval(parse(text = M1_success))
  success_prob_M2 <- eval(parse(text = M2_success))
  p_M1 <- sapply(0:(M1M2_threshold-1), dgeom, prob = success_prob_M1, simplify = FALSE)
  p_M1_mat <- Reduce(cbind, p_M1)
  p_M1_mat <- cbind(p_M1_mat, 1 - apply(p_M1_mat, 1, sum))
  p_M2 <- sapply(0:(M1M2_threshold-1), dgeom, prob = success_prob_M2, simplify = FALSE)
  p_M2_mat <- Reduce(cbind, p_M2)
  p_M2_mat <- cbind(p_M2_mat, 1 - apply(p_M2_mat, 1, sum))

  tmp2 <- expand.grid(M1 = 0:M1M2_threshold, M2 = 0:M1M2_threshold)
  rslt <- vector(mode = "list", length = 4)
  for(i in 1:4){
    tmp <- expand.grid(M1 = p_M1_mat[i,], M2 = p_M2_mat[i,])
    rslt[[i]] <- cbind(tmp2, tmp[,1] * tmp[,2])
    rslt[[i]]$C <- C$C1[i]
    rslt[[i]]$A <- A[i]
    colnames(rslt[[i]])[3] <- "p_M1M2"
  }
  plot_data <- Reduce(rbind, rslt)  

  p <- ggplot(plot_data, aes(M1, M2, fill= p_M1M2)) + 
        geom_tile() +
        facet_grid(C ~ A, labeller = label_bquote(rows = C[1]*" = "*C[2]*" = "*.(C),
                                                  cols = "A = "*.(A))) + 
        scale_fill_gradient(low="white", high="black", 
                            breaks = m_breaks) + 
        theme_bw() + 
        labs(fill = expression(q["A, "*M[1]*", "*M[2]]))
  return(p)
}
