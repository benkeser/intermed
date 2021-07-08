args <- commandArgs(TRUE)
simulation <- args[[1]]

renv::activate()
out <- get(load(here::here("output", paste0(simulation, ".RData"))))

source(here::here("code", "all_plot_fn.R"))

if(simulation != "robustness"){
  all_bias <- format_out(out, "get_bias")
  all_mse <- format_out(out, "get_mse")
  all_sd <- format_out(out, "get_sd")
  all_cover <- format_out(out, "get_cover")
  all_cover_orc <- format_out(out, "check_oracle_ci")
  all_density_orc <- format_out(out, "get_density_orc")
  all_density_est <- format_out(out, "get_density_est")
}else{
  all_bias <- format_out_mr(out, "get_bias")
  all_mse <- format_out_mr(out, "get_mse")
  all_sd <- format_out_mr(out, "get_sd")
  all_cover <- format_out_mr(out, "get_cover")
  all_cover_orc <- format_out_mr(out, "check_oracle_ci")
  all_density_orc <- format_out_mr(out, "get_density_orc")
  all_density_est <- format_out_mr(out, "get_density_est")
}

if(simulation == "robustness"){
  # plot results for point estimates
  pdf(here::here("figs", paste0("estimation_", simulation, ".pdf")),
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
}else{
  # plot results for point estimates
  pdf(here::here("figs", paste0("estimation_", simulation, ".pdf")),
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

  pdf(here::here("figs", paste0("inference_tmle_", simulation, ".pdf")),
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

  pdf(here::here("figs", paste0("inference_aiptw_", simulation, ".pdf")),  
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

}

