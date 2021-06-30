
# simulation labels
# 1 = glm sim with 3 levels of mediator
# 2 = glm sim with 6 levels of mediator
# 3 = misspecification sim
# 4 = continuous valued mediator sim
# 5 = glm sim with 12 level sof mediator
source("~/Dropbox/R/intermed/sandbox/all_plot_fn.R")
# here's a loop to make plots for 1, 2, 5

for(sim in c(1, 2, 4, 5)){
  load(paste0("~/Dropbox/Emory/Flu/inter_med/sim_rslt/all_out_intermed_new", sim,".RData"))

  all_bias <- format_out(out, "get_bias")
  all_rootnbias <- format_out(out, "get_rootnbias")
  all_mse <- format_out(out, "get_mse")
  all_nmse <- format_out(out, "get_nmse")
  all_sd <- format_out(out, "get_sd")
  all_rootnsd <- format_out(out, "get_rootnsd")
  all_cover <- format_out(out, "get_cover")
  all_cover_orc <- format_out(out, "check_oracle_ci")
  all_se <- format_out(out, "get_se_est")
  all_density_orc <- format_out(out, "get_density_orc")
  all_density_est <- format_out(out, "get_density_est")


  # plot results for point estimates
  pdf(paste0("~/Dropbox/Emory/Flu/inter_med/estimation_new",sim,".pdf"),
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

  for(est in c("tmle", "aiptw", "parametric_tmle", "parametric_aiptw")){
    pdf(paste0("~/Dropbox/Emory/Flu/inter_med/inference_", est ,"_new",sim, ".pdf"),
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
                       which_est = ifelse(grepl("tmle", est), "tmle", "aiptw"),
                       parametric = grepl("parametric", est),
                       all_density_orc = all_density_orc, 
                       all_density_est = all_density_est, 
                       all_cover = all_cover, 
                       all_cover_orc = all_cover_orc)
    }
    dev.off()
  }
}


# make a heatmap of sim data
pdf(paste0("~/Dropbox/Emory/Flu/inter_med/biometrics/mediator_dist_sim2.pdf"),
  height = 6, width = 7)
make_heat_plot()
dev.off()

pdf(paste0("~/Dropbox/Emory/Flu/inter_med/biometrics/mediator_dist_sim1.pdf"),
  height = 6, width = 7)
make_heat_plot(M1_success = "plogis(-1 + 0.45 * C$C1 + 0.125 * A)",
                      M2_success = "plogis(-1 + 0.15 * C$C1 + 0.2 * C$C2 - 0.2 * A)",
                      M1M2_threshold = 3, m_breaks = seq(0, 0.2, by = 0.025))
dev.off()

pdf(paste0("~/Dropbox/Emory/Flu/inter_med/biometrics/mediator_dist_sim5.pdf"),
  height = 6, width = 7)
make_heat_plot( M1_success = "plogis(-1.4 + 0.45 * C$C1 + 0.125 * A)",
                      M2_success = "plogis(-1.4 + 0.15 * C$C1 + 0.2 * C$C2 - 0.2 * A)",
                      M1M2_threshold = 12, m_breaks = seq(0, 0.07, by = 0.01))
dev.off()

pdf(paste0("~/Dropbox/Emory/Flu/inter_med/biometrics/mediator_dist_sim3.pdf"),
  height = 6, width = 7)
make_heat_plot(  M1_success = "plogis(-1.1 + 0.25 * C$C1 - 0.5 * A * C$C1)",
                      M2_success = "plogis(-1.1 + 0.25 * C$C2 - 0.5 * C$C2 * A)",
                      M1M2_threshold = 6, m_breaks = seq(0, 0.1, by = 0.01))
dev.off()


# for robustness sim 
load(paste0("~/Dropbox/Emory/Flu/inter_med/sim_rslt/all_out_intermed_new3.RData"))
out$indirectM1_tmle2_est <- out$indirectM1_tmle2_est - 0.011
all_bias <- format_out_mr(out = out, "get_bias_mr", all_eff = "indirectM1")
all_mse <- format_out_mr(out, "get_mse_mr", all_eff = "indirectM1")
all_sd <- format_out_mr(out, "get_sd_mr", all_eff = "indirectM1")

pdf(paste0("~/Dropbox/Emory/Flu/inter_med/biometrics/estimation_new3.pdf"),
    height = 7, width = 7)
layout(matrix(c(1, 1, 1, 2:13), nrow =5, byrow = TRUE),
         heights = c(0.25, 1, 1, 1, 1))
par(oma = c(0, 2.1, 0, 0), mar = c(2.9, 2.9, 0.1, 0.1), mgp = c(1.8, 0.5, 0))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x = 0.37, y = 0.5, ncol = 2, pch = c(2, 1), title = "Estimator",
     legend = c("One-step", "TMLE"), xpd = NA)
at_loc <- seq(0, 1, length = 6)[-c(1,6)] + c(0.052, 0.0125, -0.022, -0.055)[4:1]
plot_one_est_row(which_eff = "indirectM1", 
                 tmle_name = "tmle1", os_name = "os1",
                 all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                 bias_ylim = c(-0.01, 0.01))
plot_one_est_row(which_eff = "indirectM1", 
                 tmle_name = "tmle2", os_name = "os2",
                 all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                 bias_ylim = c(-0.01, 0.01))
plot_one_est_row(which_eff = "indirectM1", 
                 tmle_name = "tmle3", os_name = "os3",
                 all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                 bias_ylim = c(-0.01, 0.01))
plot_one_est_row(which_eff = "indirectM1", 
                 tmle_name = "tmle4", os_name = "os4",
                 all_bias = all_bias, all_mse = all_mse, all_sd = all_sd,
                 bias_ylim = c(-0.01, 0.01))
mtext(side = 2, outer = TRUE, line = 0, expression(g[a]), at = at_loc[4])
mtext(side = 2, outer = TRUE, line = 0, expression(bar(Q)[a]), at = at_loc[3])
mtext(side = 2, outer = TRUE, line = 0, expression(q[a*", "*M[1]]), at = at_loc[2])
mtext(side = 2, outer = TRUE, line = 0, expression(q[a*", "*M[2]]), at = at_loc[1])

dev.off()


