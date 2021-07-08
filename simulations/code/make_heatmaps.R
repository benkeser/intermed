renv::activate()

source(here::here("code", "all_plot_fn.R"))

pdf(here::here("figs", "mediator_dist_discrete2.pdf"), height = 6, width = 7)
make_heat_plot(M1_success = "plogis(-1 + 0.45 * C$C1 + 0.125 * A)",
               M2_success = "plogis(-1 + 0.15 * C$C1 + 0.2 * C$C2 - 0.2 * A)",
               M1M2_threshold = 3, m_breaks = seq(0, 0.2, by = 0.025))
dev.off()

pdf(here::here("figs", "mediator_dist_discrete3.pdf"), height = 6, width = 7)
make_heat_plot(M1_success = "plogis(-1 + 0.45 * C$C1 + 0.125 * A)",
               M2_success = "plogis(-1 + 0.15 * C$C1 + 0.2 * C$C2 - 0.2 * A)",
               M1M2_threshold = 12, m_breaks = seq(0, 0.07, by = 0.01))
dev.off()

pdf(here::here("figs", "mediator_dist_discrete.pdf"), height = 6, width = 7)
make_heat_plot(M1_success = "plogis(-1 + 0.45 * C$C1 + 0.125 * A)",
               M2_success = "plogis(-1 + 0.15 * C$C1 + 0.2 * C$C2 - 0.2 * A)",
               M1M2_threshold = 6, m_breaks = seq(0, 0.2, by = 0.025))
dev.off()

pdf(here::here("figs", "mediator_dist_robustness.pdf"), height = 6, width = 7)
make_heat_plot(M1_success = "plogis(-1.1 + 0.25 * C$C1 - 0.5 * A * C$C1)",
               M2_success = "plogis(-1.1 + 0.25 * C$C2 - 0.5 * C$C2 * A)",
               M1M2_threshold = 6, m_breaks = seq(0, 0.1, by = 0.01))
dev.off()

