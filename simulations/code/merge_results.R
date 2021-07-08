renv::activate()

load(here::here("output", "parm.RData"))
n_sim <- nrow(parm)

rslt_discrete <- matrix(NA, nrow = sum(parm$simulation == "discrete"), ncol = 97)
rslt_discrete2 <- matrix(NA, nrow = sum(parm$simulation == "discrete2"), ncol = 97)
rslt_discrete3 <- matrix(NA, nrow = sum(parm$simulation == "discrete3"), ncol = 97)
rslt_robustness <- matrix(NA, nrow = sum(parm$simulation == "discrete3"), ncol = 39)
rslt_continuous <- matrix(NA, nrow = sum(parm$simulation == "discrete3"), ncol = 97)

ct_discrete <- ct_discrete2 <- ct_discrete3 <- ct_continuous <- ct_robustness <- 0
for(i in seq_len(n_sim)){
	load(here::here("output", paste0("sim_rslt_", i, ".RData")))
	if(parm$simulation[i] == "discrete"){
		ct_discrete <- ct_discrete + 1
		rlst_discrete[ct_discrete, ] <- out
	}else if(parm$simulation[i] == "discrete2"){
		ct_discrete2 <- ct_discrete2 + 1
		rlst_discrete2[ct_discrete2, ] <- out
	}else if(parm$simulation[i] == "discrete3"){
		ct_discrete3 <- ct_discrete3 + 1
		rlst_discrete3[ct_discrete3, ] <- out
	}else if(parm$simulation[i] == "robustness"){
		ct_robustness <- ct_robustness + 1
		rslt_robustness[ct_robustness, ] <- out
	}else if(parm$simulation[i] == "continuous"){
		ct_continuous <- ct_continuous + 1
		rslt_continuous[ct_continuous, ] <- out
	}
}

source(here::here("code", "output_column_names.R"))

out_discrete <- data.frame(rslt_discrete); colnames(out_discrete) <- colnames_discrete
out_discrete2 <- data.frame(rslt_discrete2); colnames(out_discrete2) <- colnames_discrete
out_discrete3 <- data.frame(rslt_discrete3); colnames(out_discrete3) <- colnames_discrete
# same column names as discrete sim
out_continuous <- data.frame(rslt_continuous); colnames(out_continuous) <- colnames_discrete
out_robustness <- data.frame(rslt_robustness); colnames(out_robustness) <- colnames_robustness

save(out_discrete, file = here::here("output", "out_discrete.RData"))
save(out_discrete2, file = here::here("output", "out_discrete2.RData"))
save(out_discrete3, file = here::here("output", "out_discrete3.RData"))
save(out_continuous, file = here::here("output", "out_continuous.RData"))
save(out_robustness, file = here::here("output", "out_robustness.RData"))