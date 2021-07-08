output/parm.RData: code/make_parameters.R
	Rscript code/make_parameters.R

simulation_results: code/run_simulation.R output/parm.RData
	NUM_SIMS=$(Rscript -e "renv::activate(); load(here::here('output','parm.RData')); cat(nrow(parm))")
	for i in $(seq 1 $NUM_SIMS); do \
 	    sbatch code/run_simulation.R $i	\
	done

merge_results: code/merge_results.R
	Rscript code/merge_results.R

heatmaps: code/make_heatmaps.R
	Rscript code/make_heatmaps.R

figures:
	chmod +x code/make_figures.R
	for SIM in discrete discrete2 discrete3 continuous robustness; do \
		code/make_figures.R ${SIM} \
	done