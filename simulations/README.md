
# R/`intermed`/`simulations

> Simulation studies for manuscript "Nonparametric estimators of mediation
effects with multiple mediators"

**Author:** [David Benkeser](https://www.davidbphd.com/) and Jialu Ran

------------------------------------------------------------------------

## Description

This repository houses simulation code for replicating simulation results. There
are five simulations. Each simulation is run at four sample sizes on one
thousand simulated data sets.

The code is organized using a `Makefile` with package management via the `renv`
package. 


The general work flow for simulations is as follows.

1. load project environment

Install the `renv` package in `R` (e.g., via `install.packages`). From the
`simulations` folder, synchronize the package environment as follows.

```bash
Rscript -e "renv::restore()"
```

This command will install packages from the project lockfile.

2. `make` the simulation parameters

```bash
make output/parm.RData
```

This creates the `.RData` file that defines the parameters for simulation jobs.

3. `make` the simulation results

```bash
make simulation_results
```

This command will only work on a `slurm` cluster -- it submits 20k jobs using
`sbatch`. Some modification to the `Make` recipe may be needed depending on
cluster configuration. 

If a cluster is not available, then it should be straightforward to follow the
logic of the code to reproduce the results using some other form of high
throughput computing.

4. merge the simulation results

The `simulation_results` recipe generates 20k output files that need to be
merged.

```bash 
make merge_results
```

5. `make` the figures

```bash 
make figures
```

This makes all figures presented in the main text and supplement.

6. `make` the heatmaps illustrating mediator distributions

```bash
make heatmaps
```

------------------------------------------------------------------------

## G-computation simulation

The code for executing the G-computation simulation included in the supplement
is included in [a separate repository](https://github.com/Scarlett422301/intermed_sims).
