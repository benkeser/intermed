#! /usr/bin/env Rscript

renv::activate("..")

# get command lines arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Not enough arguments. Please use args 'listsize', 'prepare', 'run <itemsize>' or 'merge'")
}

# load parameters
load(here::here("parm.RData"))

# save directory
save_dir <- here::here("output")
code_dir <- here::here("code")

source(paste0(code_dir, "make_data.R"))
source(paste0(code_dir, "run_simulation.R"))
source(paste0(code_dir, "all_plot_fn.R"))

# load packages
library(SuperLearner)
library(intermed)

i <- as.numeric(args[1])
j <- NULL
if(parm$simulation[i] == "discrete2") j <- 2
if(parm$simulation[i] == "discrete3") j <- 3

out <- do.call(paste0('do_one_', parm$simulation[i]), 
               args = list(i = i, parm = parm, j = j))

save(out, file = here::here("output", paste0("sim_rslt_", i, ".RData")))